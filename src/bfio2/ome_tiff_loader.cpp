#include "ome_tiff_loader.h"


OmeTiffLoader::OmeTiffLoader(const std::string &fNameWithPath) : 
	xml_metadata_ptr(nullptr),
	gsTiffTileLoader(nullptr),
	fName(fNameWithPath),
	nThreads(1)
{
    if (checkTileStatus())
		{
			gsTiffTileLoader = std::make_unique<GrayscaleTiffTileLoader<uint32_t>>(nThreads, fName);
		} 
		else 
		{
			// since the file is not tiled, we provide the tile dimensions
			auto [tw, th, td]  = calculateTileDimensions(); //vector of (tw, th, td)
            gsTiffTileLoader = std::make_unique<GrayscaleTiffStripLoader<uint32_t>>(nThreads, fName, tw, th, td);
		}
};

OmeTiffLoader::~OmeTiffLoader(){
	gsTiffTileLoader = nullptr;
	xml_metadata_ptr = nullptr;	
};

size_t OmeTiffLoader::getRowTileCount() const {return gsTiffTileLoader->numberTileHeight();}
size_t OmeTiffLoader::getColumnTileCount() const {return gsTiffTileLoader->numberTileWidth();}
size_t OmeTiffLoader::getImageHeight() const {return gsTiffTileLoader->fullHeight(0);}
size_t OmeTiffLoader::getImageWidth() const {return gsTiffTileLoader->fullWidth(0);}
size_t OmeTiffLoader::getTileHeight() const {return gsTiffTileLoader->tileHeight(0);}
size_t OmeTiffLoader::getTileWidth() const {return gsTiffTileLoader->tileWidth(0);}

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::getTileData(size_t const indexRowGlobalTile, size_t const indexColGlobalTile)
{
    auto tw = gsTiffTileLoader->tileWidth(0);
    auto th = gsTiffTileLoader->tileHeight(0);
    std::shared_ptr<std::vector<uint32_t>> tileData = std::make_shared<std::vector<uint32_t>>(tw * th);
    gsTiffTileLoader->loadTileFromFile(tileData, indexRowGlobalTile, indexColGlobalTile, 0, 0);
    return tileData;
}

std::tuple<uint32_t, uint32_t, uint32_t>  OmeTiffLoader::getImageDimensions  () const
{
	TIFF *tiff_ = TIFFOpen(fName.c_str(), "r");
	if (tiff_ != nullptr) 
	{
		uint32_t w, l, ndirs;
		TIFFGetField(tiff_, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tiff_, TIFFTAG_IMAGELENGTH, &l);
		ndirs = TIFFNumberOfDirectories(tiff_);
		TIFFClose(tiff_);
		return {w, l, ndirs};	
	} 
	else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }
}

std::tuple<uint32_t, uint32_t, uint32_t>  OmeTiffLoader::calculateTileDimensions() const
{
	auto [w, h, d] = getImageDimensions();
	uint32_t defaultWidthSize = 1024;
	uint32_t defaultHeightSize = 1024;
	uint32_t defaultDepthSize = 1;
	w = std::min ({ w, defaultWidthSize });
	h = std::min ({ h, defaultHeightSize });
	d = std::min ({ d, defaultDepthSize });
	return {w, h, d};
}

bool OmeTiffLoader::checkTileStatus() const
{
	TIFF *tiff_ = TIFFOpen(fName.c_str(), "r");
	if (tiff_ != nullptr) 
	{
		if (TIFFIsTiled(tiff_) == 0) 
		{ 
			TIFFClose(tiff_);
			return false;
			} else 
			{
			TIFFClose(tiff_);
			return true;
			}
	} else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }
}

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::getTileData(size_t const indexGlobalTile)
{
	size_t columnTileCount = getColumnTileCount();	
	size_t indexRowGlobalTile = indexGlobalTile/columnTileCount;
	size_t indexColGlobalTile = indexGlobalTile%columnTileCount;
    auto tw = gsTiffTileLoader->tileWidth(0);
    auto th = gsTiffTileLoader->tileHeight(0);
    std::shared_ptr<std::vector<uint32_t>> tileData = std::make_shared<std::vector<uint32_t>>(tw * th);
    gsTiffTileLoader->loadTileFromFile(tileData, indexRowGlobalTile, indexColGlobalTile, 0, 0);
    return tileData;
}

std::pair<size_t, size_t> OmeTiffLoader::getTileContainingPixel(size_t const indexRowPixel, size_t const indexColPixel) const
{
	size_t th = getTileHeight();	
	size_t tw = getTileWidth();
	size_t indexRowGlobalTile = indexRowPixel/th;
	size_t indexColGlobalTile = indexColPixel/tw;
	return std::make_pair(indexRowGlobalTile, indexColGlobalTile);
}

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::getBoundingBoxVirtualTileData(size_t const indexRowPixelMin, size_t const indexRowPixelMax,
                                                                    size_t const indexColPixelMin, size_t const indexColPixelMax)
{

	// Convention 
	// rows are X coordinate (increasing from top to bottom)
	// cols are Y coordinate (increasing from left to right)
	// we need to transform from Local Tile Coordinate to Global Pixel Coordiate to Virtual Tile Coordinate

	auto ih = gsTiffTileLoader->fullHeight(0);
	auto iw = gsTiffTileLoader->fullWidth(0);
	auto indexTrueRowPixelMax = indexRowPixelMax > ih ? ih : indexRowPixelMax;
	auto indexTrueColPixelMax = indexColPixelMax > iw ? iw : indexColPixelMax;
	auto topLeftTile = getTileContainingPixel(indexRowPixelMin, indexColPixelMin);
	auto bottomRightTile = getTileContainingPixel(indexTrueRowPixelMax, indexTrueColPixelMax);
	auto minRowIndex = topLeftTile.first;
	auto minColIndex = topLeftTile.second;
	auto maxRowIndex = bottomRightTile.first;
	auto maxColIndex = bottomRightTile.second;

	// now loop through each tile, get tile data, fill the virtual tile vector

	auto tw = gsTiffTileLoader->tileWidth(0);
	auto th = gsTiffTileLoader->tileHeight(0);
	std::shared_ptr<std::vector<uint32_t>> tileData = std::make_shared<std::vector<uint32_t>>(tw * th);

	auto vtw = indexTrueColPixelMax-indexColPixelMin+1;
	auto vth = indexTrueRowPixelMax-indexRowPixelMin+1;
	std::shared_ptr<std::vector<uint32_t>> virtualTileData = std::make_shared<std::vector<uint32_t>>(vtw * vth);

	for (int i = minRowIndex; i <= maxRowIndex; ++i)
	{
		for (int j = minColIndex; j <= maxColIndex; ++j)
		{
			gsTiffTileLoader->loadTileFromFile(tileData, i, j, 0, 0);	
			// take row slice from local tile and place it in virtual tile
			// globalX = i*th + localX;
			// virtualX = globalX - indexRowPixelMin;
			// initialGLobalY = j*tw + initialLocalY;
			// initialVirtualY = initialGLobalY - indexColPixelMin;
			size_t initialLocalX = indexRowPixelMin > i*th ? indexRowPixelMin-i*th : 0;	
			size_t endLocalX = indexTrueRowPixelMax < (i+1)*th ? indexTrueRowPixelMax-i*th: th-1;
			size_t initialLocalY = indexColPixelMin > j*tw ? indexColPixelMin-j*tw : 0;
			size_t endLocalY = indexTrueColPixelMax < (j+1)*tw ? indexTrueColPixelMax-j*tw : tw-1;
			size_t initialVirtualY = j*tw + initialLocalY - indexColPixelMin;
#pragma omp parallel
#pragma omp for
			for (size_t localX=initialLocalX; localX<=endLocalX; ++localX){
				size_t virtualX = i*th + localX - indexRowPixelMin;
				std::copy(tileData->begin()+localX*tw+initialLocalY, tileData->begin()+localX*tw+endLocalY+1,virtualTileData->begin()+virtualX*vtw+initialVirtualY);					
			}	
		}
	}
	return virtualTileData;
}

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::getBoundingBoxVirtualTileDataStrideVersion(size_t const indexRowPixelMin, size_t const indexRowPixelMax,
                                                                    size_t rowStride, size_t const indexColPixelMin, size_t const indexColPixelMax, 
                                                                    size_t colStride)
{

	// Convention 
	// rows are X coordinate (increasing from top to bottom)
	// cols are Y coordinate (increasing from left to right)
	// we need to transform from Local Tile Coordinate to Global Pixel Coordiate to Virtual Tile Coordinate

	auto ih = gsTiffTileLoader->fullHeight(0);
	auto iw = gsTiffTileLoader->fullWidth(0);
	auto indexTrueRowPixelMax = indexRowPixelMax > ih ? ih : indexRowPixelMax;
	auto indexTrueColPixelMax = indexColPixelMax > iw ? iw : indexColPixelMax;
	auto topLeftTile = getTileContainingPixel(indexRowPixelMin, indexColPixelMin);
	auto bottomRightTile = getTileContainingPixel(indexTrueRowPixelMax, indexTrueColPixelMax);
	auto minRowIndex = topLeftTile.first;
	auto minColIndex = topLeftTile.second;
	auto maxRowIndex = bottomRightTile.first;
	auto maxColIndex = bottomRightTile.second;

	// now loop through each tile, get tile data, fill the virtual tile vector

	auto tw = gsTiffTileLoader->tileWidth(0);
	auto th = gsTiffTileLoader->tileHeight(0);
	std::shared_ptr<std::vector<uint32_t>> tileDataPtr = std::make_shared<std::vector<uint32_t>>(tw * th);
	auto tileData = tileDataPtr->data();

	auto vtw = (indexTrueColPixelMax-indexColPixelMin)/colStride+1;
	auto vth = (indexTrueRowPixelMax-indexRowPixelMin)/rowStride+1;
	std::shared_ptr<std::vector<uint32_t>> virtualTileDataPtr = std::make_shared<std::vector<uint32_t>>(vtw * vth);
	auto virtualTileData = virtualTileDataPtr->data();

	for (int i = minRowIndex; i <= maxRowIndex; ++i)
	{
		for (int j = minColIndex; j <= maxColIndex; ++j)
		{
			gsTiffTileLoader->loadTileFromFile(tileDataPtr, i, j, 0, 0);	
			// take row slice from local tile and place it in virtual tile
			// globalX = i*th + localX;
			// virtualX = globalX - indexRowPixelMin;
			// initialGLobalY = j*tw + initialLocalY;
			// initialVirtualY = initialGLobalY - indexColPixelMin;
			size_t initialLocalX = indexRowPixelMin > i*th ? indexRowPixelMin-i*th : 0;	
			// adjust for row stride
			size_t initialGlobalX = i*th + initialLocalX;
			initialGlobalX = adjustStride(indexRowPixelMin, initialGlobalX, rowStride);
			initialLocalX = initialGlobalX - i*th;
			size_t endLocalX = indexTrueRowPixelMax < (i+1)*th ? indexTrueRowPixelMax-i*th: th-1;

			size_t initialLocalY = indexColPixelMin > j*tw ? indexColPixelMin-j*tw : 0;
			// adjust for col stride
			size_t initialGlobalY = j*tw + initialLocalY;
			initialGlobalY = adjustStride(indexColPixelMin, initialGlobalY, colStride);
			initialLocalY = initialGlobalY - j*tw;			
			size_t endLocalY = indexTrueColPixelMax < (j+1)*tw ? indexTrueColPixelMax-j*tw : tw-1;
			size_t initialVirtualY = j*tw + initialLocalY - indexColPixelMin;
#pragma omp parallel
#pragma omp for

			for (size_t localX=initialLocalX; localX<=endLocalX; localX=localX+rowStride)
			{
				size_t virtualX = (i*th + localX - indexRowPixelMin)/rowStride;
				if (colStride == 1) 
				{

					std::copy(tileDataPtr->begin()+localX*tw+initialLocalY, tileDataPtr->begin()+localX*tw+endLocalY+1,virtualTileDataPtr->begin()+virtualX*vtw+initialVirtualY);
				}
				else 
				{
					

					for (size_t localY=initialLocalY; localY<=endLocalY; localY=localY+colStride)
					{
						size_t virtualY = (j*tw + localY - indexColPixelMin)/colStride;
						size_t localDataIndex = localX*tw+localY;
						size_t virtualDataIndex = virtualX*vtw+virtualY;
						virtualTileData[virtualDataIndex] = tileData[localDataIndex];						
					}
				}
			}	
		}
	}
	return virtualTileDataPtr;
}

std::string OmeTiffLoader::get_metadata_value(const std::string &metadata_key) const
{
	std::string value = "" ;
	if (xml_metadata_ptr == nullptr){
		parse_metadata();		
	}
	try {
		value = xml_metadata_ptr->at(metadata_key);
	}
	catch (const std::exception& e) {
		std::cout<<"Requested metadata key not found"<<std::endl;
	}
	return value;
}


void OmeTiffLoader::parse_metadata() const
{	
	TIFF *tiff_ = TIFFOpen(fName.c_str(), "r");
	if (tiff_ != nullptr) 
	{
		char *infobuf;
		TIFFGetField(tiff_, TIFFTAG_IMAGEDESCRIPTION , &infobuf);
	   	TIFFClose(tiff_);

		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_string(infobuf);;
		xml_metadata_ptr = std::make_shared<std::map<std::string, std::string>>();
		
		if (result){
			std::vector<pugi::xml_node> node_list;
			node_list.push_back(doc.child("OME").child("Image").child("Pixels"));
			node_list.push_back(doc.child("OME").child("Image").child("Pixels").child("TiffData"));

			for (const auto &node: node_list){
				for (const pugi::xml_attribute &attr: node.attributes()){
					xml_metadata_ptr->emplace(attr.name(), attr.value());
				}
			}
			// get channel info

			// read structured annotaion
			pugi::xml_node annotion_list = doc.child("OME").child("StructuredAnnotations");
			for(const pugi::xml_node &annotation : annotion_list){
				auto key = annotation.child("Value").child("OriginalMetadata").child("Key").child_value();
				auto value = annotation.child("Value").child("OriginalMetadata").child("Value").child_value();
				xml_metadata_ptr->emplace(key,value);
			}


		}
	
	} else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }	
}

size_t OmeTiffLoader::adjustStride (size_t startPos, size_t currentPos, size_t strideVal) const
{
	if (strideVal == 0) return currentPos; // guard against div by 0

	size_t tmp = currentPos-startPos;
	if (tmp%strideVal == 0) 
	{
		return currentPos; // no adjustment needed
	} else 
	{
		return ((tmp/strideVal)+1)*strideVal; // move to the next eligible position
	}
}

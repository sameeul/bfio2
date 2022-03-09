#include "ome_tiff_loader.h"

size_t adjustStride (size_t startPos, size_t currentPos, size_t strideVal){
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

OmeTiffLoader::OmeTiffLoader(const std::string &fNameWithPath){
	parse_metadata(fNameWithPath);

    if (checkTileStatus(fNameWithPath))
		{
			gsTiffTileLoader = std::make_unique<GrayscaleTiffTileLoader<uint32_t>>(nThreads, fNameWithPath);
		} 
		else 
		{
			// since the file is not tiled, we provide the tile dimensions
			auto tileDims = calculateTileDimensions(fNameWithPath); //vector of (tw, th, td)
            gsTiffTileLoader = std::make_unique<GrayscaleTiffStripLoader<uint32_t>>(nThreads, fNameWithPath, tileDims->at(0), 
																						tileDims->at(1), tileDims->at(2));
		}
};

OmeTiffLoader::~OmeTiffLoader(){
	gsTiffTileLoader = nullptr;	
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

std::unique_ptr<std::vector<size_t>>  OmeTiffLoader::getImageDimensions  (const std::string& filePath) const
{
	TIFF *tiff_ = TIFFOpen(filePath.c_str(), "r");
	if (tiff_ != nullptr) 
	{
		std::unique_ptr<std::vector<size_t>> imageDims = std::make_unique<std::vector<size_t>>(0);
		size_t tmp = 0;
		TIFFGetField(tiff_, TIFFTAG_IMAGEWIDTH, &tmp);
		imageDims->push_back(tmp);
      	TIFFGetField(tiff_, TIFFTAG_IMAGELENGTH, &tmp);
		imageDims->push_back(tmp);
		imageDims->push_back(TIFFNumberOfDirectories(tiff_));
	   	TIFFClose(tiff_);
	  	return std::move(imageDims);	
	} else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }
}

std::unique_ptr<std::vector<size_t>>  OmeTiffLoader::calculateTileDimensions(const std::string& filePath) const
{
	auto imageDims = getImageDimensions(filePath);
	size_t defaultWidthSize = 1024;
	size_t defaultHeightSize = 1024;
	size_t defaultDepthSize = 1;
	imageDims->at(0) = std::min({imageDims->at(0), defaultWidthSize});
	imageDims->at(1) = std::min({imageDims->at(1), defaultHeightSize});
	imageDims->at(2) = std::min({imageDims->at(2), defaultDepthSize});
	return std::move(imageDims);
}

bool OmeTiffLoader::checkTileStatus(const std::string& filePath) const
{
	TIFF *tiff_ = TIFFOpen(filePath.c_str(), "r");
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

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::getTileDataContainingPixel(size_t const indexRowPixel, size_t const indexColPixel)
{
	size_t th = getTileHeight();	
	size_t tw = getTileWidth();
	size_t indexRowGlobalTile = indexRowPixel/th;
	size_t indexColGlobalTile = indexColPixel/tw;
    std::shared_ptr<std::vector<uint32_t>> tileData = std::make_shared<std::vector<uint32_t>>(tw * th);
    gsTiffTileLoader->loadTileFromFile(tileData, indexRowGlobalTile, indexColGlobalTile, 0, 0);
    return tileData;
}


std::pair<size_t, size_t> OmeTiffLoader::getTileContainingPixel(size_t const indexRowPixel, size_t const indexColPixel)
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
	std::shared_ptr<std::vector<uint32_t>> tileData = std::make_shared<std::vector<uint32_t>>(tw * th);

	auto vtw = (indexTrueColPixelMax-indexColPixelMin)/colStride+1;
	auto vth = (indexTrueRowPixelMax-indexRowPixelMin)/rowStride+1;
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

					std::copy(tileData->begin()+localX*tw+initialLocalY, tileData->begin()+localX*tw+endLocalY+1,virtualTileData->begin()+virtualX*vtw+initialVirtualY);
				}
				else 
				{
					for (size_t localY=initialLocalY; localY<=endLocalY; localY=localY+colStride)
					{
						size_t virtualY = (j*tw + localY - indexColPixelMin)/colStride;
						size_t localDataIndex = localX*tw+localY;
						size_t virtualDataIndex = virtualX*vtw+virtualY;
						virtualTileData->at(virtualDataIndex) = tileData->at(localDataIndex);						
					}
				}
			}	
		}
	}
	return virtualTileData;
}

void OmeTiffLoader::parse_metadata(const std::string &fName){
	

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

std::shared_ptr<std::map<std::string, std::string>> OmeTiffLoader::get_xml_metadata()
{
	return xml_metadata_ptr;
}
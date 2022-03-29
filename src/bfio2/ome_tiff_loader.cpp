#include "ome_tiff_loader.h"


OmeTiffLoader::OmeTiffLoader(const std::string &fname_with_path) : 
	tile_loader_(nullptr),
	xml_metadata_ptr_(nullptr),
	fast_loader_(nullptr),
	n_threads_(8),
	fname_(fname_with_path)
{
	if (CheckTileStatus())
	{
		tile_loader_ = std::make_shared<OmeTiffGrayScaleTileLoader<uint32_t>>(n_threads_, fname_);
	}
	else 
	{
		// since the file is not tiled, we provide the tile dimensions
		auto [tw, th, td]  = CalculateTileDimensions();		
		tile_loader_ = std::make_shared<OmeTiffGrayScaleStripLoader<uint32_t>>(n_threads_, fname_, tw, th, td);
	}

    auto options = std::make_unique<fl::FastLoaderConfiguration<fl::DefaultView<uint32_t>>>(tile_loader_);
    // Set the configuration
    uint32_t radiusDepth = 0;
    uint32_t radiusHeight = 1;
    uint32_t radiusWidth = 1;

    options->radius(radiusDepth, radiusHeight, radiusWidth);
    options->ordered(true);
    options->borderCreatorConstant(0);
    options->cacheCapacity(0,100);
    options->viewAvailable(0,1);
    // Create the Fast Loader Graph	
    fast_loader_ = std::make_shared<fl::FastLoaderGraph<fl::DefaultView<uint32_t>>>(std::move(options));
    // Execute the graph
    fast_loader_->executeGraph();
};

OmeTiffLoader::~OmeTiffLoader(){

	xml_metadata_ptr_ = nullptr;	
    fast_loader_->finishRequestingViews(); 
	fast_loader_->waitForTermination();
	tile_loader_ = nullptr;
};

size_t OmeTiffLoader::GetRowTileCount() const {return tile_loader_->numberTileHeight();}
size_t OmeTiffLoader::GetColumnTileCount() const {return tile_loader_->numberTileWidth();}
size_t OmeTiffLoader::GetImageHeight() const {return tile_loader_->fullHeight(0);}
size_t OmeTiffLoader::GetImageWidth() const {return tile_loader_->fullWidth(0);}
size_t OmeTiffLoader::GetTileHeight() const {return tile_loader_->tileHeight(0);}
size_t OmeTiffLoader::GetTileWidth() const {return tile_loader_->tileWidth(0);}

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::GetTileDataByRowCol(size_t const row, size_t const col)
{
    auto tw = tile_loader_->tileWidth(0);
    auto th = tile_loader_->tileHeight(0);
    std::shared_ptr<std::vector<uint32_t>> tile_data = std::make_shared<std::vector<uint32_t>>(tw * th);
	fast_loader_->requestView(row, col, 0, 0);
	const auto &view = fast_loader_->getBlockingResult();
	auto vw = view->viewWidth();
	auto vrw = view->radiusWidth();
	auto vrh = view->radiusHeight();
	auto view_ptr = view->viewOrigin() + vrh*vw; 
	auto tile_ptr = tile_data->begin();
	for (size_t i = 0; i < th; ++i) 
	{	
		std::copy(view_ptr+i*vw+vrw, view_ptr+i*vw+vrw+tw, tile_ptr+i*tw);
	}
	view->returnToMemoryManager();
    return tile_data;
}

std::tuple<uint32_t, uint32_t, uint32_t>  OmeTiffLoader::GetImageDimensions  () const
{
	TIFF *tiff_ = TIFFOpen(fname_.c_str(), "r");
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

std::tuple<uint32_t, uint32_t, uint32_t>  OmeTiffLoader::CalculateTileDimensions() const
{
	auto [w, h, d] = GetImageDimensions();
	uint32_t default_width = 1024;
	uint32_t default_height = 1024;
	uint32_t default_depth = 1;
	w = std::min ({ w, default_width });
	h = std::min ({ h, default_height });
	d = std::min ({ d, default_depth });
	return {w, h, d};
}

bool OmeTiffLoader::CheckTileStatus() const
{
	TIFF *tiff_ = TIFFOpen(fname_.c_str(), "r");
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

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::GetTileDataByIndex(size_t const tile_index)
{
	size_t num_col_tiles = GetColumnTileCount();	
	size_t row = tile_index/num_col_tiles;
	size_t col = tile_index%num_col_tiles;
    auto tw = tile_loader_->tileWidth(0);
    auto th = tile_loader_->tileHeight(0);
    std::shared_ptr<std::vector<uint32_t>> tile_data = std::make_shared<std::vector<uint32_t>>(tw*th);
	fast_loader_->requestView(row, col, 0, 0);
	const auto &view = fast_loader_->getBlockingResult();
	auto vw = view->viewWidth();
	auto vrw = view->radiusWidth();
	auto vrh = view->radiusHeight();
	auto view_ptr = view->viewOrigin() + vrh*vw; 
	auto tile_ptr = tile_data->begin();
	for (size_t i = 0; i < th; ++i) 
	{	
		std::copy(view_ptr+i*vw+vrw, view_ptr+i*vw+vrw+tw, tile_ptr+i*tw);
	}
	view->returnToMemoryManager();

    return tile_data;
}

std::pair<size_t, size_t> OmeTiffLoader::GetTileContainingPixel(size_t const row_pixel_index, size_t const col_pixel_index) const
{
	size_t th = GetTileHeight();	
	size_t tw = GetTileWidth();
	size_t row = row_pixel_index/th;
	size_t col = col_pixel_index/tw;
	return std::make_pair(row, col);
}

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::GetBoundingBoxVirtualTileData(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const index_col_pixel_min, size_t const index_col_pixel_max)
{

	// Convention 
	// rows are X coordinate (increasing from top to bottom)
	// cols are Y coordinate (increasing from left to right)
	// we need to transform from Local Tile Coordinate to Global Pixel Coordiate to Virtual Tile Coordinate

	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);
	auto index_true_row_pixel_max = index_row_pixel_max > ih ? ih : index_row_pixel_max;
	auto index_true_col_pixel_max = index_col_pixel_max > iw ? iw : index_col_pixel_max;
	auto top_left_tile = GetTileContainingPixel(index_row_pixel_min, index_col_pixel_min);
	auto bottom_right_tile = GetTileContainingPixel(index_true_row_pixel_max, index_true_col_pixel_max);
	auto min_row_index = top_left_tile.first;
	auto min_col_index = top_left_tile.second;
	auto max_row_index = bottom_right_tile.first;
	auto max_col_index = bottom_right_tile.second;

	// now loop through each tile, get tile data, fill the virtual tile vector

	auto tw = tile_loader_->tileWidth(0);
	auto th = tile_loader_->tileHeight(0);

	auto vtw = index_true_col_pixel_max-index_col_pixel_min+1;
	auto vth = index_true_row_pixel_max-index_row_pixel_min+1;
	std::shared_ptr<std::vector<uint32_t>> virtual_tile_data_ptr = std::make_shared<std::vector<uint32_t>>(vtw * vth);

	for (int i = min_row_index; i <= max_row_index; ++i)
	{
		for (int j = min_col_index; j <= max_col_index; ++j)
		{
			fast_loader_->requestView(i, j, 0, 0);
		}
	}

	for (int i = min_row_index; i <= max_row_index; ++i)
	{
		for (int j = min_col_index; j <= max_col_index; ++j)
		{
	//		fast_loader_->requestView(i, j, 0, 0);
			const auto &view = fast_loader_->getBlockingResult();	
			// take row slice from local tile and place it in virtual tile
			// global_x = i*th + local_x;
			// virtual_x = global_x - index_row_pixel_min;
			// initial_global_y = j*tw + initial_local_y;
			// initial_virtual_y = initial_global_y - index_col_pixel_min;
			size_t initial_local_x = index_row_pixel_min > i*th ? index_row_pixel_min-i*th : 0;	
			size_t end_local_x = index_true_row_pixel_max < (i+1)*th ? index_true_row_pixel_max-i*th: th-1;
			size_t initial_local_y = index_col_pixel_min > j*tw ? index_col_pixel_min-j*tw : 0;
			size_t end_local_y = index_true_col_pixel_max < (j+1)*tw ? index_true_col_pixel_max-j*tw : tw-1;
			size_t initial_virtual_y = j*tw + initial_local_y - index_col_pixel_min;

			auto vw = view->viewWidth();
			auto vrw = view->radiusWidth();
			auto vrh = view->radiusHeight();
			auto view_ptr = view->viewOrigin() + vrh*vw; 
			auto virtual_tile_ptr = virtual_tile_data_ptr->begin();
			for (size_t local_x=initial_local_x; local_x<=end_local_x; ++local_x){
				size_t virtual_x = i*th + local_x - index_row_pixel_min;
				std::copy(view_ptr+local_x*vw+vrw+initial_local_y, view_ptr+local_x*vw+vrw+end_local_y+1,virtual_tile_ptr+virtual_x*vtw+initial_virtual_y);					
			}
			view->returnToMemoryManager();	
		}
	}
	


	return virtual_tile_data_ptr;
}

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::GetBoundingBoxVirtualTileDataStrideVersion(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const row_stride, size_t const index_col_pixel_min, size_t const index_col_pixel_max, 
                                                                    size_t const col_stride)
{

	// Convention 
	// rows are X coordinate (increasing from top to bottom)
	// cols are Y coordinate (increasing from left to right)
	// we need to transform from Local Tile Coordinate to Global Pixel Coordiate to Virtual Tile Coordinate

	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);
	auto index_true_row_pixel_max = index_row_pixel_max > ih ? ih : index_row_pixel_max;
	auto index_true_col_pixel_max = index_col_pixel_max > iw ? iw : index_col_pixel_max;
	auto top_left_tile = GetTileContainingPixel(index_row_pixel_min, index_col_pixel_min);
	auto bottom_right_tile = GetTileContainingPixel(index_true_row_pixel_max, index_true_col_pixel_max);
	auto min_row_index = top_left_tile.first;
	auto min_col_index = top_left_tile.second;
	auto max_row_index = bottom_right_tile.first;
	auto max_col_index = bottom_right_tile.second;

	// now loop through each tile, get tile data, fill the virtual tile vector

	auto tw = tile_loader_->tileWidth(0);
	auto th = tile_loader_->tileHeight(0);

	auto vtw = (index_true_col_pixel_max-index_col_pixel_min)/col_stride+1;
	auto vth = (index_true_row_pixel_max-index_row_pixel_min)/row_stride+1;
	std::shared_ptr<std::vector<uint32_t>> virtual_tile_data = std::make_shared<std::vector<uint32_t>>(vtw * vth);
	auto virtual_tile_data_ptr = virtual_tile_data->data();
	auto virtual_tile_data_begin = virtual_tile_data->begin();
	for (int i = min_row_index; i <= max_row_index; ++i)
	{
		for (int j = min_col_index; j <= max_col_index; ++j)
		{
			fast_loader_->requestView(i, j, 0, 0);
		}
	}

	for (int i = min_row_index; i <= max_row_index; ++i)
	{
		for (int j = min_col_index; j <= max_col_index; ++j)
		{
			const auto &view = fast_loader_->getBlockingResult();
			// take row slice from local tile and place it in virtual tile
			// global_x = i*th + local_x;
			// virtual_x = global_x - index_row_pixel_min;
			// initial_global_y = j*tw + initial_local_y;
			// initial_virtual_y = initial_global_y - index_col_pixel_min;
			size_t initial_local_x = index_row_pixel_min > i*th ? index_row_pixel_min-i*th : 0;	
			// adjust for row stride
			size_t initial_global_x = i*th + initial_local_x;
			initial_global_x = AdjustStride(index_row_pixel_min, initial_global_x, row_stride);
			initial_local_x = initial_global_x - i*th;
			size_t end_local_x = index_true_row_pixel_max < (i+1)*th ? index_true_row_pixel_max-i*th: th-1;

			size_t initial_local_y = index_col_pixel_min > j*tw ? index_col_pixel_min-j*tw : 0;
			// adjust for col stride
			size_t initial_global_y = j*tw + initial_local_y;
			initial_global_y = AdjustStride(index_col_pixel_min, initial_global_y, col_stride);
			initial_local_y = initial_global_y - j*tw;			
			size_t end_local_y = index_true_col_pixel_max < (j+1)*tw ? index_true_col_pixel_max-j*tw : tw-1;
			size_t initial_virtual_y = j*tw + initial_local_y - index_col_pixel_min;


			auto vw = view->viewWidth();
			auto vrw = view->radiusWidth();
			auto vrh = view->radiusHeight();
			auto view_ptr = view->viewOrigin() + vrh*vw; 

			for (size_t local_x=initial_local_x; local_x<=end_local_x; local_x=local_x+row_stride)
			{
				size_t virtual_x = (i*th + local_x - index_row_pixel_min)/row_stride;
				if (col_stride == 1) 
				{

					std::copy(view_ptr+local_x*vw+vrw+initial_local_y, view_ptr+local_x*vw+vrw+end_local_y+1,virtual_tile_data_begin+virtual_x*vtw+initial_virtual_y);
				}
				else 
				{
					

					for (size_t local_y=initial_local_y; local_y<=end_local_y; local_y=local_y+col_stride)
					{
						size_t virtual_y = (j*tw + local_y - index_col_pixel_min)/col_stride;
						// local_data_index = view_ptr+local_x*vw+vrw+local_y;
						// virtual_data_index = virtual_x*vtw+virtual_y;
						virtual_tile_data_ptr[virtual_x*vtw+virtual_y] = *(view_ptr+local_x*vw+vrw+local_y);						
					}
				}
			}
			view->returnToMemoryManager();	
		}
	}
	return virtual_tile_data;
}

std::string OmeTiffLoader::GetMetaDataValue(const std::string &metadata_key) const
{
	std::string value = "" ;
	if (xml_metadata_ptr_ == nullptr){
		ParseMetadata();		
	}
	try {
		value = xml_metadata_ptr_->at(metadata_key);
	}
	catch (const std::exception& e) {
		std::cout<<"Requested metadata key not found"<<std::endl;
	}
	return value;
}


void OmeTiffLoader::ParseMetadata() const
{	
	TIFF *tiff_ = TIFFOpen(fname_.c_str(), "r");
	if (tiff_ != nullptr) 
	{
		char *infobuf;
		TIFFGetField(tiff_, TIFFTAG_IMAGEDESCRIPTION , &infobuf);
	   	TIFFClose(tiff_);

		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_string(infobuf);;
		xml_metadata_ptr_ = std::make_shared<std::map<std::string, std::string>>();
		
		if (result){
			std::vector<pugi::xml_node> node_list;
			node_list.push_back(doc.child("OME").child("Image").child("Pixels"));
			node_list.push_back(doc.child("OME").child("Image").child("Pixels").child("TiffData"));

			for (const auto &node: node_list){
				for (const pugi::xml_attribute &attr: node.attributes()){
					xml_metadata_ptr_->emplace(attr.name(), attr.value());
				}
			}
			// get channel info

			// read structured annotaion
			pugi::xml_node annotion_list = doc.child("OME").child("StructuredAnnotations");
			for(const pugi::xml_node &annotation : annotion_list){
				auto key = annotation.child("Value").child("OriginalMetadata").child("Key").child_value();
				auto value = annotation.child("Value").child("OriginalMetadata").child("Value").child_value();
				xml_metadata_ptr_->emplace(key,value);
			}


		}
	
	} else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }	
}

size_t OmeTiffLoader::AdjustStride (size_t start_pos, size_t current_pos, size_t stride_val) const
{
	if (stride_val == 0) return current_pos; // guard against div by 0

	size_t tmp = current_pos-start_pos;
	if (tmp%stride_val == 0) 
	{
		return current_pos; // no adjustment needed
	} else 
	{
		return ((tmp/stride_val)+1)*stride_val; // move to the next eligible position
	}
}

short OmeTiffLoader::GetBitsPerSamples() const
{
	return tile_loader_->bitsPerSample();
}

short OmeTiffLoader::GetSamplesPerPixel() const
{
	return tile_loader_->samplePerPixel();
}

void OmeTiffLoader::SetViewRequests(size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride)
{

	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);

	for(size_t x=0; x<ih; x+=row_stride){
		size_t r_min = x;
		size_t r_max = x+row_stride-1;
		for(size_t y=0; y<iw; y+=col_stride){
			size_t c_min = y;
			size_t c_max = y+col_stride-1;
			auto [min_row_index, min_col_index] = GetTileContainingPixel(r_min, c_min);
			auto [max_row_index, max_col_index] = GetTileContainingPixel(r_max, c_max);
			for (auto i = min_row_index; i <= max_row_index; ++i)
			{
				for (auto j = min_col_index; j <= max_col_index; ++j)
				{
					fast_loader_->requestView(i, j, 0, 0);
				}
			}


		}
	}


}

std::shared_ptr<std::vector<uint32_t>> OmeTiffLoader::GetViewRequests(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const index_col_pixel_min, size_t const index_col_pixel_max)
{

	// Convention 
	// rows are X coordinate (increasing from top to bottom)
	// cols are Y coordinate (increasing from left to right)
	// we need to transform from Local Tile Coordinate to Global Pixel Coordiate to Virtual Tile Coordinate

	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);
	auto index_true_row_pixel_max = index_row_pixel_max > ih ? ih : index_row_pixel_max;
	auto index_true_col_pixel_max = index_col_pixel_max > iw ? iw : index_col_pixel_max;
	auto top_left_tile = GetTileContainingPixel(index_row_pixel_min, index_col_pixel_min);
	auto bottom_right_tile = GetTileContainingPixel(index_true_row_pixel_max, index_true_col_pixel_max);
	auto min_row_index = top_left_tile.first;
	auto min_col_index = top_left_tile.second;
	auto max_row_index = bottom_right_tile.first;
	auto max_col_index = bottom_right_tile.second;

	// now loop through each tile, get tile data, fill the virtual tile vector

	auto tw = tile_loader_->tileWidth(0);
	auto th = tile_loader_->tileHeight(0);

	auto vtw = index_true_col_pixel_max-index_col_pixel_min+1;
	auto vth = index_true_row_pixel_max-index_row_pixel_min+1;
	std::shared_ptr<std::vector<uint32_t>> virtual_tile_data_ptr = std::make_shared<std::vector<uint32_t>>(vtw * vth);

	for (int i = min_row_index; i <= max_row_index; ++i)
	{
		for (int j = min_col_index; j <= max_col_index; ++j)
		{
			const auto &view = fast_loader_->getBlockingResult();	
			// take row slice from local tile and place it in virtual tile
			// global_x = i*th + local_x;
			// virtual_x = global_x - index_row_pixel_min;
			// initial_global_y = j*tw + initial_local_y;
			// initial_virtual_y = initial_global_y - index_col_pixel_min;
			size_t initial_local_x = index_row_pixel_min > i*th ? index_row_pixel_min-i*th : 0;	
			size_t end_local_x = index_true_row_pixel_max < (i+1)*th ? index_true_row_pixel_max-i*th: th-1;
			size_t initial_local_y = index_col_pixel_min > j*tw ? index_col_pixel_min-j*tw : 0;
			size_t end_local_y = index_true_col_pixel_max < (j+1)*tw ? index_true_col_pixel_max-j*tw : tw-1;
			size_t initial_virtual_y = j*tw + initial_local_y - index_col_pixel_min;

			auto vw = view->viewWidth();
			auto vrw = view->radiusWidth();
			auto vrh = view->radiusHeight();
			auto view_ptr = view->viewOrigin() + vrh*vw; 
			auto virtual_tile_ptr = virtual_tile_data_ptr->begin();
			for (size_t local_x=initial_local_x; local_x<=end_local_x; ++local_x){
				size_t virtual_x = i*th + local_x - index_row_pixel_min;
				std::copy(view_ptr+local_x*vw+vrw+initial_local_y, view_ptr+local_x*vw+vrw+end_local_y+1,virtual_tile_ptr+virtual_x*vtw+initial_virtual_y);					
			}
			view->returnToMemoryManager();	
		}
	}

	return virtual_tile_data_ptr;

}
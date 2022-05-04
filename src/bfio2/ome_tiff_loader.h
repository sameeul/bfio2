#include <future>
#include <cstring>
#include <string>
#include <tuple>
#include <memory>
#include <vector>

#include <fast_loader/fast_loader.h>
#include <omp.h>
#include <pugixml.hpp>
#ifdef WITH_PYTHON_H
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

#include "bfio_tile_loader.h"
#include "ome_tiff_gs_strip_loader.h"
#include "ome_tiff_gs_tile_loader.h"
#include "sequence.h"

template <class SampleType>
class OmeTiffLoader{
    private:
        std::shared_ptr<BfioTileLoader<SampleType>> tile_loader_;
        mutable std::shared_ptr<std::map<std::tuple<size_t, size_t, size_t>, size_t>> ifd_data_ptr_;

        mutable std::shared_ptr<std::map<std::string, std::string>> xml_metadata_ptr_;
        std::shared_ptr<fl::FastLoaderGraph<fl::DefaultView<SampleType>>> fast_loader_;
        size_t n_threads_, nc_, nt_, nz_, ifd_offset_;
        short dim_order_;
        std::string fname_;
        
        std::tuple<uint32_t, uint32_t>  GetImageDimensions() const;
        std::tuple<uint32_t, uint32_t>  CalculateTileDimensions() const;
        bool CheckTileStatus() const;

        void ParseMetadata() const;
        size_t AdjustStride (size_t start_pos, size_t current_pos, size_t stride_val) const;

        void SetZCT();
        void CopyToVirtualTile(const Seq& rows, const Seq& cols, std::shared_ptr<std::vector<SampleType>> virtual_tile_data, const std::map<size_t, size_t>& ifd_offset_lookup);

    public:
        OmeTiffLoader(const std::string &fNameWithPath, const int num_threads=1);
        ~OmeTiffLoader();
 
        std::shared_ptr<std::vector<SampleType>> GetTileData(size_t const row, size_t const col, size_t const layer=0, size_t const channel=0, size_t const tstep=0);
        std::shared_ptr<std::vector<SampleType>> GetTileDataByIndex(size_t const tile_index, size_t const channel=0, size_t const tstep=0);
        std::shared_ptr<std::vector<SampleType>> GetVirtualTileData(const Seq& rows, const Seq& cols, const Seq& layers = Seq(0,0), const Seq& channels = Seq(0,0), const Seq& tsteps = Seq(0,0));
        std::shared_ptr<std::vector<SampleType>> GetVirtualTileDataStrided(const Seq& rows, const Seq& cols, const Seq& layers = Seq(0,0), const Seq& channels = Seq(0,0), const Seq& tsteps = Seq(0,0));
        size_t GetRowTileCount () const;
        size_t GetColumnTileCount () const;
        size_t GetImageHeight() const ;
        size_t GetImageWidth () const ;
        size_t GetImageDepth () const ;
        size_t GetTileHeight () const ;
        size_t GetTileWidth () const ;
        size_t GetTileDepth () const ;
        short GetBitsPerSamples () const ;
        short GetChannelCount () const;
        size_t GetTstepCount () const;
        std::pair<size_t, size_t> GetTileContainingPixel(size_t const row_pixel_index, size_t const col_pixel_index) const;
        std::string GetMetaDataValue(const std::string &metadata_key) const;
        void SetViewRequests(size_t const tile_width, size_t const tile_height, size_t const row_stride, size_t const col_stride);
        std::shared_ptr<std::vector<SampleType>> GetViewRequests(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const index_col_pixel_min, size_t const index_col_pixel_max);
        std::vector< std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> tile_coordinate_list_;
        size_t CalcIFDIndex (size_t Z, size_t C, size_t T) const; 
};

template <class SampleType>
OmeTiffLoader<SampleType>::OmeTiffLoader(const std::string &fname_with_path, const int num_threads) : 
	tile_loader_(nullptr),
	ifd_data_ptr_(nullptr),
	xml_metadata_ptr_(nullptr),
	fast_loader_(nullptr),
	n_threads_(num_threads),
	nc_(1),
	nt_(1),
	nz_(1),
	ifd_offset_(0),
	dim_order_(1),
	fname_(fname_with_path),
	tile_coordinate_list_(0)
{
	ParseMetadata();
	SetZCT();
	if (CheckTileStatus())
	{
		tile_loader_ = std::make_shared<OmeTiffGrayScaleTileLoader<SampleType>>(n_threads_, fname_);
	}
	else 
	{
		// since the file is not tiled, we provide the tile dimensions
		auto [tw, th]  = CalculateTileDimensions();		
		tile_loader_ = std::make_shared<OmeTiffGrayScaleStripLoader<SampleType>>(n_threads_, fname_, tw, th, 1);
	}

    auto options = std::make_unique<fl::FastLoaderConfiguration<fl::DefaultView<SampleType>>>(tile_loader_);
    // Set the configuration
    uint32_t radiusDepth = 0;
    uint32_t radiusHeight = 0;
    uint32_t radiusWidth = 1;

    options->radius(radiusDepth, radiusHeight, radiusWidth);
    options->ordered(true);
    options->borderCreatorConstant(0);
    options->cacheCapacity(0,120);
    options->viewAvailable(0,10);
    // Create the Fast Loader Graph	
    fast_loader_ = std::make_shared<fl::FastLoaderGraph<fl::DefaultView<SampleType>>>(std::move(options));
    // Execute the graph
    fast_loader_->executeGraph();
};

template <class SampleType>
OmeTiffLoader<SampleType>::~OmeTiffLoader(){
	xml_metadata_ptr_ = nullptr;
	ifd_data_ptr_ = nullptr;	
    fast_loader_->finishRequestingViews(); 
	fast_loader_->waitForTermination();
	tile_loader_ = nullptr;
};

template <class SampleType> size_t OmeTiffLoader<SampleType>::GetRowTileCount() const {return tile_loader_->numberTileHeight();}
template <class SampleType> size_t OmeTiffLoader<SampleType>::GetColumnTileCount() const {return tile_loader_->numberTileWidth();}
template <class SampleType> size_t OmeTiffLoader<SampleType>::GetImageHeight() const {return tile_loader_->fullHeight(0);}
template <class SampleType> size_t OmeTiffLoader<SampleType>::GetImageWidth() const {return tile_loader_->fullWidth(0);}
template <class SampleType> size_t OmeTiffLoader<SampleType>::GetImageDepth() const {return nz_;}
template <class SampleType> size_t OmeTiffLoader<SampleType>::GetTileHeight() const {return tile_loader_->tileHeight(0);}
template <class SampleType> size_t OmeTiffLoader<SampleType>::GetTileWidth() const {return tile_loader_->tileWidth(0);}
template <class SampleType> size_t OmeTiffLoader<SampleType>::GetTileDepth() const {return tile_loader_->tileDepth(0);}

template <class SampleType>
std::shared_ptr<std::vector<SampleType>> OmeTiffLoader<SampleType>::GetTileData(size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep)
{
    auto tw = tile_loader_->tileWidth(0);
    auto th = tile_loader_->tileHeight(0);
	auto td = tile_loader_->tileDepth(0);
	auto iw = tile_loader_->fullWidth(0);
	auto ih = tile_loader_->fullHeight(0);
	auto actual_tw = iw > (col+1)*tw -1 ? tw : iw - col*tw;
	auto actual_th = ih > (row+1)*th -1 ? th : ih - row*th;

    std::shared_ptr<std::vector<SampleType>> tile_data = std::make_shared<std::vector<SampleType>>(actual_tw * actual_th * td);
	auto ifd_index = CalcIFDIndex(layer, channel, tstep);
	fast_loader_->requestView(row, col, ifd_index, 0);
	const auto &view = fast_loader_->getBlockingResult();
	auto vw = view->viewWidth();
	auto vh = view->viewHeight();
	auto vrw = view->radiusWidth();
	auto vrh = view->radiusHeight();
	auto vrd = view->radiusDepth();
	auto view_ptr = view->viewOrigin() + vrd*vw*vh + vrh*vw;
	auto tile_ptr = tile_data->data();
	for (size_t i = 0; i < actual_th; ++i) 
	{	
		for(size_t j = 0; j < actual_tw; ++j)
		{
			tile_ptr[i*actual_tw+j] = *(view_ptr+i*vw+vrw+j);
		}
	}
	
	view->returnToMemoryManager();
    return tile_data;
}

template <class SampleType>
std::tuple<uint32_t, uint32_t>  OmeTiffLoader<SampleType>::GetImageDimensions  () const
{
	TIFF *tiff_ = TIFFOpen(fname_.c_str(), "r");
	if (tiff_ != nullptr) 
	{
		uint32_t w, l, ndirs;
		TIFFGetField(tiff_, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tiff_, TIFFTAG_IMAGELENGTH, &l);
		TIFFClose(tiff_);
		return {w, l};	
	} 
	else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }
}

template <class SampleType>
std::tuple<uint32_t, uint32_t>  OmeTiffLoader<SampleType>::CalculateTileDimensions() const
{
	auto [w, h] = GetImageDimensions();
	uint32_t default_width = 1024;
	uint32_t default_height = 1024;
	w = std::min ({ w, default_width });
	h = std::min ({ h, default_height });
	return {w, h};
}

template <class SampleType>
bool OmeTiffLoader<SampleType>::CheckTileStatus() const
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

template <class SampleType>
std::shared_ptr<std::vector<SampleType>> OmeTiffLoader<SampleType>::GetTileDataByIndex(size_t const tile_index, size_t const channel, size_t const tstep)
{
	size_t num_col_tiles = GetColumnTileCount();	
	size_t num_row_tiles = GetRowTileCount();
	size_t layer = tile_index/(num_col_tiles*num_row_tiles);
	size_t tile_index_2d = tile_index%(num_col_tiles*num_row_tiles); 
	size_t row = tile_index_2d/num_col_tiles;
	size_t col = tile_index_2d%num_col_tiles;
	auto tile_data = GetTileData(row, col, layer, channel, tstep);
    return tile_data;
}

template <class SampleType>
std::pair<size_t, size_t> OmeTiffLoader<SampleType>::GetTileContainingPixel(size_t const row_pixel_index, size_t const col_pixel_index) const
{
	size_t th = GetTileHeight();	
	size_t tw = GetTileWidth();
	size_t ih = GetImageHeight();
	size_t iw = GetImageWidth();
	size_t row = row_pixel_index/th;

	if (row_pixel_index >= ih)
	{
		row = (ih-1)/th;
	} 
	
	size_t col = col_pixel_index/tw;	
	if (col_pixel_index >= iw)
	{
		col = (iw-1)/tw;
	} 

	return std::make_pair(row, col);
}

template <class SampleType>
std::string OmeTiffLoader<SampleType>::GetMetaDataValue(const std::string &metadata_key) const
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

template <class SampleType>
void OmeTiffLoader<SampleType>::ParseMetadata() const
{	
	TIFF *tiff_ = TIFFOpen(fname_.c_str(), "r");
	if (tiff_ != nullptr) 
	{
		char* infobuf;
		TIFFGetField(tiff_, TIFFTAG_IMAGEDESCRIPTION , &infobuf);
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_string(infobuf);;
		xml_metadata_ptr_ = std::make_shared<std::map<std::string, std::string>>();
		ifd_data_ptr_ = std::make_shared<std::map<std::tuple<size_t, size_t, size_t>, size_t>>();
		if (result){
			pugi::xml_node pixel = doc.child("OME").child("Image").child("Pixels");

			for (const pugi::xml_attribute &attr: pixel.attributes()){
				xml_metadata_ptr_->emplace(attr.name(), attr.value());
			}
		

			// get TiffData info
        	for (pugi::xml_node tiff_data: pixel.children("TiffData"))
			{
				size_t c=0, t=0, z=0, ifd;
				for (pugi::xml_attribute attr: tiff_data.attributes())
				{
					if (strcmp(attr.name(),"FirstC") == 0) {c = atoi(attr.value());}
					else if (strcmp(attr.name(),"FirstZ") == 0) {z = atoi(attr.value());}
					else if (strcmp(attr.name(),"FirstT") == 0) {t = atoi(attr.value());}
					else if (strcmp(attr.name(),"IFD") == 0) {ifd = atoi(attr.value());}
					else {continue;}
				} 
				ifd_data_ptr_->emplace(std::make_pair(std::make_tuple(z,c,t),ifd));
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
		TIFFClose(tiff_);
	
	} else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }	
}

template <class SampleType> 
size_t OmeTiffLoader<SampleType>::AdjustStride (size_t start_pos, size_t current_pos, size_t stride_val) const
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
template <class SampleType>
short OmeTiffLoader<SampleType>::GetBitsPerSamples() const
{
	return tile_loader_->bitsPerSample();
}

template <class SampleType>
short OmeTiffLoader<SampleType>::GetChannelCount() const
{
	return nc_;
}

template <class SampleType> 
size_t OmeTiffLoader<SampleType>::GetTstepCount() const
{
	return nt_;
}

template <class SampleType>
void OmeTiffLoader<SampleType>::SetViewRequests(size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride)
{

	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);
	tile_coordinate_list_.clear();

	for (size_t t =0; t<nt_; ++t)
	{
		for (size_t c =0; c<nc_; ++c)
		{
			for (size_t z =0; z<nz_; ++z)
			{
				auto ifd_index = CalcIFDIndex(z,c,t);
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
								fast_loader_->requestView(i, j, ifd_index, 0);
								tile_coordinate_list_.push_back(std::make_tuple(r_min,r_max,c_min,c_max));
							}
						}
					}
				}
			}
		}
	}



}

template <class SampleType>
std::shared_ptr<std::vector<SampleType>> OmeTiffLoader<SampleType>::GetViewRequests(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const index_col_pixel_min, size_t const index_col_pixel_max)
{

	// Convention 
	// rows are X coordinate (increasing from top to bottom)
	// cols are Y coordinate (increasing from left to right)
	// we need to transform from Local Tile Coordinate to Global Pixel Coordiate to Virtual Tile Coordinate

	Seq rows = Seq(index_row_pixel_min, index_row_pixel_max);
	Seq cols = Seq(index_col_pixel_min, index_col_pixel_max);

	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);
	auto index_true_row_pixel_max = index_row_pixel_max > ih ? ih-1 : index_row_pixel_max;
	auto index_true_col_pixel_max = index_col_pixel_max > iw ? iw-1 : index_col_pixel_max;
	auto top_left_tile = GetTileContainingPixel(index_row_pixel_min, index_col_pixel_min);
	auto bottom_right_tile = GetTileContainingPixel(index_true_row_pixel_max, index_true_col_pixel_max);
	auto min_row_index = top_left_tile.first;
	auto min_col_index = top_left_tile.second;
	auto max_row_index = bottom_right_tile.first;
	auto max_col_index = bottom_right_tile.second;

	auto vtw = index_true_col_pixel_max-index_col_pixel_min+1;
	auto vth = index_true_row_pixel_max-index_row_pixel_min+1;
	std::shared_ptr<std::vector<SampleType>> virtual_tile_data = std::make_shared<std::vector<SampleType>>(vtw * vth);

	// now loop through each tile, get tile data, fill the virtual tile vector
	std::map<size_t, size_t> ifd_offset_lookup;
	for (int i = min_row_index; i <= max_row_index; ++i)
	{
		for (int j = min_col_index; j <= max_col_index; ++j)
		{
			CopyToVirtualTile(rows, cols, virtual_tile_data, ifd_offset_lookup);	
		}
	}

	return virtual_tile_data;
}

template <class SampleType>
std::shared_ptr<std::vector<SampleType>> OmeTiffLoader<SampleType>::GetVirtualTileData(const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)
{
	// Convention 
	// rows are X coordinate (increasing from top to bottom)
	// cols are Y coordinate (increasing from left to right)
	// we need to transform from Local Tile Coordinate to Global Pixel Coordiate to Virtual Tile Coordinate

	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);
	auto id = GetImageDepth();
	auto index_true_row_pixel_max = rows.Stop() > ih ? ih-1 : rows.Stop();
	auto index_true_col_pixel_max = cols.Stop() > iw ? iw-1 : cols.Stop();
	auto top_left_tile = GetTileContainingPixel(rows.Start(), cols.Start());
	auto bottom_right_tile = GetTileContainingPixel(index_true_row_pixel_max, index_true_col_pixel_max);
	auto min_row_index = top_left_tile.first;
	auto min_col_index = top_left_tile.second;
	auto max_row_index = bottom_right_tile.first;
	auto max_col_index = bottom_right_tile.second;

	auto index_true_min_layer = layers.Start() > 0? layers.Start() : 0;
	auto index_true_max_layer = layers.Stop() > nz_-1? nz_-1 : layers.Stop();
	auto index_true_min_channel = channels.Start() > 0? channels.Start() : 0;
	auto index_true_max_channel = channels.Stop() > nc_-1? nc_-1 : channels.Stop();
	auto index_true_min_tstep = tsteps.Start() > 0? tsteps.Start() : 0;
	auto index_true_max_tstep = tsteps.Stop() > nt_-1? nt_-1 : tsteps.Stop();

	auto vtw = index_true_col_pixel_max-cols.Start()+1;
	auto vth = index_true_row_pixel_max-rows.Start()+1;
	auto vtd = (index_true_max_layer-index_true_min_layer)/layers.Step()+1;
	auto num_channels = (index_true_max_channel-index_true_min_channel)/channels.Step()+1;
	auto num_tsteps = (index_true_max_tstep-index_true_min_tstep)/tsteps.Step()+1;
	std::shared_ptr<std::vector<SampleType>> virtual_tile_data = std::make_shared<std::vector<SampleType>>(vtw * vth * vtd * num_channels * num_tsteps);

#ifdef WITH_PYTHON_H
 	py::gil_scoped_acquire acquire;
#endif
	size_t virtual_tstep = 0;
	size_t total_views = 0;
	std::map<size_t, size_t> ifd_offset_lookup;

	for (auto m = index_true_min_tstep; m<=index_true_max_tstep; m=m+tsteps.Step())
	{
		size_t t_offset = virtual_tstep*vtw*vth*vtd*num_channels;
		size_t virtual_ch = 0;
		for (auto l = index_true_min_channel; l<=index_true_max_channel; l=l+channels.Step())
		{	
			size_t ch_offset = virtual_ch*vtw*vth*vtd;
			size_t virtual_z = 0;
			for (auto k = index_true_min_layer; k<=index_true_max_layer; k=k+layers.Step())
			{
				size_t z_offset = virtual_z*vtw*vth;
				auto offset = t_offset+ch_offset+z_offset;
				for (auto i = min_row_index; i <= max_row_index; ++i)
				{
				//	#pragma omp parallel for
					for (auto j = min_col_index; j <= max_col_index; ++j)
					{	
						auto layer = CalcIFDIndex(k,l,m);
						fast_loader_->requestView(i, j, layer, 0);
						total_views++;
						ifd_offset_lookup.emplace(layer, offset);
					}
				}
				++virtual_z;
			}
			++virtual_ch;
		}
		++virtual_tstep;
	}

//	#pragma omp parallel for
	for(auto i=0; i<total_views; i++){
		CopyToVirtualTile(rows, cols, virtual_tile_data, ifd_offset_lookup);
	}

#ifdef WITH_PYTHON_H
 	py::gil_scoped_release release;
#endif
	return virtual_tile_data;
}


template <class SampleType>
std::shared_ptr<std::vector<SampleType>> OmeTiffLoader<SampleType>::GetVirtualTileDataStrided(const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)
{

	// Convention 
	// rows are X coordinate (increasing from top to bottom)
	// cols are Y coordinate (increasing from left to right)
	// we need to transform from Local Tile Coordinate to Global Pixel Coordiate to Virtual Tile Coordinate

	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);
	auto id = GetImageDepth();
	auto index_true_row_pixel_max = rows.Stop() > ih ? ih-1 : rows.Stop();
	auto index_true_col_pixel_max = cols.Stop() > iw ? iw-1 : cols.Stop();
	auto top_left_tile = GetTileContainingPixel(rows.Start(), cols.Start());
	auto bottom_right_tile = GetTileContainingPixel(index_true_row_pixel_max, index_true_col_pixel_max);
	auto min_row_index = top_left_tile.first;
	auto min_col_index = top_left_tile.second;
	auto max_row_index = bottom_right_tile.first;
	auto max_col_index = bottom_right_tile.second;

	auto index_true_min_layer = layers.Start() > 0? layers.Start() : 0;
	auto index_true_max_layer = layers.Stop() > nz_-1? nz_-1 : layers.Stop();
	auto index_true_min_channel = channels.Start() > 0? channels.Start() : 0;
	auto index_true_max_channel = channels.Stop() > nc_-1? nc_-1 : channels.Stop();
	auto index_true_min_tstep = tsteps.Start() > 0? tsteps.Start() : 0;
	auto index_true_max_tstep = tsteps.Stop() > nt_-1? nt_-1 : tsteps.Stop();

	auto vtw = index_true_col_pixel_max-cols.Start()+1;
	auto vth = index_true_row_pixel_max-rows.Start()+1;
	auto vtd = (index_true_max_layer-index_true_min_layer)/layers.Step()+1;
	auto num_channels = (index_true_max_channel-index_true_min_channel)/channels.Step()+1;
	auto num_tsteps = (index_true_max_tstep-index_true_min_tstep)/tsteps.Step()+1;
	std::shared_ptr<std::vector<SampleType>> virtual_tile_data = std::make_shared<std::vector<SampleType>>(vtw * vth * vtd * num_channels * num_tsteps) ;

	size_t virtual_tstep = 0;
	size_t total_views = 0;
	std::map<size_t, size_t> ifd_offset_lookup;
	for (auto m = index_true_min_tstep; m<=index_true_max_tstep; m=m+tsteps.Step())
	{
		size_t t_offset = virtual_tstep*vtw*vth*vtd*num_channels;
		size_t virtual_ch = 0;
		for (auto l = index_true_min_channel; l<=index_true_max_channel; l=l+channels.Step())
		{	
			size_t ch_offset = virtual_ch*vtw*vth*vtd;
			size_t virtual_z = 0;
			for (auto k = index_true_min_layer; k<=index_true_max_layer; k=k+layers.Step())
			{
				size_t z_offset = virtual_z*vtw*vth;
				
				for (auto i = min_row_index; i <= max_row_index; ++i){
					for (auto j = min_col_index; j <= max_col_index; ++j)
					{
						auto offset = t_offset+ch_offset+z_offset;
						auto layer = CalcIFDIndex(k,l,m);
						fast_loader_->requestView(i, j, layer, 0);
						total_views++;
						ifd_offset_lookup.emplace(layer, offset);
					}
				}
				++virtual_z;
			}
			++virtual_ch;
		}
		++virtual_tstep;
	}
	#pragma omp parallel for
	for(auto i=0; i<total_views; i++){
		CopyToVirtualTile(rows, cols, virtual_tile_data, ifd_offset_lookup);
	}

	return virtual_tile_data;
}

template <class SampleType>
void OmeTiffLoader<SampleType>::SetZCT()
{
	auto it = xml_metadata_ptr_->find("SizeC");
	if (it != xml_metadata_ptr_->end()) nc_ = stoi(it->second);
	
	it = xml_metadata_ptr_->find("SizeZ");
	if (it != xml_metadata_ptr_->end()) nz_ = stoi(it->second);

	it = xml_metadata_ptr_->find("SizeT");
	if (it != xml_metadata_ptr_->end()) nt_ = stoi(it->second);
	
	auto ifd_it = ifd_data_ptr_->find(std::make_tuple(size_t(0), size_t(0), size_t(0)));
	if (ifd_it != ifd_data_ptr_->end()) ifd_offset_ = ifd_it->second;

	it = xml_metadata_ptr_->find("DimensionOrder");
	if (it != xml_metadata_ptr_->end())
	{
		auto dim_order_str = it->second;
		if (dim_order_str == "XYZTC") { dim_order_ = 1;}
		else if (dim_order_str == "XYZCT") { dim_order_ = 2;}
		else if (dim_order_str == "XYTCZ") { dim_order_ = 4;}
		else if (dim_order_str == "XYTZC") { dim_order_ = 8;}
		else if (dim_order_str == "XYCTZ") { dim_order_ = 16;}
		else if (dim_order_str == "XYCZT") { dim_order_ = 32;}
		else { dim_order_ = 1;}
	}
}

template <class SampleType> 
size_t OmeTiffLoader<SampleType>::CalcIFDIndex (size_t z, size_t c, size_t t) const
{
	size_t ifd_dir = 0;
	if (ifd_data_ptr_ != nullptr){
		auto ifd_it = ifd_data_ptr_->find(std::make_tuple(z,c,t));
		if (ifd_it != ifd_data_ptr_->end()){
			ifd_dir = ifd_it->second;
		} 
	}
	else 
	{
		if (dim_order_ == 1) {ifd_dir = nz_*nt_*c + nz_*t + z + ifd_offset_;}
		else if (dim_order_ == 2) {ifd_dir = nz_*nc_*t + nz_*c + z + ifd_offset_;}
		else if (dim_order_ == 4) {ifd_dir = nt_*nc_*z + nt_*c + t + ifd_offset_;}
		else if (dim_order_ == 8) {ifd_dir = nt_*nz_*c + nt_*z + t + ifd_offset_;}
		else if (dim_order_ == 16) {ifd_dir = nc_*nt_*z + nc_*t + c + ifd_offset_;}
		else if (dim_order_ == 32) {ifd_dir = nc_*nz_*t + nc_*z + c + ifd_offset_;}
		else {ifd_dir = nz_*nt_*c + nz_*t + z + ifd_offset_;}
	}

	return ifd_dir;
}

template <class SampleType>
void OmeTiffLoader<SampleType>::CopyToVirtualTile(const Seq& rows, const Seq& cols, std::shared_ptr<std::vector<SampleType>> virtual_tile_data, const std::map<size_t, size_t>& ifd_offset_lookup)
{
	auto ih = tile_loader_->fullHeight(0);
	auto iw = tile_loader_->fullWidth(0);
	auto tw = tile_loader_->tileWidth(0);
	auto th = tile_loader_->tileHeight(0);

	auto index_true_row_pixel_max = rows.Stop() > ih ? ih-1 : rows.Stop();
	auto index_true_col_pixel_max = cols.Stop() > iw ? iw-1 : cols.Stop();
	auto vtw = index_true_col_pixel_max-cols.Start()+1;
	auto vth = index_true_row_pixel_max-rows.Start()+1;

	const auto &view = fast_loader_->getBlockingResult();
	auto i = view->tileRowIndex();
	auto j = view->tileColIndex();
	auto offset = ifd_offset_lookup.at(view->tileLayerIndex());	
	// take row slice from local tile and place it in virtual tile
	// global_x = i*th + local_x;
	// virtual_x = global_x - index_row_pixel_min;
	// initial_global_y = j*tw + initial_local_y;
	// initial_virtual_y = initial_global_y - index_col_pixel_min;
	size_t initial_local_x = rows.Start() > i*th ? rows.Start()-i*th : 0;	
	// adjust for row stride
	if (rows.Step()!=1)
	{
		size_t initial_global_x = i*th + initial_local_x;
		initial_global_x = AdjustStride(rows.Start(), initial_global_x, rows.Step());
		initial_local_x = initial_global_x - i*th;
	}
	size_t end_local_x = index_true_row_pixel_max < (i+1)*th ? index_true_row_pixel_max-i*th: th-1;
	size_t initial_local_y = cols.Start() > j*tw ? cols.Start()-j*tw : 0;
	// adjust for col stride
	if (cols.Step()!=1)
	{
		size_t initial_global_y = j*tw + initial_local_y;
		initial_global_y = AdjustStride(cols.Start(), initial_global_y, cols.Step());
		initial_local_y = initial_global_y - j*tw;		
	}

	size_t end_local_y = index_true_col_pixel_max < (j+1)*tw ? index_true_col_pixel_max-j*tw : tw-1;
	size_t initial_virtual_y = j*tw + initial_local_y - cols.Start();

	auto vw = view->viewWidth();
	auto vh = view->viewHeight();
	auto vrw = view->radiusWidth();
	auto vrh = view->radiusHeight();
	auto vrd = view->radiusDepth();
	auto view_ptr = view->viewOrigin() + vrd*vw*vh + vrh*vw;
	auto virtual_tile_data_ptr = virtual_tile_data->data();
	auto virtual_tile_data_begin = virtual_tile_data->begin();


	if (cols.Step() == 1){
		for (size_t local_x=initial_local_x; local_x<=end_local_x; ++local_x){
			size_t virtual_x = (i*th + local_x - rows.Start())/rows.Step();
			std::copy(view_ptr+local_x*vw+vrw+initial_local_y, view_ptr+local_x*vw+vrw+end_local_y+1,virtual_tile_data_begin+offset+virtual_x*vtw+initial_virtual_y);					
		}
	}
	else 
	{
		for (size_t local_x=initial_local_x; local_x<=end_local_x; ++local_x){
			size_t virtual_x = (i*th + local_x - rows.Start())/rows.Step();
			for (size_t local_y=initial_local_y; local_y<=end_local_y; local_y=local_y+cols.Step())
			{
				size_t virtual_y = (j*tw + local_y - cols.Start())/cols.Step();
				virtual_tile_data_ptr[offset+virtual_x*vtw+virtual_y] = *(view_ptr+local_x*vw+vrw+local_y);
			}
		}
	}
	view->returnToMemoryManager();
}
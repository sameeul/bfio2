#include <future>
#include <cstring>
#include <string>
#include <tuple>
#include <memory>
#include <regex>
#include <vector>
#include <execution>
#include <thread>
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
#include "omezarr.h"
#include "sequence.h"
#include "utilities.h"
#include "thread_pool.hpp"
#ifndef FOLLY_UMH_H
#define FOLLY_UMH_H
#include "UninitializedMemoryHacks.h"
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(uint16_t)
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(uint32_t)
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(uint64_t)
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(int8_t)
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(int16_t)
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(int32_t)
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(int64_t)
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(float)
FOLLY_DECLARE_VECTOR_RESIZE_WITHOUT_INIT(double)
#endif

template <class SampleType>
class BioReader{
    private:
        std::shared_ptr<BfioTileLoader<SampleType>> tile_loader_;
        mutable std::unique_ptr<std::vector<size_t>> ifd_data_ptr_;

        mutable std::shared_ptr<std::map<std::string, std::string>> xml_metadata_ptr_;
        std::unique_ptr<fl::FastLoaderGraph<fl::DefaultView<SampleType>>> fast_loader_;
        size_t n_threads_, nc_, nt_, nz_, ifd_offset_;
        short dim_order_;
		std::string fname_;
        short fl_cut_off = 49;
        bool is_zarr_ = false;

        void ParseMetadata();
        void ParseTiffMetadata();
        void ParseZarrMetadata();
        size_t AdjustStride (size_t start_pos, size_t current_pos, size_t stride_val) const;
        void CopyToVirtualTile(const Seq& rows, const Seq& cols, const std::shared_ptr<std::vector<SampleType>> virtual_tile_data, const std::map<size_t, size_t>& ifd_offset_lookup);

    public:
        BioReader(const std::string &fNameWithPath, const int num_threads=1);
        ~BioReader();
 
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
BioReader<SampleType>::BioReader(const std::string &fname_with_path, const int num_threads) : 
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
	tile_coordinate_list_(0),
    is_zarr_(false)
{
    if 	(std::filesystem::path(fname_).extension() == ".zarr") {is_zarr_ = true;}

    if (is_zarr_)
    {
        tile_loader_ = std::make_shared<OmeZarrTileLoader<SampleType>>(n_threads_, fname_);
        nc_ = tile_loader_->numberChannels();
        nz_ = static_cast<size_t>(ceil(1.0*tile_loader_->fullDepth(0)/ tile_loader_->tileDepth(0)));
        nt_ = tile_loader_->numTsteps();
    }
    else 
    {
        if (CheckTileStatus(fname_))
        {
            tile_loader_ = std::make_shared<OmeTiffGrayScaleTileLoader<SampleType>>(n_threads_, fname_);
        }
        else 
        {
            tile_loader_ = std::make_shared<OmeTiffGrayScaleStripLoader<SampleType>>(n_threads_, fname_);
        }
    }

	ParseMetadata();

    auto options = std::make_unique<fl::FastLoaderConfiguration<fl::DefaultView<SampleType>>>(tile_loader_);
    // Set the configuration
    uint32_t radiusDepth = 0;
    uint32_t radiusHeight = 0;
    uint32_t radiusWidth = 0;

    options->radius(radiusDepth, radiusHeight, radiusWidth);
    options->ordered(true);
    options->borderCreatorConstant(0);
	auto processor_count = std::thread::hardware_concurrency();
	if (processor_count == 0){processor_count = 1;}
	if (nc_*nz_*nt_ < fl_cut_off) 
	{
		options->cacheCapacity(0,0);
		options->viewAvailable(0,processor_count*2);
	}
	else
	{
		options->cacheCapacity(0,processor_count*24);
		options->viewAvailable(0,processor_count*3);
	}
    // Create the Fast Loader Graph	
    fast_loader_ = std::make_unique<fl::FastLoaderGraph<fl::DefaultView<SampleType>>>(std::move(options));
    // Execute the graph
    fast_loader_->executeGraph();
};

template <class SampleType>
BioReader<SampleType>::~BioReader(){
	xml_metadata_ptr_ = nullptr;
	ifd_data_ptr_ = nullptr;	
    fast_loader_->finishRequestingViews(); 
	fast_loader_->waitForTermination();
	tile_loader_ = nullptr;
};

template <class SampleType> size_t BioReader<SampleType>::GetRowTileCount() const {return tile_loader_->numberTileHeight();}
template <class SampleType> size_t BioReader<SampleType>::GetColumnTileCount() const {return tile_loader_->numberTileWidth();}
template <class SampleType> size_t BioReader<SampleType>::GetImageHeight() const {return tile_loader_->fullHeight(0);}
template <class SampleType> size_t BioReader<SampleType>::GetImageWidth() const {return tile_loader_->fullWidth(0);}
template <class SampleType> size_t BioReader<SampleType>::GetImageDepth() const {return nz_;}
template <class SampleType> size_t BioReader<SampleType>::GetTileHeight() const {return tile_loader_->tileHeight(0);}
template <class SampleType> size_t BioReader<SampleType>::GetTileWidth() const {return tile_loader_->tileWidth(0);}
template <class SampleType> size_t BioReader<SampleType>::GetTileDepth() const {return tile_loader_->tileDepth(0);}

template <class SampleType>
std::shared_ptr<std::vector<SampleType>> BioReader<SampleType>::GetTileData(size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep)
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
std::shared_ptr<std::vector<SampleType>> BioReader<SampleType>::GetTileDataByIndex(size_t const tile_index, size_t const channel, size_t const tstep)
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
std::pair<size_t, size_t> BioReader<SampleType>::GetTileContainingPixel(size_t const row_pixel_index, size_t const col_pixel_index) const
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
std::string BioReader<SampleType>::GetMetaDataValue(const std::string &metadata_key) const
{
	std::string value = "" ;
	// if (xml_metadata_ptr_ == nullptr){
	// 	ParseMetadata();		
	// }
	try {
		value = xml_metadata_ptr_->at(metadata_key);
	}
	catch (const std::exception& e) {
		std::cout<<"Requested metadata key not found"<<std::endl;
	}
	return value;
}

template <class SampleType>
void BioReader<SampleType>::ParseMetadata() 
{
    if(is_zarr_) { ParseZarrMetadata();}
    else {ParseTiffMetadata();}
}

template <class SampleType>
void BioReader<SampleType>::ParseZarrMetadata() 
{	
	std::unique_ptr zarr_ptr = std::make_unique<z5::filesystem::handle::File>(fname_.c_str());
	nlohmann::json attributes;
	z5::readAttributes(*zarr_ptr, attributes);
	std::string metadata = attributes["metadata"].dump();

	std::regex ome("ome:");
	metadata = std::regex_replace(metadata, ome, "");

	std::regex quote("\\\\\"");
	metadata = std::regex_replace(metadata, quote, "\"");

	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_string(metadata.c_str());;
	xml_metadata_ptr_ = std::make_shared<std::map<std::string, std::string>>();
	if (result){
		pugi::xml_node pixel = doc.child("OME").child("Image").child("Pixels");

		for (const pugi::xml_attribute &attr: pixel.attributes()){
			xml_metadata_ptr_->emplace(attr.name(), attr.value());
		}
	// get channel info -ToDo

	// read structured annotaion
		pugi::xml_node annotion_list = doc.child("OME").child("StructuredAnnotations");
		for(const pugi::xml_node &annotation : annotion_list){
			auto key = annotation.child("Value").child("OriginalMetadata").child("Key").child_value();
			auto value = annotation.child("Value").child("OriginalMetadata").child("Value").child_value();
			xml_metadata_ptr_->emplace(key,value);
		}
	}
}

template <class SampleType>
void BioReader<SampleType>::ParseTiffMetadata() 
{	
	TIFF *tiff_ = TIFFOpen(fname_.c_str(), "r");
	if (tiff_ != nullptr) 
	{
		char* infobuf;
		TIFFGetField(tiff_, TIFFTAG_IMAGEDESCRIPTION , &infobuf);
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_string(infobuf);;
		xml_metadata_ptr_ = std::make_shared<std::map<std::string, std::string>>();
		ifd_data_ptr_ = std::make_unique<std::vector<size_t>>(1);
		if (result){
			pugi::xml_node pixel = doc.child("OME").child("Image").child("Pixels");

			for (const pugi::xml_attribute &attr: pixel.attributes()){
				xml_metadata_ptr_->emplace(attr.name(), attr.value());
			}
		

			auto it = xml_metadata_ptr_->find("SizeC");
			if (it != xml_metadata_ptr_->end()) nc_ = stoi(it->second);

			it = xml_metadata_ptr_->find("SizeZ");
			if (it != xml_metadata_ptr_->end()) nz_ = stoi(it->second);

			it = xml_metadata_ptr_->find("SizeT");
			if (it != xml_metadata_ptr_->end()) nt_ = stoi(it->second);

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

			//set up look up vector for ifp_data index
			auto ifd_nums = nc_*nz_*nt_;
			ifd_data_ptr_->at(0) = 99999999; //hack for now
			for(auto i=1; i<ifd_nums;++i){ifd_data_ptr_->emplace_back(99999999);}

			// get TiffData info
        	for (pugi::xml_node tiff_data: pixel.children("TiffData"))
			{
				size_t c=0, t=0, z=0, ifd=0;
				for (pugi::xml_attribute attr: tiff_data.attributes())
				{
					if (strcmp(attr.name(),"FirstC") == 0) {c = atoi(attr.value());}
					else if (strcmp(attr.name(),"FirstZ") == 0) {z = atoi(attr.value());}
					else if (strcmp(attr.name(),"FirstT") == 0) {t = atoi(attr.value());}
					else if (strcmp(attr.name(),"IFD") == 0) {ifd = atoi(attr.value());}
					else {continue;}
				} 

				ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z) = ifd;
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
size_t BioReader<SampleType>::AdjustStride (size_t start_pos, size_t current_pos, size_t stride_val) const
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
short BioReader<SampleType>::GetBitsPerSamples() const
{
	return tile_loader_->bitsPerSample();
}

template <class SampleType>
short BioReader<SampleType>::GetChannelCount() const
{
	return nc_;
}

template <class SampleType> 
size_t BioReader<SampleType>::GetTstepCount() const
{
	return nt_;
}

template <class SampleType>
void BioReader<SampleType>::SetViewRequests(size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride)
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
std::shared_ptr<std::vector<SampleType>> BioReader<SampleType>::GetViewRequests(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
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
std::shared_ptr<std::vector<SampleType>> BioReader<SampleType>::GetVirtualTileData(const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)
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
	std::shared_ptr<std::vector<SampleType>> virtual_tile_data = std::make_shared<std::vector<SampleType>>(0) ;
	folly::resizeWithoutInitialization(*virtual_tile_data, vtw * vth * vtd * num_channels * num_tsteps);

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


#ifdef WITH_PYTHON_H
 	py::gil_scoped_release release;
#endif
	auto processor_count = std::thread::hardware_concurrency();
	if (processor_count == 0){processor_count = 1;}
	short pool_worker;
	if (total_views < fl_cut_off) {pool_worker = processor_count;}
	else {pool_worker = 3*processor_count;}
	thread_pool pool(pool_worker);
	pool.parallelize_loop(0, total_views, 
							[&rows, &cols, virtual_tile_data, &ifd_offset_lookup, this](const size_t &a, const size_t &b)
							{
								for (size_t i = a; i < b; i++)
								this->CopyToVirtualTile(rows, cols, virtual_tile_data, ifd_offset_lookup);
							}
							);
	// #pragma omp parallel for
	// for(auto i=0; i<total_views; i++){
	// 	CopyToVirtualTile(rows, cols, virtual_tile_data, ifd_offset_lookup);
	// }

#ifdef WITH_PYTHON_H
 	py::gil_scoped_acquire acquire;
#endif
	return virtual_tile_data;
}


template <class SampleType>
std::shared_ptr<std::vector<SampleType>> BioReader<SampleType>::GetVirtualTileDataStrided(const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)
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
	std::shared_ptr<std::vector<SampleType>> virtual_tile_data = std::make_shared<std::vector<SampleType>>(0) ;
	folly::resizeWithoutInitialization(*virtual_tile_data, vtw * vth * vtd * num_channels * num_tsteps);

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
//	#pragma omp parallel for
	for(auto i=0; i<total_views; i++){
		CopyToVirtualTile(rows, cols, virtual_tile_data, ifd_offset_lookup);
	}

	return virtual_tile_data;
}


template <class SampleType> 
size_t BioReader<SampleType>::CalcIFDIndex (size_t z, size_t c, size_t t) const
{
    if (is_zarr_) {return t*nz_*nc_ + c*nz_ + z;}

    // come here for Tiff
	size_t ifd_dir = 0;
	ifd_dir = ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z);
	if (ifd_dir == 99999999){
		switch (dim_order_)
		{
			case 1:
				ifd_dir = nz_*nt_*c + nz_*t + z + ifd_offset_;
				ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z) = ifd_dir;
				break;
			case 2:
				ifd_dir = nz_*nc_*t + nz_*c + z + ifd_offset_;
				ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z) = ifd_dir;
				break;
			case 4:
				ifd_dir = nt_*nc_*z + nt_*c + t + ifd_offset_;
				ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z) = ifd_dir;
				break;
			case 8:
				ifd_dir = nt_*nz_*c + nt_*z + t + ifd_offset_;
				ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z) = ifd_dir;
				break;
			case 16:
				ifd_dir = nc_*nt_*z + nc_*t + c + ifd_offset_;
				ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z) = ifd_dir;
				break;
			case 32:
				ifd_dir = nc_*nz_*t + nc_*z + c + ifd_offset_;
				ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z) = ifd_dir;
				break;
			
			default:
				ifd_dir = nz_*nt_*c + nz_*t + z + ifd_offset_;
				ifd_data_ptr_->at(t*(nz_*nc_) + c*nz_ + z) = ifd_dir;
				break;
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
void BioReader<SampleType>::CopyToVirtualTile(const Seq& rows, const Seq& cols, const std::shared_ptr<std::vector<SampleType>> virtual_tile_data, const std::map<size_t, size_t>& ifd_offset_lookup)
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
	if (view != nullptr){
		auto i = view->tileRowIndex();
		auto j = view->tileColIndex();
		size_t offset = 0;
		if (ifd_offset_lookup.size()!=0) {
			offset = ifd_offset_lookup.at(view->tileLayerIndex());
		}
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
		   // #pragma omp parallel for
			for (size_t local_x=initial_local_x; local_x<end_local_x+1; ++local_x){
				size_t virtual_x = (i*th + local_x - rows.Start())/rows.Step();
				std::copy(std::execution::par_unseq, view_ptr+local_x*vw+vrw+initial_local_y, view_ptr+local_x*vw+vrw+end_local_y+1,virtual_tile_data_begin+offset+virtual_x*vtw+initial_virtual_y);					
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

}
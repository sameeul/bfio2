#include <cstring>
#include <string>
#include <tuple>
#include <memory>
#include <vector>
#include <fast_loader/fast_loader.h>
#include <omp.h>
#include <pugixml.hpp>
#include "bfio_tile_loader.h"
#include "ome_tiff_gs_strip_loader.h"
#include "ome_tiff_gs_tile_loader.h"
#include "sequence.h"


class OmeTiffLoader{

    private:
        std::shared_ptr<BfioTileLoader<uint32_t>> tile_loader_;
        mutable std::shared_ptr<std::map<std::tuple<size_t, size_t, size_t>, size_t>> ifd_data_ptr_;

        mutable std::shared_ptr<std::map<std::string, std::string>> xml_metadata_ptr_;
        std::shared_ptr<fl::FastLoaderGraph<fl::DefaultView<uint32_t>>> fast_loader_;
        size_t n_threads_, nc_, nt_, nz_, ifd_offset_;
        short dim_order_;
        std::string fname_;
        
        std::tuple<uint32_t, uint32_t>  GetImageDimensions() const;
        std::tuple<uint32_t, uint32_t>  CalculateTileDimensions() const;
        bool CheckTileStatus() const;

        void ParseMetadata() const;
        size_t AdjustStride (size_t start_pos, size_t current_pos, size_t stride_val) const;

        void SetZCT();

    public:
        OmeTiffLoader(const std::string &fNameWithPath, const int num_threads=1);
        ~OmeTiffLoader();
 
        std::shared_ptr<std::vector<uint32_t>> GetTileData(size_t const row, size_t const col, size_t const layer=0, size_t const channel=0, size_t const tstep=0);
        std::shared_ptr<std::vector<uint32_t>> GetTileDataByIndex(size_t const tile_index, size_t const channel=0, size_t const tstep=0);
        std::shared_ptr<std::vector<uint32_t>> GetVirtualTileData(const Seq& rows, const Seq& cols, const Seq& layers = Seq(0,0), const Seq& channels = Seq(0,0), const Seq& tsteps = Seq(0,0));
        std::shared_ptr<std::vector<uint32_t>> GetVirtualTileDataStrided(const Seq& rows, const Seq& cols, const Seq& layers = Seq(0,0), const Seq& channels = Seq(0,0), const Seq& tsteps = Seq(0,0));
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
        std::shared_ptr<std::vector<uint32_t>> GetViewRequests(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const index_col_pixel_min, size_t const index_col_pixel_max);
        std::vector< std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> tile_coordinate_list_;
        size_t CalcIFDIndex (size_t Z, size_t C, size_t T) const;
};
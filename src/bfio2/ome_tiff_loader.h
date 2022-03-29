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

class OmeTiffLoader{

    private:
        std::shared_ptr<BfioTileLoader<uint32_t>> tile_loader_;
        mutable std::shared_ptr<std::map<std::string, std::string>> xml_metadata_ptr_;
        std::shared_ptr<fl::FastLoaderGraph<fl::DefaultView<uint32_t>>> fast_loader_;
        size_t n_threads_;
        std::string fname_;
        
        std::tuple<uint32_t, uint32_t, uint32_t>  GetImageDimensions() const;
        std::tuple<uint32_t, uint32_t, uint32_t>  CalculateTileDimensions() const;
        bool CheckTileStatus() const;
        std::pair<size_t, size_t> GetTileContainingPixel(size_t const row_pixel_index, size_t const col_pixel_index) const;
        void ParseMetadata() const;
        size_t AdjustStride (size_t start_pos, size_t current_pos, size_t stride_val) const;

    public:
        OmeTiffLoader(const std::string &fNameWithPath);
        ~OmeTiffLoader();
 
        std::shared_ptr<std::vector<uint32_t>> GetTileDataByRowCol(size_t const row, size_t const col);
        std::shared_ptr<std::vector<uint32_t>> GetTileDataByIndex(size_t const tile_index);
        std::shared_ptr<std::vector<uint32_t>> GetBoundingBoxVirtualTileData(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const index_col_pixel_min, size_t const index_col_pixel_max);
        std::shared_ptr<std::vector<uint32_t>> GetBoundingBoxVirtualTileDataStrideVersion(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const row_stride, size_t const index_col_pixel_min, size_t const index_col_pixel_max, 
                                                                    size_t const col_stride);
        size_t GetRowTileCount () const;
        size_t GetColumnTileCount () const;
        size_t GetImageHeight() const ;
        size_t GetImageWidth () const ;
        size_t GetTileHeight () const ;
        size_t GetTileWidth () const ;
        short GetSamplesPerPixel () const ;
        short GetBitsPerSamples () const ;

        std::string GetMetaDataValue(const std::string &metadata_key) const;
        void SetViewRequests(size_t const tile_width, size_t const tile_height, size_t const row_stride, size_t const col_stride);
        std::shared_ptr<std::vector<uint32_t>> GetViewRequests(size_t const index_row_pixel_min, size_t const index_row_pixel_max,
                                                                    size_t const index_col_pixel_min, size_t const index_col_pixel_max);
};
#include <string>
#include <tuple>
#include <memory>
#include <vector>
#include <fast_loader/fast_loader.h>
#include <omp.h>
#include <pugixml.hpp>
#include "ome_tiff_tile_loader.h"

class OmeTiffLoader{

    private:
        std::shared_ptr<OmeTiffTileLoader<uint32_t>> tile_loader_;
        mutable std::shared_ptr<std::map<std::string, std::string>> xml_metadata_ptr_;
        std::shared_ptr<fl::FastLoaderGraph<fl::DefaultView<uint32_t>>> fast_loader_;
        size_t n_threads_;
        std::string fname_;
        
        std::tuple<uint32_t, uint32_t, uint32_t>  GetImageDimensions() const;
        std::tuple<uint32_t, uint32_t, uint32_t>  CalculateTileDimensions() const;
        bool CheckTileStatus() const;
        std::pair<size_t, size_t> GetTileContainingPixel(size_t const indexRowPixel, size_t const indexColPixel) const;
        void ParseMetadata() const;
        size_t AdjustStride (size_t startPos, size_t currentPos, size_t strideVal) const;

    public:
        OmeTiffLoader(const std::string &fNameWithPath);
        ~OmeTiffLoader();
 
        std::shared_ptr<std::vector<uint32_t>> GetTileDataByRowCol(size_t const indexRowGlobalTile, size_t const indexColGlobalTile);
        std::shared_ptr<std::vector<uint32_t>> GetTileDataByIndex(size_t const indexGlobalTile);
        std::shared_ptr<std::vector<uint32_t>> GetBoundingBoxVirtualTileData(size_t const indexRowPixelMin, size_t const indexRowPixelMax,
                                                                    size_t const indexColPixelMin, size_t const indexColPixelMax);
        std::shared_ptr<std::vector<uint32_t>> GetBoundingBoxVirtualTileDataStrideVersion(size_t const indexRowPixelMin, size_t const indexRowPixelMax,
                                                                    size_t rowStride, size_t const indexColPixelMin, size_t const indexColPixelMax, 
                                                                    size_t colStride);
        size_t GetRowTileCount () const;
        size_t GetColumnTileCount () const;
        size_t GetImageHeight() const ;
        size_t GetImageWidth () const ;
        size_t GetTileHeight () const ;
        size_t GetTileWidth () const ;

        std::string GetMetaDataValue(const std::string &metadata_key) const;

};
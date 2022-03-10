#include <string>
#include <tuple>
#include <memory>
#include <vector>
#include <fast_loader/specialised_tile_loader/grayscale_tiff_strip_loader.h>
#include <fast_loader/specialised_tile_loader/grayscale_tiff_tile_loader.h>
#include <omp.h>
#include <pugixml.hpp>

class OmeTiffLoader{

    private:
        std::unique_ptr<fl::AbstractTileLoader<fl::DefaultView<uint32_t>>> gsTiffTileLoader;
        std::shared_ptr<std::map<std::string, std::string>> xml_metadata_ptr;
        size_t nThreads;
        std::string fName;
        
        std::tuple<uint32_t, uint32_t, uint32_t>  getImageDimensions() const;
        std::tuple<uint32_t, uint32_t, uint32_t>  calculateTileDimensions() const;
        bool checkTileStatus() const;
        std::pair<size_t, size_t> getTileContainingPixel(size_t const indexRowPixel, size_t const indexColPixel) const;
        void parse_metadata();
        size_t adjustStride (size_t startPos, size_t currentPos, size_t strideVal) const;

    public:
        OmeTiffLoader(const std::string &fNameWithPath);
        ~OmeTiffLoader();
 
        std::shared_ptr<std::vector<uint32_t>> getTileData(size_t const indexRowGlobalTile, size_t const indexColGlobalTile);
        std::shared_ptr<std::vector<uint32_t>> getTileData(size_t const indexGlobalTile);
        std::shared_ptr<std::vector<uint32_t>> getBoundingBoxVirtualTileData(size_t const indexRowPixelMin, size_t const indexRowPixelMax,
                                                                    size_t const indexColPixelMin, size_t const indexColPixelMax);
        std::shared_ptr<std::vector<uint32_t>> getBoundingBoxVirtualTileDataStrideVersion(size_t const indexRowPixelMin, size_t const indexRowPixelMax,
                                                                    size_t rowStride, size_t const indexColPixelMin, size_t const indexColPixelMax, 
                                                                    size_t colStride);
        size_t getRowTileCount () const;
        size_t getColumnTileCount () const;
        size_t getImageHeight() const ;
        size_t getImageWidth () const ;
        size_t getTileHeight () const ;
        size_t getTileWidth () const ;
        
        std::shared_ptr<std::map<std::string, std::string>> get_xml_metadata();

};
// followed the example of grayscale_strip_tiff_loader from fastloader repo


#ifndef BFIO_TILE_LOADER_H
#define BFIO_TILE_LOADER_H
#include <fast_loader/fast_loader.h>

#ifdef __APPLE__
#define uint64 uint64_hack_
#define int64 int64_hack_
#include <tiffio.h>
#undef uint64
#undef int64
#else
#include <tiffio.h>
#endif

/// @brief Tile Loader for 2D OMETiff Grayscale tiff strip files
/// @tparam DataType AbstractView's internal type
template<class DataType>
class BfioTileLoader : public fl::AbstractTileLoader<fl::DefaultView<DataType>> {
public:
  BfioTileLoader(std::string_view const &name, size_t numberThreads, std::string const &filePath)
      : fl::AbstractTileLoader<fl::DefaultView<DataType>>(name, numberThreads, filePath) {}
};

#endif //BFIO_LOADER_H
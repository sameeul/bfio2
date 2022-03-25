// followed the example of grayscale_tile_tiff_loader from fastloader repo

#ifndef OMETIFF_TIFF_TILE_LOADER_H
#define OMETIFF_TIFF_TILE_LOADER_H
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

/// @brief Tile Loader for 2D OMETiff Grayscale tiff files
/// @tparam DataType AbstractView's internal type
template<class DataType>
 class OmeTiffTileLoader : public fl::AbstractTileLoader<fl::DefaultView<DataType>> {
  TIFF *
      tiff_ = nullptr;             ///< Tiff file pointer

  size_t
      full_height_ = 0,           ///< Full height in pixel
      full_width_ = 0,            ///< Full width in pixel
      tile_height_ = 0,            ///< Tile height
      tile_width_ = 0;             ///< Tile width

  short
      sample_format_ = 0,          ///< Sample format as defined by libtiff
      bits_per_sample_ = 0;         ///< Bit Per Sample as defined by libtiff
 public:

  /// @brief OmeTiffTileLoader unique constructor
  /// @param numberThreads Number of threads associated
  /// @param filePath Path of tiff file
  OmeTiffTileLoader(size_t numberThreads, std::string const &filePath)
      : fl::AbstractTileLoader<fl::DefaultView<DataType>>("OmeTiffTileLoader", numberThreads, filePath) {
    short samples_per_pixel = 0;

    // Open the file
    tiff_ = TIFFOpen(filePath.c_str(), "rm");
    if (tiff_ != nullptr) {
      if (TIFFIsTiled(tiff_) == 0) { throw (std::runtime_error("Tile Loader ERROR: The file is not tiled.")); }
      // Load/parse header
      // the tags are uint32_t, so passing a size_t does not work. So, we need to play this game
      uint32_t tmp;
      TIFFGetField(tiff_, TIFFTAG_IMAGEWIDTH, &tmp);
      this->full_width_ = size_t(tmp);
      TIFFGetField(tiff_, TIFFTAG_IMAGELENGTH, &tmp);
      this->full_height_ = size_t(tmp);
      TIFFGetField(tiff_, TIFFTAG_TILEWIDTH, &tmp);
      this->tile_width_ = size_t(tmp);
      TIFFGetField(tiff_, TIFFTAG_TILELENGTH, &tmp);
      this->tile_height_ = size_t(tmp);
      TIFFGetField(tiff_, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
      TIFFGetField(tiff_, TIFFTAG_BITSPERSAMPLE, &(this->bits_per_sample_));
      TIFFGetField(tiff_, TIFFTAG_SAMPLEFORMAT, &(this->sample_format_));

      // Test if the file is greyscale
      if (samples_per_pixel != 1) {
        std::stringstream message;
        message << "Tile Loader ERROR: The file is not greyscale: SamplesPerPixel = " << samples_per_pixel << ".";
        throw (std::runtime_error(message.str()));
      }
      // Interpret undefined data format as unsigned integer data
      if (sample_format_ < 1 || sample_format_ > 3) { sample_format_ = 1; }
    } else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }
  }

  /// @brief OmeTiffTileLoader destructor
  ~OmeTiffTileLoader() override {
    if (tiff_) {
      TIFFClose(tiff_);
      tiff_ = nullptr;
    }
  }

  /// @brief Load a tiff tile from a view
  /// @param tile Tile to copy into
  /// @param index_row_global_tile Tile row index
  /// @param index_col_global_tile Tile column index
  /// @param level Tile's level
  void loadTileFromFile(std::shared_ptr<std::vector<DataType>> tile,
                        size_t index_row_global_tile, size_t index_col_global_tile, [[maybe_unused]] size_t index_layer_global_tile,
                        [[maybe_unused]] size_t level) override {
    tdata_t tiff_tile = nullptr;
    tiff_tile = _TIFFmalloc(TIFFTileSize(tiff_));
    TIFFReadTile(tiff_, tiff_tile, index_col_global_tile * tile_width_, index_row_global_tile * tile_height_, 0, 0);
    std::stringstream message;
    switch (sample_format_) {
      case 1 :
        switch (bits_per_sample_) {
          case 8:loadTile<uint8_t>(tiff_tile, tile);
            break;
          case 16:loadTile<uint16_t>(tiff_tile, tile);
            break;
          case 32:loadTile<uint32_t>(tiff_tile, tile);
            break;
          case 64:loadTile<uint64_t>(tiff_tile, tile);
            break;
          default:
            message
                << "Tile Loader ERROR: The data format is not supported for unsigned integer, number bits per pixel = "
                << bits_per_sample_;
            throw (std::runtime_error(message.str()));
        }
        break;
      case 2:
        switch (bits_per_sample_) {
          case 8:loadTile<int8_t>(tiff_tile, tile);
            break;
          case 16:loadTile<int16_t>(tiff_tile, tile);
            break;
          case 32:loadTile<int32_t>(tiff_tile, tile);
            break;
          case 64:loadTile<int64_t>(tiff_tile, tile);
            break;
          default:
            message
                << "Tile Loader ERROR: The data format is not supported for signed integer, number bits per pixel = "
                << bits_per_sample_;
            throw (std::runtime_error(message.str()));
        }
        break;
      case 3:
        switch (bits_per_sample_) {
          case 8:
          case 16:
          case 32:
            loadTile<float>(tiff_tile, tile);
            break;
          case 64:
            loadTile<double>(tiff_tile, tile);
            break;
          default:
            message
                << "Tile Loader ERROR: The data format is not supported for float, number bits per pixel = "
                << bits_per_sample_;
            throw (std::runtime_error(message.str()));
        }
        break;
      default:message << "Tile Loader ERROR: The data format is not supported, sample format = " << sample_format_;
        throw (std::runtime_error(message.str()));
    }

    _TIFFfree(tiff_tile);
  }

  /// @brief Copy Method for the OmeTiffTileLoader
  /// @return Return a copy of the current OmeTiffTileLoader
  std::shared_ptr<fl::AbstractTileLoader<fl::DefaultView<DataType>>> copyTileLoader() override {
    return std::make_shared<OmeTiffTileLoader<DataType>>(this->numberThreads(), this->filePath());
  }

  /// @brief Tiff file height
  /// @param level Tiff level [not used]
  /// @return Full height
  [[nodiscard]] size_t fullHeight([[maybe_unused]] size_t level) const override { return full_height_; }
  /// @brief Tiff full width
  /// @param level Tiff level [not used]
  /// @return Full width
  [[nodiscard]] size_t fullWidth([[maybe_unused]] size_t level) const override { return full_width_; }
  /// @brief Tiff tile width
  /// @param level Tiff level [not used]
  /// @return Tile width
  [[nodiscard]] size_t tileWidth([[maybe_unused]] size_t level) const override { return tile_width_; }
  /// @brief Tiff tile height
  /// @param level Tiff level [not used]
  /// @return Tile height
  [[nodiscard]] size_t tileHeight([[maybe_unused]] size_t level) const override { return tile_height_; }
  /// @brief Tiff bits per sample
  /// @return Size of a sample in bits
  [[nodiscard]] short bitsPerSample() const override { return bits_per_sample_; }
  /// @brief Level accessor
  /// @return 1
  [[nodiscard]] size_t numberPyramidLevels() const override { return 1; }

 private:
  /// @brief Private function to copy and cast the values
  /// @tparam FileType Type inside the file
  /// @param src Piece of memory coming from libtiff
  /// @param dest Piece of memory to fill
  template<typename FileType>
  void loadTile(tdata_t src, std::shared_ptr<std::vector<DataType>> &dest) {
    for (size_t i = 0; i < tile_height_ * tile_width_; ++i) { dest->data()[i] = (DataType) ((FileType *) (src))[i]; }
  }

};

#endif //OMETIFF_TIFF_TILE_LOADER_H

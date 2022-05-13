// followed the example of grayscale_tile_tiff_loader from fastloader repo

#ifndef OMETIFF_GS_TIFF_TILE_LOADER_H
#define OMETIFF_GS_TIFF_TILE_LOADER_H
#include "bfio_tile_loader.h"
#include <omp.h>

/// @brief Tile Loader for 2D OMETiff Grayscale tiff tile files
/// @tparam DataType AbstractView's internal type
template<class DataType>
 class OmeTiffGrayScaleTileLoader : public BfioTileLoader<DataType> {
  TIFF *
      tiff_ = nullptr;             ///< Tiff file pointer

  size_t
      full_height_ = 0,           ///< Full height in pixel
      full_width_ = 0,            ///< Full width in pixel
      full_depth_ = 0,           ///< Full depth in pixel
      tile_height_ = 0,            ///< Tile height
      tile_width_ = 0,           ///< Tile width
      tile_depth_ = 0;          // Tile Depth

  short
      sample_format_ = 0,          ///< Sample format as defined by libtiff
      bits_per_sample_ = 0,         ///< Bit Per Sample as defined by libtiff
      samples_per_pixel_ = 0;
 public:

  /// @brief OmeTiffGrayScaleTileLoader unique constructor
  /// @param number_threads Number of threads associated
  /// @param file_path Path of tiff file
  OmeTiffGrayScaleTileLoader(size_t number_threads, std::string const &file_path)
      : BfioTileLoader<DataType>("OmeTiffGrayScaleTileLoader", number_threads, file_path) {

    // Open the file
    tiff_ = TIFFOpen(file_path.c_str(), "rm");
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
      TIFFGetField(tiff_, TIFFTAG_SAMPLESPERPIXEL, &(this->samples_per_pixel_));
      TIFFGetField(tiff_, TIFFTAG_BITSPERSAMPLE, &(this->bits_per_sample_));
      TIFFGetField(tiff_, TIFFTAG_SAMPLEFORMAT, &(this->sample_format_));
      full_depth_ = TIFFNumberOfDirectories(tiff_); 
      tile_depth_ = 1;
      // Test if the file is greyscale
      if (samples_per_pixel_ != 1) {
        std::stringstream message;
        message << "Tile Loader ERROR: The file is not greyscale: SamplesPerPixel = " << samples_per_pixel_ << ".";
        throw (std::runtime_error(message.str()));
      }
      // Interpret undefined data format as unsigned integer data
      if (sample_format_ < 1 || sample_format_ > 3) { sample_format_ = 1; }
    } else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }
  }

  /// @brief OmeTiffGrayScaleTileLoader destructor
  ~OmeTiffGrayScaleTileLoader() override {
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
                        size_t index_row_global_tile, size_t index_col_global_tile, size_t index_layer_global_tile,
                        [[maybe_unused]] size_t level) override {
    TIFFSetDirectory(tiff_, index_layer_global_tile);
    TIFFReadTile(tiff_, (void*)(tile->data()), index_col_global_tile * tile_width_, index_row_global_tile * tile_height_, 0, 0);
  }

  /// @brief Copy Method for the OmeTiffGrayScaleTileLoader
  /// @return Return a copy of the current OmeTiffGrayScaleTileLoader
  std::shared_ptr<fl::AbstractTileLoader<fl::DefaultView<DataType>>> copyTileLoader() override {
    return std::make_shared<OmeTiffGrayScaleTileLoader<DataType>>(this->numberThreads(), this->filePath());
  }

  /// @brief Tiff file height
  /// @param level Tiff level [not used]
  /// @return Full height
  [[nodiscard]] size_t fullHeight([[maybe_unused]] size_t level) const override { return full_height_; }
  /// @brief Tiff full width
  /// @param level Tiff level [not used]
  /// @return Full width
  [[nodiscard]] size_t fullWidth([[maybe_unused]] size_t level) const override { return full_width_; }
  [[nodiscard]] size_t fullDepth([[maybe_unused]] size_t level) const override { return full_depth_; }
  /// @brief Tiff tile width
  /// @param level Tiff level [not used]
  /// @return Tile width
  [[nodiscard]] size_t tileWidth([[maybe_unused]] size_t level) const override { return tile_width_; }
  /// @brief Tiff tile height
  /// @param level Tiff level [not used]
  /// @return Tile height
  [[nodiscard]] size_t tileHeight([[maybe_unused]] size_t level) const override { return tile_height_; }
  [[nodiscard]] size_t tileDepth([[maybe_unused]] size_t level) const override { return tile_depth_; }
  /// @brief Tiff bits per sample
  /// @return Size of a sample in bits
  [[nodiscard]] short bitsPerSample() const override { return bits_per_sample_; }
  /// @brief Level accessor
  /// @return 1
  [[nodiscard]] size_t numberPyramidLevels() const override { return 1; }
  /// \brief Getter to the number of channels (default 1)
  /// \return Number of pixel's channels
  [[nodiscard]] virtual size_t numberChannels() const override {
    return samples_per_pixel_;
  }

};

#endif //OMETIFF_GS_TIFF_TILE_LOADER_H

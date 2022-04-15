// followed the example of grayscale_strip_tiff_loader from fastloader repo


#ifndef OMETIFF_GS_TIFF_STRIP_LOADER_H
#define OMETIFF_GS_TIFF_STRIP_LOADER_H
#include "bfio_tile_loader.h"

/// @brief Tile Loader for 2D OMETiff Grayscale tiff strip files
/// @tparam DataType AbstractView's internal type
template<class DataType>
class OmeTiffGrayScaleStripLoader : public BfioTileLoader<DataType> {
  TIFF *
      tiff_ = nullptr;             ///< Tiff file pointer

  size_t
      full_height_ = 0,          ///< Full height in pixel
  full_width_ = 0,           ///< Full width in pixel
  full_depth_ = 0,           ///< Full depth in pixel
  tile_width_ = 0,           ///< Tile width
  tile_height_ = 0,          ///< Tile height
  tile_depth_ = 0;           ///< Tile depth

  short
      sample_format_ = 0,        ///< Sample format as defined by libtiff
      bits_per_sample_ = 0,       ///< Bit Per Sample as defined by libtiff
      samples_per_pixel_ = 0;
 public:
  /// @brief OmeTiffGrayScaleStripLoader constructor
  /// @param number_threads Number of threads associated
  /// @param file_path Path of tiff file
  /// @param tile_width Tile width requested
  /// @param tile_height Tile height requested
  /// @param tile_depth Tile depth requested
  OmeTiffGrayScaleStripLoader(
      size_t number_threads,
      std::string const &file_path,
      size_t tile_width, size_t tile_height, size_t tile_depth)
      : BfioTileLoader<DataType>("OmeTiffGrayScaleStripLoader", number_threads, file_path),
        tile_width_(tile_width), tile_height_(tile_height), tile_depth_(tile_depth) {

    // Open the file
    tiff_ = TIFFOpen(file_path.c_str(), "r");
    if (tiff_ != nullptr) {
      // Load/parse header
      // the tags are uint32_t, so passing a size_t does not work. So, we need to play this game
      uint32_t tmp;
      TIFFGetField(tiff_, TIFFTAG_IMAGEWIDTH, &tmp);
      this->full_width_ = size_t(tmp);
      TIFFGetField(tiff_, TIFFTAG_IMAGELENGTH, &tmp);
      this->full_height_ = size_t(tmp);
      TIFFGetField(tiff_, TIFFTAG_SAMPLESPERPIXEL, &(this->samples_per_pixel_));
      TIFFGetField(tiff_, TIFFTAG_BITSPERSAMPLE, &(this->bits_per_sample_));
      TIFFGetField(tiff_, TIFFTAG_SAMPLEFORMAT, &(this->sample_format_));

      full_depth_ = TIFFNumberOfDirectories(tiff_);

      // Test if the file is grayscale
      // if (samples_per_pixel_ != 1) {
      //   std::stringstream message;
      //   message << "Tile Loader ERROR: The file is not grayscale: SamplesPerPixel = " << samples_per_pixel_ << ".";
      //   throw (std::runtime_error(message.str()));
      // }
      // Interpret undefined data format as unsigned integer data
      if (sample_format_ < 1 || sample_format_ > 3) {
        sample_format_ = 1;
      }
    } else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }
  }

  /// @brief OmeTiffGrayScaleStripLoader destructor
  ~OmeTiffGrayScaleStripLoader() override {
    if (tiff_) {
      TIFFClose(tiff_);
      tiff_ = nullptr;
    }
  }

  /// @brief Load a tiff tile from a view
  /// @param tile Tile to copy into
  /// @param index_row_global_tile Tile row index
  /// @param index_col_global_tile Tile column index
  /// @param index_layer_global_tile Tile layer index
  /// @param level Tile's level
  void loadTileFromFile(std::shared_ptr<std::vector<DataType>> tile,
                        size_t index_row_global_tile,
                        size_t index_col_global_tile,
                        size_t index_layer_global_tile,
                        [[maybe_unused]] size_t level) override {

    tdata_t buf;
    size_t row, layer;

    buf = _TIFFmalloc(TIFFScanlineSize(tiff_));

    size_t
        start_layer = index_layer_global_tile * tile_depth_,
        end_layer = std::min((index_layer_global_tile + 1) * tile_depth_, full_depth_),
        start_row = index_row_global_tile * tile_height_,
        end_row = std::min((index_row_global_tile + 1) * tile_height_, full_height_),
        start_col = index_col_global_tile * tile_width_,
        end_col = std::min((index_col_global_tile + 1) * tile_width_, full_width_);

    for (layer = start_layer; layer < end_layer; ++layer) {
      TIFFSetDirectory(tiff_, layer);
      for (row = start_row; row < end_row; row++) {
        TIFFReadScanline(tiff_, buf, row);
        std::stringstream message;
        switch (sample_format_) {
          case 1 :
            switch (bits_per_sample_) {
              case 8:copyRow<uint8_t>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
                break;
              case 16:copyRow<uint16_t>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
                break;
              case 32:copyRow<size_t>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
                break;
              case 64:copyRow<uint64_t>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
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
              case 8:copyRow<int8_t>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
                break;
              case 16:copyRow<int16_t>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
                break;
              case 32:copyRow<int32_t>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
                break;
              case 64:copyRow<int64_t>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
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
              case 32:copyRow<float>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
                break;
              case 64:copyRow<double>(buf, tile, layer - start_layer, row - start_row, start_col, end_col);
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
      }
    }
    _TIFFfree(buf);
  }

  /// @brief Copy Method for the OmeTiffGrayScaleStripLoader
  /// @return Return a copy of the current OmeTiffGrayScaleStripLoader
  std::shared_ptr<fl::AbstractTileLoader<fl::DefaultView<DataType>>> copyTileLoader() override {
    return std::make_shared<OmeTiffGrayScaleStripLoader<DataType>>(this->numberThreads(),
                                                                this->filePath(),
                                                                this->tile_width_,
                                                                this->tile_height_,
                                                                this->tile_depth_);
  }

  /// @brief Tiff file height
  /// @param level Tiff level [not used]
  /// @return Full height
  [[nodiscard]] size_t fullHeight([[maybe_unused]] size_t level) const override { return full_height_; }
  /// @brief Tiff full width
  /// @param level Tiff level [not used]
  /// @return Full width
  [[nodiscard]] size_t fullWidth([[maybe_unused]] size_t level) const override { return full_width_; }
  /// @brief Tiff full depth
  /// @param level Tiff level [not used]
  /// @return Full Depth
  [[nodiscard]] size_t fullDepth([[maybe_unused]] size_t level) const override { return full_depth_; }

  /// @brief Tiff tile width
  /// @param level Tiff level [not used]
  /// @return Tile width
  [[nodiscard]] size_t tileWidth([[maybe_unused]] size_t level) const override { return tile_width_; }
  /// @brief Tiff tile height
  /// @param level Tiff level [not used]
  /// @return Tile height
  [[nodiscard]] size_t tileHeight([[maybe_unused]] size_t level) const override { return tile_height_; }
  /// @brief Tiff tile depth
  /// @param level Tiff level [not used]
  /// @return Tile depth
  [[nodiscard]] size_t tileDepth([[maybe_unused]] size_t level) const override { return tile_depth_; }

  /// @brief Tiff bits per sample
  /// @return Size of a sample in bits
  [[nodiscard]] short bitsPerSample() const override { return bits_per_sample_; }
  /// @brief Level accessor
  /// @return 1
  [[nodiscard]] size_t numberPyramidLevels() const override { return 1; }
  /// \brief Getter to the number of channels (default 1)
  /// \return Number of pixel's channels
  [[nodiscard]] virtual size_t numberChannels() const {
    return samples_per_pixel_;
  }
 private:
  /// @brief Private function to copy and cast the values
  /// @tparam FileType Type inside the file
  /// @param src Piece of memory coming from libtiff
  /// @param dest Piece of memory to fill
  /// @param layer Destination layer
  /// @param row Destination row
  /// @param start_col Starting column tile to copy
  /// @param end_col End column tile to copy
  template<typename FileType>
  void copyRow(tdata_t src,
               std::shared_ptr<std::vector<DataType>> &dest,
               size_t layer,
               size_t row,
               size_t start_col,
               size_t end_col) {
    for (size_t col = start_col; col < end_col; col++) {
      dest->data()[
          tile_width_ * tile_height_ * layer
              + tile_width_ * row
              + col - start_col] = (DataType) ((FileType *) (src))[col];
    }
  }

};

#endif //OMETIFF_GS_TIFF_STRIP_LOADER_H

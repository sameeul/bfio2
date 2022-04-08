#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "../ome_tiff_loader.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<uint32_t>);

template <typename Sequence>
inline py::array_t<typename Sequence::value_type> as_pyarray_shared(std::shared_ptr<Sequence> seq_ptr) {
    auto size = seq_ptr->size();
    auto data = seq_ptr->data();
    auto capsule = py::capsule(new auto (seq_ptr), [](void *p) {delete reinterpret_cast<decltype(seq_ptr)*>(p);});
    return py::array(size, data, capsule);
 
}

template <typename Sequence>
inline py::array_t<typename Sequence::value_type> as_pyarray_shared_2d(std::shared_ptr<Sequence> seq_ptr, size_t num_rows, size_t num_cols) {
    auto size = seq_ptr->size();
    auto data = seq_ptr->data();
    auto capsule = py::capsule(new auto (seq_ptr), [](void *p) {delete reinterpret_cast<decltype(seq_ptr)*>(p);});
    return py::array(size, data, capsule).reshape({num_rows, num_cols});
 
}

template <typename Sequence>
inline py::array_t<typename Sequence::value_type> as_pyarray_shared_3d(std::shared_ptr<Sequence> seq_ptr, size_t num_rows, size_t num_cols, size_t num_layers) {
    auto size = seq_ptr->size();
    auto data = seq_ptr->data();
    auto capsule = py::capsule(new auto (seq_ptr), [](void *p) {delete reinterpret_cast<decltype(seq_ptr)*>(p);});
    return py::array(size, data, capsule).reshape({num_layers, num_rows, num_cols});
 
}

PYBIND11_MODULE(libbfio2, m) {
  py::class_<OmeTiffLoader, std::shared_ptr<OmeTiffLoader>>(m, "OmeTiffLoader")
    .def(py::init<const std::string &>())
    .def("get_image_height", &OmeTiffLoader::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader::GetColumnTileCount)

    .def("get_bits_per_sample", &OmeTiffLoader::GetBitsPerSamples)

    .def("get_sample_per_pixel", &OmeTiffLoader::GetSamplesPerPixel)

    .def("get_tile_data_2d_by_index",
        [](OmeTiffLoader& tl, size_t const index_global_tile) -> py::array_t<uint32_t> {
            auto tmp = tl.GetTileDataByIndex(index_global_tile);
            auto iw = tl.GetImageWidth();
            auto ih = tl.GetImageHeight();
            auto tw = tl.GetTileWidth();
            auto th = tl.GetTileHeight();
            size_t num_col_tiles = tl.GetColumnTileCount();	
	        size_t row = index_global_tile/num_col_tiles;
	        size_t col = index_global_tile%num_col_tiles;
            auto actual_tw = iw > (col+1)*tw -1 ? tw : iw - col*tw;
	        auto actual_th = ih > (row+1)*th -1 ? th : ih - row*th;
            return as_pyarray_shared_2d(tmp, actual_th, actual_tw) ;;
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer",
        [](OmeTiffLoader& tl, size_t const index_row_global_tile, size_t const index_col_global_tile, size_t const index_layer_global_tile) -> py::array_t<uint32_t> {
            auto tmp = tl.GetTileDataByRowColLayer(index_row_global_tile, index_col_global_tile, index_layer_global_tile);
            auto iw = tl.GetImageWidth();
            auto ih = tl.GetImageHeight();
            auto tw = tl.GetTileWidth();
            auto th = tl.GetTileHeight();
            auto actual_tw = iw > (index_col_global_tile+1)*tw -1 ? tw : iw - index_col_global_tile*tw;
	        auto actual_th = ih > (index_row_global_tile+1)*th -1 ? th : ih - index_row_global_tile*th;
            return as_pyarray_shared_2d(tmp, actual_th, actual_tw) ;;
        }, py::return_value_policy::reference)

        // .def("get_virtual_tile_data_bounding_box_2d",
        // [](OmeTiffLoader& tl, size_t const index_row_min_pixel, size_t const index_row_max_pixel, size_t const index_col_min_pixel, size_t const index_col_max_pixel) -> py::array_t<uint32_t> {
        //     auto tmp = tl.GetBoundingBoxVirtualTileData(index_row_min_pixel, index_row_max_pixel, index_col_min_pixel, index_col_max_pixel);
        //     auto ih = tl.GetImageHeight();
	    //     auto iw = tl.GetImageWidth();
	    //     auto index_true_row_pixel_max = index_row_max_pixel > ih ? ih-1 : index_row_max_pixel;
	    //     auto index_true_col_pixel_max = index_col_max_pixel > iw ? iw-1 : index_col_max_pixel;
        //     size_t num_rows = index_true_row_pixel_max - index_row_min_pixel + 1;
        //     size_t num_cols = index_true_col_pixel_max - index_col_min_pixel + 1;
        //     return as_pyarray_shared_2d(tmp, num_rows, num_cols) ;
        // }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_bounding_box_3d",
        [](OmeTiffLoader& tl, size_t const index_row_min_pixel, size_t const index_row_max_pixel, size_t const index_col_min_pixel, size_t const index_col_max_pixel, size_t const index_layer_min, size_t const index_layer_max, size_t const layer_stride) -> py::array_t<uint32_t> {
            auto tmp = tl.GetBoundingBoxVirtualTileData(index_row_min_pixel, index_row_max_pixel, index_col_min_pixel, index_col_max_pixel, index_layer_min, index_layer_max, layer_stride);
            auto ih = tl.GetImageHeight();
	        auto iw = tl.GetImageWidth();
            auto id = tl.GetImageDepth();
            auto index_true_min_layer = index_layer_min > 0? index_layer_min : 0;
	        auto index_true_max_layer = index_layer_max > id-1? id-1 : index_layer_max;
	        auto index_true_row_pixel_max = index_row_max_pixel > ih ? ih-1 : index_row_max_pixel;
	        auto index_true_col_pixel_max = index_col_max_pixel > iw ? iw-1 : index_col_max_pixel;
            size_t num_rows = index_true_row_pixel_max - index_row_min_pixel + 1;
            size_t num_cols = index_true_col_pixel_max - index_col_min_pixel + 1;
            size_t num_layers = (index_true_max_layer - index_true_min_layer)/layer_stride + 1;
            return as_pyarray_shared_3d(tmp, num_rows, num_cols, num_layers) ;
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_bounding_box_3d_strided",
        [](OmeTiffLoader& tl, size_t const index_row_min_pixel, size_t const index_row_max_pixel, size_t const row_stride,  size_t const index_col_min_pixel, size_t const index_col_max_pixel, size_t const col_stride, size_t const index_layer_min, size_t const index_layer_max, size_t const layer_stride) -> py::array_t<uint32_t> {
            auto tmp = tl.GetBoundingBoxVirtualTileDataStrideVersion(index_row_min_pixel, index_row_max_pixel, row_stride, index_col_min_pixel, index_col_max_pixel, col_stride, index_layer_min, index_layer_max, layer_stride);
            auto ih = tl.GetImageHeight();
	        auto iw = tl.GetImageWidth();
            auto id = tl.GetImageDepth();
	        auto index_true_row_pixel_max = index_row_max_pixel > ih ? ih-1 : index_row_max_pixel;
	        auto index_true_col_pixel_max = index_col_max_pixel > iw ? iw-1 : index_col_max_pixel;
            auto index_true_min_layer = index_layer_min > 0? index_layer_min : 0;
	        auto index_true_max_layer = index_layer_max > id-1? id-1 : index_layer_max;
            size_t num_rows = (index_true_row_pixel_max - index_row_min_pixel)/row_stride+1;
            size_t num_cols = (index_true_col_pixel_max - index_col_min_pixel)/col_stride+1;
            size_t num_layers = (index_true_max_layer - index_true_min_layer)/layer_stride+1;
            return as_pyarray_shared_3d(tmp, num_rows, num_cols, num_layers) ;
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader& tl, const std::string& metadata_key ) -> py::str {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) -> void {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) -> py::array_t<uint32_t> {
            auto tmp = tl.GetViewRequests(index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
            auto ih = tl.GetImageHeight();
	        auto iw = tl.GetImageWidth();
	        auto index_true_row_pixel_max = index_row_pixel_max > ih ? ih-1 : index_row_pixel_max;
	        auto index_true_col_pixel_max = index_col_pixel_max > iw ? iw-1 : index_col_pixel_max;
            size_t num_rows = index_true_row_pixel_max - index_row_pixel_min + 1;
            size_t num_cols = index_true_col_pixel_max - index_col_pixel_min + 1;
            return as_pyarray_shared_2d(tmp, num_rows, num_cols) ;
        }, py::return_value_policy::reference)
        .def("__iter__", [](OmeTiffLoader& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;
}
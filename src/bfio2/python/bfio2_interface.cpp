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

PYBIND11_MODULE(libbfio2, m) {
  py::class_<OmeTiffLoader, std::shared_ptr<OmeTiffLoader>>(m, "OmeTiffLoader")
    .def(py::init<const std::string &>())
    .def("get_image_height", &OmeTiffLoader::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader::GetImageWidth)

    .def("get_tile_height", &OmeTiffLoader::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader::GetTileWidth)

    .def("get_row_tile_count", &OmeTiffLoader::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader::GetColumnTileCount)

    .def("get_tile_data_2d_by_index",
        [](OmeTiffLoader& tl, size_t const indexGlobalTile) -> py::array_t<uint32_t> {
            auto tmp = tl.GetTileDataByIndex(indexGlobalTile);
            return as_pyarray_shared_2d(tmp, tl.GetTileHeight(), tl.GetTileWidth()) ;;
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col",
        [](OmeTiffLoader& tl, size_t const indexRowGlobalTile, size_t const indexColGlobalTile) -> py::array_t<uint32_t> {
            auto tmp = tl.GetTileDataByRowCol(indexRowGlobalTile, indexColGlobalTile);
            return as_pyarray_shared_2d(tmp, tl.GetTileHeight(), tl.GetTileWidth()) ;;
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_bounding_box_2d",
        [](OmeTiffLoader& tl, size_t const indexRowMinPixel, size_t const indexRowMaxPixel, size_t const indexColMinPixel, size_t const indexColMaxPixel) -> py::array_t<uint32_t> {
            auto tmp = tl.GetBoundingBoxVirtualTileData(indexRowMinPixel, indexRowMaxPixel, indexColMinPixel, indexColMaxPixel);
            auto ih = tl.GetImageHeight();
	        auto iw = tl.GetImageWidth();
	        auto indexTrueRowPixelMax = indexRowMaxPixel > ih ? ih : indexRowMaxPixel;
	        auto indexTrueColPixelMax = indexColMaxPixel > iw ? iw : indexColMaxPixel;
            size_t num_rows = indexTrueRowPixelMax - indexRowMinPixel + 1;
            size_t num_cols = indexTrueColPixelMax - indexColMinPixel + 1;
            return as_pyarray_shared_2d(tmp, num_rows, num_cols) ;
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_bounding_box_2d_strided",
        [](OmeTiffLoader& tl, size_t const indexRowMinPixel, size_t const indexRowMaxPixel, size_t const rowStride,  size_t const indexColMinPixel, size_t const indexColMaxPixel, size_t const colStride) -> py::array_t<uint32_t> {
            auto tmp = tl.GetBoundingBoxVirtualTileDataStrideVersion(indexRowMinPixel, indexRowMaxPixel, rowStride, indexColMinPixel, indexColMaxPixel, colStride);
            auto ih = tl.GetImageHeight();
	        auto iw = tl.GetImageWidth();
	        auto indexTrueRowPixelMax = indexRowMaxPixel > ih ? ih : indexRowMaxPixel;
	        auto indexTrueColPixelMax = indexColMaxPixel > iw ? iw : indexColMaxPixel;
            size_t num_rows = (indexTrueRowPixelMax - indexRowMinPixel)/rowStride+1;
            size_t num_cols = (indexTrueColPixelMax - indexColMinPixel)/colStride+1;
            return as_pyarray_shared_2d(tmp, num_rows, num_cols) ;
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader& tl, const std::string& metadata_key ) -> py::str {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal);
}
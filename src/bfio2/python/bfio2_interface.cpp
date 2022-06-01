#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "../ome_tiff_loader.h"
#include "../utilities.h"

namespace py = pybind11;
template <typename SampleType>
inline py::array_t<SampleType> as_pyarray_shared_5d(std::shared_ptr<std::vector<SampleType>> seq_ptr, size_t num_rows, size_t num_cols, size_t num_layers=1, size_t num_channels=1, size_t num_tsteps=1 ) {
    auto size = seq_ptr->size();
    auto data = seq_ptr->data();
    auto capsule = py::capsule(new auto (seq_ptr), [](void *p) {delete reinterpret_cast<decltype(seq_ptr)*>(p);});
    return py::array(size, data, capsule).reshape({num_tsteps, num_channels, num_layers, num_rows, num_cols}).squeeze();
 
}


template <typename SampleType>
py::array_t<SampleType>  get_virtual_tile_data(OmeTiffLoader<SampleType>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
    auto tmp = tl.GetVirtualTileData(rows, cols, layers, channels, tsteps);
    auto ih = tl.GetImageHeight();
    auto iw = tl.GetImageWidth();
    auto id = tl.GetImageDepth();
    auto nc = tl.GetChannelCount();
    auto nt = tl.GetTstepCount();
    auto index_true_row_pixel_max = rows.Stop() > ih ? ih-1 : rows.Stop();
    auto index_true_col_pixel_max = cols.Stop() > iw ? iw-1 : cols.Stop();
    auto index_true_min_layer = layers.Start() > 0? layers.Start() : 0;
    auto index_true_max_layer = layers.Stop() > id-1? id-1 : layers.Stop();
    auto index_true_min_channel = channels.Start() > 0? channels.Start() : 0;
    auto index_true_max_channel = channels.Stop() > nc-1? nc-1 : channels.Stop();
    auto index_true_min_tstep = tsteps.Start() > 0? tsteps.Start() : 0;
    auto index_true_max_tstep = tsteps.Stop() > nt-1? nt-1 : tsteps.Stop();

    size_t num_rows = (index_true_row_pixel_max - rows.Start())/rows.Step()+1;
    size_t num_cols = (index_true_col_pixel_max - cols.Start())/cols.Step()+1;
    size_t num_layers = (index_true_max_layer - index_true_min_layer)/layers.Step()+1;
    size_t num_channels = (index_true_max_channel - index_true_min_channel)/channels.Step()+1;
    size_t num_tsteps = (index_true_max_tstep-index_true_min_tstep)/tsteps.Step()+1;
    return as_pyarray_shared_5d(tmp, num_rows, num_cols, num_layers, num_channels, num_tsteps) ;
}

template <typename SampleType>
py::array_t<SampleType>  get_virtual_tile_data_strided(OmeTiffLoader<SampleType>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
    auto tmp = tl.GetVirtualTileDataStrided(rows, cols, layers, channels, tsteps);
    auto ih = tl.GetImageHeight();
    auto iw = tl.GetImageWidth();
    auto id = tl.GetImageDepth();
    auto nc = tl.GetChannelCount();
    auto nt = tl.GetTstepCount();
    auto index_true_row_pixel_max = rows.Stop() > ih ? ih-1 : rows.Stop();
    auto index_true_col_pixel_max = cols.Stop() > iw ? iw-1 : cols.Stop();
    auto index_true_min_layer = layers.Start() > 0? layers.Start() : 0;
    auto index_true_max_layer = layers.Stop() > id-1? id-1 : layers.Stop();
    auto index_true_min_channel = channels.Start() > 0? channels.Start() : 0;
    auto index_true_max_channel = channels.Stop() > nc-1? nc-1 : channels.Stop();
    auto index_true_min_tstep = tsteps.Start() > 0? tsteps.Start() : 0;
    auto index_true_max_tstep = tsteps.Stop() > nt-1? nt-1 : tsteps.Stop();

    size_t num_rows = (index_true_row_pixel_max - rows.Start())/rows.Step()+1;
    size_t num_cols = (index_true_col_pixel_max - cols.Start())/cols.Step()+1;
    size_t num_layers = (index_true_max_layer - index_true_min_layer)/layers.Step()+1;
    size_t num_channels = (index_true_max_channel - index_true_min_channel)/channels.Step()+1;
    size_t num_tsteps = (index_true_max_tstep-index_true_min_tstep)/tsteps.Step()+1;
    return as_pyarray_shared_5d(tmp, num_rows, num_cols, num_layers, num_channels, num_tsteps) ;
}


template <typename SampleType>
py::array_t<SampleType> get_tile_data_2d_by_index_channel(OmeTiffLoader<SampleType>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
    auto tmp = tl.GetTileDataByIndex(index_global_tile, channel, tstep);
    auto iw = tl.GetImageWidth();
    auto ih = tl.GetImageHeight();
    auto tw = tl.GetTileWidth();
    auto th = tl.GetTileHeight();
    size_t num_col_tiles = tl.GetColumnTileCount();	
    size_t row = index_global_tile/num_col_tiles;
    size_t col = index_global_tile%num_col_tiles;
    auto actual_tw = iw > (col+1)*tw -1 ? tw : iw - col*tw;
    auto actual_th = ih > (row+1)*th -1 ? th : ih - row*th;
    return as_pyarray_shared_5d(tmp, actual_th, actual_tw, 1, 1, 1);
}
template <typename SampleType>
py::array_t<SampleType> get_tile_data_2d_by_row_col_layer_channel(OmeTiffLoader<SampleType>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
    auto tmp = tl.GetTileData(row, col, layer, channel, tstep);
    auto iw = tl.GetImageWidth();
    auto ih = tl.GetImageHeight();
    auto tw = tl.GetTileWidth();
    auto th = tl.GetTileHeight();
    auto actual_tw = iw > (col+1)*tw -1 ? tw : iw - col*tw;
    auto actual_th = ih > (row+1)*th -1 ? th : ih - row*th;
    return as_pyarray_shared_5d(tmp, actual_th, actual_tw, 1, 1, 1);
}
template <typename SampleType>
py::array_t<SampleType> get_iterator_requested_tile_data(OmeTiffLoader<SampleType>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
    auto tmp = tl.GetViewRequests(index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
    auto ih = tl.GetImageHeight();
    auto iw = tl.GetImageWidth();
    auto index_true_row_pixel_max = index_row_pixel_max > ih ? ih-1 : index_row_pixel_max;
    auto index_true_col_pixel_max = index_col_pixel_max > iw ? iw-1 : index_col_pixel_max;
    size_t num_rows = index_true_row_pixel_max - index_row_pixel_min + 1;
    size_t num_cols = index_true_col_pixel_max - index_col_pixel_min + 1;
    return as_pyarray_shared_5d(tmp, num_rows, num_cols, 1, 1, 1) ;
}
PYBIND11_MODULE(libbfio2, m) {

  m.def("get_tiff_type", &GetTiffType);

  py::class_<Seq, std::shared_ptr<Seq>>(m, "Seq")  
    .def(py::init<const size_t, const size_t, const size_t>());

  py::class_<OmeTiffLoader<uint16_t>, std::shared_ptr<OmeTiffLoader<uint16_t>>>(m, "OmeTiffLoaderUint16")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<uint16_t>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<uint16_t>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<uint16_t>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<uint16_t>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<uint16_t>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<uint16_t>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<uint16_t>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<uint16_t>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<uint16_t>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<uint16_t>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<uint16_t>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<uint16_t>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<uint16_t>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<uint16_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<uint16_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<uint16_t>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<uint16_t>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<uint16_t>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<uint16_t>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<uint32_t>, std::shared_ptr<OmeTiffLoader<uint32_t>>>(m, "OmeTiffLoaderUint32")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<uint32_t>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<uint32_t>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<uint32_t>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<uint32_t>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<uint32_t>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<uint32_t>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<uint32_t>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<uint32_t>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<uint32_t>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<uint32_t>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<uint32_t>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<uint32_t>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<uint32_t>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<uint32_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<uint32_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<uint32_t>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<uint32_t>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<uint32_t>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<uint32_t>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<uint8_t>, std::shared_ptr<OmeTiffLoader<uint8_t>>>(m, "OmeTiffLoaderUint8")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<uint8_t>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<uint8_t>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<uint8_t>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<uint8_t>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<uint8_t>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<uint8_t>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<uint8_t>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<uint8_t>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<uint8_t>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<uint8_t>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<uint8_t>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<uint8_t>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<uint8_t>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<uint8_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<uint8_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<uint8_t>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<uint8_t>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<uint8_t>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<uint8_t>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<uint64_t>, std::shared_ptr<OmeTiffLoader<uint64_t>>>(m, "OmeTiffLoaderUint64")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<uint64_t>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<uint64_t>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<uint64_t>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<uint64_t>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<uint64_t>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<uint64_t>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<uint64_t>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<uint64_t>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<uint64_t>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<uint64_t>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<uint64_t>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<uint64_t>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<uint64_t>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<uint64_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<uint64_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<uint64_t>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<uint64_t>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<uint64_t>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<uint64_t>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<int8_t>, std::shared_ptr<OmeTiffLoader<int8_t>>>(m, "OmeTiffLoaderInt8")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<int8_t>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<int8_t>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<int8_t>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<int8_t>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<int8_t>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<int8_t>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<int8_t>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<int8_t>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<int8_t>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<int8_t>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<int8_t>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<int8_t>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<int8_t>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<int8_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<int8_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<int8_t>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<int8_t>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<int8_t>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<int8_t>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<int16_t>, std::shared_ptr<OmeTiffLoader<int16_t>>>(m, "OmeTiffLoaderInt16")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<int16_t>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<int16_t>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<int16_t>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<int16_t>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<int16_t>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<int16_t>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<int16_t>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<int16_t>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<int16_t>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<int16_t>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<int16_t>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<int16_t>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<int16_t>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<int16_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<int16_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<int16_t>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<int16_t>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<int16_t>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<int16_t>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<int32_t>, std::shared_ptr<OmeTiffLoader<int32_t>>>(m, "OmeTiffLoaderInt32")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<int32_t>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<int32_t>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<int32_t>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<int32_t>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<int32_t>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<int32_t>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<int32_t>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<int32_t>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<int32_t>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<int32_t>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<int32_t>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<int32_t>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<int32_t>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<int32_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<int32_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<int32_t>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<int32_t>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<int32_t>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<int32_t>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<int64_t>, std::shared_ptr<OmeTiffLoader<int64_t>>>(m, "OmeTiffLoaderInt64")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<int64_t>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<int64_t>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<int64_t>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<int64_t>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<int64_t>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<int64_t>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<int64_t>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<int64_t>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<int64_t>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<int64_t>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<int64_t>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<int64_t>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<int64_t>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<int64_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<int64_t>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<int64_t>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<int64_t>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<int64_t>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<int64_t>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<float>, std::shared_ptr<OmeTiffLoader<float>>>(m, "OmeTiffLoaderFloat")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<float>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<float>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<float>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<float>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<float>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<float>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<float>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<float>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<float>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<float>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<float>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<float>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<float>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<float>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<float>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<float>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<float>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<float>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<float>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;

  py::class_<OmeTiffLoader<double>, std::shared_ptr<OmeTiffLoader<double>>>(m, "OmeTiffLoaderDouble")
    .def(py::init<const std::string &, const int>())
    .def("get_image_height", &OmeTiffLoader<double>::GetImageHeight)

    .def("get_image_width", &OmeTiffLoader<double>::GetImageWidth)

    .def("get_image_depth", &OmeTiffLoader<double>::GetImageDepth)

    .def("get_tile_height", &OmeTiffLoader<double>::GetTileHeight)

    .def("get_tile_width", &OmeTiffLoader<double>::GetTileWidth)

    .def("get_tile_depth", &OmeTiffLoader<double>::GetTileDepth)

    .def("get_row_tile_count", &OmeTiffLoader<double>::GetRowTileCount)

    .def("get_column_tile_count", &OmeTiffLoader<double>::GetColumnTileCount)

    .def("get_channel_count", &OmeTiffLoader<double>::GetChannelCount)

    .def("get_tstep_count", &OmeTiffLoader<double>::GetTstepCount)

    .def("get_bits_per_sample", &OmeTiffLoader<double>::GetBitsPerSamples)

    .def("get_tile_data_2d_by_index_channel",
        [](OmeTiffLoader<double>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep);
        }, py::return_value_policy::reference)

    .def("get_tile_data_2d_by_row_col_layer_channel",
        [](OmeTiffLoader<double>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep);
        }, py::return_value_policy::reference)

        .def("get_virtual_tile_data_5d",
        [](OmeTiffLoader<double>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  {
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

       .def("get_virtual_tile_data_5d_strided",
        [](OmeTiffLoader<double>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps);
        }, py::return_value_policy::reference)

        .def("get_metadata_value",
        [](OmeTiffLoader<double>& tl, const std::string& metadata_key ) {
            return tl.GetMetaDataValue(metadata_key);
        }, py::return_value_policy::reference_internal)

        .def("send_iterator_view_requests",
        [](OmeTiffLoader<double>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) {
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride);
        })

        .def("get_iterator_requested_tile_data",
        [](OmeTiffLoader<double>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max);
        }, py::return_value_policy::reference)

        .def("__iter__", [](OmeTiffLoader<double>& tl){ 
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end());
            }, 
            py::keep_alive<0, 1>()); // Keep vector alive while iterator is used 
        ;
}
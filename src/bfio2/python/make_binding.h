#define makeBinding(SampleType, ReaderType, Name)	\
  py::class_<ReaderType<SampleType>, std::shared_ptr<ReaderType<SampleType>>>(m, Name) \
    .def(py::init<const std::string &, const int>()) \
    .def("get_image_height", &ReaderType<SampleType>::GetImageHeight) \
    .def("get_image_width", &ReaderType<SampleType>::GetImageWidth) \
    .def("get_image_depth", &ReaderType<SampleType>::GetImageDepth) \
    .def("get_tile_height", &ReaderType<SampleType>::GetTileHeight) \
    .def("get_tile_width", &ReaderType<SampleType>::GetTileWidth) \
    .def("get_tile_depth", &ReaderType<SampleType>::GetTileDepth) \
    .def("get_row_tile_count", &ReaderType<SampleType>::GetRowTileCount) \
    .def("get_column_tile_count", &ReaderType<SampleType>::GetColumnTileCount) \
    .def("get_channel_count", &ReaderType<SampleType>::GetChannelCount) \
    .def("get_tstep_count", &ReaderType<SampleType>::GetTstepCount) \
    .def("get_bits_per_sample", &ReaderType<SampleType>::GetBitsPerSamples) \
    .def("get_tile_data_2d_by_index_channel",  \
        [](ReaderType<SampleType>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) { \
            return get_tile_data_2d_by_index_channel(tl, index_global_tile, channel, tstep); \
        }, py::return_value_policy::reference) \
    .def("get_tile_data_2d_by_row_col_layer_channel", \
        [](ReaderType<SampleType>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) { \
            return get_tile_data_2d_by_row_col_layer_channel(tl, row, col, layer, channel, tstep); \
        }, py::return_value_policy::reference) \
        .def("get_virtual_tile_data_5d", \
        [](ReaderType<SampleType>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps)  { \
            return get_virtual_tile_data(tl, rows, cols, layers, channels, tsteps); \
        }, py::return_value_policy::reference) \
       .def("get_virtual_tile_data_5d_strided", \
        [](ReaderType<SampleType>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) { \
            return get_virtual_tile_data_strided(tl, rows, cols, layers, channels, tsteps); \
        }, py::return_value_policy::reference) \
        .def("get_metadata_value", \
        [](ReaderType<SampleType>& tl, const std::string& metadata_key ) { \
            return tl.GetMetaDataValue(metadata_key); \
        }, py::return_value_policy::reference_internal) \
        .def("send_iterator_view_requests", \
        [](ReaderType<SampleType>& tl, size_t const tile_height, size_t const tile_width, size_t const row_stride, size_t const col_stride) { \
            tl.SetViewRequests(tile_height, tile_width, row_stride, col_stride); \
        }) \
        .def("get_iterator_requested_tile_data", \
        [](ReaderType<SampleType>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) { \
            return get_iterator_requested_tile_data(tl, index_row_pixel_min, index_row_pixel_max, index_col_pixel_min, index_col_pixel_max); \
        }, py::return_value_policy::reference) \
        .def("__iter__", [](ReaderType<SampleType>& tl){ \
            return py::make_iterator(tl.tile_coordinate_list_.begin(), tl.tile_coordinate_list_.end()); \
            }, py::keep_alive<0, 1>()); 
   
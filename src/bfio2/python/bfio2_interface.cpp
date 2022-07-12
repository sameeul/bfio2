#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "../ome_tiff_loader.h"
#include "../ome_zarr_loader.h"
#include "../utilities.h"
#include "make_binding.h"

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

template <typename SampleType>
py::array_t<SampleType>  get_virtual_tile_data(OmeZarrLoader<SampleType>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
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
py::array_t<SampleType>  get_virtual_tile_data_strided(OmeZarrLoader<SampleType>& tl, const Seq& rows, const Seq& cols, const Seq& layers, const Seq& channels, const Seq& tsteps) {
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
py::array_t<SampleType> get_tile_data_2d_by_index_channel(OmeZarrLoader<SampleType>& tl, size_t const index_global_tile, size_t const channel, size_t const tstep) {
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
py::array_t<SampleType> get_tile_data_2d_by_row_col_layer_channel(OmeZarrLoader<SampleType>& tl, size_t const row, size_t const col, size_t const layer, size_t const channel, size_t const tstep) {
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
py::array_t<SampleType> get_iterator_requested_tile_data(OmeZarrLoader<SampleType>& tl, size_t const index_row_pixel_min, size_t const index_row_pixel_max, size_t const index_col_pixel_min, size_t const index_col_pixel_max) {
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
  m.def("get_zarr_type", &GetZarrType);

  py::class_<Seq, std::shared_ptr<Seq>>(m, "Seq")  
    .def(py::init<const size_t, const size_t, const size_t>());

    //OmeTiff
    makeBinding(uint8_t, OmeTiffLoader, "OmeTiffLoaderUint8");
    makeBinding(uint16_t, OmeTiffLoader, "OmeTiffLoaderUint16");
    makeBinding(uint32_t, OmeTiffLoader, "OmeTiffLoaderUint32");
    makeBinding(uint64_t, OmeTiffLoader, "OmeTiffLoaderUint64");
    makeBinding(int8_t, OmeTiffLoader, "OmeTiffLoaderInt8");
    makeBinding(int16_t, OmeTiffLoader, "OmeTiffLoaderInt16");
    makeBinding(int32_t, OmeTiffLoader, "OmeTiffLoaderInt32");
    makeBinding(int64_t, OmeTiffLoader, "OmeTiffLoaderInt64");
    makeBinding(float, OmeTiffLoader, "OmeTiffLoaderFloat");
    makeBinding(double, OmeTiffLoader, "OmeTiffLoaderDouble");
    //OmeZarr
    makeBinding(uint8_t, OmeZarrLoader, "OmeZarrLoaderUint8");
    makeBinding(uint16_t, OmeZarrLoader, "OmeZarrLoaderUint16");
    makeBinding(uint32_t, OmeZarrLoader, "OmeZarrLoaderUint32");
    makeBinding(uint64_t, OmeZarrLoader, "OmeZarrLoaderUint64");
    makeBinding(int8_t, OmeZarrLoader, "OmeZarrLoaderInt8");
    makeBinding(int16_t, OmeZarrLoader, "OmeZarrLoaderInt16");
    makeBinding(int32_t, OmeZarrLoader, "OmeZarrLoaderInt32");
    makeBinding(int64_t, OmeZarrLoader, "OmeZarrLoaderInt64");
    makeBinding(float, OmeZarrLoader, "OmeZarrLoaderFloat");
    makeBinding(double, OmeZarrLoader, "OmeZarrLoaderDouble");

}
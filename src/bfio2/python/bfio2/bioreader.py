from . import libbfio2 as bfio2
from .libbfio2 import Seq, get_tiff_type,\
                    BioReaderUint8, \
                    BioReaderUint16, \
                    BioReaderUint32, \
                    BioReaderUint64, \
                    BioReaderInt8, \
                    BioReaderInt16, \
                    BioReaderInt32, \
                    BioReaderInt64, \
                    BioReaderFloat, \
                    BioReaderDouble

import numpy


class BioReader:

    READ_ONLY_MESSAGE = "{} is read-only."
    


    def __init__(self, file_name, num_threads=1):

        img_cls_dict_ = {
            "uint8_t":BioReaderUint8,
            "uint16_t":BioReaderUint16,
            "uint32_t":BioReaderUint32,
            "uint64_t":BioReaderUint64,
            "int8_t":BioReaderInt8,
            "int16_t":BioReaderInt16,
            "int32_t":BioReaderInt32,
            "int64_t":BioReaderInt64,
            "float":BioReaderFloat,
            "double":BioReaderDouble,
        }

        self._file_name = file_name
        self._file_type = ""
        self._DIMS = {}
        if file_name.endswith('.ome.tif') or file_name.endswith('.ome.tiff'):
            self._file_type = "ome_tiff"
        elif file_name.endswith('.ome.zarr'):
            self._file_type = "ome_zarr"
        else:
            raise TypeError("Only OMETiff and OMEZarr file formats are supported")

        data_type = "uint16_t" # default
        if self._file_type == "ome_zarr":
            data_type = bfio2.get_zarr_type(file_name)
        else: 
            data_type = bfio2.get_tiff_type(file_name)
        self._image_reader = img_cls_dict_[data_type](file_name, num_threads)
        
        self._Y = self._image_reader.get_image_height()
        self._X = self._image_reader.get_image_width()
        self._Z = self._image_reader.get_image_depth()
        self._C = self._image_reader.get_channel_count()
        self._T = self._image_reader.get_tstep_count()
        self._DIMS['X'] = self._X
        self._DIMS['Y'] = self._Y
        self._DIMS['Z'] = self._Z
        self._DIMS['C'] = self._C
        self._DIMS['T'] = self._T

        
        self._physical_size_x = None
        self._physical_size_y = None
        self._physical_size_z = None
        self._physical_size_x_unit = None
        self._physical_size_y_unit = None
        self._physical_size_z_unit = None
        self._samples_per_pixel = None
        self._bits_per_sample = None
        self._tile_height = None
        self._tile_width = None
        self._row_tile_count = None
        self._column_tile_count = None
    
    def data(self, row=None, col=None, layer=0, channel=0, tstep=0):
        if col != None:
            return self._image_reader.get_tile_data_2d_by_row_col_layer_channel(row, col, layer, channel, tstep)
        else:
            return self._image_reader.get_tile_data_2d_by_index_channel(row, channel, tstep)

    @property
    def Y(self):
        return self._Y
    @Y.setter
    def Y(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def X(self):
        return self._X

    @X.setter
    def X(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))


    @property
    def Z(self):
        return self._Z

    @Z.setter
    def Z(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def C(self):
        return self._C

    @C.setter
    def C(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def shape(self):
        return (self._X, self._Y, self._Z, self._C, self._T)

    @shape.setter
    def shape(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    def image_size(self):
        return (self._X, self._Y)

    @property    
    def tile_height(self):
        if self._tile_height == None:
            self._tile_height = self._image_reader.get_tile_height()
        return self._tile_height

    @tile_height.setter
    def tile_height(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property    
    def tile_width(self):
        if self._tile_width == None:
            self._tile_width = self._image_reader.get_tile_width()
        return self._tile_width


    @tile_width.setter
    def tile_width(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property  
    def row_tile_count(self):
        if self._row_tile_count == None:
            self._row_tile_count = self._image_reader.get_row_tile_count()
        return self._row_tile_count

    @row_tile_count.setter
    def row_tile_count(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def column_tile_count(self):
        if self._column_tile_count == None:
            self._column_tile_count = self._image_reader.get_column_tile_count()
        return self._column_tile_count

    @column_tile_count.setter
    def column_tile_count(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))


    def channel_names(self) :
        pass

    @property
    def physical_size_x(self):
        if self._physical_size_x == None:
            self._physical_size_x = self._image_reader.get_metadata_value("PhysicalSizeX")

        if self._physical_size_x_unit == None:
            self._physical_size_x_unit = self._image_reader.get_metadata_value("PhysicalSizeXUnit")

        return (self._physical_size_x, self._physical_size_x_unit)

    @physical_size_x.setter
    def physical_size_x(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def physical_size_y(self):
        if self._physical_size_y == None:
            self._physical_size_y = self._image_reader.get_metadata_value("PhysicalSizeY")

        if self._physical_size_y_unit == None:
            self._physical_size_y_unit = self._image_reader.get_metadata_value("PhysicalSizeYUnit")

        return (self._physical_size_y, self._physical_size_y_unit)

    @physical_size_y.setter
    def physical_size_y(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def physical_size_z(self):
        if self._physical_size_z == None:
            self._physical_size_z = self._image_reader.get_metadata_value("PhysicalSizeZ")

        if self._physical_size_z_unit == None:
            self._physical_size_z_unit = self._image_reader.get_metadata_value("PhysicalSizeZUnit")

        return (self._physical_size_z, self._physical_size_z_unit)

    @physical_size_z.setter
    def physical_size_z(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))
    
    @property
    def samples_per_pixel(self):
        if self._samples_per_pixel == None:
            self._samples_per_pixel = self._image_reader.get_channel_count() 

        return int(self._samples_per_pixel)

    @samples_per_pixel.setter
    def samples_per_pixel(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def bits_per_sample(self):
        if self._bits_per_sample == None:
            self._bits_per_sample = self._image_reader.get_bits_per_sample() 

        return int(self._bits_per_sample)

    @bits_per_sample.setter
    def bits_per_sample(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    @property
    def bytes_per_pixel(self):
        return self.bits_per_sample*self.samples_per_pixel/8

    @bytes_per_pixel.setter
    def bytes_per_pixel(self):
        raise AttributeError(self._READ_ONLY_MESSAGE.format("read_only"))

    def __getitem__(self, keys):
        slice_items = (self._parse_slice(keys))
        X = self._val_dims(slice_items['X'], 'X')
        Y = self._val_dims(slice_items['Y'], 'Y')
        Z = self._val_dims(slice_items['Z'], 'Z')
        C = self._val_dims(slice_items['C'], 'C')
        T = self._val_dims(slice_items['T'], 'T')

        row_min = Y[0]
        row_max = Y[1]-1
        if row_max < 0:
           row_max += self._DIMS['Y']

        if len(Y) == 3:
            row_step = Y[2]
        else:
            row_step = 1

        rows = Seq(row_min, row_max, row_step)

        col_min = X[0]
        col_max = X[1]-1
        if col_max < 0:
            col_max += self._DIMS['X']
 
        if len(X) == 3:
            col_step = X[2]
        else:
            col_step = 1
        cols = Seq(col_min, col_max, col_step)

        layer_min = Z[0]
        layer_max = Z[1]-1
        if layer_max < 0:
            layer_max += self._DIMS['Z']
        if len(Z) == 3:
            layer_step = Z[2]
        else:
            layer_step = 1
        layers = Seq(layer_min, layer_max, layer_step)

        channel_min = C[0]
        channel_max = C[1]-1
        if channel_max < 0:
            channel_max += self._DIMS['C']
        if len(C) == 3:
            channel_step = C[2]
        else:
            channel_step = 1
        channels = Seq(channel_min, channel_max, channel_step)

        time_min = T[0]
        time_max = T[1]-1
        if time_max < 0:
            time_max += self._DIMS['T']
        if len(T) == 3:
            time_step = C[2]
        else:
            time_step = 1
        times = Seq(time_min, time_max, time_step)

        if row_step != 1 or col_step != 1:
            return self._image_reader.get_virtual_tile_data_5d_strided(rows, cols, layers, channels, times)
        else:
            return self._image_reader.get_virtual_tile_data_5d(rows, cols, layers, channels, times)

    def _parse_slice(self,keys):
        # Dimension ordering and index initialization
        dims = 'YXZCT'
        ind = {d:None for d in dims}
        
        # If an empty slice, load the whole image
        # may be aslo setup a failsafe / warning for too large image
        if not isinstance(keys,tuple):
            if isinstance(keys,slice) and keys.start == None and keys.stop == None and keys.step==None:
                pass
            else:
                raise ValueError('If only a single index is supplied, it must be an empty slice ([:]).')
            
        # If not an empty slice, parse the key tuple
        else:
            
            # At most, 5 indices can be indicated
            if len(keys) > 5:
                raise ValueError('Found {} indices, but at most 3 indices may be supplied.'.format(len(keys)))
            
            # If the first key is an ellipsis, read backwards
            if keys[0] == Ellipsis:
                
                # If the last key is an ellipsis, throw an error
                if keys[-1]==Ellipsis:
                    raise ValueError('Ellipsis (...) may be used in either the first or last index, not both.')
                
                dims = ''.join([d for d in reversed(dims)])
                keys = [k for k in reversed(keys)]
            
            # Get key values
            for dim,key in zip(dims,keys):
                if isinstance(key,slice):
                    start = 0 if key.start == None else key.start
                    stop = self._DIMS.get(dim, 0) if key.stop == None else key.stop
                    
                    # For CT dimensions, generate a list from slices
                    if dim in 'XYZCT':
                        step = 1 if key.step is None else key.step
                        ind[dim] = [start,stop, step]

                elif key==Ellipsis:
                    if dims.find(dim)+1 < len(keys):
                        raise ValueError('Ellipsis may only be used in the first or last index.')
                    
                elif isinstance(key,(int,tuple,list)) or numpy.issubdtype(key,numpy.integer):
                    # Only the last three dimensions can use int, tuple, or list indexing
                    if dim in 'XYZCT':
                        if isinstance(key,int) or numpy.issubdtype(key,numpy.integer):
                            ind[dim] = [int(key),int(key)+1]
                    else:
                        raise ValueError('The index in position {} must be a slice type.'.format(dims.find(dim)))
                else:
                    raise ValueError('Did not recognize indexing value of type: {}'.format(type(key)))

        return ind

    def __call__(self, tile_size, tile_stride=None) :
        # Iterate through tiles of an image
        self._iter_tile_size = tile_size
        self._iter_tile_stride = tile_stride        
        return self

    def __iter__(self):
        tile_size = self._iter_tile_size
        tile_stride = self._iter_tile_stride


        if tile_size is None:
            raise SyntaxError(
                "Cannot directly iterate over a BioReader object."
                + "Call it (i.e. for i in bioreader(256,256))"
            )

        # input error checking
        assert len(tile_size) == 2, "tile_size must be a list with 2 elements"
        if tile_stride is not None:
            assert len(tile_stride) == 2, "stride must be a list with 2 elements"
        else:
            tile_stride = tile_size

        self._iter_tile_size = None
        self._iter_tile_stride = None
        # request views
        self._image_reader.send_iterator_view_requests(tile_size[0], tile_size[1], tile_stride[0], tile_stride[1])
        # collect views
        tile_count = 0
        col_per_row = self._image_reader.get_column_tile_count()
        for data in self._image_reader.__iter__():
            row_index = int(tile_count/col_per_row)
            col_index = tile_count%col_per_row
            yield (row_index, col_index), self._image_reader.get_iterator_requested_tile_data(data[0], data[1], data[2], data[3])
            tile_count+=1


    def test_iter(self):
        for data in self._image_reader.__iter__():
            yield data

    """ -------------------- """
    """ -Validation methods- """
    """ -------------------- """
    
    def _val_dims(self, xyz, axis):
        assert axis in 'XYZCT'
        if xyz == None:
            xyz = [0,self._DIMS[axis]]
        else:
            assert len(xyz) == 2 or len(xyz) == 3, \
                '{} must be a list or tuple of length 2 or 3.'.format(axis)
            assert xyz[0] >= 0, \
                '{}[0] must be greater than or equal to 0.'.format(axis)
            assert xyz[1] <= self._DIMS[axis], \
                '{}[1] cannot be greater than the maximum of the dimension ({}).'.format(axis,self._DIMS[axis])
                
        return xyz

    def close(self):
        # need to figure out what to do with
        # OMETiff object
        pass

    def __enter__(self):
        return self

    def __del__(self):
        self.close()
        

    def __exit__(self, type_class, value, traceback):
        self.close()

    def __setitem__(self,keys):
        raise NotImplementedError('Cannot set values for {} class.'.format(self.__class__.__name__))
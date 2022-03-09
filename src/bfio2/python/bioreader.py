import libbfio2
import numpy

class BioReader:
    def __init__(self, file_name):
        self._file_name = file_name
        self._DIMS = {}
        if file_name.endswith('.ome.tif'):
            self._image_reader = libbfio2.OmeTiffLoader(file_name)
            self._meta_data = self.get_xml_metadata()
            self._image_height = self._image_reader.get_image_height()
            self._image_width = self._image_reader.get_image_width()
            self._DIMS['Y'] = self._image_height
            self._DIMS['X'] = self._image_width
        else:
            raise TypeError("Only OMETiff file format is supported")

    def get_xml_metadata(self):
        return self._image_reader.get_xml_metadata()
    
    def get_tile_data(self, row=None, col=None):
        if col != None:
            return self._image_reader.get_tile_data_2d_by_row_col(row, col)
        else:
            return self._image_reader.get_tile_data_2d_by_index(row)

    def get_image_height(self):
        return self._image_height

    def get_image_width(self):
        return self._image_width

    def get_tile_height(self):
        return self._image_reader.get_tile_height()

    def get_tile_width(self):
        return self._image_reader.get_tile_width()

    def get_row_tile_count(self):
        return self._image_reader.get_row_tile_count()

    def get_column_tile_count(self):
        return self._image_reader.get_column_tile_count()

    def channel_names(self) :
        pass

    def get_physical_size_x(self):
        return self._meta_data.get("PhysicalSizeX", "") + self._meta_data.get("PhysicalSizeXUnit", "")

    def get_physical_size_y(self):
        return self._meta_data.get("PhysicalSizeY", "") + self._meta_data.get("PhysicalSizeYUnit", "")

    def get_physical_size_z(self):
        pass
    
    def samples_per_pixel(self):
        #check how to do it for multi channel 
        return int(self._meta_data.get("SamplesPerPixel", "0"))

    def bits_per_sample(self):
        return int(self._meta_data.get("BitsPerSample", "0"))
    
    def bytes_per_pixel(self):
        return self.bits_per_sample()*self.samples_per_pixel/8

    def __getitem__(self, keys):
        slice_items = (self._parse_slice(keys))
        print(keys)
        print(slice_items)
        X = self._val_xyz(slice_items['X'], 'X')
        Y = self._val_xyz(slice_items['Y'], 'Y')

        row_min = Y[0]
        row_max = Y[1]
        if len(Y) == 3:
            row_step = Y[2]
        else:
            row_step = 1
        col_min = X[0]
        col_max = X[1]
        if len(X) == 3:
            col_step = X[2]
        else:
            col_step = 1
        
        if row_step != 1 or col_step != 1:
            return self._image_reader.get_virtual_tile_data_bounding_box_2d_strided(row_min, row_max, row_step, col_min, col_max, col_step)             
        else:
            return self._image_reader.get_virtual_tile_data_bounding_box_2d(row_min, row_max, col_min, col_max)

    def get_vtile_data_strided(self, row_min, row_max, row_stride, col_min, col_max, col_stride):
        return self._image_reader.get_virtual_tile_data_bounding_box_2d_strided(row_min, row_max, row_stride, col_min, col_max, col_stride)

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
                raise ValueError('Found {} indices, but at most 5 indices may be supplied.'.format(len(keys)))
            
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
                    stop = getattr(self,dim) if key.stop == None else key.stop
                    
                    # For CT dimensions, generate a list from slices
                    if dim in 'CT':
                        step = 1 if key.step is None else key.step
                        ind[dim] = list(range(start,stop,step))
                    
                    # For XYZ dimensions, get start and stop of slice, ignore step
                    else:
                        step = 1 if key.step is None else key.step
                        stop = min(self._DIMS[dim]-1, key.stop)
                        ind[dim] = [start,stop, step]
                        
                elif key==Ellipsis:
                    if dims.find(dim)+1 < len(keys):
                        raise ValueError('Ellipsis may only be used in the first or last index.')
                    
                elif isinstance(key,(int,tuple,list)) or numpy.issubdtype(key,numpy.integer):
                    # Only the last three dimensions can use int, tuple, or list indexing
                    if dim in 'CT':
                        if isinstance(key,int) or numpy.issubdtype(key,numpy.integer):
                            ind[dim] = [key]
                            if dim == 'Z':
                                ind[dim].append(ind[dim][0]+1)
                        else:
                            ind[dim] = key
                    elif isinstance(key,int) or numpy.issubdtype(key,numpy.integer):
                        ind[dim] = [int(key),int(key)+1]
                    else:
                        raise ValueError('The index in position {} must be a slice type.'.format(dims.find(dim)))
                else:
                    raise ValueError('Did not recognize indexing value of type: {}'.format(type(key)))
        return ind

    def __call__(self) :
        # Iterate through tiles of an image
        return self

    def __iter__(self):
        total_tiles = self.get_column_tile_count()*self.get_row_tile_count()    
        # start looping through batches
        for index in range(total_tiles):
            image = self.get_tile_data(index)
            yield image, index


    """ -------------------- """
    """ -Validation methods- """
    """ -------------------- """
    
    def _val_xyz(self, xyz, axis):
        assert axis in 'XY'
        print(xyz)
        if xyz == None:
            xyz = [0,self._DIMS[axis]-1]
        else:
            if axis=='Z' and isinstance(xyz,int):
                xyz = [xyz,xyz+1]
            assert len(xyz) == 2 or len(xyz) == 3, \
                '{} must be a list or tuple of length 2 or 3.'.format(axis)
            assert xyz[0] >= 0, \
                '{}[0] must be greater than or equal to 0.'.format(axis)
            assert xyz[1] <= self._DIMS[axis], \
                '{}[1] cannot be greater than the maximum of the dimension ({}).'.format(axis,self._DIMS[axis])
                
        return xyz

import abc, threading, logging
import numpy
from concurrent.futures import ThreadPoolExecutor
from queue import Queue
import multiprocessing, typing
from pathlib import Path
import bfio2.python.bfio2 as bfio2

class BioBase(object,metaclass=abc.ABCMeta) :
    """ Abstract class for reading/writing OME tiled tiff images
    
    Attributes:
        dtype: Gets/sets the pixel type (e.g. uint8)
        channel_names: Gets/sets the names of each channel
        samples_per_pixel: Numbers of numbers per pixel location
        bytes_per_pixel: Number of bytes per pixel
        x: Get/set number of pixels in the x-dimension (width)
        y: Get/set number of pixels in the y-dimension (height)
        z: Get/set number of pixels in the z-dimension (depth)
        c: Get/set number of channels in the image
        t: Get/set number of timepoints in the image
        physical_size_x: Get/set the physical size of the x-dimension
        physical_size_y: Get/set the physical size of the y-dimension
        physical_size_z: Get/set the physical size of the z-dimension
        metadata: OmeXml object for the image
        cnames: Same as channel_names
        spp: Same as samples_per_pixel
        bpp: Same as bytes_per_pixel
        X: Same as x attribute
        Y: Same as y attribute
        Z: Same as z attribute
        C: Same as c attribute
        T: same as t attribute
        ps_x: Same as physical_size_x
        ps_y: Same as physical_size_y
        ps_z: Same as physical_size_z

    """
    
    # protected attribute to hold metadata
    _metadata = None
    _backend = None
    
    
    def __init__(self,
                 file_path: typing.Union[str,Path]):
        """__init__ Initialize BioBase object

         Args:
            file_path: Path to output file
        """
        
        # Internally, keep the file_path as a Path object
        if isinstance(file_path,str):
            file_path = Path(file_path)
        self._file_path = file_path

        
    
    def __setitem__(self,keys: typing.Union[list,tuple],values: numpy.ndarray):
        raise NotImplementedError('Cannot set values for {} class.'.format(self.__class__.__name__))
    
    def __getitem__(self,keys):
        raise NotImplementedError('Cannot get values for {} class.'.format(self.__class__.__name__))
    
    def _parse_slice(self,keys):
        
        # Dimension ordering and index initialization
        dims = 'YXZCT'
        ind = {d:None for d in dims}
        
        # If an empty slice, load the whole image
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
                        ind[dim] = [start,stop]
                        
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
                
    """ ------------------------------ """
    """ -Get/Set Dimension Properties- """
    """ ------------------------------ """
    @property
    def channel_names(self) -> typing.List[str]:
        """Get the channel names for the image"""
        image = self._metadata.image()
        return [image.Pixels.Channel(i).Name for i in range(0, self.C)]
        

            
    def __physical_size(self,dimension,psize,units):
        if psize != None and units != None:
            assert not self._read_only, self._READ_ONLY_MESSAGE.format("physical_size_{}".format(dimension.lower()))
            setattr(self._metadata.image(0).Pixels,'PhysicalSize{}'.format(dimension.upper()),psize)
            setattr(self._metadata.image(0).Pixels,'PhysicalSize{}Unit'.format(dimension.upper()),units)

    @property
    def physical_size_x(self) -> typing.Tuple[float,str]:
        """Physical size of pixels in x-dimension

        Returns:
            Units per pixel, Units (i.e. "cm" or "mm")
        """
        return (self._metadata.image(0).Pixels.PhysicalSizeX, self._metadata.image(0).Pixels.PhysicalSizeXUnit)
        
    @property
    def physical_size_y(self) -> typing.Tuple[float,str]:
        """Physical size of pixels in y-dimension

        Returns:
            Units per pixel, Units (i.e. "cm" or "mm")
        """
        return (self._metadata.image(0).Pixels.PhysicalSizeY, self._metadata.image(0).Pixels.PhysicalSizeYUnit)

        
    @property
    def physical_size_z(self) -> typing.Tuple[float,str]:
        """Physical size of pixels in z-dimension

        Returns:
            Units per pixel, Units (i.e. "cm" or "mm")
        """
        return (self._metadata.image(0).Pixels.PhysicalSizeZ, self._metadata.image(0).Pixels.PhysicalSizeZUnit)

    """ -------------------- """
    """ -Validation methods- """
    """ -------------------- """
    
    def _val_xyz(self, xyz: typing.Union[int,list,tuple], axis: str) -> typing.List[int]:
        """_val_xyz Utility function for validating image dimensions

        Args:
            xyz: Pixel value of x, y, or z dimension.
                If None, returns the maximum range of the dimension
            axis: Must be 'x', 'y', or 'z'

        Returns:
            list of ints indicating the first and last index in the dimension
        """
        assert axis in 'XYZ'
        
        if xyz == None:
            xyz = [0,self._DIMS[axis]]
        else:
            if axis=='Z' and isinstance(xyz,int):
                xyz = [xyz,xyz+1]
            assert len(xyz) == 2, \
                '{} must be a list or tuple of length 2.'.format(axis)
            assert xyz[0] >= 0, \
                '{}[0] must be greater than or equal to 0.'.format(axis)
            assert xyz[1] <= self._DIMS[axis], \
                '{}[1] cannot be greater than the maximum of the dimension ({}).'.format(axis,self._DIMS[axis])
                
        return xyz

    def _val_ct(self, ct: typing.Union[int,list], axis: str) -> typing.List[int]:
        """_val_ct Utility function for validating image dimensions

        Args:
            ct: List of ints indicating the channels or timepoints to load
                If None, returns a list of ints
            axis: Must be 'c', 't'

        Returns:
            list of ints indicating the first and last index in the dimension
        """

        assert axis in 'CT'
        
        if ct == None:
            # number of timepoints
            ct = list(range(0, self._DIMS[axis]))
        else:
            assert numpy.any(numpy.greater(self._DIMS[axis], ct)), \
            'At least one of the {}-indices was larger than largest index ({}).'.format(axis, self._DIMS[axis] - 1)
            assert numpy.any(numpy.less_equal(0, ct)), \
            'At least one of the {}-indices was less than 0.'.format(axis)
            assert len(ct) != 0, \
            'At least one {}-index must be selected.'.format(axis)
            
        return ct

    """ ------------------- """
    """ -Pixel information- """
    """ ------------------- """
    
        
    @property
    def samples_per_pixel(self) -> int:
        """Number of samples per pixel"""
        return self._metadata.image().Pixels.Channel().SamplesPerPixel
    
    @property
    def bytes_per_pixel(self) -> int:
        """Number of bytes per pixel"""
        return self._BPP[self._metadata.image().Pixels.get_PixelType()]
    
 
    """ -------------------------- """
    """ -Other Methods/Properties- """
    """ -------------------------- """
    @property
    def metadata(self) -> bfio2.OmeXml.OMEXML:
        """Get the metadata for the image

        Returns:
            Dictionary containing metadata for the image
        """        
        return self._metadata
    def close(self):
        """Close the image"""
        if self._backend is not None:
            self._backend.close()
        
    def __enter__(self):
        """Handle entrance to a context manager
        
        This code is called when a `with` statement is used. This allows a
        BioBase object to be used like this:
        
        with bfio.BioReader('Path/To/File.ome.tif') as reader:
            ...
            
        with bfio.BioWriter('Path/To/File.ome.tif') as writer:
            ...
        """
        return self
    
    def __del__(self):
        """Handle file deletion

        This code runs when an object is deleted..
        """
        self.close()
        
    
    def __exit__(self, type_class, value, traceback):
        """Handle exit from the context manager

        This code runs when exiting a `with` statement.
        """
        self.close()
        
class AbstractBackend(object,metaclass=abc.ABCMeta):
    
    @abc.abstractmethod
    def __init__(self,frontend):
        self.frontend = frontend
        self._lock = threading.Lock()

    def _image_io(self,X,Y,Z,C,T,image):
            
        # Define tile bounds
        ts = self.frontend._TILE_SIZE
        X_tile_shape = X[1] - X[0]
        Y_tile_shape = Y[1] - Y[0]
        Z_tile_shape = Z[1] - Z[0]
        
        # Set the output for asynchronous reading
        self._image = image

        # Set up the tile indices
        self._tile_indices = []
        for t in range(len(T)):
            for c in range(len(C)):
                for z in range(Z_tile_shape):
                    for y in range(0,Y_tile_shape,ts):
                        for x in range(0,X_tile_shape,ts):
                            self._tile_indices.append(((x,X[0]+x),
                                                       (y,Y[0]+y),
                                                       (z,Z[0]+z),
                                                       (c,C[c]),
                                                       (t,T[t])))
        
        self.logger.debug('_image_io(): _tile_indices = {}'.format(self._tile_indices))
        
    @abc.abstractmethod
    def close(self):
        pass

class AbstractReader(AbstractBackend):
    
    @abc.abstractmethod
    def __init__(self,frontend):
        super().__init__(frontend)
    
    @abc.abstractmethod
    def read_metadata(self):
        pass
    
    def read_image(self,*args):
        with self._lock:
            self._image_io(*args)
            self._read_image(*args)

    @abc.abstractmethod
    def _read_image(self,X,Y,Z,C,T,output):
        pass

class AbstractWriter(AbstractBackend):
    
    _writer = None
    
    @abc.abstractmethod
    def __init__(self,frontend):
        super().__init__(frontend)
        self.initialized = False
    
    def write_image(self,*args):
        with self._lock:
            if not self.initialized:
                self._init_writer()
                self.frontend.__read_only = True
                self.initialized = True
            
            self._image_io(*args)
            self._write_image(*args)
        
    @abc.abstractmethod
    def _init_writer(self):
        pass
    
    @abc.abstractmethod
    def _write_image(*args):
        pass
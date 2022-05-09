from mimetypes import init
from bfio2 import BioReader as NewBioReader
import time
from bfio import BioReader
import os

def bfio_loop_all_tiles(file_name, num_threads=1):
    test_name = "Bfio_all_tiles_access"
    start_time = time.time()
    tile_grid_size = 1
    tile_size = tile_grid_size * 1024

    br = BioReader(file_name, max_workers=num_threads)
    init_time = time.time()
    for t in range(br.T):
        for c in range(br.C):
            for z in range(br.Z):
                for y in range(0, br.Y, tile_size):
                    y_max = min([br.Y,y+tile_size])
                    for x in range(0, br.X, tile_size):
                        x_max = min([br.X,x+tile_size])
                        data = br[y:y_max,x:x_max,z:z+1,c,t]
                        sum = data.sum()



    end_time = time.time()

    print(f"Test Name : {test_name} Init time : {init_time - start_time} IO time : {end_time-init_time} ")
    return init_time - start_time, end_time-init_time

def bfio2_loop_all_tiles_unplanned(file_name, num_threads=1):
    test_name = "Bfio2_all_tiles_unplanned_access"
    start_time = time.time()
    tile_grid_size = 1
    tile_size = tile_grid_size * 1024

    br = NewBioReader(file_name, num_threads)
    init_time = time.time()
    for t in range(br.T):
        for c in range(br.C):
            for z in range(br.Z):
                for y in range(0, br.Y, tile_size):
                    y_max = min([br.Y,y+tile_size])
                    for x in range(0, br.X, tile_size):
                        x_max = min([br.X,x+tile_size])
                        data = br[y:y_max,x:x_max,z:z+1,c,t]
                        sum = data.sum()


    end_time = time.time()
    print(f"Test Name : {test_name} Init time : {init_time - start_time} IO time : {end_time-init_time} ")
    return init_time - start_time, end_time-init_time


def bfio2_loop_all_tiles_planned(file_name, num_threads=1):
    test_name = "Bfio2_all_tiles_planned_access"
    start_time = time.time()
    br = NewBioReader(file_name, num_threads)
    init_time = time.time()
    for ind, tile in br(tile_size=[1024,1024],tile_stride=[1024,1024]):
        sum = tile.sum()

    end_time = time.time()
    print(f"Test Name : {test_name} Init time : {init_time - start_time} IO time : {end_time-init_time} ")
    return init_time - start_time, end_time-init_time


def bfio2_access_full_image(file_name, num_threads=1):
    test_name = "Bfio2_access_full_image"
    start_time = time.time()
    br = NewBioReader(file_name, num_threads)
    init_time = time.time()
    data = br[:,:,:,:,:]
    end_time = time.time()
    print(data.sum())
    print(f"Test Name : {test_name} Init time : {init_time - start_time} IO time : {end_time-init_time} ")
    return init_time - start_time, end_time-init_time

def bfio_access_full_image(file_name, num_threads=1):
    test_name = "Bfio_access_full_image"
    start_time = time.time()
    br = BioReader(file_name, num_threads)
    init_time = time.time()
    data = br[:,:,:,:,:]
    end_time = time.time()
    print(data.sum())
    print(f"Test Name : {test_name} Init time : {init_time - start_time} IO time : {end_time-init_time} ")
    return init_time - start_time, end_time-init_time

if __name__ == "__main__":
    image_list = [
        '/home/ec2-user/data/s_1_t_10_c_3_z_1_tiff_tiles.ome.tif',
        '/home/ec2-user/data/3d-cell-viewer_ome_tiff_tiles.ome.tif',
       # '/home/ec2-user/data/r001_c001_z000.ome.tif',
       # '/home/ec2-user/data/s_1_t_10_c_3_z_1_tiff_tiles_uncompressed.ome.tif',
        '/home/ec2-user/data/3d-cell-viewer_ome_tiff_tiles_uncompressed.ome.tif',
        '/home/ec2-user/data/r001_c001_z000_uncompressed.ome.tif',
    ]

    thread_cfg = [4,]

    test_dict = {}
    test_dict["access full image"] = {}
    test_dict["loop through image"] = {}

    test_dict["access full image"]["bfio"] = bfio_access_full_image
    test_dict["access full image"]["bfio2"] = bfio2_access_full_image
    test_dict["loop through image"]["bfio"] = bfio_loop_all_tiles
    test_dict["loop through image"]["bfio2 planned"] = bfio2_loop_all_tiles_planned
    test_dict["loop through image"]["bfio2 unplanned"] = bfio2_loop_all_tiles_unplanned    

    num_trials = 5
    with open("perf_data_"+str(int(time.time()))+".txt", "w") as fp:
        for image in image_list:
            for thread in thread_cfg:
                for test_type in test_dict:
                    test_data = {}
                    for test_name in test_dict[test_type]:
                        function_name = test_dict[test_type][test_name]
                        sum_init = 0
                        sum_io = 0
                        for i in range(num_trials):
                            init_time, io_time = function_name(image, thread)
                            sum_init += init_time
                            sum_io += io_time
                            time.sleep(3)
                        test_data[test_name] = (sum_init/num_trials, sum_io/num_trials)
                
                    tmp, image_name = os.path.split(image)
                    init_result = f"{test_type}, {image_name},Init Time , {thread}, "
                    io_result = f"{test_type}, {image_name}, IO Time, {thread}, "

                    for test_name in test_data:
                        result = test_data[test_name]
                        init_result = init_result + ", " + test_name + ", " + str(result[0])
                        io_result = io_result + ", " + test_name + ", " + str(result[1])

                    init_result = init_result+"\n"
                    io_result = io_result+"\n"
                    print(init_result)
                    print(io_result)
                    fp.write(init_result)
                    fp.write(io_result)
                    fp.flush()                



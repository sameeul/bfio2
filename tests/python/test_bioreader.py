from bfio2 import BioReader

#loop through tiles
print("Checking Pythonic loop through tiles")
with BioReader("r01_x10_y05_z08.ome.tif") as br:
    for index, tile_data in br:
        print(f"Tile {index} - All pixel sum {tile_data.sum()}")


# Access using index - result should match access by (row,col)
print("Checking access by slice index")
img = BioReader("r01_x10_y05_z08.ome.tif")
img_by_slice_index = img[0:1024,0:1024]
img_by_row_col = img.get_tile_data(0,0)

if (img_by_slice_index-img_by_row_col).sum() != 0:
    print("Results do not match")
else:
    print("Resuls match!")


#print various metadata
print("Printing metadata")
print(f"Image height: {img.get_image_height()}")
print(f"Image width: {img.get_image_width()}")
print(f"Tile height: {img.get_tile_height()}")
print(f"Tile width: {img.get_tile_width()}")
print(f"Number of tile rows: {img.get_row_tile_count()}")
print(f"Number of tile column: {img.get_column_tile_count()}")
print(f"Physical size x: {img.get_physical_size_x()}")
print(f"Physical size y: {img.get_physical_size_y()}")

# test various row stride
print("Checking virtual tile access with stride")
whole_img_by_silce_index = img[0:1080,0:1080]
match = True
for i in range(1, 10):
    for j in range(1, 10):
        row_min = 0
        row_max = 1079
        row_stride = i

        col_min = 0
        col_max = 1079
        col_stride = j
        test2 = img[row_min:row_max+1:row_stride,col_min:col_max+1:col_stride]

        test3 = whole_img_by_silce_index[0:row_max+1:row_stride, 0:col_max+1:col_stride]
        if (test3-test2).sum() != 0:
            print("Alas!")
            match = False
            break
if match:
    print("Virtual tile with non-unit stride works!")
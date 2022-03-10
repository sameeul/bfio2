from bfio2 import BioReader

img = BioReader("r01_x10_y05_z08.ome.tif")
test1 = img[0:1079,0:1079]

for i in range(1, 10):
    for j in range(1, 10):
        row_min = 0
        row_max = 1079
        row_stride = i

        col_min = 0
        col_max = 1079
        col_stride = j
        test2 = img[row_min:row_max:row_stride,col_min:col_max:col_stride]
        #print(test2.shape)

        test3 = test1[0:1080:row_stride, 0:1080:col_stride]
        if (test3-test2).sum() != 0:
            print("Alas!")
            break
# print(test3)
# print(test2)
from cgi import test
from bioreader import BioReader
import matplotlib.pyplot as plt
import numpy as np

img = BioReader("r01_x10_y05_z08.ome.tif")

print(img.get_xml_metadata())
meta_data = dict(img.get_xml_metadata())
print(meta_data)
test = img.get_tile_data(1)
print(test)
for tile_data, index in img:
    print(tile_data.sum(), index)

test = img[0:1079,0:1079]
print(np.shape(test))
print(test)
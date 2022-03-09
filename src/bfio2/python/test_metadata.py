from cgi import test
from bfio2.python.bioreader import OmeTiffReader
import matplotlib.pyplot as plt
img = OmeTiffReader("/home/samee/axle_dev/fl_test/r01_x10_y05_z08.ome.tif")
print(img.get_xml_metadata())
meta_data = dict(img.get_xml_metadata())
print(meta_data)
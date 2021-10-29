import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from pyhdf.SD import SD, SDC
import sys

data_fname = '/mnt/lfs/d30/vegeo/SAT/DATA/MODIS/MODIS-MCD43D51/tiles/EPS_based_dates/2015/02/25/MCD43D51.A2015056.006.2016166224249.hdf'
qflag_fname = '/mnt/lfs/d30/vegeo/SAT/DATA/MODIS/MODIS-MCD43D31/tiles/EPS_based_dates/2015/02/25/MCD43D31.A2015056.006.2016166224247.hdf'

data_f = SD(data_fname, SDC.READ)
qflag_f = SD(qflag_fname, SDC.READ)


for i in data_f.datasets().items():
    print(i)
print(data_f.attributes()["ArchiveMetadata.0"])
print(data_f.attributes()["StructMetadata.0"])
sys.exit()
for i in data_f.attributes().items():
    print(i)
for i in qflag_f.datasets().items():
    print(i)

if 0:
    data = data_f.select('BRDF_Albedo_BSA_Shortwave')[5150:5650,22000:22500] # Alps
    qflag = qflag_f.select('BRDF_Quality')[5150:5650,22000:22500] # Alps
else:
    stride = 10
    print(data_f.select('BRDF_Albedo_BSA_Shortwave').dimensions())
    sys.exit()
    data = data_f.select('BRDF_Albedo_BSA_Shortwave')[::stride,::stride] # World
    qflag = qflag_f.select('BRDF_Quality')[::stride,::stride] # World

good_quality_mask = qflag==0

data[data==32767] = -1
#data[~good_quality_mask] = -1

#ims = plt.imshow(data, cmap='jet', vmin=0, vmax=600)
ims = plt.imshow(data, cmap='jet', vmin=0)

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "5%", pad="3%")
plt.colorbar(ims, cax=cax)
plt.show()



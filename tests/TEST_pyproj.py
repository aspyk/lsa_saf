from pyproj import CRS, Transformer
import numpy as np


crs0 = CRS.from_proj4("+proj=eqc +units=km")
crs1 = CRS.from_proj4("+proj=sinu +units=km")
crs2 = CRS.from_proj4("+proj=geos +h=35774290 +units=km")
crs3 = CRS.from_proj4("+proj=lonlat +units=km")

if 0:
    t0 = Transformer.from_crs(crs0, crs3)

    x = np.array([-20000, 0, 20000])
    y = np.array([-10000, 0, 10000])
    
    res = t0.transform(x,y)

if 1:    
    t1 = Transformer.from_crs(crs3, crs1)
    
    x1 = np.array([-180, 0, 180])
    y1 = np.array([90, 0, -90])
    
    x1, y1 = np.meshgrid(x1,y1)
    
    res = t1.transform(x1,y1)
    
    print(res[0])
    print(res[1])

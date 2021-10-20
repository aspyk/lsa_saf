import h5py
import numpy as np


pix_size = 1/120
lat = np.linspace(90-pix_size/2, -90+pix_size/2, 21600)
lon = np.linspace(-180+pix_size/2, 180-pix_size/2, 2*21600)
lon,lat = np.meshgrid(lon,lat)

with h5py.File('lat_modis.h5','w') as h5f:
    #dataset = h5f.create_dataset('lat', chunks=True, compression='gzip', fletcher32=True, shape=lat.shape, dtype=int)
    dataset = h5f.create_dataset('lat', chunks=True, compression='gzip', fletcher32=True, shape=lat.shape)
    dataset[:,:] = lat

with h5py.File('lon_modis.h5','w') as h5f:
    #dataset = h5f.create_dataset('lat', chunks=True, compression='gzip', fletcher32=True, shape=lat.shape, dtype=int)
    dataset = h5f.create_dataset('lon', chunks=True, compression='gzip', fletcher32=True, shape=lat.shape)
    dataset[:,:] = lon



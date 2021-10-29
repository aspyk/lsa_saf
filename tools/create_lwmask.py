import h5py

with h5py.File('/mnt/lfs/d30/vegeo/SAT/DATA/EPS/Backprocessed/ETAL/globe/2016/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_201607150000', 'r') as h5f:
    lw = h5f['Q-Flag'][:]
lw = lw & 0b11

with h5py.File('etal_lwmask.h5', 'w') as h5f:
    dataset = h5f.create_dataset('lwmask', chunks=True, compression='gzip', shape=lw.shape, dtype=int)
    dataset[:,:] = lw.astype('i4')

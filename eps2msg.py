"""Project Global EPS results onto MSG disk

This script takes a global EPS-ETAL dataset and projects it onto MSG-disk pixel locations.
"""


import h5py
import numpy as np
from datetime import datetime as dtime
import matplotlib.pyplot as plt
import copy
from time import process_time

from tools import SimpleTimer


# probably need a different interpolation method for dealing with NaN's, maybe: https://stackoverflow.com/questions/34408293/2-d-interpolation-ignoring-nan-values
# alternative approach: for each MDAL pixel form ETAL estimate by taking mean in certain radius? <---- THIS IS IMPLEMENTED AT THE MOMENT

ti = SimpleTimer()

albVar = 'AL-BB-BH'
date_ = '2020-DEC-05'
dateObj = dtime.strptime(date_, '%Y-%b-%d')
dateForFile = dtime.strftime(dateObj, '%Y%m%d')

# ETAL
with h5py.File(f'/cnrm/vegeo/juncud/NO_SAVE/ETAL/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_{dateForFile}0000','r+') as F:
    E = F[albVar][:]/ F[albVar].attrs['SCALING_FACTOR']
    
# MSG
# Latitude
with h5py.File('/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LAT_MSG-Disk_4bytesPrecision', 'r') as F:
    msgLat = F['LAT'] / F['LAT'].attrs['SCALING_FACTOR']
    msgLat[msgLat == 91] = np.nan
# Longitude
with h5py.File('/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LON_MSG-Disk_4bytesPrecision', 'r') as F:
    msgLon = F['LON'] / F['LON'].attrs['SCALING_FACTOR']
    msgLon[msgLon == 91] = np.nan
# actual MSG data
with h5py.File(f'/cnrm/vegeo/juncud/NO_SAVE/aod-test-2/v1/HDF5_LSASAF_MSG_ALBEDO_MSG-Disk_{dateForFile}0000', 'r') as F:
    mdalv1 = F[albVar] / F[albVar].attrs['SCALING_FACTOR']
    mdalv1[mdalv1<0.] = np.nan

# only use every nth pixel
stride = 4
stride = 10

# lon and lat limits
# (already discarding part of ETAL data that is for sure not in the MSG disk)
iLatStartEps = 850
iLatEndEps   = 17150 # 6501
iLonStartEps = 9000 # 22000
iLonEndEps   = 27000 # 25501
#iLonEndEps   = 15000 # 25501

# averaging parameters
radius_ = 20e3 # 10e3
latBoxSize = radius_ / 110574 # approximative length for 1 deg lat.: 110.574km

# albedo values for given lon/lat range
eps = E[iLatStartEps:iLatEndEps:stride,iLonStartEps:iLonEndEps:stride]
del E
del F

eps[eps<0.] = np.nan
nEpsLat = np.shape(eps)[0]
nEpsLon = np.shape(eps)[1]

# split the scene into 4 quadrants, for performance reasons
eps_q1 = eps[0:nEpsLat//2, 0:nEpsLon//2]
eps_q2 = eps[nEpsLat//2:nEpsLat, 0:nEpsLon//2]
eps_q3 = eps[nEpsLat//2:nEpsLat, nEpsLon//2:nEpsLon]
eps_q4 = eps[0:nEpsLat//2, nEpsLon//2:nEpsLon]
del eps

# reshape into 1d arrays
eps_q1_1d = np.reshape(eps_q1, (np.size(eps_q1), 1)).astype('f4')
eps_q2_1d = np.reshape(eps_q2, (np.size(eps_q2), 1)).astype('f4')
eps_q3_1d = np.reshape(eps_q3, (np.size(eps_q3), 1)).astype('f4')
eps_q4_1d = np.reshape(eps_q4, (np.size(eps_q4), 1)).astype('f4')
del eps_q1
del eps_q2
del eps_q3
del eps_q4
#nEps = np.size(eps_1d)

# create EPS coordinate arrays
epsLat1d = np.linspace(90., -90., 18001)
epsLon1d = np.linspace(-180., 179.99, 36000)
epsLat1d = epsLat1d[iLatStartEps:iLatEndEps:stride]
epsLon1d = epsLon1d[iLonStartEps:iLonEndEps:stride]
epsLat, epsLon = np.meshgrid(epsLat1d, epsLon1d)
X = copy.deepcopy(epsLon.T)

# to get the actual EPS longitude
epsLon /= np.cos(np.deg2rad(epsLat))

epsLon = epsLon.T
epsLat = epsLat.T

# lon and lat meshgrids for each quadrant
epsLon_q1 = epsLon[0:nEpsLat//2, 0:nEpsLon//2]
epsLon_q2 = epsLon[nEpsLat//2:nEpsLat, 0:nEpsLon//2]
epsLon_q3 = epsLon[nEpsLat//2:nEpsLat, nEpsLon//2:nEpsLon]
epsLon_q4 = epsLon[0:nEpsLat//2, nEpsLon//2:nEpsLon]
epsLat_q1 = epsLat[0:nEpsLat//2, 0:nEpsLon//2]
epsLat_q2 = epsLat[nEpsLat//2:nEpsLat, 0:nEpsLon//2]
epsLat_q3 = epsLat[nEpsLat//2:nEpsLat, nEpsLon//2:nEpsLon]
epsLat_q4 = epsLat[0:nEpsLat//2, nEpsLon//2:nEpsLon]
del(epsLon)
del(epsLat)

# lon and lat 1d arrays for each quadrant
epsLon_q1_1d = np.reshape(epsLon_q1, (np.size(epsLon_q1), 1))
epsLon_q2_1d = np.reshape(epsLon_q2, (np.size(epsLon_q2), 1))
epsLon_q3_1d = np.reshape(epsLon_q3, (np.size(epsLon_q3), 1))
epsLon_q4_1d = np.reshape(epsLon_q4, (np.size(epsLon_q4), 1))
epsLat_q1_1d = np.reshape(epsLat_q1, (np.size(epsLat_q1), 1))
epsLat_q2_1d = np.reshape(epsLat_q2, (np.size(epsLat_q2), 1))
epsLat_q3_1d = np.reshape(epsLat_q3, (np.size(epsLat_q3), 1))
epsLat_q4_1d = np.reshape(epsLat_q4, (np.size(epsLat_q4), 1))
del(epsLon_q1)
del(epsLon_q2)
del(epsLon_q3)
del(epsLon_q4)
del(epsLat_q1)
del(epsLat_q2)
del(epsLat_q3)
del(epsLat_q4)

# mask for only those pixels that are relevant
# (within MSG-disk, no NAN's)
epsMask_q1 = np.logical_and(epsLon_q1_1d > -90., epsLon_q1_1d < 90.) & \
          np.logical_and(epsLat_q1_1d > np.nanmin(msgLat), epsLat_q1_1d < np.nanmax(msgLat)) & \
          ~np.isnan(eps_q1_1d)
epsMask_q2 = np.logical_and(epsLon_q2_1d > -90., epsLon_q2_1d < 90.) & \
          np.logical_and(epsLat_q2_1d > np.nanmin(msgLat), epsLat_q2_1d < np.nanmax(msgLat)) & \
          ~np.isnan(eps_q2_1d)
epsMask_q3 = np.logical_and(epsLon_q3_1d > -90., epsLon_q3_1d < 90.) & \
          np.logical_and(epsLat_q3_1d > np.nanmin(msgLat), epsLat_q3_1d < np.nanmax(msgLat)) & \
          ~np.isnan(eps_q3_1d)
epsMask_q4 = np.logical_and(epsLon_q4_1d > -90., epsLon_q4_1d < 90.) & \
          np.logical_and(epsLat_q4_1d > np.nanmin(msgLat), epsLat_q4_1d < np.nanmax(msgLat)) & \
          ~np.isnan(eps_q4_1d)

eps_q1_1d = eps_q1_1d[epsMask_q1]
epsLon_q1_1d = epsLon_q1_1d[epsMask_q1]
epsLat_q1_1d = epsLat_q1_1d[epsMask_q1]
eps_q2_1d = eps_q2_1d[epsMask_q2]
epsLon_q2_1d = epsLon_q2_1d[epsMask_q2]
epsLat_q2_1d = epsLat_q2_1d[epsMask_q2]
eps_q3_1d = eps_q3_1d[epsMask_q3]
epsLon_q3_1d = epsLon_q3_1d[epsMask_q3]
epsLat_q3_1d = epsLat_q3_1d[epsMask_q3]
eps_q4_1d = eps_q4_1d[epsMask_q4]
epsLon_q4_1d = epsLon_q4_1d[epsMask_q4]
epsLat_q4_1d = epsLat_q4_1d[epsMask_q4]

nLon = np.shape(msgLon)[1]
nLat = np.shape(msgLat)[0]
nAll = nLon * nLat
loopIndex = 1
milestone = 1
indexArray_q1 = np.indices(epsLon_q1_1d.shape)[0]
indexArray_q2 = np.indices(epsLon_q2_1d.shape)[0]
indexArray_q3 = np.indices(epsLon_q3_1d.shape)[0]
indexArray_q4 = np.indices(epsLon_q4_1d.shape)[0]
#minDist_1d = np.zeros((np.size(epsLon_1d), ), dtype=bool)
eps_geos = np.full((nLat, nLon), np.nan, dtype='f4')
actualILat = 0


ti('preproc')


#TIME RESULTS
#Find quadrant: 4.773999989993172e-06 s
#Create mask: 0.009312423999972452 s
#Average & store: 0.0004993200000171782 s
#inner loop: 0.00991558300000861 s

# loop over MSG pixels. 
# At each pixel find neighbouring EPS pixels and get mean value
for iLat in range(nLat):
    for iLon in range(nLon):
        #t_inner = process_time() 

        # continue, if outside of MSG disk
        if np.isnan(msgLon[iLat, iLon]) or np.isnan(mdalv1[iLat, iLon]):
            #print('out of disk, continue')
            loopIndex += 1
            continue

        #t0 = process_time() 
        # check in which quadrant we are in (using quadrants makes it faster creating the mask below)
        if msgLon[iLat, iLon] <= 0. and msgLat[iLat, iLon] >= 0.:
            eps_Q = eps_q1_1d
            epsLon_Q = epsLon_q1_1d
            epsLat_Q = epsLat_q1_1d
            #ind = indexArray_q1
        elif msgLon[iLat, iLon] <= 0. and msgLat[iLat, iLon] < 0.:
            eps_Q = eps_q2_1d
            epsLon_Q = epsLon_q2_1d
            epsLat_Q = epsLat_q2_1d
            #ind = indexArray_q2
        elif msgLon[iLat, iLon] > 0. and msgLat[iLat, iLon] < 0.:
            eps_Q = eps_q3_1d
            epsLon_Q = epsLon_q3_1d
            epsLat_Q = epsLat_q3_1d
            #ind = indexArray_q3
        elif msgLon[iLat, iLon] > 0. and msgLat[iLat, iLon] >= 0.:
            eps_Q = eps_q4_1d
            epsLon_Q = epsLon_q4_1d
            epsLat_Q = epsLat_q4_1d
            #ind = indexArray_q4
        #print(f'Find quadrant: {process_time() - t0} s') 

        # approximative length for 1 deg lon.: 111.320*cos(lat.) km
        lonBoxSize = radius_ / 111320 * np.cos(np.deg2rad(msgLat[iLat,iLon]))
        
        #t0 = process_time()

        # mask for pixels within lonbox & latbox around MSG pixel
        neighborMask = ((epsLon_Q > msgLon[iLat,iLon]-lonBoxSize) & \
                       (epsLon_Q < msgLon[iLat,iLon]+lonBoxSize)) & \
                       ((epsLat_Q > msgLat[iLat,iLon]-latBoxSize) & \
                       (epsLat_Q < msgLat[iLat,iLon]+latBoxSize))
        
        #t0 = process_time() 
        
        # the re-projected value is the mean of EPS values near (within given radius) a MSG coordinate
        eps_geos[iLat, iLon] = np.nanmean(eps_Q[neighborMask])
        #eps_geos[iLat, iLon] = np.nanmean(eps_Q[ind])
        
        #print(f'Average & store: {process_time() - t0} s') 
        loopProgress = loopIndex / nAll * 100.
        #print(f'{int(loopProgress)} %')
        if loopProgress >= milestone:
            print(f'{int(loopProgress)} %')
            print(f'latitude index: {iLat}')
            print(f'Lon / Lat: {msgLon[iLat,iLon]} / {msgLat[iLat,iLon]}')
            print(f'last no. of pixels in box: {np.sum(neighborMask)}')
            milestone += 1
            ti(f'loopProgress{round(loopProgress)}%,iLat{iLat}')
        loopIndex += 1
        #print(f'inner loop: {process_time() - t_inner} s') 

ti('endofproc')

ti.show()

with h5py.File(f'eps_geos_{dateForFile}0000', 'w') as f:
    dataset = f.create_dataset('AL-BB-BH', chunks=True, compression='gzip', fletcher32=True, shape=(nLat, nLon), dtype=float)
    dataset[:,:] = eps_geos


# some plotting

#clrmap= plt.get_cmap('viridis')
#clrmap.set_over('black')
#clrmap.set_under('black')

#fig, (ax1,ax2,ax3) = plt.subplots(1,3)
#p1 = ax1.pcolormesh(X,epsLat,epsLon, vmin = -120, vmax=120, cmap=clrmap)
#p2 = ax2.pcolormesh(X,epsLat,EPS)
#p3 = ax3.pcolormesh(epsLon,epsLat,epsLon-X, vmin=-50, vmax=50)
#fig.colorbar(p1, ax=ax1)
#fig.colorbar(p2, ax=ax2)
#fig.colorbar(p3, ax=ax3)

#fig, ax = plt.subplots()
#ax.pcolormesh(msgLon, msgLat, eps_geos, cmap=clrmap)
#
#plt.show()





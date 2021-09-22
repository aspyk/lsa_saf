import h5py
import numpy as np
#import vtk
#from vtk.util import numpy_support as ns

import os,sys
from tools import SimpleTimer

import matplotlib.pyplot as plt

import datetime as dt

import traceback
import itertools

from vtk_tools import vtk_interpolation

import pandas as pd
import pathlib



"""
TODO
----
- make a class SatelliteData with .lon, .lat and .var attribute. And after a MSG() and EPS() subclass.
"""


def load_msg_lonlat(stride=10, **param):
    msg_lat_file = '/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LAT_MSG-Disk_4bytesPrecision'
    msg_lon_file = '/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LON_MSG-Disk_4bytesPrecision'
    #'hdf5_lsasaf_msg_lat_msg-disk_4bytesprecision'
    #'hdf5_lsasaf_msg_lon_msg-disk_4bytesprecision'
    #'hdf5_lsasaf_usgs-igbp_lwmask_msg-disk'
   
    if 1: # full disc
        lat1, lat2 = None, None
        lon1, lon2 = None, None
    else: # zoom
        lat1 = 900
        lat2 = 980
        lon1 = 2320
        lon2 = 2400

    with h5py.File(msg_lat_file, 'r') as flat:
        lat = flat['LAT'][lat1:lat2:stride,lon1:lon2:stride]/10000.
    with h5py.File(msg_lon_file, 'r') as flon:
        lon = flon['LON'][lat1:lat2:stride,lon1:lon2:stride]/10000.
    with h5py.File('hdf5_lsasaf_usgs-igbp_lwmask_msg-disk', 'r') as flwmask:
        lwmask = flwmask['LWMASK'][lat1:lat2:stride,lon1:lon2:stride]
    
    out_shape = lwmask.shape

    if 0:
        ## Set nan to non valid points
        lat[lat==91.] = np.nan
        print('min/max lat:', np.nanmin(lat), np.nanmax(lat))
        lon[lon==91.] = np.nan
        print('min/max lon:', np.nanmin(lon), np.nanmax(lon))
    else:
        ## Remove non valid points
        valid_mask = np.where(lwmask==1)
        lat = lat[valid_mask]
        lon = lon[valid_mask]
        print('min/max lat:', lat.min(), lat.max())
        print('min/max lon:', lon.min(), lon.max())

    return lon, lat, valid_mask, out_shape


def load_etal_lonlat_var(etal_file, stride=100, var=None, date=None, **param):

    t0 = SimpleTimer()

    # already discarding part of ETAL data that is for sure not in the MSG disk
    if 1:
        iLatStartEps = 850
        iLatEndEps   = 17150 
        iLonStartEps = 9000 
        iLonEndEps   = 27000 
    else: # zoom
        #iLatStartEps = 6220
        #iLatEndEps   = 6480 
        #iLonStartEps = 19290 
        #iLonEndEps   = 19530 
        iLatStartEps = 5900
        iLatEndEps   = 6900 
        iLonStartEps = 18900 
        iLonEndEps   = 19900 
    

    print(f'--- Read {etal_file}...')
    with h5py.File(etal_file,'r') as h5f:
        etal_shape  = h5f[var].shape
        etal = h5f[var][iLatStartEps:iLatEndEps:stride,iLonStartEps:iLonEndEps:stride]
    #etal = h5py.File('.\\input_data\\HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_202012250000', 'r')[var]

    ## Debug: set checkerboard
    if 0:
        etal = 6000*(np.indices(etal.shape).sum(axis=0) % 2)
        print(etal)

    ## Debug: plot input etal to image
    if 0:
        var_dic = {'data':etal, 'name':'etal'}
        export_to_image(var_dic, vmin=0, vmax=6000)

    t0('h5_slice')

    ## Filter #1: discard non valid etal
    ## Version without np.where is a bit quicker, you can test it with the following if:
    if 1:
        mask_nonvalid = etal!=-1. # Discard points without valid albedo value (outside projection and in the ocean)
        t0('mask_direct')
    else:
        mask_nonvalid = np.where(etal!=-1.) # Discard points without valid albedo value (outside projection and in the ocean)
        t0('mask_where')

    print(mask_nonvalid.shape)
    s0 = etal.size
    etal = etal[mask_nonvalid]
    s1= etal.size
    print("--- {:.2f} % data masked".format(100*(s0-s1)/s0))
    
    t0('mask_nonvalid_etal')
    
    ## Get lonlat from file (quicker than creating them manually when no downsampling is used)
    if 1:
        metop_lonlat_file = '/mnt/lfs/d30/vegeo/SAT/DATA/EPS/metop_lonlat.nc'
        with h5py.File(metop_lonlat_file, 'r') as flonlat:
            lon = flonlat['lon'][iLatStartEps:iLatEndEps:stride,iLonStartEps:iLonEndEps:stride][mask_nonvalid]
            lat = flonlat['lat'][iLatStartEps:iLatEndEps:stride,iLonStartEps:iLonEndEps:stride][mask_nonvalid]
    
    ## Create lonlat coords manually for etal sinusoidal projection
    else:
        lon = np.linspace(-180, 180, etal_shape[1])[iLonStartEps:iLonEndEps:stride]
        lat = np.linspace(90, -90, etal_shape[0])[iLatStartEps:iLatEndEps:stride]

        lon, lat = np.meshgrid(lon,lat)
        lon = lon/np.cos(np.deg2rad(lat))

        lon = lon[mask_nonvalid] 
        lat = lat[mask_nonvalid]
    
    print('min/max lat:', lat.min(), lat.max())
    print('min/max lon:', lon.min(), lon.max())

    t0('mask_nonvalid_lonlat')

    mask_fov = np.logical_and(np.abs(lon)<81, np.abs(lat)<81)  
    lon = lon[mask_fov] 
    lat = lat[mask_fov]
    etal = etal[mask_fov]
    s1= lon.size
    print("--- {:.2f} % data masked".format(100*(s0-s1)/s0))
    print('min/max lat:', lat.min(), lat.max())
    print('min/max lon:', lon.min(), lon.max())

    t0('mask_fov')

    dic_var = {var:etal}

    t0.show()

    return lon, lat, dic_var


def param_to_string(param, sep1='_', sep2='-'):
    """
    Normalize parameter dict

    Keep only alphanumeric char in parameters and create a unique long string with custom separators.
    Separators are use like this:
    {p1:v1, p2:v2} -> p1<sep2>v1<sep1>p2<sep2>v2
    """
    ## Remove some unwanted parameters if any
    if 'string_exclude' in param.keys():
        for k in param['string_exclude']:
            param.pop(k, None)
    param.pop('string_exclude', None)
    delchars = ''.join(c for c in map(chr, range(256)) if not c.isalnum())
    param_str = sep1.join([f"{k}{sep2}{str(v).translate(str.maketrans('','',delchars))}" for k,v in param.items()])
    return param_str


class Plots:

    def __init__(self, zoom=False):
        if zoom:
            lat1 = 380
            lat2 = 480
            lon1 = 1950
            lon2 = 2050
        else:
            lat1 = None
            lat2 = None
            lon1 = None
            lon2 = None
        self.lat = slice(lat1,lat2)
        self.lon = slice(lon1,lon2)

    def initialize(self):
        plt.clf()
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
    
    def imshow(self, data, plot_param={}, **param):
        self.type = 'imshow'
        self.initialize()

        data = data[self.lat,self.lon]
        ims = self.ax.imshow(data, **plot_param)
        plt.tight_layout()
        plt.colorbar(ims)
        
        self.finalize(**param)

    def scatter(self, data1, data2, plot_param={}, **param):
        self.type = 'scatter'
        self.initialize()

        data1 = data1[self.lat,self.lon]
        data2 = data2[self.lat,self.lon]
        al_max = 9000
        counts, xedges, yedges, im = self.ax.hist2d(data1.ravel(), data2.ravel(), bins=[np.linspace(0, al_max, 200)]*2, **plot_param)
        self.ax.plot([0,al_max], [0,al_max], c='k', lw=1)
        #plt.tight_layout()
        
        self.ax.set_xlim(0,al_max)
        self.ax.set_ylim(0,al_max)
        self.ax.set_xlabel(param['s1'])
        self.ax.set_ylabel(param['s2'])

        plt.axis('square')

        self.finalize(**param)

    def finalize(self, **param):
        from textwrap import wrap
        param_str = '\n'.join(wrap(param_to_string(param, sep1=' | ', sep2=':'), 60))
        self.ax.set_title(param_str, fontsize=10 )
        #self.ax.set_title(param_to_string(param, sep1=' | ', sep2=':'), fontsize=9 )
        
        #self.fig.subplots_adjust(top=0.9)
        plt.tight_layout()
        
        param_str = param_to_string(param)
        im_name = f"res_proj2vtk_{self.type}_{param_to_string(param)}.png"
        plt.savefig(im_name, dpi=200)
        print(f"--- Output image saved to: {im_name}")
    

def export_to_h5(data, **param):
    h5_name = f"res_proj2vtk_output_{param_to_string(param)}.h5"
    with h5py.File(h5_name, 'w') as h5_file:
        h5_file[param['var']] = data
    print(f"--- Output h5 saved to: {h5_name}")

def load_h5_var(h5_path, var, slicing=None):
    if h5_path.is_file():
        with h5py.File(h5_path, 'r') as h5_file:
            if slicing is None:
                data = h5_file[var][:]
            else:
                data = h5_file[var][slicing]
        print(f"--- Read h5 from: {h5_path}")
        return data
    else:
        print(f"--- File: {h5_path} not found.")
        return None
        
def compare_with_real_msg(msg_interp, **param):
    date = param['date']
    var = param['var']
    
    ## MTAL (MTAL-R not available at this date)
    if 1: 
        msg_ref_file = f'/mnt/lfs/d30/vegeo/fransenr/CODES/DATA/HDF5_LSASAF_MSG_ALBEDO-D10_MSG-Disk_{date}0000'
        source_ref = 'RefMTAL'
    
    ## MDAL
    else: 
        msg_ref_file = f'/cnrm/vegeo/juncud/NO_SAVE/aod-test-2/v1/HDF5_LSASAF_MSG_ALBEDO_MSG-Disk_{date}0000'
        source_ref = 'RefMDAL'
    
    with h5py.File(msg_ref_file) as h5ref:
        msg_ref = h5ref[var][:]

    msg_daniel_file = f'/cnrm/vegeo/juncud/scripts/eps_geos_{date}0000'
    with h5py.File(msg_daniel_file) as h5ref:
        msg_daniel = 1e4*h5ref[var][:]

    with h5py.File('hdf5_lsasaf_usgs-igbp_lwmask_msg-disk', 'r') as flwmask:
        lwmask = flwmask['LWMASK'][:]

    nonvalid_mask = np.where(lwmask!=1)
    msg_ref[nonvalid_mask] = -1
    msg_daniel[nonvalid_mask] = -1
    
    albedo = {'cmap':'jet', 'vmin':0, 'vmax':6000}
    albedo_diff = {'cmap':'seismic', 'vmin':-2000, 'vmax':2000}
    albedo_scatter = {'cmin':1, 'cmap':'jet'}

    source_vtk = 'Vtk'
    source_daniel = 'Daniel'

    p = Plots(zoom=0)

    p.imshow(msg_interp,            plot_param=albedo,      **param,                              source=source_vtk)
    p.imshow(msg_ref,               plot_param=albedo,      var=param['var'], date=param['date'], source=source_ref)
    p.imshow(msg_daniel,            plot_param=albedo,      var=param['var'], date=param['date'], source=source_daniel)

    p.imshow(msg_ref-msg_interp,    plot_param=albedo_diff, **param,                              source='diff'+source_ref+source_vtk)
    p.imshow(msg_ref-msg_daniel,    plot_param=albedo_diff, var=param['var'], date=param['date'], source='diff'+source_ref+source_daniel)
    p.imshow(msg_daniel-msg_interp, plot_param=albedo_diff, **param,                              source='diff'+source_daniel+source_vtk)

    p.scatter(msg_ref, msg_interp,    plot_param=albedo_scatter, **param, s1=source_ref, s2=source_vtk)
    p.scatter(msg_ref, msg_daniel,    plot_param=albedo_scatter, **param, s1=source_ref, s2=source_daniel)
    p.scatter(msg_daniel, msg_interp, plot_param=albedo_scatter, **param, s1=source_daniel, s2=source_vtk)

def get_time_range(start, end, days=[], **param):
    dseries = pd.date_range(start, end, freq='D')
    
    if len(days)>0:
        #self.start_01 = self.start.replace(day=1)
        #self.end_31 = self.end.replace(day=pd.Period(self.end, freq='D').days_in_month) # freq='xxx' is here just to avoid a bug
        dseries = dseries[[d in days for d in dseries.day]]

    return dseries

def interpolate_etal_to_msg(etal_file, **param):
        ti = SimpleTimer()

        print('### MSG  extraction (target)')
        msg_lon, msg_lat, msg_valid_mask, msg_shape = load_msg_lonlat(stride=1)
        ti('read MSG')
        
        print('### ETAL extraction (source)')
        etal_lon, etal_lat, etal_dic_var = load_etal_lonlat_var(etal_file, stride=1, **param)
        ti('read ETAL')

        print('### Interpolation')
        interp = vtk_interpolation(**param) 
        ti('interp init')
        interp.set_source(etal_lon, etal_lat, etal_dic_var)
        ti('interp set_source')
        interp.set_target(msg_lon, msg_lat)
        ti('interp set_destination')
        interp.run()
        ti('interp run')
        interp_var = interp.get_output()
        ti('interp get_output')

        data = np.zeros(msg_shape)-1.
        data[msg_valid_mask] = interp_var[param['var']]
        
        export_to_h5(data, **param)
        
        return data

def get_etal_on_msg(etal_path, slicing=None, **param):
    """If cache file exists, load it, otherwise interpolate"""
    for k,v in param.items():
        print(f'{k}:{v}')
    
    ## Cache files are  stored in working directory for now
    h5_path = pathlib.Path(f"res_proj2vtk_output_{param_to_string(param)}.h5")
    data = load_h5_var(h5_path, param['var'], slicing)
    if data is not None:
        return data

    else:
        return interpolate_etal_to_msg(etal_path, **param)

def load_mtalr(mtalr_path, slicing=None, **param):
    return load_h5_var(mtalr_path, param['var'], slicing)


def process_etal_series(**param):

    mtalr_root = pathlib.Path('/cnrm/vegeo/SAT/DATA/MSG/Reprocessed-on-2017/MTAL/2015')
    etal_root = pathlib.Path('/cnrm/vegeo/SAT/DATA/EPS/NRT-Operational/ETAL/2015')
    
    
    ## Input dates
    days = [5,15,25]
    dseries = get_time_range(**param, days=days)
    print(dseries)
    print(len(dseries))
    print(param_to_string(param))
    
    ## Init main dataframe with desired dates
    empty_list = [None]*len(dseries)
    df = pd.DataFrame({'date': dseries, 'etal_path':empty_list, 'mtalr_path':empty_list})
    df = df.set_index('date')
    print(df)

    ## Fill the dataframe with paths
    for p in mtalr_root.glob('**/*MSG_ALBEDO*'):
        file_date = dt.datetime.strptime(p.parts[-1].split('_')[-1], '%Y%m%d%H%M')
        if file_date in df.index:
            print(f'{file_date} IN dataframe')
            df.at[file_date, 'mtalr_path'] = p
        else:
            print(f'{file_date} NOT IN dataframe')

    for p in etal_root.glob('**/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE*'):
        # [-2] is to remove "03" at the end of the filename (is it seconds ? it does not exist in the current etal files from lsa saf website)
        file_date = dt.datetime.strptime(p.parts[-1].split('_')[-1][:-2], '%Y%m%d%H%M')
        if file_date in df.index:
            print(f'{file_date} IN dataframe')
            df.at[file_date, 'etal_path'] = p
        else:
            print(f'{file_date} NOT IN dataframe')

    print(df)

    ## Get ETAL on MSG grid 
    interp_param = {
            'var' : 'AL-BB-BH',
            'date' : '20201205',
            'kernel' : 'inverse_distance',
            #'kernel' : ['mean','inverse_distance','gaussian'],
            #'radius' : [5,10],
            'radius' : 5,
            'null_points' : 'closest',
            #'null_points' : -1.,
           }

    ## Data accumulation

    ## DEBUG: subsampling
    if 1:
        slicing = (slice(None, None, 4), slice(None, None, 4))
        interp_param['slicing'] = slicing
        interp_param['string_exclude'] = ['slicing']

    with h5py.File('hdf5_lsasaf_usgs-igbp_lwmask_msg-disk', 'r') as flwmask:
        try:
            lwmask = flwmask['LWMASK'][slicing]
        except:
            lwmask = flwmask['LWMASK'][:]

    nonvalid_mask = np.where(lwmask!=1)
    res_bias = np.zeros_like(lwmask, dtype=np.float)
    nbias = 0


    ## Plot params
    albedo = {'cmap':'jet', 'vmin':0, 'vmax':6000}
    albedo_diff = {'cmap':'seismic', 'vmin':-2000, 'vmax':2000}
    source_vtk = 'VtkETAL'
    source_ref = 'RefMTALR'
    p = Plots(zoom=0)
    for index,row in df.iterrows():
        if (row['etal_path'] is None) or (row['mtalr_path'] is None):
            continue

        interp_param['date'] = index.strftime('%Y%m%d')
        etal_on_msg = get_etal_on_msg(row['etal_path'], **interp_param)
        mtalr_on_msg = load_mtalr(row['mtalr_path'],  **interp_param)
        mtalr_on_msg[nonvalid_mask] = -1
        bias = mtalr_on_msg-etal_on_msg
        res_bias += bias
        nbias += 1
        print(pd.DataFrame(etal_on_msg.ravel()).describe().applymap('{:.2f}'.format))
        #print(pd.DataFrame(mtalr_on_msg.ravel()).describe().applymap('{:.2f}'.format))

        p.imshow(etal_on_msg, plot_param=albedo, **interp_param, source=source_vtk)
        p.imshow(mtalr_on_msg, plot_param=albedo, var=interp_param['var'], date=interp_param['date'], source=source_ref)
        p.imshow(bias, plot_param=albedo_diff, **interp_param, source='diff'+source_ref+source_vtk)

    res_bias /= nbias
    param2 = {k:v for k,v in interp_param.items()}
    param2.pop('date', None)
    param2['start'] = param['start']
    param2['end'] = param['end']
    p.imshow(res_bias, plot_param=albedo_diff, **param2, source='diff'+source_ref+source_vtk)

    
    ## Init empty data containers to store results
    # temporal domain (10x10 km area)
    # spatial domain (global bias map)

def main(param):

    ti = SimpleTimer()



    if 1:
        process_etal_series(**param)        
    
    ## Load cache file
    else:
        data = load_h5(**param)

    #print('### Compare results')
    #compare_with_real_msg(data, **param)
    #ti('comparison')


    ti.show()


def main_vtk(param):

    ti = SimpleTimer() 


    if 1:

        print('### MSG  extraction (target)')
        tlon, tlat, valid_mask, msg_shape = msg_to_vtk(stride=1, lonlat_only=True)
        ti('read MSG')
        
        print('### ETAL extraction (source)')
        slon, slat, dic_var = etal_to_vtk(stride=1, lonlat_only=True, **param)
        ti('read ETAL')

        print('### Interpolation')
        interp = vtk_interpolation(**param) 
        ti('interp init')
        interp.set_source(slon, slat, dic_var)
        ti('interp set_source')
        interp.set_target(tlon, tlat)
        ti('interp set_destination')
        interp.run()
        ti('interp run')
        interp_var = interp.get_output()
        ti('interp get_output')

        data = np.zeros(msg_shape)-1.
        data[valid_mask] = interp_var[param['var']]
        
    ## Load cache file
    else:
        data = load_h5(**param)

    print('### Compare results')
    compare_with_real_msg(data, **param)
    ti('comparison')

    ti.show()


def combine_param(param_dic):
    """
    Iterator returning a dict made of combination of all list values in input dict, using the same keys.
    An id is also add in the output dict.

    Example:
    input dict: {'a':[1,2], 'b'=[x,y]}
    output dicts:
    #1: {'id':'00', 'a':1, 'b':x}
    #2: {'id':'01', 'a':1, 'b':y}
    #3: {'id':'02', 'a':2, 'b':x}
    #4: {'id':'03', 'a':2, 'b':y}
    """
    ii = 0
    n_combi = len(list(itertools.product(*list(param_dic.values()))))

    for v in itertools.product(*list(param_dic.values())):
        print(f'--- Parameter combination {ii+1}/{n_combi}')
        tmp = {}
        tmp['id'] = str(ii).zfill(len(str(n_combi))+1)
        for ik,k in enumerate(param_dic.keys()):
            tmp[k] = v[ik]
        yield tmp
        
        ii+=1


if __name__=='__main__':
    
    test = {
            'var' : ['AL-BB-BH'],
            'start' : ['2015-01-05'],
            'end' : ['2015-12-25'],
            'size' : [10],
           }

    for param in combine_param(test):
        for k,v in param.items():
            print(k, ':', v)
        try:
            main(param)
        except Exception as e:
            print('--- ERROR IN PROCESSING ---')
            print(traceback.format_exc())
            #logging.error(traceback.format_exc())

        print('---------------')
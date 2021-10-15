import h5py
import numpy as np

from pyhdf.SD import SD, SDC # To read HDF-EOS (HDF4 type) MODIS format.

import os,sys
from tools import SimpleTimer

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    
    lat_slice = slice(iLatStartEps, iLatEndEps, stride)
    lon_slice = slice(iLonStartEps, iLonEndEps, stride)

    print(f'--- Read {etal_file}...')
    with h5py.File(etal_file,'r') as h5f:
        etal_shape  = h5f[var].shape
        etal = h5f[var][lat_slice, lon_slice]

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
        mask_valid = etal!=-1. # Discard points without valid albedo value (outside projection and in the ocean)
        t0('mask_direct')
    else:
        mask_valid = np.where(etal!=-1.) # Discard points without valid albedo value (outside projection and in the ocean)
        t0('mask_where')

    print(mask_valid.shape)
    s0 = etal.size
    etal = etal[mask_valid]
    s1= etal.size
    print("--- {:.2f} % data masked".format(100*(s0-s1)/s0))
    
    t0('mask_valid_etal')
    
    ## Get lonlat from file (quicker than creating them manually when no downsampling is used)
    if 1:
        metop_lonlat_file = '/mnt/lfs/d30/vegeo/SAT/DATA/EPS/metop_lonlat.nc'
        with h5py.File(metop_lonlat_file, 'r') as flonlat:
            lon = flonlat['lon'][lat_slice, lon_slice][mask_valid]
            lat = flonlat['lat'][lat_slice, lon_slice][mask_valid]
    
        #lon, lat = load_etal_lonlat(lat_slice, lon_slice, mask=mask_valid)
    
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


def param_to_string(param, sep1='_', sep2='-', shorten=False):
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
    ## Shorten key if required
    if shorten:
        param = {k[:3]:v for k,v in param.items()}
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
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
    


    def imshow(self, data, plot_param={}, noaxis=False, **param):
        self.type = 'imshow'
        self.initialize()

        data = data[self.lat,self.lon]
        ims = self.ax.imshow(data, **plot_param)
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cb = plt.colorbar(ims, cax=cax)
        
        self.finalize(noaxis, **param)

    def scatter(self, data1, data2, plot_param={}, noaxis=False, **param):
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

        self.finalize(noaxis, **param)

    def finalize(self, noaxis=False, **param):
        from textwrap import wrap
        param_str = '\n'.join(wrap(param_to_string(param, sep1=' | ', sep2=':'), 60))
        self.ax.set_title(param_str, fontsize=10 )
        plt.tight_layout()
        
        im_name = f"res_proj2vtk_{self.type}_{param_to_string(param, shorten=True)}.png"
        if noaxis:
            self.save_no_whitespace(im_name)
        else:
            plt.savefig(im_name, bbox_inches='tight')
        print(f"--- Output image saved to: {im_name}")
        plt.close(self.fig)
   
    def save_no_whitespace(self, filepath):
        '''Save the current image with no whitespace'''
        plt.subplots_adjust(0,0,1,1,0,0)
        if len(self.fig.axes)==2:
            plt.delaxes(self.fig.axes[1]) # axes[0] correspond to the main plot.
        self.ax.set_title("")
        self.ax.axis('off')
        self.ax.margins(0,0)
        self.ax.xaxis.set_major_locator(plt.NullLocator())
        self.ax.yaxis.set_major_locator(plt.NullLocator())
        self.fig.savefig(filepath, pad_inches=0, bbox_inches='tight', dpi=200)

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
        
        export_to_h5(data, 'etal2msg', **param)
        
        return data

class SatelliteTools:
    def __init__(self, product, var, slicing=slice(None)):
        self.product = product
        self.var = var
        self.slicing= slicing

        self._full_shape = None
        self._shape = None

        self._ground_mask = True 
        self._mask = True

        self._lat = None
        self._lon = None

    @property
    def full_shape(self):
        """Get global shape before slicing"""
        if self._full_shape is None:
            with h5py.File(self.ground_mask_conf['file'], 'r') as flw:
                self._full_shape = flw[self.ground_mask_conf['var']].shape
        return self._full_shape

    @property
    def shape(self):
        """Compute shape after slicing without loading data"""
        if self._shape is None:
            out_shape = [] 
            for s,l in zip(self.slicing, self.full_shape):
                out_shape.append(len(range(*s.indices(l))))
            self._shape = tuple(out_shape)
        return self._shape
    
    @property
    def mask(self):
        if self._mask is True:
            self._mask = self.ground_mask
        return self._mask

    @mask.setter
    def mask(self, value):
        self._mask = value

    @property
    def ground_mask(self):
        """
        Load ground mask if asked 
        """
        self.ti()
        if self._ground_mask is True:
            if self.ground_mask_conf['type'] is not None:
                print(f"--- Read {self.product.upper()} mask in {self.ground_mask_conf['file']} ...")
                with h5py.File(self.ground_mask_conf['file'], 'r') as flw:
                    lwmask = flw[self.ground_mask_conf['var']][self.slicing]
                if self.ground_mask_conf['type']=='land': mask_value = 1
                self._ground_mask = lwmask==mask_value
            ## Update main mask
            self.mask = self._mask & self._ground_mask
        return self._ground_mask
        self.ti('get_ground_mask')

    @property
    def lat(self):
        if self._lat is None:
            print(f"--- Read {self.product.upper()} lat in {self.lat_conf['file']} ...")
            with h5py.File(self.lat_conf['file'], 'r') as fl:
                self._lat = fl[self.lat_conf['var']][self.slicing]*self.lat_conf['scaling']
        return self._lat[self.mask]

    @property
    def lon(self):
        if self._lon is None:
            print(f"--- Read {self.product.upper()} lon in {self.lon_conf['file']} ...")
            with h5py.File(self.lon_conf['file'], 'r') as fl:
                self._lon = fl[self.lon_conf['var']][self.slicing]*self.lon_conf['scaling']
        return self._lon[self.mask]

    def interpolate_on(self, target, source_file, use_cache, **param):
        """If cache file exists, load it, otherwise interpolate"""
        
        for k,v in param.items():
            print(f'{k}:{v}')

        interp_type = f'{self.product}2{type(target).__name__}'

        interpolation_required = ~use_cache
        if use_cache:
            ## Cache files are  stored in working directory for now
            h5_path = pathlib.Path(f"./cache/cache_{interp_type}_{param_to_string(param)}.h5")
            #data = self.load_h5_cache(h5_path, param['var'], target.slicing)
            data = self.load_h5_cache(h5_path, param['var']) # return cache as it without slicing (cache may have already been sliced)
            if data is not None:
                return data
            else:
                interpolation_required = True

        if interpolation_required:

            self.get_data(source_file)

            ti = SimpleTimer()

            ## Interpolation
            print('### Interpolation')
            interp = vtk_interpolation(**param) 
            ti('interp init')
            interp.set_source(self.lon, self.lat, {self.var:self.data})
            ti('interp set_source')
            interp.set_target(target.lon, target.lat)
            ti('interp set_destination')
            interp.run()
            ti('interp run')
            interp_var = interp.get_output()
            ti('interp get_output')

            ## Extract and save data
            data = np.full(target.shape, -1.)
            data[target.ground_mask] = interp_var[param['var']]
            self.export_to_h5(data, interp_type, **param)
            
            return data

    def export_to_h5(self, data, name, **param):
        h5_name = f"./cache/cache_{name}_{param_to_string(param)}.h5"
        with h5py.File(h5_name, 'w') as h5_file:
            ## Proper way to save data as integer
            if 1:
                dataset = h5_file.create_dataset(
                    param['var'], chunks=True, compression='gzip',
                    fletcher32=True, shape=data.shape, dtype=int)
                dataset.attrs.create('SCALING_FACTOR', 1e4)
                dataset[:,:] = data.astype('i4')
            ## else quick write keeping the type of data
            else:
                h5_file[param['var']] = data
        print(f"--- Output h5 saved to: {h5_name}")

    def load_h5_cache(self, h5_path, var, slicing=slice(None)):
        if h5_path.is_file():
            with h5py.File(h5_path, 'r') as h5_file:
                data = h5_file[var][slicing]
            print(f"--- Read h5 from: {h5_path}")
            return data
        else:
            print(f"--- File: {h5_path} not found.")
            return None
            
    def describe(self):
        print(self.data.shape)
        print(self.lon.shape)

        print_stats(self.data)
        print(np.count_nonzero(self.data==-1))
    

class MODIS(SatelliteTools):
    def __init__(self, product, var, slicing=slice(None), mask_type='land'):
        
        super().__init__(product, var, slicing)

        self.ti = SimpleTimer()

        self.ground_mask_conf = {'file': '/mnt/lfs/d30/vegeo/fransenr/CODES/DATA/NO_SAVE/ETAL/etal_lwmask.h5',
                                 'var': 'lwmask',
                                 'type': mask_type}
        self.lat_conf = {'file': '/mnt/lfs/d30/vegeo/SAT/DATA/EPS/metop_lonlat.nc',
                         'var': 'lat',
                         'scaling': 1.}
        self.lon_conf = {'file': '/mnt/lfs/d30/vegeo/SAT/DATA/EPS/metop_lonlat.nc',
                         'var': 'lon',
                         'scaling': 1.}

        self.data = None

    def get_data(self, data_file, mask_value=None):
        """
        Read a data file and compute a variable mask
        
        :mask_value: None: return a flatten array with only valid data
                     <number> : return an array with the shape of the read data and <number> used as non valid data.
        """
        self.ti()
        print(f'--- Read {self.product.upper()}/{self.var} in {data_file} ...')
        with h5py.File(data_file,'r') as h5f:
            self.data0 = h5f[self.var][self.slicing]
        self.mask = self.ground_mask & (self.data0!=-1)
        if mask_value is None:
            self.data = self.data0[self.mask]
        else:
            self.data = np.zeros(self.shape) + mask_value
            self.data[self.mask] = self.data0[self.mask]
        self.ti('get_data')

    def load_hdf_var(self, data_file, var, slicing):
        with h5py.File(data_file,'r') as h5f:
            self.data0 = h5f[var][slicing]

class EPS(SatelliteTools):
    def __init__(self, product, var, slicing=slice(None), mask_type='land'):
        
        super().__init__(product, var, slicing)

        self.ti = SimpleTimer()

        self.ground_mask_conf = {'file': '/mnt/lfs/d30/vegeo/fransenr/CODES/DATA/NO_SAVE/ETAL/etal_lwmask.h5',
                                 'var': 'lwmask',
                                 'type': mask_type}
        self.lat_conf = {'file': '/mnt/lfs/d30/vegeo/SAT/DATA/EPS/metop_lonlat.nc',
                         'var': 'lat',
                         'scaling': 1.}
        self.lon_conf = {'file': '/mnt/lfs/d30/vegeo/SAT/DATA/EPS/metop_lonlat.nc',
                         'var': 'lon',
                         'scaling': 1.}

        self.data = None

    def get_data(self, data_file, mask_value=None):
        """
        Read a data file and compute a variable mask
        
        :mask_value: None: return a flatten array with only valid data
                     <number> : return an array with the shape of the read data and <number> used as non valid data.
        """
        self.ti()
        print(f'--- Read {self.product.upper()}/{self.var} in {data_file} ...')
        with h5py.File(data_file,'r') as h5f:
            self.data0 = h5f[self.var][self.slicing]
        self.mask = self.ground_mask & (self.data0!=-1)
        if mask_value is None:
            self.data = self.data0[self.mask]
        else:
            self.data = np.zeros(self.shape) + mask_value
            self.data[self.mask] = self.data0[self.mask]
        self.ti('get_data')

    def load_hdf_var(self, data_file, var, slicing):
        with h5py.File(data_file,'r') as h5f:
            self.data0 = h5f[var][slicing]

class MSG(SatelliteTools):
    def __init__(self, product, var, slicing=slice(None), mask_type='land'):

        super().__init__(product, var, slicing)

        self.ti = SimpleTimer()

        self.ground_mask_conf = {'file': 'hdf5_lsasaf_usgs-igbp_lwmask_msg-disk',
                                 'var': 'LWMASK',
                                 'type': mask_type}
        self.lat_conf = {'file': '/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LAT_MSG-Disk_4bytesPrecision',
                         'var': 'LAT',
                         'scaling': 0.0001}
        self.lon_conf = {'file': '/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LON_MSG-Disk_4bytesPrecision',
                         'var': 'LON',
                         'scaling': 0.0001}

        self.data = None
        

def get_modis_on_etal(modis_path, slicing=None, **param):
    """If cache file exists, load it, otherwise interpolate"""
    for k,v in param.items():
        print(f'{k}:{v}')
    
    ## Cache files are  stored in working directory for now
    h5_path = pathlib.Path(f"cache_modis2etal_{param_to_string(param)}.h5")
    data = load_h5_var(h5_path, param['var'], slicing)
    if data is not None:
        return data

    else:
        return interpolate_modis_to_etal(modis_path, **param)

def get_all_filenames(**param):
    path_dic = {}
    path_dic['mtalr'] = {'root':pathlib.Path('/cnrm/vegeo/SAT/DATA/MSG/Reprocessed-on-2017/MTAL'),
                              'format':'%Y/%m/%d/HDF5_LSASAF_MSG_ALBEDO-D10_MSG-Disk_%Y%m%d0000'}
    #path_dic['etal']  = {'root':pathlib.Path('/mnt/lfs/d30/vegeo/fransenr/CODES/DATA/NO_SAVE/ETAL'),
    #                          'format':'%Y/%m/%d/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_%Y%m%d0000'}
    path_dic['etalr'] = {'root':pathlib.Path('/cnrm/vegeo/fransenr/CODES/DATA/NO_SAVE/EPS_Reprocess/ETAL'),
                              'format':'%Y/%m/%d/HDF5_LSASAF_M02-AVHR_ETAL_GLOBE_%Y%m%d0000'}
    path_dic['modis'] = {'root':pathlib.Path('/mnt/lfs/d30/vegeo/SAT/DATA/MODIS/MODIS-MCD43D51/tiles'),
                              'format':'%Y/MCD43D51.A%Y%j.006.?????????????.hdf'}
    
    
    ## Input dates
    days = [5,15,25]
    dseries = get_time_range(**param, days=days)
    print(dseries)
    print(len(dseries))
    print(param_to_string(param))
    
    ## Init main dataframe with desired dates
    empty_list = [None]*len(dseries)
    #df = pd.DataFrame({'date': dseries, 'etal_path':empty_list, 'mtalr_path':empty_list})
    df = pd.DataFrame({'date': dseries, **{f'{k}_path':empty_list for k in path_dic.keys()}})
    df = df.set_index('date')
    print(df)

    ## Check paths and fill the dataframe with valids
    for prod in path_dic.keys():
        for d in df.index:
            root = path_dic[prod]['root']
            path_format = path_dic[prod]['format']
            
            if prod=='modis': # Unpredictable part in MODIS file name, need to glob...
                fpath = list(root.glob(d.strftime(path_format)))
                if len(fpath)==1:
                    fpath = fpath[0]
                    print(f"{prod.upper()} {d:%Y%m%d} found.")
                    df.at[d, f'{prod}_path'] = fpath
                elif len(fpath)==0:
                    print(f"{prod.upper()} {d:'%Y%m%d'} NOT found.")
                else:
                    print(f"{prod.upper()} {d:'%Y%m%d'} ERROR: several files found. Exiting.")
                    sys.exit()

            else:
                fpath = root/d.strftime(path_format)
                if fpath.exists():
                    print(f"{prod.upper()} {d:%Y%m%d} found.")
                    df.at[d, f'{prod}_path'] = fpath
                else:
                    print(f"{prod.upper()} {d:'%Y%m%d'} NOT found.")
    
    print(df)

    return df


def process_etal_series(**param):

    df = get_all_filenames(**param)

    ## Get ETAL on MSG grid 
    interp_param = {
            'var' : 'AL-BB-BH',
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
        #slicing_msg = (slice(None, None, 2), slice(None, None, 2))
        slicing_msg = (slice(None), slice(None))
        slicing_eps = (slice(None, None, 10), slice(None, None, 10))
        #slicing_eps = (slice(850, 17150, 1), slice(9000, 27000, 1)) # Remove points not on MSG disc
        #slicing = (slice(380, 480), slice(1950, 2050)) # France
        #slicing = (slice(50, 700), slice(1550, 3250)) # Euro
        #slicing = (slice(700, 1850), slice(1240, 3450)) # NAfr
        #slicing = (slice(1850, 3040), slice(2140, 3350)) # SAfr
        #slicing = (slice(1460, 2970), slice(40, 740)) # SAme
        interp_param['slicing'] = slicing_eps
        interp_param['string_exclude'] = ['slicing']

    ## Create product objects
    mtalr = MSG(product='mtalr', var=interp_param['var'], slicing=slicing_msg, mask_type='land')   
    etalr = EPS(product='etalr', var=interp_param['var'], slicing=slicing_eps, mask_type='land')   
    etal = EPS(product='etal', var=interp_param['var'], slicing=slicing_eps, mask_type='land')   
    etal = EPS(product='etal', var=interp_param['var'], slicing=slicing_eps, mask_type='land')   

    res_bias = np.zeros_like(mtalr.shape, dtype=float)
    res_bias2 = np.zeros_like(mtalr.shape, dtype=float)
    nbias = 0

    ## Plot params
    albedo = {'cmap':'jet', 'vmin':0, 'vmax':0.6}
    blank = {'cmap':'binary', 'vmin':0, 'vmax':0.6}
    albedo_diff = {'cmap':'seismic', 'vmin':-0.15, 'vmax':0.15}
    albedo_biasstd = {'cmap':'jet', 'vmin':0.0, 'vmax':0.08}
    source_vtk = 'VtkETAL'
    source_ref = 'RefMTALR'
    p = Plots(zoom=0)

    # Recursive plot
    rec_plot = 1
    if rec_plot:
        #recfig = plt.figure()
        #recax1 = recfig.add_subplot(211)
        #recax2 = recfig.add_subplot(212)
        recfig, (recax1, recax2) = plt.subplots(2, 1, sharex=True)

        recax1.axhline(20.0, ls='--', c='r')
        recax1.axhline(10.0, ls='--', c='y')
        recax1.axhline(5.0, ls='--', c='g')
        recax1.axhline(0.0, ls='-', c='k')
        recax1.axhline(-5.0, ls='--', c='g')
        recax1.axhline(-10.0, ls='--', c='y')
        recax1.axhline(-20.0, ls='--', c='r')
        
        recax2.axhline(0.015, ls='--', c='r')
        recax2.axhline(0.0, ls='-', c='k')
        recax2.axhline(-0.015, ls='--', c='r')

    means = [[],[]]
    stds = [[],[]]
    dates = []

    tglob = SimpleTimer()

    for i, (index,row) in enumerate(df.iterrows()):
        print(f'\n================ {index:%Y-%m-%d} ================')

        params = {}
        params['id'] = f'{i:03d}'
        params['date'] = f'{index:%Y%m%d}'

        ## DEBUG ON: show input data (disable to avoid loading input data if cache file already exists)
        if 1:
            if row['etalr_path'] is None:
            #if row['etal_path'] is None:
                p.imshow(np.zeros(etalr.shape), plot_param=blank, noaxis=False, **params, source='None')
                #p.imshow(np.zeros(etal.shape), plot_param=blank, noaxis=False, **params, source='None')
                continue
            
            etalr.get_data(row['etalr_path'], mask_value=-1)
            p.imshow(etalr.data*0.0001, plot_param=albedo, noaxis=False, **params, source='etalr')
            #etal.get_data(row['etal_path'], mask_value=-1)
            #p.imshow(etal.data*0.0001, plot_param=albedo, noaxis=False, **params, source='etal')
        ## DEBUG OFF
        else:
            if row['etalr_path'] is None:
                print('--- No data file found for this date.')
                tglob(params['id'])
                continue

        #etalr_on_msg = etalr.interpolate_on(mtalr, source_file=row['etalr_path'], use_cache=True, date=params['date'], **interp_param)

        if 0:
            p.imshow(etalr_on_msg*0.0001, plot_param=albedo, noaxis=False, **params, **interp_param, source=source_vtk)

        tglob(params['id'])

        if (i==6) & 0:
            print('--- BREAK LOOP')
            break

    tglob.show()
    sys.exit()

    for index,row in df.iterrows():
        if (row['etalr_path'] is None) or (row['mtalr_path'] is None):
            continue

        etal = ETAL(slicing, interp_param['var'])
        etal.get_data(row['etal_path'])
        etal.get_latlon()
        etal.describe()
        sys.exit()

        interp_param['date'] = index.strftime('%Y%m%d')
        etal_on_msg = get_etal_on_msg(row['etal_path'], **interp_param)
        mtalr_on_msg = load_mtalr(row['mtalr_path'],  **interp_param)
        
        # Filtering
        mtalr_on_msg[nonvalid_mask] = -1
        mtalr_on_msg[etal_on_msg==-1] = -1
        etal_on_msg[mtalr_on_msg==-1] = -1
        mtalr_on_msg = mtalr_on_msg*0.0001
        etal_on_msg = etal_on_msg*0.0001

        print('mtalr:',np.count_nonzero(mtalr_on_msg<0), 'etal:', np.count_nonzero(etal_on_msg<0))

        # Compute bias
        bias = mtalr_on_msg-etal_on_msg
        res_bias += bias
        res_bias2 += bias*bias
        nbias += 1
        print_stats(mtalr_on_msg, etal_on_msg, label=['mtalr', 'etal'])

        if 1:
            p.imshow(etal_on_msg, plot_param=albedo, noaxis=True, **interp_param, source=source_vtk)
            p.imshow(mtalr_on_msg, plot_param=albedo, noaxis=True, var=interp_param['var'], date=interp_param['date'], source=source_ref)
            p.imshow(bias, plot_param=albedo_diff, noaxis=True, **interp_param, source='diff'+source_ref+source_vtk)

        # Fill recursive plot
        if rec_plot:
            #recax.violinplot([bias.ravel()], positions=[nbias], showmeans=False, showmedians=True, showextrema=True, widths=0.9)
            
            mask_015 = mtalr_on_msg.ravel()>0.15
            bias_above = 100*bias.ravel()[mask_015]/mtalr_on_msg.ravel()[mask_015]
            bias_below = bias.ravel()[~mask_015]

            print_stats(bias_above, bias_below, label=['bias_above_0.15', 'bias_below_0.15'], fmt='{:.3f}')

            means[0].append(bias_above.mean())
            means[1].append(bias_below.mean())
            stds[0].append(bias_above.std())
            stds[1].append(bias_below.std())
            dates.append(index.strftime('%Y-%m-%d'))

            #recax1.boxplot([bias_above], positions=[nbias], labels=[index.strftime('%Y-%m-%d')], widths=0.75, whis='range')
            #recax2.boxplot([bias_below], positions=[nbias], labels=[index.strftime('%Y-%m-%d')], widths=0.75, whis='range')
            
            recax1.vlines(nbias-1, means[0][-1]-stds[0][-1], means[0][-1]+stds[0][-1])
            recax2.vlines(nbias-1, means[1][-1]-stds[1][-1], means[1][-1]+stds[1][-1])
            
            #x = np.random.normal(nbias, 0.04, size=len(bias.ravel()))
            #recax.scatter(x, bias.ravel(), alpha=0.2)

            if (nbias == 4) and 1:
                break

    if rec_plot:
        recax1.plot(means[0], 'o:', c='k', markeredgecolor='k', markerfacecolor='r')
        recax2.plot(means[1], 'o:', c='k', markeredgecolor='k', markerfacecolor='r')

        recax1.set_ylabel('Relative bias(%)')
        recax2.set_ylabel('MBE and std')
        recax1.set_xticks(range(nbias))
        recax1.set_xticklabels(dates)
        recax1.xaxis.set_tick_params(rotation=90)
        recax2.set_xticks(range(nbias))
        recax2.set_xticklabels(dates)
        recax2.xaxis.set_tick_params(rotation=90)
        plt.tight_layout() 
        plt.show()
        #sys.exit()
    
    res_bias /= nbias
    res_bias2 /= nbias
    param2 = {k:v for k,v in interp_param.items()}
    param2.pop('date', None)
    param2['start'] = param['start']
    param2['end'] = param['end']
    p.imshow(res_bias, plot_param=albedo_diff, **param2, source='bias'+source_ref+source_vtk)
    p.imshow(np.sqrt(res_bias2-res_bias**2), plot_param=albedo_biasstd, **param2, source='biasSTD'+source_ref+source_vtk)

    
    ## Init empty data containers to store results
    # temporal domain (10x10 km area)
    # spatial domain (global bias map)


def process_etal_series_OLD(**param):

    df = get_all_filenames(**param)

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
        slicing = (slice(None, None, 100), slice(None, None, 100))
        #slicing = (slice(380, 480), slice(1950, 2050)) # France
        #slicing = (slice(50, 700), slice(1550, 3250)) # Euro
        #slicing = (slice(700, 1850), slice(1240, 3450)) # NAfr
        #slicing = (slice(1850, 3040), slice(2140, 3350)) # SAfr
        #slicing = (slice(1460, 2970), slice(40, 740)) # SAme
        interp_param['slicing'] = slicing
        interp_param['string_exclude'] = ['slicing']

    with h5py.File('hdf5_lsasaf_usgs-igbp_lwmask_msg-disk', 'r') as flwmask:
        try:
            lwmask = flwmask['LWMASK'][slicing]
        except:
            lwmask = flwmask['LWMASK'][:]

    nonvalid_mask = np.where(lwmask!=1)
    res_bias = np.zeros_like(lwmask, dtype=float)
    res_bias2 = np.zeros_like(lwmask, dtype=float)
    nbias = 0

    ## Plot params
    albedo = {'cmap':'jet', 'vmin':0, 'vmax':0.6}
    albedo_diff = {'cmap':'seismic', 'vmin':-0.15, 'vmax':0.15}
    albedo_biasstd = {'cmap':'jet', 'vmin':0.0, 'vmax':0.08}
    source_vtk = 'VtkETAL'
    source_ref = 'RefMTALR'
    p = Plots(zoom=0)

    # Recursive plot
    rec_plot = 1
    if rec_plot:
        #recfig = plt.figure()
        #recax1 = recfig.add_subplot(211)
        #recax2 = recfig.add_subplot(212)
        recfig, (recax1, recax2) = plt.subplots(2, 1, sharex=True)

        recax1.axhline(20.0, ls='--', c='r')
        recax1.axhline(10.0, ls='--', c='y')
        recax1.axhline(5.0, ls='--', c='g')
        recax1.axhline(0.0, ls='-', c='k')
        recax1.axhline(-5.0, ls='--', c='g')
        recax1.axhline(-10.0, ls='--', c='y')
        recax1.axhline(-20.0, ls='--', c='r')
        
        recax2.axhline(0.015, ls='--', c='r')
        recax2.axhline(0.0, ls='-', c='k')
        recax2.axhline(-0.015, ls='--', c='r')

    means = [[],[]]
    stds = [[],[]]
    dates = []

    for index,row in df.iterrows():
        if (row['etal_path'] is None) or (row['mtalr_path'] is None):
            continue

        etal = ETAL(slicing, interp_param['var'])
        etal.get_data(row['etal_path'])
        etal.get_latlon()
        etal.describe()
        sys.exit()

        interp_param['date'] = index.strftime('%Y%m%d')
        etal_on_msg = get_etal_on_msg(row['etal_path'], **interp_param)
        mtalr_on_msg = load_mtalr(row['mtalr_path'],  **interp_param)
        
        # Filtering
        mtalr_on_msg[nonvalid_mask] = -1
        mtalr_on_msg[etal_on_msg==-1] = -1
        etal_on_msg[mtalr_on_msg==-1] = -1
        mtalr_on_msg = mtalr_on_msg*0.0001
        etal_on_msg = etal_on_msg*0.0001

        print('mtalr:',np.count_nonzero(mtalr_on_msg<0), 'etal:', np.count_nonzero(etal_on_msg<0))

        # Compute bias
        bias = mtalr_on_msg-etal_on_msg
        res_bias += bias
        res_bias2 += bias*bias
        nbias += 1
        print_stats(mtalr_on_msg, etal_on_msg, label=['mtalr', 'etal'])

        if 1:
            p.imshow(etal_on_msg, plot_param=albedo, noaxis=True, **interp_param, source=source_vtk)
            p.imshow(mtalr_on_msg, plot_param=albedo, noaxis=True, var=interp_param['var'], date=interp_param['date'], source=source_ref)
            p.imshow(bias, plot_param=albedo_diff, noaxis=True, **interp_param, source='diff'+source_ref+source_vtk)

        # Fill recursive plot
        if rec_plot:
            #recax.violinplot([bias.ravel()], positions=[nbias], showmeans=False, showmedians=True, showextrema=True, widths=0.9)
            
            mask_015 = mtalr_on_msg.ravel()>0.15
            bias_above = 100*bias.ravel()[mask_015]/mtalr_on_msg.ravel()[mask_015]
            bias_below = bias.ravel()[~mask_015]

            print_stats(bias_above, bias_below, label=['bias_above_0.15', 'bias_below_0.15'], fmt='{:.3f}')

            means[0].append(bias_above.mean())
            means[1].append(bias_below.mean())
            stds[0].append(bias_above.std())
            stds[1].append(bias_below.std())
            dates.append(index.strftime('%Y-%m-%d'))

            #recax1.boxplot([bias_above], positions=[nbias], labels=[index.strftime('%Y-%m-%d')], widths=0.75, whis='range')
            #recax2.boxplot([bias_below], positions=[nbias], labels=[index.strftime('%Y-%m-%d')], widths=0.75, whis='range')
            
            recax1.vlines(nbias-1, means[0][-1]-stds[0][-1], means[0][-1]+stds[0][-1])
            recax2.vlines(nbias-1, means[1][-1]-stds[1][-1], means[1][-1]+stds[1][-1])
            
            #x = np.random.normal(nbias, 0.04, size=len(bias.ravel()))
            #recax.scatter(x, bias.ravel(), alpha=0.2)

            if (nbias == 4) and 1:
                break

    if rec_plot:
        recax1.plot(means[0], 'o:', c='k', markeredgecolor='k', markerfacecolor='r')
        recax2.plot(means[1], 'o:', c='k', markeredgecolor='k', markerfacecolor='r')

        recax1.set_ylabel('Relative bias(%)')
        recax2.set_ylabel('MBE and std')
        recax1.set_xticks(range(nbias))
        recax1.set_xticklabels(dates)
        recax1.xaxis.set_tick_params(rotation=90)
        recax2.set_xticks(range(nbias))
        recax2.set_xticklabels(dates)
        recax2.xaxis.set_tick_params(rotation=90)
        plt.tight_layout() 
        plt.show()
        #sys.exit()
    
    res_bias /= nbias
    res_bias2 /= nbias
    param2 = {k:v for k,v in interp_param.items()}
    param2.pop('date', None)
    param2['start'] = param['start']
    param2['end'] = param['end']
    p.imshow(res_bias, plot_param=albedo_diff, **param2, source='bias'+source_ref+source_vtk)
    p.imshow(np.sqrt(res_bias2-res_bias**2), plot_param=albedo_biasstd, **param2, source='biasSTD'+source_ref+source_vtk)

    
    ## Init empty data containers to store results
    # temporal domain (10x10 km area)
    # spatial domain (global bias map)

def print_stats(*arrays, label=None, fmt='{:.2f}'):
    res = []
    for ida,a in enumerate(arrays):
        if label is None:
            new = str(ida)
        else:
            new = label[ida]
        res.append(pd.DataFrame(a.ravel()).describe().applymap(fmt.format).rename(columns={0:new}))
    res = pd.concat(res, axis=1)
    print(res)

def main(param):

    ti = SimpleTimer()

    if 1:
        process_etal_series(**param)        
    
    ## Load cache file
    else:
        data = load_h5(**param)

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
            #'start' : ['2015-01-05'],
            'start' : ['2007-01-05'],
            #'end' : ['2017-12-25'],
            'end' : ['2009-12-25'],
            #'size' : [10],
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

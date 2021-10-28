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



#------------ CONFIG FILE TO GLOBAL --------------------
import yaml

## Read config yaml file
conf_fname = './config_cnrm.yml'
print('--- Read config file {conf_fname} ...')
with open(conf_fname, 'r') as stream:
    try:
        conf = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()


#------------ INIT AND HELPER CLASSES AND FUNCTIONS --------------------

def get_time_range(start, end, days=[], **param):
    dseries = pd.date_range(start, end, freq='D')
    
    if len(days)>0:
        #self.start_01 = self.start.replace(day=1)
        #self.end_31 = self.end.replace(day=pd.Period(self.end, freq='D').days_in_month) # freq='xxx' is here just to avoid a bug
        dseries = dseries[[d in days for d in dseries.day]]

    return dseries

def get_all_filenames(**param):
    
    path_dic = conf['data_paths']
    
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
            root = pathlib.Path(path_dic[prod]['root'])
            path_format = path_dic[prod]['format']
            
            #if prod=='modis': # Unpredictable part in MODIS file name, need to glob...
            if 'MODIS' in root.as_posix(): # Unpredictable part in MODIS file name, need to glob...
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

            else: # EPS and MSG file names are predictable.
                fpath = root/d.strftime(path_format)
                if fpath.exists():
                    print(f"{prod.upper()} {d:%Y%m%d} found.")
                    df.at[d, f'{prod}_path'] = fpath
                else:
                    print(f"{prod.upper()} {d:'%Y%m%d'} NOT found.")
    
    print(df)

    return df

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

    def imshow(self, data, plot_param={}, show_axis=True, dpi='figure', **param):
        self.type = 'imshow' 
        self.initialize()

        data = data[self.lat,self.lon]
        ims = self.ax.imshow(data, **plot_param)
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cb = plt.colorbar(ims, cax=cax)
        
        self.finalize(show_axis, dpi, **param)

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

    def finalize(self, show_axis=True, dpi='figure', **param):
        from textwrap import wrap
        param_str = '\n'.join(wrap(param_to_string(param, sep1=' | ', sep2=':'), 60))
        self.ax.set_title(param_str, fontsize=10 )
        plt.tight_layout()
        
        im_name = f"res_proj2vtk_{self.type}_{param_to_string(param, shorten=True)}.png"
        if show_axis:
            plt.savefig(im_name, bbox_inches='tight', dpi=dpi)
        else:
            self.save_no_whitespace(im_name, dpi)
        print(f"--- Output image saved to: {im_name}")
        plt.close(self.fig)
   
    def save_no_whitespace(self, filepath, dpi='figure'):
        '''Save the current image with no whitespace'''
        plt.subplots_adjust(0,0,1,1,0,0)
        if len(self.fig.axes)==2:
            plt.delaxes(self.fig.axes[1]) # axes[0] correspond to the main plot.
        self.ax.set_title("")
        self.ax.axis('off')
        self.ax.margins(0,0)
        self.ax.xaxis.set_major_locator(plt.NullLocator())
        self.ax.yaxis.set_major_locator(plt.NullLocator())
        self.fig.savefig(filepath, pad_inches=0, bbox_inches='tight', dpi=dpi)

#------------ SATELLITE CLASSES FOR DATA PROCESSING ----------------

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

        self.get_config()

    @property
    def full_shape(self):
        """Get global shape before slicing"""
        if self._full_shape is None:
            with h5py.File(self.ground_mask_conf['file'], 'r') as flw:
                self._full_shape = flw[self.ground_mask_conf['var']].shape
        return self._full_shape

    @full_shape.setter
    def full_shape(self, value):
        self._full_shape = value

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

    def interpolate_on(self, target, from_source, use_cache, cache_slicing=slice(None), **param):
        """If cache file exists, load it, otherwise interpolate"""
        
        print('Interpolation parameters:')
        for k,v in param.items():
            print(f'    {k}:{v}')

        interp_type = f'{self.product}2{type(target).__name__}'

        interpolation_required = ~use_cache
        if use_cache:
            ## Cache files are  stored in working directory for now
            h5_path = pathlib.Path(f"./cache/cache_{interp_type}_{param_to_string(param)}.h5")
            data = self.load_h5_cache(h5_path, param['var'], cache_slicing)
            if data is not None:
                self.data_interp = data
                return data
            else:
                interpolation_required = True

        if interpolation_required:

            self.get_data(from_source=from_source)

            ti = SimpleTimer()

            ## Interpolation
            print('### Interpolation')
            interp = vtk_interpolation(**param) 
            ti('interp init')
            interp.set_source(self.lon, self.lat, {self.var:self.data})
            ti('interp set_source')
            interp.set_target(target.lon, target.lat)
            ti('interp set_destination')
            print(f'--- Interpolate {self.lon.size} source pts on {target.lon.size} target pts')
            interp.run()
            ti('interp run')
            interp_var = interp.get_output()
            ti('interp get_output')

            ## Extract and save data
            data = np.full(target.shape, -1.)
            ti('interp np_full')
            data[target.mask] = interp_var[param['var']]
            ti('interp fill_with_mask')
            self.export_to_h5(data, interp_type, **param)
            ti('interp export_to_h5')
            
            self.data_interp = data
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
            print(f"--- Read h5 from: {h5_path} ...")
            with h5py.File(h5_path, 'r') as h5_file:
                data = h5_file[var][slicing]
            return data
        else:
            print(f"--- File: {h5_path} not found.")
            return None
            
    def describe(self):
        print(self.data.shape)
        print(self.lon.shape)

        print_stats(self.data)
        print(np.count_nonzero(self.data==-1))

    def whoami(self):
        return type(self).__name__

    def get_config(self):
        for key in conf[self.whoami()]:
            setattr(self, key, conf[self.whoami()][key])

class MODIS(SatelliteTools):
    def __init__(self, product, var, slicing=slice(None)):
        
        self.ti = SimpleTimer()

        super().__init__(product, var, slicing)

        self.data_scaling = 0.001
        self.data = None

    def get_data(self, from_source, mask_value=None):
        """
        Read a data file and compute a variable mask
        
        :mask_value: None: return a flatten array with only valid data
                     <number> : return an array with the shape of the read data and <number> used as non valid data.
        """
        self.ti()
        print(f"--- Read {self.product.upper()}/{self.var} in {from_source['data']} ...")
        h4f = SD(from_source['data'].as_posix(), SDC.READ)
        self.full_shape = (h4f.select(self.var).dimensions()['YDim:Grid_Parameter'], h4f.select(self.var).dimensions()['XDim:Grid_Parameter'])
        self.data0 = h4f.select(self.var)[self.slicing] 
        self.data0[self.data0==32767] = -1
        self.mask = self.ground_mask & (self.data0!=-1)
        if mask_value is None:
            self.data = self.data0[self.mask]
        else:
            self.data = np.zeros(self.shape) + mask_value
            self.data[self.mask] = self.data0[self.mask]
        self.ti('get_data')

class EPS(SatelliteTools):
    def __init__(self, product, var, slicing=slice(None)):

        self.ti = SimpleTimer()
        
        super().__init__(product, var, slicing)

        self.data_scaling = 0.0001
        self.data = None

    def get_data(self, from_source, mask_value=None):
        """
        Read a data file and compute a variable mask
        
        :mask_value: None: return a flatten array with only valid data
                     <number> : return an array with the shape of the read data and <number> used as non valid data.
        """
        self.ti()
        print(f"--- Read {self.product.upper()}/{self.var} in {from_source['data']} ...")
        with h5py.File(from_source['data'],'r') as h5f:
            self.data0 = h5f[self.var][self.slicing]
        self.mask = self.ground_mask & (self.data0!=-1)
        if mask_value is None:
            self.data = self.data0[self.mask]
        else:
            self.data = np.zeros(self.shape) + mask_value
            self.data[self.mask] = self.data0[self.mask]
        self.ti('get_data')
        return self.data

class MSG(SatelliteTools):
    def __init__(self, product, var, slicing=slice(None)):

        self.ti = SimpleTimer()

        super().__init__(product, var, slicing)

        self.data_scaling = 0.0001
        self.data = None
        
    def get_data(self, from_source, mask_value=None):
        """
        Read a data file and compute a variable mask
        
        :mask_value: None: return a flatten array with only valid data
                     <number> : return an array with the shape of the read data and <number> used as non valid data.
        """
        self.ti()
        print(f"--- Read {self.product.upper()}/{self.var} in {from_source['data']} ...")
        with h5py.File(from_source['data'],'r') as h5f:
            self.data0 = h5f[self.var][self.slicing]
        self.mask = self.ground_mask & (self.data0!=-1)
        if mask_value is None:
            self.data = self.data0[self.mask]
        else:
            self.data = np.zeros(self.shape) + mask_value
            self.data[self.mask] = self.data0[self.mask]
        self.ti('get_data')
        return self.data


#------------ MAIN FUNCTIONS -------------

def process_etal_series(**param):

    df = get_all_filenames(**param)

    interp_param = conf['interp_param'] 

    ## Data accumulation

    ## DEBUG: subsampling
    if 1:
        slicing_msg = (slice(None, None, 2), slice(None, None, 2))
        #slicing_msg = (slice(None), slice(None))
        #slicing_eps = (slice(None, None, 50), slice(None, None, 50))
        slicing_eps = (slice(850, 17150, 10), slice(9000, 27000, 10)) # Remove points not on MSG disc
        slicing_modis = (slice(None, None, 50), slice(None, None, 50))
        #slicing = (slice(380, 480), slice(1950, 2050)) # France
        #slicing = (slice(50, 700), slice(1550, 3250)) # Euro
        #slicing = (slice(700, 1850), slice(1240, 3450)) # NAfr
        #slicing = (slice(1850, 3040), slice(2140, 3350)) # SAfr
        #slicing = (slice(1460, 2970), slice(40, 740)) # SAme
        interp_param['slicing'] = slicing_eps
        interp_param['string_exclude'] = ['slicing']

    ## Create product objects
    mtalr = MSG(product='mtalr', var='AL-BB-BH', slicing=slicing_msg)   
    etalr = EPS(product='etalr', var='AL-BB-BH', slicing=slicing_eps)   
    # bsa-sw = black sky albedo shortwave
    modis = MODIS(product='bsa-sw', var='BRDF_Albedo_BSA_Shortwave', slicing=slicing_modis)   

    #compare_two(new=etalr, ref=modis, grid='new', df_paths=df, **interp_param)
    compare_two(new=etalr, ref=mtalr, grid='ref', df_paths=df, **interp_param)


def compare_two(new, ref, grid, df_paths, **interp_param):

    if grid=='new':
        source = ref
        target = new
    elif grid=='ref':
        source = new
        target = ref
    else:
        print("--- ERROR: wrong grid parameter ('new' or 'ref')")
        sys.exit()

    ## Create a path organizer 
    path_list = {'source_paths':{}, 'target_paths':{}}
    path_list['source_col'] = [col for col in [*df_paths] if col.startswith(source.product)]
    path_list['target_col'] = [col for col in [*df_paths] if col.startswith(target.product)]
    for obj in ['source', 'target']:
        for c in path_list[f'{obj}_col']:
            if ':' in c.rstrip('_path'):
                path_list[f'{obj}_paths'][c.rstrip('_path').split(':')[-1]] = c
            else:
                path_list[f'{obj}_paths']['data'] = c

    ## Create cache di if not already exists
    pathlib.Path('./cache').mkdir(parents=True, exist_ok=True)

    ## Init statistical arrays
    res_bias = np.zeros(target.shape, dtype=float)
    res_bias2 = np.zeros(target.shape, dtype=float)
    nbias = 0

    ## Plot params
    albedo = {'cmap':'jet', 'vmin':0, 'vmax':0.6}
    albedo_modis = {'cmap':'jet'}
    blank = {'cmap':'binary', 'vmin':0, 'vmax':0.6}
    albedo_diff = {'cmap':'seismic', 'vmin':-0.15, 'vmax':0.15}
    albedo_biasstd = {'cmap':'jet', 'vmin':0.0, 'vmax':0.08}
    source_name = 'Vtk' + source.product.upper()
    target_name = target.product.upper()
    p = Plots(zoom=0)

    ## Setup recursive plot
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

    for i, (index,row) in enumerate(df_paths.iterrows()):
        print(f'\n================ {index:%Y-%m-%d} ================')

        params = {}
        params['id'] = f'{i:03d}'
        params['date'] = f'{index:%Y%m%d}'

        if (row[f'{source.product}_path'] is None) or (row[f'{target.product}_path'] is None):
            print('--- Data file(s) not found for this date.')
            tglob(params['id'])
            continue

        ## DEBUG ON: show input data (disable to avoid loading input data if cache file already exists)
        if 0:
            if row['modis_path'] is None:
            #if row['etalr_path'] is None:
            #if row['etal_path'] is None:
                #p.imshow(np.zeros(etalr.shape), plot_param=blank, show_axis=True, **params, source='None')
                #p.imshow(np.zeros(etal.shape), plot_param=blank, show_axis=True, **params, source='None')
                p.imshow(np.zeros(modis.shape), plot_param=blank, show_axis=True, **params, source='None')
                continue
            
            #etalr.get_data(row['etalr_path'], mask_value=-1)
            #p.imshow(etalr.data*0.0001, plot_param=albedo, show_axis=True, **params, source='etalr')

            modis.get_data(row['modis_path'], mask_value=-1)
            p.imshow(modis.data*0.001, plot_param=albedo, show_axis=True, **params, source='modis')

            #etal.get_data(row['etal_path'], mask_value=-1)
            #p.imshow(etal.data*0.0001, plot_param=albedo, show_axis=True, **params, source='etal')

        source_paths = {k:row[v] for k,v in path_list['source_paths'].items()}
        target_paths = {k:row[v] for k,v in path_list['target_paths'].items()}
        
        source_on_target_grid = source.interpolate_on(target, from_source=source_paths,
                                                    use_cache=False, cache_slicing=target.slicing,
                                                    date=params['date'], **interp_param)
        target_on_target_grid = target.get_data(from_source=target_paths, mask_value=-1)

        print_stats(source_on_target_grid, target_on_target_grid, label=['source.product', 'target.product'])

        ## Filtering
        source_on_target_grid[target_on_target_grid==-1] = -1
        target_on_target_grid[source_on_target_grid==-1] = -1
        source_on_target_grid = source_on_target_grid*source.data_scaling
        target_on_target_grid = target_on_target_grid*target.data_scaling
        # Note: var *= 0.001 does not work if not the same type (ex var has int and 0.001 is float)

        print(f'{source.product}:',np.count_nonzero(source_on_target_grid<0),
              f'{target.product}:', np.count_nonzero(target_on_target_grid<0))

        if 0:
            pass
            #p.imshow(etalr_on_msg*0.0001, plot_param=albedo, show_axis=True, **params, **interp_param, source=source_vtk)
            #p.imshow(modis_on_eps*0.001, plot_param=albedo, show_axis=True, dpi=200, **params, **interp_param, source='modis-vtk')


        ## Compute bias
        bias = source_on_target_grid - target_on_target_grid
        print(bias.shape)
        print(res_bias.shape)
        res_bias += bias
        res_bias2 += bias*bias
        nbias += 1

        ## Plot global maps
        if 1:
            p.imshow(target_on_target_grid, plot_param=albedo, show_axis=True, **params, source=target_name)
            p.imshow(source_on_target_grid, plot_param=albedo, show_axis=True, **params, **interp_param, source=source_name)
            p.imshow(bias, plot_param=albedo_diff, show_axis=True, **params, **interp_param, source='diff'+source_name+target_name)

        tglob(params['id'])

        ## Fill recursive plot
        if rec_plot:
            #recax.violinplot([bias.ravel()], positions=[nbias], showmeans=False, showmedians=True, showextrema=True, widths=0.9)
            
            mask_015 = source_on_target_grid.ravel()>0.15
            bias_above = 100*bias.ravel()[mask_015]/source_on_target_grid.ravel()[mask_015]
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

        if (i==4) & 1:
            print('--- BREAK LOOP')
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
        #plt.show()
        plt.savefig('res_stats.png')
    
    res_bias /= nbias
    res_bias2 /= nbias
    param2 = {k:v for k,v in interp_param.items()}
    param2.pop('date', None)
    param2['start'] = param['start']
    param2['end'] = param['end']
    p.imshow(res_bias, plot_param=albedo_diff, **param2, source='bias'+source_name+target_name)
    p.imshow(np.sqrt(np.abs(res_bias2-res_bias**2)), plot_param=albedo_biasstd, **param2, source='biasSTD'+source_name+target_name)

    tglob.show()

def main(param):
    ti = SimpleTimer()

    process_etal_series(**param)        
    
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
    
    for param in combine_param(conf['global_param']):
        for k,v in param.items():
            print(k, ':', v)
        try:
            main(param)
        except Exception as e:
            print('--- ERROR IN PROCESSING ---')
            print(traceback.format_exc())
            #logging.error(traceback.format_exc())

        print('---------------')

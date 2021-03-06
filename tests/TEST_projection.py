import numpy as np
import xarray as xr
#import dask
#import dask.dataframe as dd
#import dask.array as da
import sys
import h5py
import yaml
from tools import SimpleTimer

import matplotlib.pyplot as plt

import cartopy.crs as ccrs



class ETALTools:

    def __init__(self, yfile):
        self.var = 'AL-BB-DH'
        #self.var = 'Q-Flag'
        self.read_yaml_config_file(yfile)

    def load_nc_file(self, h5file):
        ds = xr.open_dataset(h5file) 

    def load_h5_file(self, h5file):
        sub = 10
        self.hf = h5py.File(h5file, 'r')[self.var][::sub,::sub]
        print(self.hf)
        print(self.hf.shape, self.hf.itemsize)

    def create_lon_lat(self):
        ''' Create lon/lat arrays for input sinusoidal data (same shape)'''
        s = self.hf.shape
        lon = np.linspace(-180, 180, s[1])
        lat = np.linspace(90, -90, s[0])
        self.lon, self.lat = np.meshgrid(lon,lat)
        self.lon = self.lon/np.cos(np.deg2rad(self.lat))
        self.mask = np.abs(self.lon)>180.
        self.lon[self.mask] = np.nan
        self.lat[self.mask] = np.nan
        print(self.lon)
        print(self.lat)

    def load_nc_file_list(self, h5list):
        ds = xr.open_mfdataset(h5list, parallel=True, chunks='auto') 
        
    def plot_sinusoidal_proj(self, array):
        ''' Plot data using cartopy processing'''
        if 1:
            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.Sinusoidal())
        elif 0:
            fig = plt.figure(figsize=(10, 10))
            #ax = fig.add_subplot(1, 1, 1, projection=ccrs.Geostationary())
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    
        # make the map global rather than have it zoom in to
        # the extents of any plotted data
        ax.set_global()
    
        #ax.stock_img()
        ax.coastlines()
        
        sinus_transfo = ccrs.Sinusoidal()
        # Get the x,y extent of the transfo to fit the image in it
        img_extent = list(sinus_transfo.x_limits) + list(sinus_transfo.y_limits)
        ax.imshow(array, origin='upper', extent=img_extent, transform=sinus_transfo)
    
        plt.tight_layout()
    
        imname = 'res_global_map_{}.png'.format(self.var)
        imname = 'res_global_map_Sinusoidal_{}.png'.format(self.var)
        plt.savefig(imname, dpi=100)
        print('--- Global map saved to: {}'.format(imname))

    def plot_2D_array(self, array):
        ''' Simply plot a 2D array (useful for debug)'''
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
   
        ax.imshow(array)

        # Write the value in each pixel
        #for (j,i),label in np.ndenumerate(array):
        for (j,i),label in np.ndenumerate(self.hf):
            ax.text(i, j, label, ha='center', va='center', fontsize=6)

        imname = 'res_2darray_{}.png'.format(self.var)
        plt.savefig(imname)
        print('--- Array image saved to: {}'.format(imname))

    def read_yaml_config_file(self, yfile):
        with open(yfile, 'r') as stream:
            try:
                self.config = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
                sys.exit()
        print('--- Read config file OK.')
 
    def transform_with_pyproj(self):
        from  pyproj import CRS, Transformer
        
        crs0 = CRS.from_proj4("+proj=eqc +units=km")
        crs3 = CRS.from_proj4("+proj=lonlat +units=km")
        crs1 = CRS.from_proj4("+proj=sinu +units=km")
        crs2 = CRS.from_proj4("+proj=geos +h=35774290 +units=km")
        
        t0 = Transformer.from_crs(crs0, crs3)
        t1 = Transformer.from_crs(crs3, crs0)
        
        x = np.array([-20000, 0, 20000])
        y = np.array([-10000, 0, 10000])
        
        x1 = np.array([-180, 0, 180])
        y1 = np.array([-90, 0, 90])
        
        res = t0.transform(x,y)
        res = t1.transform(x1,y1)
        
        print(res)

    def process(self):
        if 0:
            h5list = '/cnrm/vegeo/juncud/NO_SAVE/ETAL/HDF5_LSASAF*'
            self.load_nc_file_list(h5list)

        if 1:
            ti = SimpleTimer()
            #h5file = '/cnrm/vegeo/juncud/NO_SAVE/ETAL/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_202012250000'
            h5file = self.config['input_path']['single_h5_file']
            self.load_h5_file(h5file)
            self.create_lon_lat()
            ti('load')
            #self.plot_2D_array(np.where(self.mask, np.nan, self.hf))
            ti('imshow')
            self.plot_sinusoidal_proj(self.hf)
            ti('project')
            ti.show()

        if 0:
            self.plot_sinusoidal_proj()


if __name__=='__main__':
     
    et = ETALTools('config_cnrm.yml')
    et.process()



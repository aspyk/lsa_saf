import h5py
import numpy as np
import vtk
from vtk.util import numpy_support as ns

import os,sys
from tools import SimpleTimer

import matplotlib.pyplot as plt

from datetime import datetime as dtime

import traceback
import itertools

from vtk_tools import vtk_interpolation


def numpy_to_vtkpoints(pts_array):
    '''
    pts_array must be a (N,3) numpy array.
    Return a vtkPoints
    '''
    points = vtk.vtkPoints()
    if 0:
        points.SetNumberOfPoints(pts_array.shape[0])
        for ip,pt in enumerate(pts_array):
            #points.InsertNextPoint(pt)
            points.SetPoint(ip,pt)
        points.Modified()
    else:
        points.SetData(ns.numpy_to_vtk(pts_array, deep=True))
    return points 

def array_to_vtk(x, y, z, var):
    n_pts_dim = 6


    if n_pts_dim==2:
        pass
    elif n_pts_dim==3:
        pass
    elif n_pts_dim==4:
        pass
    elif n_pts_dim==5:
        pass
    elif n_pts_dim==6:
        n_pts_cell = 21
        to_vtk_order = np.array([0,5,20,1,2,3,4,10,14,17,19,18,15,11,6,7,9,16,8,13,12])

    np_pts = np.array((x,y,z)).T
    vtk_pts = numpy_to_vtkpoints(np_pts) 

    # Create vtk point data
    var = ns.numpy_to_vtk(num_array=var.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    var.SetNumberOfComponents(1)
    var.SetName('rho')

    # Create vtkCellArray
    n_cell = int(len(x)/n_pts_cell)
    connectivity = np.array([[n_pts_cell]+list(to_vtk_order+i*n_pts_cell) for i in range(n_cell)]).ravel()
    vtk_connectivity = ns.numpy_to_vtk(num_array=connectivity, deep=True, array_type=vtk.VTK_ID_TYPE)
    cellarray = vtk.vtkCellArray()
    cellarray.SetCells(n_cell, vtk_connectivity)

    # Merge everything into a vtkUnstructuredGrid
    lag_tri = vtk.vtkUnstructuredGrid()
    lag_tri.SetPoints(vtk_pts)
    lag_tri.GetPointData().AddArray(var)
    lag_tri.SetCells(vtk.VTK_LAGRANGE_TRIANGLE, cellarray)

    # Write the mesh
    if 0:
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName("test.vtu")
        writer.SetInputData(lag_tri)
        #writer->SetDataModeToAscii() #Optional for debug - set the mode. The default is binary.
        writer.Write()

def array_to_vtk_scatter(x, y, z, dic_var=None, fname='res_test'):
    x = x.ravel()
    y = y.ravel()
    z = z.ravel()

    np_pts = np.array((x,y,z)).T

    t1 = SimpleTimer()
    vtk_pts = numpy_to_vtkpoints(np_pts) 
    t1('vtkpoints')
    t1.show()

    # Create vtk point data
    vars = []
    if dic_var is not None:
        for vname, vdata in dic_var.items():
            vars.append(ns.numpy_to_vtk(num_array=vdata.ravel(), deep=True, array_type=vtk.VTK_FLOAT))
            vars[-1].SetNumberOfComponents(1)
            vars[-1].SetName(vname)

    # Merge everything into a vtkUnstructuredGrid
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(vtk_pts)
    for var in vars:
        grid.GetPointData().AddArray(var)

    if 0:
        # Write the mesh
        writer = vtk.vtkXMLUnstructuredGridWriter()
        fname = "{}.vtu".format(fname)
        writer.SetFileName(fname)
        writer.SetInputData(grid)
        #writer->SetDataModeToAscii() #Optional for debug - set the mode. The default is binary.
        writer.Write()
        print("--- VTK file saved to: {}".format(fname))

    return grid

def lonlat_to_xyz(lon, lat, R=1.0, unit='rad'):
    if unit=='deg':
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)
    return x,y,z

def msg_to_vtk(stride=10, lonlat_only=False, **param):
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

    if lonlat_only:
        return lon, lat, valid_mask, out_shape

    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    x,y,z = lonlat_to_xyz(lon, lat, unit='rad')

    print(x.shape)

    msg_grid = array_to_vtk_scatter(x, y, z, fname='res_msg')
    return msg_grid, valid_mask, out_shape

def etal_to_vtk(stride=100, lonlat_only=False, var=None, date=None, **param):

    t0 = SimpleTimer()

    if date is None:
        date_ = '2020-DEC-05'
        dateObj = dtime.strptime(date_, '%Y-%b-%d')
        dateForFile = dtime.strftime(dateObj, '%Y%m%d')
    else:
        dateForFile = date
    
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
        ## 500x500 Alps
        if 1:
            iLatStartEps = 4200
            iLatEndEps   = 4700 
            iLonStartEps = 18150 
            iLonEndEps   = 18650 
        ## 100x100 Alps
        if 0:
            iLatStartEps = 4350
            iLatEndEps   = 4450 
            iLonStartEps = 18430 
            iLonEndEps   = 18530 
    
    etal_file = f'/cnrm/vegeo/juncud/NO_SAVE/ETAL/2020/12/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_{dateForFile}0000'
    #etal_file = '/cnrm/vegeo/juncud/NO_SAVE/ETAL/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_202012050000'

    print(f'--- Read {etal_file}...')
    with h5py.File(etal_file,'r') as h5f:
        etal_shape  = h5f[var].shape
        etal = h5f[var][iLatStartEps:iLatEndEps:stride,iLonStartEps:iLonEndEps:stride]
        qflag = h5f['Q-Flag'][iLatStartEps:iLatEndEps:stride,iLonStartEps:iLonEndEps:stride]
    #etal = h5py.File('.\\input_data\\HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_202012250000', 'r')[var]

    ## DEBUG: set checkerboard
    if 0:
        etal = 6000*(np.indices(etal.shape).sum(axis=0) % 2)
        print(etal)

    ## DEBUG: plot input etal to image
    if 0:
        albedo = {'cmap':'jet', 'vmin':0, 'vmax':6000}
        qflag_p = {'cmap':'jet'}
        p = Plots(zoom=0)
        qflag = (qflag & 3)==1
        etal[qflag==0] = -1
        p.imshow(etal, plot_param=albedo, **param, source='ETAL')
        #p.imshow(qflag, plot_param=qflag_p, **param, source='ETAL')
        sys.exit()

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

    if lonlat_only:
        return lon, lat, dic_var

    x,y,z = lonlat_to_xyz(lon, lat, unit='deg')

    t0('to_xyz')

    print(x.shape)

    etal_grid = array_to_vtk_scatter(x, y, z, dic_var, fname='res_etal')
    
    t0('to_vtk')

    t0.show()

    return etal_grid

def etal_to_msg(source, target, kernel='mean', radius=0.001, **param):

    radius /= 6371 # radius is given in km while interpolation is done on unit sphere

    ## Setup the locator
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(source)
    locator.BuildLocator()

    if kernel=='mean':
        vtk_kernel = vtk.vtkLinearKernel()
    elif kernel=='gaussian':
        vtk_kernel = vtk.vtkGaussianKernel()
        vtk_kernel.SetSharpness(4)
    elif kernel=='inverse_distance':
        vtk_kernel = vtk.vtkShepardKernel()
        #vtk_kernel.SetPowerParameter(2) # Default power is 2

    vtk_kernel.SetRadius(radius)


    ## Setup interpolator
    interpolator = vtk.vtkPointInterpolator()
    interpolator.SetInputData(target)
    interpolator.SetSourceData(source)
    interpolator.SetKernel(vtk_kernel)
    interpolator.SetLocator(locator)
    if 1:
        interpolator.SetNullPointsStrategyToClosestPoint()
    else:
        interpolator.SetNullPointsStrategyToNullValue()
        interpolator.SetNullValue(-1.)

    ## Run the interpolation
    interpolator.Update()

    # Write the result
    if 0:
        writer = vtk.vtkXMLUnstructuredGridWriter()
        fname = 'res_msg_interpolated'
        fname = "{}.vtu".format(fname)
        writer.SetFileName(fname)
        writer.SetInputData(interpolator.GetOutput())
        #writer->SetDataModeToAscii() #Optional for debug - set the mode. The default is binary.
        writer.Write()
        print("--- VTK file saved to: {}".format(fname))

    return interpolator

def msg_var_vtk_to_numpy(interp, valid_mask, msg_data_shape):
    for a in range(interp.GetOutput().GetPointData().GetNumberOfArrays()):
        vname = interp.GetOutput().GetPointData().GetArrayName(a)
        numpy_var = ns.vtk_to_numpy(interp.GetOutput().GetPointData().GetArray(vname))

        eps_on_msg = np.zeros(msg_data_shape)-1.
        print(eps_on_msg.shape)
        eps_on_msg[valid_mask] = numpy_var
        yield vname,eps_on_msg

def params_to_string(param, sep1='_', sep2='-'):
    """
    Normalize parameter dict

    Keep only alphanumeric char in parameters and create a unique long string with custom separators.
    Separators are use like this:
    {p1:v1, p2:v2} -> p1<sep2>v1<sep1>p2<sep2>v2
    """
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
        param_str = '\n'.join(wrap(params_to_string(param, sep1=' | ', sep2=':'), 60))
        self.ax.set_title(param_str, fontsize=10 )
        #self.ax.set_title(params_to_string(param, sep1=' | ', sep2=':'), fontsize=9 )
        
        #self.fig.subplots_adjust(top=0.9)
        plt.tight_layout()
        
        param_str = params_to_string(param)
        im_name = f"res_proj2vtk_{self.type}_{params_to_string(param)}.png"
        plt.savefig(im_name, dpi=200)
        print(f"--- Output image saved to: {im_name}")
    

def export_to_h5(data, **param):
    h5_name = f"res_proj2vtk_output_{params_to_string(param)}.h5"
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

def load_h5(**param):
    h5_name = f"res_proj2vtk_output_{params_to_string(param)}.h5"
    with h5py.File(h5_name, 'r') as h5_file:
        data = h5_file[param['var']][:]
    print(f"--- Read h5 from: {h5_name}")
    return data

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
        msg_daniel = h5ref[var][:]

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

    p = Plots(zoom=1)

    p.imshow(msg_interp,            plot_param=albedo,      **param,                              source=source_vtk)
    p.imshow(msg_ref,               plot_param=albedo,      var=param['var'], date=param['date'], source=source_ref)
    p.imshow(msg_daniel,            plot_param=albedo,      var=param['var'], date=param['date'], source=source_daniel)

    p.imshow(msg_ref-msg_interp,    plot_param=albedo_diff, **param,                              source='diff'+source_ref+source_vtk)
    p.imshow(msg_ref-msg_daniel,    plot_param=albedo_diff, var=param['var'], date=param['date'], source='diff'+source_ref+source_daniel)
    p.imshow(msg_daniel-msg_interp, plot_param=albedo_diff, **param,                              source='diff'+source_daniel+source_vtk)

    p.scatter(msg_ref, msg_interp,    plot_param=albedo_scatter, **param, s1=source_ref, s2=source_vtk)
    p.scatter(msg_ref, msg_daniel,    plot_param=albedo_scatter, **param, s1=source_ref, s2=source_daniel)
    p.scatter(msg_daniel, msg_interp, plot_param=albedo_scatter, **param, s1=source_daniel, s2=source_vtk)

def DEPRECATED_main(param):

    ti = SimpleTimer()

    ## Perform interpolation
    if 0:
        print('### MSG extraction')
        msg_grid, valid_mask, msg_shape = msg_to_vtk(stride=1)
        ti('MSG')
        
        print('### ETAL extraction')
        etal_grid = etal_to_vtk(stride=1, **param)
        ti('ETAL')
        

        print('### Interpolation')
        kernel = 'mean'
        #kernel = 'inverse_distance'
        #kernel = 'gaussian'
        radius = 3
        #interpolation = etal_to_msg(etal_grid, msg_grid, kernel=kernel, radius=radius)
        interpolation = etal_to_msg(etal_grid, msg_grid, **param)
        ti('interpolation')
        
        print('### Export interpolated data')
        for name,data in msg_var_vtk_to_numpy(interpolation, valid_mask, msg_shape):
            ## Debug
            if 0:
                var_dic[param['var']] = 6000*(np.indices(var_dic[param['var']].shape).sum(axis=0) % 2)
            export_to_image(data, **param)
            export_to_h5(data, **param)
            ti('export')    
    
    ## Load cache file
    else:
        data = load_h5(**param)

    print('### Compare results')
    compare_with_real_msg(data, **param)
    ti('comparison')


    ti.show()


def main_vtk(param):
    
    ti = SimpleTimer()

    if 1:

        print('### MSG extraction (target)')
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
            #'date' : ['20201205','20201215','20201225'],
            'date' : ['20201225'],
            #'kernel' : ['inverse_distance','gaussian'],
            #'kernel' : ['mean','inverse_distance','gaussian'],
            'kernel' : ['mean','inverse_distance'],
            #'radius' : [5,10],
            'radius' : [3,5,10],
            'null_points' : ['closest'],
            #'null_points' : [-1.],
           }

    for param in combine_param(test):
        for k,v in param.items():
            print(k, ':', v)
        try:
            #interp = main(param)
            interp = main_vtk(param)
            pass
        except Exception as e:
            print('--- ERROR IN PROCESSING ---')
            print(traceback.format_exc())
            #logging.error(traceback.format_exc())

        #sys.exit()
        print('---------------')

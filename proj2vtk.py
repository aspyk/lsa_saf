import h5py
import numpy as np
import vtk
from vtk.util import numpy_support as ns

import os,sys
from tools import SimpleTimer

import matplotlib.pyplot as plt

from datetime import datetime as dtime


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

def msg_to_vtk(stride=10, **param):
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

    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    x,y,z = lonlat_to_xyz(lon, lat, unit='rad')

    print(x.shape)

    msg_grid = array_to_vtk_scatter(x, y, z, fname='res_msg')
    return msg_grid, valid_mask, out_shape

def etal_to_vtk(stride=100, var=None, date=None, **param):

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
        iLatStartEps = 5900
        iLatEndEps   = 6900 
        iLonStartEps = 18900 
        iLonEndEps   = 19900 
    
    etal_file = f'/cnrm/vegeo/juncud/NO_SAVE/ETAL/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_{dateForFile}0000'
    #etal_file = '/cnrm/vegeo/juncud/NO_SAVE/ETAL/HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_202012050000'

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

    ## Create lonlat coords for etal sinusoidal projection
    lon = np.linspace(-180, 180, etal_shape[1])
    lat = np.linspace(90, -90, etal_shape[0])

    t0('lon_lat')

    lat = lat[iLatStartEps:iLatEndEps:stride]
    lon = lon[iLonStartEps:iLonEndEps:stride]

    t0('lon_lat_slice')

    lon, lat = np.meshgrid(lon,lat)
    lon = lon/np.cos(np.deg2rad(lat))

    ## Version without np.where is a bit quicker, you can test it with the following if:
    if 1:
        mask_nonvalid = etal!=-1. # Discard points without valid albedo value (outside projection and in the ocean)
        t0('mask_direct')
    else:
        mask_nonvalid = np.where(etal!=-1.) # Discard points without valid albedo value (outside projection and in the ocean)
        t0('mask_where')
    s0 = lon.size
    lon = lon[mask_nonvalid] 
    lat = lat[mask_nonvalid]
    etal = etal[mask_nonvalid]
    s1= lon.size
    print("--- {:.2f} % data masked".format(100*(s0-s1)/s0))
    print('min/max lat:', lat.min(), lat.max())
    print('min/max lon:', lon.min(), lon.max())

    t0('mask_nonvalid')

    mask_fov = np.logical_and(np.abs(lon)<81, np.abs(lat)<81)  
    lon = lon[mask_fov] 
    lat = lat[mask_fov]
    etal = etal[mask_fov]
    s1= lon.size
    print("--- {:.2f} % data masked".format(100*(s0-s1)/s0))
    print('min/max lat:', lat.min(), lat.max())
    print('min/max lon:', lon.min(), lon.max())

    t0('mask_fov')

    x,y,z = lonlat_to_xyz(lon, lat, unit='deg')

    t0('to_xyz')

    print(x.shape)

    dic_var = {var:etal}
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

def export_to_image(data, plot_param={}, **param):
    plt.clf()
    #ims = plt.imshow(data, cmap='jet', vmin=0, vmax=6000)
    ims = plt.imshow(data, **plot_param)
    plt.tight_layout()
    plt.colorbar(ims)
    plt.title(params_to_string(param, sep1=' | ', sep2=':'), fontsize=10)
    param_str = params_to_string(param)
    im_name = f"res_proj2vtk_output_{params_to_string(param)}.png"
    plt.savefig(im_name, dpi=200)
    print(f"--- Output image saved to: {im_name}")

def export_to_h5(data, **param):
    h5_name = f"res_proj2vtk_output_{params_to_string(param)}.h5"
    with h5py.File(h5_name, 'w') as h5_file:
        h5_file[param['var']] = data
    print(f"--- Output h5 saved to: {h5_name}")

def compare_with_real_msg(msg_interp, **param):
    date = param['date']
    var = param['var']
    msg_ref_file = f'/cnrm/vegeo/juncud/NO_SAVE/aod-test-2/v1/HDF5_LSASAF_MSG_ALBEDO_MSG-Disk_{date}0000'
    with h5py.File(msg_ref_file) as h5ref:
        msg_ref = h5ref[var][:]

    msg_daniel_file = f'/cnrm/vegeo/juncud/scripts/eps_geos_{date}0000'
    with h5py.File(msg_daniel_file) as h5ref:
        msg_daniel = 1e4*h5ref[var][:]

    albedo = {'cmap':'jet', 'vmin':0, 'vmax':6000}
    albedo_diff = {'cmap':'seismic', 'vmin':-2000, 'vmax':2000}

    export_to_image(msg_interp, source='interp', plot_param=albedo, **param)
    export_to_image(msg_ref, source='ref', plot_param=albedo, date=param['date'])
    export_to_image(msg_ref-msg_interp, source='diffRefVtk', plot_param=albedo_diff, **param)
    export_to_image(msg_ref-msg_daniel, source='diffRefDaniel', plot_param=albedo_diff, var=param['var'], date=param['date'])
    export_to_image(msg_daniel-msg_interp, source='diffDanielVtk', plot_param=albedo_diff, **param)


def main():

    ti = SimpleTimer()

    param = {
            #'var' : 'AL-BB-DH',
            'var' : 'AL-BB-BH',
            'date' : '20201215',
            'kernel' : 'mean',
            #'kernel' : 'inverse_distance',
            #'kernel' : 'gaussian',
            'radius' : 3,
            }

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

    print('### Compare results')
    compare_with_real_msg(data, **param)
    ti('comparison')


    ti.show()

if __name__=='__main__':
    
    interp = main()


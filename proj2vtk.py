import h5py
import numpy as np
import vtk
from vtk.util import numpy_support as ns

import os,sys
from tools import SimpleTimer

import matplotlib.pyplot as plt


def numpy_to_vtkpoints(pts_array):
    '''
    pts_array must be a (N,3) numpy array.
    Return a vtkPoints
    '''
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(pts_array.shape[0])
    for ip,pt in enumerate(pts_array):
        #points.InsertNextPoint(pt)
        points.SetPoint(ip,pt)
    points.Modified()
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

    vtk_pts = numpy_to_vtkpoints(np_pts) 

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

def msg_to_vtk(stride=8):
    lat = h5py.File('hdf5_lsasaf_msg_lat_msg-disk_4bytesprecision', 'r')['LAT']
    lon = h5py.File('hdf5_lsasaf_msg_lon_msg-disk_4bytesprecision', 'r')['LON']
    lwmask = h5py.File('hdf5_lsasaf_usgs-igbp_lwmask_msg-disk', 'r')['LWMASK']
    
    lat = lat[::stride,::stride]/10000.
    lon = lon[::stride,::stride]/10000.
    lwmask = lwmask[::stride,::stride]
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

def etal_to_vtk(stride=100):
    var = 'AL-BB-DH'
    etal = h5py.File('.\\input_data\\HDF5_LSASAF_M01-AVHR_ETAL_GLOBE_202012250000', 'r')[var]

    ## Create lonlat coords for etal sinusoidal projection
    lon = np.linspace(-180, 180, etal.shape[1])
    lat = np.linspace(90, -90, etal.shape[0])

    # already discarding part of ETAL data that is for sure not in the MSG disk
    iLatStartEps = 850
    iLatEndEps   = 17150 # 6501
    iLonStartEps = 9000 # 22000
    iLonEndEps   = 27000 # 25501
    etal = etal[iLatStartEps:iLatEndEps:stride,iLonStartEps:iLonEndEps:stride]
    lat = lat[iLatStartEps:iLatEndEps:stride]
    lon = lon[iLonStartEps:iLonEndEps:stride]

    if 1:
        plt.imshow(etal)
        plt.tight_layout()
        plt.savefig("res_proj2vtk_input.png")

    lon, lat = np.meshgrid(lon,lat)
    lon = lon/np.cos(np.deg2rad(lat))

    mask_nonvalid = etal==-1. # Discard points without valid albedo value (outside projection and in the ocean)
    s0 = lon.size
    lon = lon[~mask_nonvalid] 
    lat = lat[~mask_nonvalid]
    etal = etal[~mask_nonvalid]
    s1= lon.size
    print("--- {:.2f} % data masked".format(100*(s0-s1)/s0))
    print('min/max lat:', lat.min(), lat.max())
    print('min/max lon:', lon.min(), lon.max())

    mask_fov = np.logical_and(np.abs(lon)<81, np.abs(lat)<81)  
    lon = lon[mask_fov] 
    lat = lat[mask_fov]
    etal = etal[mask_fov]
    s1= lon.size
    print("--- {:.2f} % data masked".format(100*(s0-s1)/s0))
    print('min/max lat:', lat.min(), lat.max())
    print('min/max lon:', lon.min(), lon.max())

    x,y,z = lonlat_to_xyz(lon, lat, unit='deg')

    print(x.shape)

    dic_var = {var:etal}

    etal_grid = array_to_vtk_scatter(x, y, z, dic_var, fname='res_etal')
    return etal_grid

def etal_to_msg(source, target):
    # Reuse the locator
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(source)
    locator.BuildLocator()

    # Use a gaussian kernel
    gaussianKernel = vtk.vtkGaussianKernel()
    gaussianKernel.SetRadius(0.1)
    gaussianKernel.SetSharpness(4)

    # Use a linear kernel
    linearKernel = vtk.vtkLinearKernel()
    linearKernel.SetRadius(0.005)

    interpolator = vtk.vtkPointInterpolator()
    interpolator.SetInputData(target)
    interpolator.SetSourceData(source)
    # interpolator.SetKernel(gaussianKernel)
    interpolator.SetKernel(linearKernel)
    interpolator.SetLocator(locator)
    interpolator.SetNullPointsStrategyToClosestPoint()

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

def export_interpolated_msg(interp, valid_mask, msg_data_shape):
    for a in range(interp.GetOutput().GetPointData().GetNumberOfArrays()):
        vname = interp.GetOutput().GetPointData().GetArrayName(a)
        numpy_var = ns.vtk_to_numpy(interp.GetOutput().GetPointData().GetArray(vname))

        eps_on_msg = np.zeros(msg_data_shape)-1.
        print(eps_on_msg.shape)
        eps_on_msg[valid_mask] = 1000*numpy_var

        plt.clf()
        plt.imshow(eps_on_msg)
        plt.tight_layout()
        plt.savefig("res_proj2vtk_output.png")

def main():

    ti = SimpleTimer()

    print('### MSG extraction')
    msg_grid, valid_mask, msg_shape = msg_to_vtk(stride=2)
    ti('MSG')
    
    print('### ETAL extraction')
    etal_grid = etal_to_vtk(stride=20)
    ti('ETAL')
    
    print('### Interpolation')
    interp = etal_to_msg(etal_grid, msg_grid)
    ti('interpolation')
    
    print('### Export interpolated data')
    export_interpolated_msg(interp, valid_mask, msg_shape)
    ti('export')    

    ti.show()

    return interp

if __name__=='__main__':
    
    interp = main()


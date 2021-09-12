#import h5py
import numpy as np
import vtk
from vtk.util import numpy_support as ns

import os,sys
from tools import SimpleTimer


class vtk_interpolation:
    def __init__(self, kernel='mean', radius=1, null_points_strategy=-1, **param):
        """
        - radius: float in [km]
            Radius of the sphere where all the neighbor points are selected.
        - kernel choices are 'mean', 'gaussian' and 'inverse'.
            Weight computation method for the mean of the neighbor points:
            - 'mean': All the weights are 1.
            - 'gaussian': The weights are computed as: exp(-(s*r/R)^2) (R being the radius and s the sharpness, set by default to 2)
            - 'inverse': The weights are computed using 1/r^a. a=2 by default but it can be modified.
        - null_points_strategy: 'closest' or a float.
            Strategy to apply if no point are found in the neihborhood sphere:
            - 'closest': the value of the closest point outside the sphere is given.
            - float : the value of the float is given.
        
        See https://vtk.org/doc/nightly/html/classvtkGeneralizedKernel.html for detailed description.
        """
        self.kernel = kernel
        self.radius = radius
        #if null_points_strategy=='closest':
        #    self.null_points_strategy = null_points_strategy
        #else:
        #    try:
        #        self.null_points_strategy = float(null_points_strategy)

    def set_source(self, lon, lat, var={}, unit='deg', fname=None):
        if len(var.keys())==0:
            print('--- ERROR: no input variable has been given.')
            sys.exit()
        x,y,z = self.lonlat_to_xyz(lon, lat, unit)
        self.source = self.array_to_vtk_scatter(x, y, z, var, fname)

    def set_target(self, lon, lat, unit='deg', fname=None):
        x,y,z = self.lonlat_to_xyz(lon, lat, unit)
        self.target = self.array_to_vtk_scatter(x, y, z, None, fname)

    def lonlat_to_xyz(self, lon, lat, unit='deg'):
        if unit=='deg':
            lon = np.deg2rad(lon)
            lat = np.deg2rad(lat)
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)
        return x,y,z

    def array_to_vtk_scatter(self, x, y, z, dic_var=None, fname=None):
        x = x.ravel()
        y = y.ravel()
        z = z.ravel()
    
        np_pts = np.array((x,y,z)).T
    
        t1 = SimpleTimer()
        #vtk_pts = self.numpy_to_vtkpoints(np_pts) 
        vtk_pts = vtk.vtkPoints()
        vtk_pts.SetData(ns.numpy_to_vtk(np_pts, deep=True))
    
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
    
        if fname is not None:
            # Write the mesh
            writer = vtk.vtkXMLUnstructuredGridWriter()
            fname = "{}.vtu".format(fname)
            writer.SetFileName(fname)
            writer.SetInputData(grid)
            #writer->SetDataModeToAscii() #Optional for debug - set the mode. The default is binary.
            writer.Write()
            print("--- VTK file saved to: {}".format(fname))
    
        return grid

    def run(self):
    
        self.radius /= 6371 # radius is given in km while interpolation is done on unit sphere
    
        ## Setup the locator
        locator = vtk.vtkStaticPointLocator()
        locator.SetDataSet(self.source)
        locator.BuildLocator()
    
        if self.kernel=='mean':
            vtk_kernel = vtk.vtkLinearKernel()
        elif self.kernel=='gaussian':
            vtk_kernel = vtk.vtkGaussianKernel()
            #vtk_kernel.SetSharpness(2) # Default is 2. As the sharpness increases the effects of distant points are reduced.
        elif self.kernel=='inverse_distance':
            vtk_kernel = vtk.vtkShepardKernel()
            #vtk_kernel.SetPowerParameter(2) # Default power is 2
    
        vtk_kernel.SetRadius(self.radius)
    
    
        ## Setup interpolator
        interpolator = vtk.vtkPointInterpolator()
        interpolator.SetInputData(self.target)
        interpolator.SetSourceData(self.source)
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
    
        #return interpolator
        self.interpolator = interpolator

    def get_output(self):
        out_var = {}
        for a in range(self.interpolator.GetOutput().GetPointData().GetNumberOfArrays()):
            vname = self.interpolator.GetOutput().GetPointData().GetArrayName(a)
            numpy_var = ns.vtk_to_numpy(self.interpolator.GetOutput().GetPointData().GetArray(vname))
            out_var[vname] = numpy_var
        return out_var


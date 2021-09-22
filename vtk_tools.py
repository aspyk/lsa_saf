import numpy as np
import vtk
from vtk.util import numpy_support as ns

import os,sys


class vtk_interpolation:
    def __init__(self, kernel='mean', radius=1, null_points=-1, **param):
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
        print(null_points)
        print(radius)
        try:
            self.null_points_strategy = float(null_points)
        except:
            if null_points=='closest':
                self.null_points_strategy = null_points
            else:
                print("--- ERROR: wrong null_points_strategy parameter ('closest' or float)")

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
    
        #vtk_pts = self.numpy_to_vtkpoints(np_pts) 
        vtk_pts = vtk.vtkPoints()
        vtk_pts.SetData(ns.numpy_to_vtk(np_pts, deep=True))
    
        ## Create vtk point data
        vars = []
        if dic_var is not None:
            for vname, vdata in dic_var.items():
                vars.append(ns.numpy_to_vtk(num_array=vdata.ravel(), deep=True, array_type=vtk.VTK_FLOAT))
                vars[-1].SetNumberOfComponents(1)
                vars[-1].SetName(vname)
    
        ## Merge everything into a vtkUnstructuredGrid
        grid = vtk.vtkUnstructuredGrid()
        grid.SetPoints(vtk_pts)
        for var in vars:
            grid.GetPointData().AddArray(var)
    
        if fname is not None:
            self.write_to_vtk_file(grid, fname)
    
        return grid

    def run(self, fname=None):
    
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
        if self.null_points_strategy=='closest':
            interpolator.SetNullPointsStrategyToClosestPoint()
        else:
            interpolator.SetNullPointsStrategyToNullValue()
            interpolator.SetNullValue(self.null_points_strategy)
    
        ## Run the interpolation
        interpolator.Update()
    
        # Write the result
        if fname is not None:
            self.write_to_vtk_file(interpolator.GetOutput(), fname)
    
        #return interpolator
        self.interpolator = interpolator

    def write_to_vtk_file(self, input_data, fname):
        ext = 'vtu'
        writer = vtk.vtkXMLUnstructuredGridWriter()
        fname = f"{fname}.{ext}"
        writer.SetFileName(fname)
        writer.SetInputData(input_data)
        #writer.SetDataModeToAscii() #Optional for debug - set the mode. The default is binary.
        writer.Write()
        print("--- VTK file saved to: {}".format(fname))


    def get_output(self):
        out_var = {}
        for a in range(self.interpolator.GetOutput().GetPointData().GetNumberOfArrays()):
            vname = self.interpolator.GetOutput().GetPointData().GetArrayName(a)
            numpy_var = ns.vtk_to_numpy(self.interpolator.GetOutput().GetPointData().GetArray(vname))
            out_var[vname] = numpy_var
        return out_var

def test(plot=False):
        
    param = {
            'var' : 'AL-BB-BH',
            'kernel' : 'inverse_distance', # choices are 'mean', 'gaussian' and 'inverse'
            'radius' : 5, # float in [km]
            'null_points' : 'closest', # choices are 'closest' or a float
            #'null_points' : -1.,
           }

    ## Create fake data
    s0 = 15
    s1 = 10

    source_lon = np.linspace(-s0, s0, 50)
    source_lat = np.linspace(-s0, s0, 50)
    source_lon, source_lat = np.meshgrid(source_lon, source_lat)
    source_dic_var = {} 
    source_dic_var['var1'] = np.sin(3*np.pi*source_lon/s0)*np.sin(3*np.pi*source_lat/s0) 

    target_lon = np.linspace(-s1, s1, 20)
    target_lat = np.linspace(-s1, s1, 20)
    target_lon, target_lat = np.meshgrid(target_lon, target_lat)


    ## Interpolate
    interp = vtk_interpolation(**param) 
    interp.set_source(source_lon, source_lat, source_dic_var, fname='source')
    interp.set_target(target_lon, target_lat, fname='target')
    interp.run(fname='result')
    interp_var = interp.get_output()

    if plot:
        import matplotlib.pyplot as plt

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, dpi=200)
        fig.suptitle('VTK interpolation example')
        
        ax1.set_title('source')
        ax1.scatter(source_lon, source_lat, c=source_dic_var['var1'], s=2)
        ax1.set_aspect('equal', 'box')
        
        ax2.set_title('target')
        ax2.scatter(target_lon, target_lat, c='w', s=6, edgecolors='k', linewidths=0.5)
        ax2.set_aspect('equal', 'box')
        
        ax3.set_title('target on source')
        ax3.scatter(source_lon, source_lat, c=source_dic_var['var1'], s=2)
        ax3.scatter(target_lon, target_lat, c='w', s=6, edgecolors='k', linewidths=0.5)
        ax3.set_aspect('equal', 'box')
        
        ax4.set_title('result')
        ax4.scatter(target_lon, target_lat, c=interp_var['var1'], s=6)
        ax4.set_aspect('equal', 'box')

        plt.tight_layout()

        plt.show()


if __name__=='__main__':

    test(plot=True)



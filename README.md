
# LSA-SAF ETAL-R validation tool

This tool can read MSG, EPS and MODIS data files and interpolate them on each other in order to compare them.


## Installation

- Create a new python env with python version >= 3.7
- Install dependencies:
    - h5py (read input files)
    - pyhdf (read input files)
    - numpy (process arrays)
    - matplotlib (output plots)
    - vtk (interpolate grid)
    - pandas (manage file paths)

This tool has been tested with the following versions:
```
python=3.7.10
h5py=3.3.0
pyhdf=0.10.3
numpy=1.20.3
matplotlib=3.3.4
vtk=8.2.0
pandas=1.2.4
```

Then you can clone it:
```
git clone https://github.com/aspyk/lsa_saf.git
```
    
## Usage/Examples

1. Choose variable and dates to process in `__main__`
1. Check and change if necessary the paths in the script:
    - path to data files in `get_all_filenames()`
    - path to longitude, latitude and land mask (if any) files in each SatelliteTools subclass
1. Choose your parameters for interpolation in `process_etal_series()`
1. Then run it:
```
python etalr_validation.py
```

TODO: create an unique external input file with all these parameters
  
## Documentation

### Interpolation parameters
From `vtk_tools.py`:
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

### Satellite data classes

```python
class SatelliteTools:

    self.product
    self.var
    self.slicing

    @property
    full_shape
    shape
    mask
    ground_mask
    lat
    lon

    interpolate_on()
    export_to_h5()
    load_h5_cache()
    describe()
```

```python
class <satellite_type>(SatelliteTools):

    self.ground_mask_conf = {'file': None,
                             'var': 'lwmask',
                             'type': mask_type}
    self.lat_conf = {'file': '/mnt/lfs/d30/vegeo/fransenr/CODES/DATA/NO_SAVE/MODIS/lat_modis.h5',
                     'var': 'lat',
                     'scaling': 1.}
    self.lon_conf = {'file': '/mnt/lfs/d30/vegeo/fransenr/CODES/DATA/NO_SAVE/MODIS/lon_modis.h5',
                     'var': 'lon',
                     'scaling': 1.}
        
    self.data = None
    self.data_scaling

    get_data()
```
that is for now:

```python
class MODIS(SatelliteTools):
    [...]

class EPS(SatelliteTools):
    [...]

class MSG(SatelliteTools):
    [...]
```

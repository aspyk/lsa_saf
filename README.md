
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

1. Choose variable and dates to process
1. Check and change if necessary the paths in the script:
    - path to data files
    - path to longitude and latitude files
    - path to land mask
1. Choose your parameters for interpolation
1. Then run it:
```
python etalr_validation.py
```

  
## Documentation

### Interpolation parameters

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

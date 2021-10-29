
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

1. In the YAML config file choose/change all the paths, parameters and variables you need.
1 In the code in `process_etal_series()` function, you still have to choose the slicing you want to apply on the data (full map or or only a part to debug or analyze a region)
1. Then run it:
```
python etalr_validation.py config_cnrm.yml
```

## Data source

### MODIS
Version 6 has been used for now:
- MCD43D51: black-sky albedo (directional hemispherical reflectance)
    - main page: https://lpdaac.usgs.gov/products/mcd43d51v006/
    - detailed specifications: https://ladsweb.modaps.eosdis.nasa.gov/filespec/MODIS/6/MCD43D51 
- MCD43D31: Albedo QA BRDF Quality:
    - main page: https://lpdaac.usgs.gov/products/mcd43d31v006/
    - detailed specifications: https://ladsweb.modaps.eosdis.nasa.gov/filespec/MODIS/6/MCD43D31

Files may be directly download from this website:
https://e4ftl01.cr.usgs.gov

Path example for a MCD43D51 product file:
https://e4ftl01.cr.usgs.gov/MOTA/MCD43D51.006/2014.01.14/MCD43D51.A2014014.006.2016146201118.hdf

#### Note about V6.1 from MODIS on 28/10/2021:

> Historic processing continues for MODIS V6.1 data products for years 2011 through 2018. It is the intent to complete processing of the archive by the end of the 2021 calendar year. The V6 land products will be available until the V6.1 processing is completed, upon which time the MODIS V6 products will be decommissioned.

V6.1 web pages:
- https://lpdaac.usgs.gov/products/mcd43d51v061/
- https://lpdaac.usgs.gov/products/mcd43d31v061/

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

A main class `SatelitteTools` has all the common functions to deal with 2D satellite data, and then subclasses are derived from it to add or override specific methods and attributes:


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
    whoami()
    get_config()
```

And the subclasses:

```python
class <satellite_type>(SatelliteTools):

    self.data = None
    self.data_scaling

    get_data()
```
Note that the only difference here are how to access the data. For now we have the following data types:

```python
class MODIS(SatelliteTools):
    [...]

class EPS(SatelliteTools):
    [...]

class MSG(SatelliteTools):
    [...]
```

### Outputs

Outputs may be of several types:
- global maps of a variable or a bias between two variables.
- temporal graph showing statistics for each date.

The class `Plots` is a helper class that make global map plotting easier using only one line to save a plot. Temporal graph are still done manually.

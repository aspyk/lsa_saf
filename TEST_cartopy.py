from urllib.request import urlopen
import ssl
from io import BytesIO

import cartopy.crs as ccrs
import matplotlib.pyplot as plt

'''
from https://scitools.org.uk/cartopy/docs/v0.16/\
            gallery/geostationary.html#sphx-glr-gallery-geostationary-py
'''
def geos_image():
    """
    Return a specific SEVIRI image by retrieving it from a github gist URL.

    Returns
    -------
    img : numpy array
        The pixels of the image in a numpy array.
    img_proj : cartopy CRS
        The rectangular coordinate system of the image.
    img_extent : tuple of floats
        The extent of the image ``(x0, y0, x1, y1)`` referenced in
        the ``img_proj`` coordinate system.
    origin : str
        The origin of the image to be passed through to matplotlib's imshow.

    """
    #img = plt.imread('EIDA50_201211061300_clip2.png')
    img = plt.imread('SMALL_EIDA50_201211061300_clip2.png')
    img_proj = ccrs.Geostationary(satellite_height=35786000)
    img_extent = [-5500000, 5500000, -5500000, 5500000]
    return img, img_proj, img_extent, 'upper'

def main():
    print('--- Retrieving image...')
    img, crs, extent, origin = geos_image()

    print('--- Orthographic projection...')
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0))
    ax.coastlines()
    ax.set_global()
    ax.imshow(img, transform=crs, extent=extent, origin=origin, cmap='gray')
    plt.savefig('res_test_ortho.png')

    plt.clf()


    print('--- Sinusoidal projection...')
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1, projection=\
                         ccrs.Sinusoidal(central_longitude=0.0, \
                            false_easting=0.0, false_northing=0.0))
    ax.coastlines()
    ax.set_global()
    print('Projecting and plotting image (this may take a while)...')
    ax.imshow(img, transform=crs, extent=extent, origin=origin, cmap='gray')

    plt.savefig('res_test_sinus.png')

if __name__=='__main__':
    main()

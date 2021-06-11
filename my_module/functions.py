"""A collection of function for doing my project."""
# Python standard-library
from urllib.parse import urlencode
from urllib.request import urlretrieve

# Third-party dependencies
from astropy import units as u
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
from astropy.utils.data import download_file
from IPython.display import Image
import astropy.coordinates as coord
import pandas as pd
import numpy as np

# Set up matplotlib and use a nicer set of plot parameters
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
plt.style.use(astropy_mpl_style)


def astroImage(name,im_size=12, im_pixels = 1024):
    """extract images within the Sloan Digital Sky Survey(SDSS) footptint.
    Parameters
    ----------
    name: string
        String to be called by SkyCoord.from_name(), used to search in SDSS catalog.
    im_size: int, default = 12
        int to controll the size of the image. Forms a 12 arcmin square by default. 1 arcmin = 1/60 degree
    im_pixels: int
        int to controll the pixels of the image.
        
    Returns
    -------
    Image(name): image
        Image extracted from SDSS
    
    """
    im_size = u.arcmin*im_size
    objCoord = SkyCoord.from_name(name) #obtain the coordinate of my object
    cutoutbaseurl = 'http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx'
    query_string = urlencode(dict(ra=objCoord.ra.deg, #R.A. means “right ascension” and Dec. means “declination”
                              dec=objCoord.dec.deg, #these are the two parts of the equatorial coordinate system. 
                              width=im_pixels, height=im_pixels,
                              scale=im_size.to(u.arcsec).value/im_pixels))
    url = cutoutbaseurl + '?' + query_string
    urlretrieve(url, name +'_cutout.jpg')
    print(name,'coordinate', objCoord) #print coordinates of the object in degree
    return Image(name +'_cutout.jpg') #output image

def fourimgshow(image, cmap1='hsv',cmap2 = 'twilight_shifted', cmap3='inferno',cmap4 = 'twilight'):
    """print four images within the fits image file downloaded previously.
    Parameters
    ----------
    image: numpy.ndarray
        numpy.ndarray data obtained from fits.getdata().
    cmap1: string, default = 'hsv'
        Cyclic colormap 'hsv'
    cmap2: string, default = 'twilight_shifted'
        Cyclic colormap 'twilight_shifted'
    cmap3: string, default = 'inferno'
        Sequential colormap 'inferno'
    cmap4: string, default = 'twilight'
        Cyclic colormap 'twilight'
        
    Returns
    -------
    Four images from matplotlib subplots
    
    """
    plt.figure()
    
    #subplot(r,c) provide the no. of rows and columns
    f, ax = plt.subplots(2,2,figsize=(10,10)) #creating 2 by 2 subplots with figure size 10 by 10
    
    # use the created array to output your multiple images. 
    ax[0,0].imshow(image, cmap = cmap1)
    ax[0,0].grid(False)
    ax[1,0].imshow(image, cmap = cmap2)
    ax[1,0].grid(False)
    ax[0,1].imshow(image, cmap = cmap3)
    ax[0,1].grid(False)
    ax[1,1].imshow(image, cmap = cmap4)
    ax[1,1].grid(False)
    
def basic_stats(fits_filelink):
    """print four images within the fits image file downloaded previously.
    Parameters
    ----------
    fits_filelink: string
        string of the .fits web link
        
    Returns
    -------
    image_min: numpy.int16
        minimum value of the image_data
    image_max: numpy.int16
        maximum value of the image_data
    image_mean: numpy.float64
        mean value of the image_data
    image_std: numpy.float64
        standard deviation value of the image_data
    
    """
    
    image_file = download_file(fits_filelink, cache=True ) #downloading fits file from link
    image_data = fits.getdata(image_file) #gets data from fits file
    
    #getting the basic stats using numpy
    image_min = np.min(image_data)
    image_max = np.max(image_data)
    image_mean = np.mean(image_data)
    image_std = np.std(image_data)
    
    #printing out basic stats
    print('Min:', image_min,'Max:', image_max,'Mean:', image_mean,'Stdev:', image_std)

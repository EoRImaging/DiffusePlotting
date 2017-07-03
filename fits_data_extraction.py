from scipy.io.idl import readsav
from netCDF4 import Dataset, num2date
from astropy.io import fits
import numpy as np
import healpy as hp
import scipy as sp
import matplotlib
from matplotlib.colors import LogNorm
import os


def fits_extractor(filename_Q, filename_U):
    if not filename_U.endswith('.fits'):
        raise ValueError('Stokes U file is not a fits file. Files must be .fits files')
    if not filename_Q.endswith('.fits'):
        raise ValueError('Stokes Q file is not a fits file. Files must be .fits files')
    # get header and data info from the fits file
    contents_Q = fits.open(filename_Q)
    pixelnum_Q = contents_Q[1].header['naxis1']
    data_Q = contents_Q[1].data
    nside_Q = contents_Q[1].header['nside']
    ordering_Q = contents_Q[1].header['ordering']

    contents_U = fits.open(filename_U)
    pixelnum_U = contents_U[1].header['naxis1']
    data_U = contents_U[1].data
    nside_U = contents_U[1].header['nside']
    ordering_U = contents_U[1].header['ordering']

    if pixelnum_U != pixelnum_Q:
        raise ValueError('files do not have same number of pixels.')

    if nside_U != nside_Q:
        raise ValueError('files do not have nside.')

    if ordering_U != ordering_Q:
        raise ValueError('files do not have nside.')

    # extract data from specified files
    pixels_Q = data_Q.field('PIXEL')
    signal_Q = data_Q.field('SIGNAL')
    pixels_U = data_U.field('PIXEL')
    signal_U = data_U.field('SIGNAL')

    if np.any(pixels_Q != pixels_U):
        raise ValueError('files do not have set of pixels.')

    return signal_Q, signal_U, nside_Q, pixels_Q, ordering_Q

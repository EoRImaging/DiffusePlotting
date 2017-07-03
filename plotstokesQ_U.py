from scipy.io.idl import readsav
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap
from astropy.io import fits
import numpy as np
import healpy as hp
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import plothealpix_map
import healpix_to_RA_dec


def plotstokesQ_U(filename_Q, filename_U, plotfile_base=None, plotdirectory=None,
                  save_show='show', projection='ortho'):
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

    # Finding x and y from Stokes parameters U and Q
    Q = signal_Q
    U = signal_U
    K = np.sqrt(U**2 + Q**2)
    U_pos = U[np.where(U >= 0)]
    U_neg = U[np.where(U < 0)]
    # theta is in radians
    theta = (.5 * np.arccos(((K + Q) / K) - 1))
    theta[np.where(U >= 0)] = theta
    theta[np.where(U < 0)] = theta + np.pi / 2
    # the x and y components of the theta-mag points.
    x_stokes = K * np.cos(theta)
    y_stokes = K * np.sin(theta)

    if plotfile_base is None:
        plotfile_base = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]

    if plotdirectory is None:
        plotdirectory = os.getcwd()
    if not os.path.isdir(plotdirectory):
        os.makedirs(plotdirectory)
    plotfile_base = os.path.join(plotdirectory, plotfile_base)
    plotfile = plotfile_base + '.png'
    # manipulate histogram for number of bins, color scheme.
    if filename_Q.endswith('.fits'):
        healpix_to_RA_dec.healpix_to_RA_dec(nside_Q, pixels_Q, ordering=ordering_Q)
    plothealpix_map.mapping(ra, dec, plotfile, theta, projection=projection,
                            save_show=save_show)
    plt.hist2d(x_stokes, y_stokes, bins=150, norm=LogNorm())
    plt.colorbar()

    # either show the graphs, or save them to a location.
    if save_show == 'show':
        plt.show()
    elif save_show == 'save':
        plt.savefig(plotfile)
        print 'saved plot to ' + plotfile
    else:
        raise ValueError('save_show needs to be equal to "save" or "show" to save or show the image.')
    # Using the function defined in plothealpix_map to graph the data on a globe.
    # Will save or show the plothealpix_map graph depending on option selected for histogram.

    # plt.savefig()
    return x_stokes, y_stokes

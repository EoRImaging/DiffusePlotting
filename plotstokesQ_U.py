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
import fits_data_extraction


def plotstokesQ_U(filename_Q, filename_U, plotfile_base=None, plotdirectory=None,
                  save_show='show', projection='ortho'):
    # Finding x and y from Stokes parameters U and Q
    signal_Q, signal_U, nside_Q, pixels_Q, ordering_Q = fits_data_extraction.fits_extractor(filename_Q, filename_U)
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

    if filename_Q.endswith('.fits'):

        ra, dec = healpix_to_RA_dec.healpix_to_RA_dec(nside_Q, pixels_Q, ordering=ordering_Q)
    plothealpix_map.mapping(ra, dec, plotfile, theta, projection=projection,
                            save_show=save_show)

    # manipulate histogram for number of bins, color scheme.
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

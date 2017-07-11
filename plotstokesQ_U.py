# Use Stokes Q and U files(linear polarizations) to plot a histogramatic plot which shows the polarization angle
# and magnitude of polarization, where color shows density of points.
# Create a map, plotted on a spherical, RA-dec grid, which plots the polarization angle
# of the data on a color scale.
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
import re
import plothealpix_map
import healpix_to_RA_dec
import fits_data_extraction


def plotstokesQ_U(filename_Q, filename_U, filename_I, histogram_plotfile_base=None, map_plotfile_base=None, plotdirectory=None,
                  save_show='show', file_extension='.png', projection='ortho'):
    # we set information from a fits file. theta is in radians.
    # x_stokes and y_stokes are the x and y components of points found using theta and magnitude.
    signal_Q, signal_U, signal_I, nside_Q, pixels_Q, ordering_Q = fits_data_extraction.fits_extractor(filename_Q, filename_U, filename_I)
    # Finding x and y from Stokes parameters U and Q
    Q = signal_Q
    U = signal_U
    I = signal_I
    K = np.sqrt(U**2 + Q**2)
    Kratio = K / I
    theta = np.arctan(U/(np.sqrt(Q**2 + U**2)+Q))
    theta[np.where(theta < 0)] = theta[np.where(theta < 0)] + np.pi
    x_stokes = np.cos(theta) * K
    y_stokes = np.sin(theta) * K
    Kratio[np.where(I < K)] = 1

    outliers_K = [K[np.all([K >= .1], axis=0)]]
    outliers_theta = [theta[np.all([K >= .1], axis=0)]]
    # setting file and path names
    filename_Q = os.path.split(filename_Q)[-1]
    filename_U = os.path.split(filename_U)[-1]
    obsID = re.findall('\d+',filename_Q)
    obsID = 'ObsID_' + ''.join(obsID)
    if histogram_plotfile_base is None:
        histogram_plotfile = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]
        histogram_plotfile_base = obsID + '_histogram'
    if map_plotfile_base is None:
        map_plotfile = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]
        map_plotfile_base = obsID
    if plotdirectory is None:
        plotdirectory = os.getcwd()
    if not os.path.isdir(plotdirectory):
        os.makedirs(plotdirectory)

    path_histogram_plotfile_base = os.path.join(plotdirectory, histogram_plotfile_base)
    path_map_plotfile_base = os.path.join(plotdirectory, map_plotfile_base)
    histogram_plotfile = path_histogram_plotfile_base + file_extension
    map_plotfile = path_map_plotfile_base + file_extension

    # convert a fits healpix-organized file to ra, dec coordinate system.
    if filename_Q.endswith('.fits'):
        ra, dec = healpix_to_RA_dec.healpix_to_RA_dec(nside_Q, pixels_Q, ordering=ordering_Q)

    # use plotting module.
    plothealpix_map.mapping(ra, dec, theta, 'theta', newplotfile_base=map_plotfile_base, projection=projection,
                                            save_show=save_show, file_extension=file_extension)
    #plothealpix_map.mapping(ra, dec, K, 'K', newplotfile_base=map_plotfile_base, projection=projection,
    #                                        save_show=save_show, file_extension=file_extension)
    # manipulate histogram for number of bins, color scheme.
    # histogram is of 'angle' of polarization between Stokes Q and Stokes U.
    plt.hist2d(x_stokes, y_stokes, bins=150, norm=LogNorm(), cmap=plt.cm.plasma)
    plt.colorbar()
    plt.title(histogram_plotfile_base)
    #plt.title(os.path.split(histogram_plotfile_base)[-1])

    # either show the histogram, or save it to a location.
    if save_show == 'show':
        plt.show(histogram_plotfile)
    elif save_show == 'save':
        plt.savefig(histogram_plotfile)
        print 'saved polarization histogram to ' + histogram_plotfile
    else:
        raise ValueError('save_show needs to be equal to "save" or "show" to save or show the image.')

    # x_stokes and y_stokes are the vector components necessary to create a drapery plot.
    return x_stokes, y_stokes, K

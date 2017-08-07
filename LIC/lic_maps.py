import numpy as np
import pylab as plt
from scipy.interpolate import griddata
import math
from numpy import interp
import lic_internal
import stokes_math
import plothealpix_map
import fits_data_extraction
import copy
import matplotlib.pyplot as plt
# automatically focuses at center


def LIC(dpi, size, length=31, full_image=True, disp_drapery='save', name_of_plot='flow-image.png'):
    filename_Q = 'data/4pol/1130789944_uniform_Residual_Q_HEALPix.fits'
    filename_U = 'data/4pol/1130789944_uniform_Residual_U_HEALPix.fits'
    filename_I = 'data/4pol/1130789944_uniform_Residual_I_HEALPix.fits'

    signal_Q, signal_U, signal_I, ra, dec =\
        fits_data_extraction.fits_extractor(filename_Q, filename_U, filename_I)
    K, theta, x_stokes, y_stokes = stokes_math.math(signal_I, signal_Q, signal_U)

    dpi = dpi
    if full_image is True:
        mean_ra = np.mean(ra)
        mean_dec = np.mean(dec)
        lower_ra = np.min(ra)
        upper_ra = np.max(ra)
        lower_dec = np.min(dec)
        upper_dec = np.max(dec)
    elif full_image is False:
        upper_ra, lower_ra, upper_dec, lower_dec =\
            fits_data_extraction.cutout_square(ra, dec)
        mean_ra = (upper_ra + lower_ra) / 2
        mean_dec = (upper_dec + lower_dec) / 2

    video = False

    x_range = upper_ra - lower_ra
    y_range = upper_dec - lower_dec
    aspect_ratio = y_range / x_range
    x_size = size
    y_size = int(round(size * aspect_ratio))

    xi = np.linspace(0, x_size, x_size)
    yi = np.linspace(0, y_size, y_size)
    # scale = divide by lenght of x, multoply by the max-min, add min ra. then check in right range
    ra_reg = (xi / max(xi))*(upper_ra - lower_ra) + lower_ra
    dec_reg = (yi / max(yi))*(upper_dec - lower_dec) + lower_dec

    # yi = ((yi / max(yi)) * (upper_range_dec - lower_range_dec)) + lower_range_dec

    indices = (np.where((ra >= lower_ra) & (ra <= upper_ra) &
                        (dec >= lower_dec) & (dec <= upper_dec)))
    #print indices
    x_stokes = x_stokes[indices]
    y_stokes = y_stokes[indices]
    ra = ra[indices]
    dec = dec[indices]
    #plothealpix_map.mapping(ra, dec, K, 'K', )
    x = ra_reg[:, np.newaxis]
    y = dec_reg[np.newaxis, :]
    x = np.repeat(x, y_size, axis=1).flatten()
    y = np.repeat(y, x_size, axis=0).flatten()
    #print ra_reg
    #print dec_reg
    # theta = np.arctan(y_stokes / x_stokes) + np.pi/2
    xval = griddata((ra, dec), x_stokes, (x, y), method='nearest')
    yval = griddata((ra, dec), y_stokes, (x, y), method='nearest')

    #this is just some graphing stuff
    xval_copy = copy.copy(xval)
    xval_copy = xval_copy.flatten()
    yval_copy = copy.copy(yval)
    yval_copy = yval_copy.flatten()
    K_copy = np.sqrt(xval_copy**2 + yval_copy**2)
    plothealpix_map.mapping(x, y, K_copy, 'K_copy', save_show=disp_drapery,
                            newplotfile_base='1130789944_uniform_Residual',
                            file_extension='.eps', projection='cyl')
    #plothealpix_map.mapping(ra, dec, K, 'K', )
    #plot = plt.tricontourf(x, y, theta_copy, 50, cmap=plt.cm.plasma)
    # plt.show()
    #plt.savefig('suzan', dpi=dpi)

    xval = np.reshape(xval, (x_size, y_size))
    yval = np.reshape(yval, (x_size, y_size))

    vectors = np.zeros((len(ra_reg), len(dec_reg), 2), dtype=np.float32)
    vectors[:, :, 0] = xval
    vectors[:, :, 1] = yval

    # lic_internal is configured for texture shape in the order (y, x)
    # instead of (x, y)
    texture = np.random.rand(y_size, x_size).astype(np.float32)
    plt.figure()

    plt.bone()
    frame = 0

    if video:
        kernellen = length
        for t in np.linspace(0, 1, 16 * 5):
            kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)*(1+np.sin(2*np.pi*5*(np.arange(kernellen)/float(kernellen)+t)))

            kernel = kernel.astype(np.float32)

            image = lic_internal.line_integral_convolution(vectors, texture, kernel)

            plt.clf()
            plt.axis('off')
            plt.figimage(image)
            plt.gcf().set_size_inches((x_size / float(dpi), y_size / float(dpi)))
            plt.savefig("flow-%04d.png" % frame, dpi=dpi)
            frame += 1
    else:
        kernellen = length
        kernel = np.sin(np.arange(kernellen) * np.pi / kernellen)
        kernel = kernel.astype(np.float32)
        image = lic_internal.line_integral_convolution(vectors, texture, kernel)

        plt.clf()
        plt.axis('off')
        plt.figimage(image)
        plt.gcf().set_size_inches((x_size / float(dpi), y_size / float(dpi)))
        if disp_drapery == 'save':
            plt.savefig(name_of_plot, dpi=dpi)
            print 'saved drapery to ' + name_of_plot
        elif disp_drapery == 'show':
            plt.show()
        elif disp_drapery == 'none':
            print "drapery plot not created"

from scipy.io.idl import readsav
from mpl_toolkits.basemap import Basemap, projection_params
from astropy.io import fits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import warnings
import matplotlib.cbook

warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


def mapping(nside, pixelnum, plotfile, data, ordering='ring', projection='ortho'):

    # file pixels may be organized as ring or nest
    print "diffuse plotting master branch version"
    if ordering.lower() == 'ring':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=False, lonlat=True)
    if ordering.lower() == 'nested':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=True, lonlat=True)

    ra[np.where(ra > 180)] -= 360
    # plot of Galactic gas with coordinate projection

    min_ra = np.min(ra)
    max_ra = np.max(ra)
    mean_ra = np.mean(ra)
    min_dec = np.min(dec)
    max_dec = np.max(dec)
    mean_dec = np.mean(dec)
    print('RA min, max, mean: ', min_ra, max_ra, mean_ra)
    print('Dec min, max, mean: ', min_dec, max_dec, mean_dec)

    print('{p} projection'.format(p=projection))
    if 'no corners' in projection_params[projection]:
        if 'bounding_lat' in projection_params[projection]:
            if projection in ['nplaea', 'npstere', 'npaeqd']:
                # northern polar projection
                if max_dec < 0:
                    print('polar projection in wrong hemisphere: data is in '
                          'the southern hemisphere')
                    boundinglat = 10
                else:
                    boundinglat = abs(min_dec)
            else:
                # southern polar projection
                if max_dec > 0:
                    print('polar projection in wrong hemisphere: data is in '
                          'the northern hemisphere')
                    boundinglat = -10
                else:
                    boundinglat = max_dec

            m = Basemap(projection=projection, lon_0=mean_ra,
                        boundinglat=boundinglat, resolution=None)
        else:
            # examples are 'moll', 'hammer', 'eck4'
            print('corners are not supported in this projection')
            # these projections can only make full-globe maps, any corner information is ignored
            m = Basemap(projection=projection, lon_0=mean_ra, lat_0=mean_dec,
                        resolution=None)
    elif 'llcrnrx' in projection_params[projection]:
        # example is 'orth'
        # can't use llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat even though the docs say you can
        # you also can't use width or height even though the docs say you can
        # you have to use llcrnrx, llcrnry, urcrnrx, urcrnry but I'm not sure how to calculate them....
        # these numbers are based on some online examples and experimentation
        deg_scale_factor = 120000
        llcrnrx = -1 * (max_ra - min_ra) * deg_scale_factor / 2
        urcrnrx = (max_ra - min_ra) * deg_scale_factor / 2
        if projection not in ['geos']:
            llcrnry = -1 * (max_dec - min_dec) * deg_scale_factor / 2
            urcrnry = (max_dec - min_dec) * deg_scale_factor / 2
            m = Basemap(projection=projection, lon_0=np.mean(ra), lat_0=np.mean(dec),
                        resolution=None, llcrnrx=llcrnrx, llcrnry=llcrnry,
                        urcrnrx=urcrnrx, urcrnry=urcrnry)
        else:
            # need a new scale factor for some reason
            deg_scale_factor = 175000
            llcrnry = min_dec * deg_scale_factor / 2
            urcrnry = max_dec * deg_scale_factor / 2
            m = Basemap(projection=projection, lon_0=np.mean(ra), resolution=None,
                        llcrnrx=llcrnrx, llcrnry=llcrnry, urcrnrx=urcrnrx,
                        urcrnry=urcrnry)
    elif 'lat_1,lat_2,lon_1,lon_2' in projection_params[projection]:
        m = Basemap(projection=projection, lon_0=mean_ra, lat_0=mean_dec,
                    lon_1=mean_ra, lat_1=min_dec, lon_2=mean_ra, lat_2=max_dec,
                    llcrnrlon=min_ra, llcrnrlat=min_dec, urcrnrlon=max_ra,
                    urcrnrlat=max_dec, resolution=None)
    elif 'lon_0' in projection_params[projection]:
        if projection in ['rotpole']:
            m = Basemap(projection=projection, lon_0=np.mean(ra), lat_0=np.mean(dec),
                        o_lat_p=0, o_lon_p=0,
                        llcrnrlon=min_ra, llcrnrlat=min_dec, urcrnrlon=max_ra,
                        urcrnrlat=max_dec, resolution=None)
        else:
            # example: 'aeqd'
            m = Basemap(projection=projection, lon_0=np.mean(ra), lat_0=np.mean(dec),
                        llcrnrlon=min_ra, llcrnrlat=min_dec, urcrnrlon=max_ra,
                        urcrnrlat=max_dec, resolution=None)
    else:
        m = Basemap(projection=projection, llcrnrlon=min_ra, llcrnrlat=min_dec,
                    urcrnrlon=max_ra, urcrnrlat=max_dec, resolution=None)

    fig = plt.figure()
    x, y = m(ra, dec)

    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    # draw parallels and meridians. Labels are 1/0 as [Top,bottom,right,left]
    # m.drawparallels(np.arange(-90., 120., 10.), labels=[1, 0, 0, 0])
    # m.drawmeridians(np.arange(0., 420., 10.), labels=[0, 0, 0, 1])
    # creates a scatter plot of the selected data on a globe.

    m.scatter(x, y, 3, marker='o', linewidths=.1, c=data, cmap=plt.cm.coolwarm)
    m.colorbar()
    plt.show()

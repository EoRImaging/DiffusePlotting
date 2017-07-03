import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import warnings
import matplotlib.cbook


# def healpix_to_RA_dec(nside, pixelnum, plotfile, data, ordering='ring', projection='ortho', save_show='show'):
def healpix_to_RA_dec(nside, pixelnum, ordering='ring'):
    # file pixels may be organized as ring or nest
    if ordering.lower() == 'ring':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=False, lonlat=True)
    if ordering.lower() == 'nested':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=True, lonlat=True)
    return ra, dec

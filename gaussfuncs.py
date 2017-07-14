# gaussfuncs.py - contains functions for plotting Gaussian functions on
# HealPIX pixel maps.
from math import exp
import healpy as hp
import numpy as np
from havnp import haversine
import itertools
from numba import autojit

@autojit
def computeGauss(amp, dist, var):
    prod =  amp*(np.exp(-np.square(dist)/(2*(var**2))))
    #if any(prod < 0):
    #    print('computeGauss returning negative number' + str(prod) + 'with inputs:')
    #    print('amp: '+str(amp))
    #    print('dist: '+str(dist))
    #    print('var: '+str(var))
    return prod

@autojit    
def makeGauss(nside, ra, dec, amp, varparam=np.sqrt(2*np.pi)): #bounds=None):#radius=6447797.0):
    print('RA length:' + str(len(ra)))
    print('dec length:' + str(len(dec)))
    print('amp length:' + str(len(amp)))
    var = varparam*1.0
    #raMin = min(ra)-2.0
    #raMax = max(ra)+2.0
    #decMin = min(dec)-2.0
    #decMax = max(dec)+2.0
    #m = np.arange(hp.nside2npix(nside))
    #vertices = hp.ang2vec([raMin,raMin,raMax,raMax],[decMin,decMax,decMax,decMin],lonlat=True)
    #m = hp.query_polygon(nside,vertices)
    #m_ra, m_dec = hp.pix2ang(nside, m, lonlat=True)
    #m_amp = np.zeros_like(m,dtype='float64')
    #ra = ra*1.0
    #dec = dec*1.0
    #numpix = m.size
    out_ra = np.array([])
    out_dec = np.array([])
    out_amp = np.array([])
    ### TO IMPLEMENT: Check ra, dec, amp are of same size
    #if isinstance(ra, np.ndarray):
    for rai, deci, ampi in itertools.izip(ra, dec, amp):
        # Used to convert from ra/dec to pixel to vec; should not be necessary
        #nearpix = hp.pix2vec(nside,hp.ang2pix(nside,rai,deci,lonlat=True))
        nearpix = hp.ang2vec(rai,deci,lonlat=True)
        # Calculate Gaussian on pixels within following radius (in degrees)
        radius = 3.0/60.0 # Hard-coded; remove eventually
        pxls = hp.query_disc(nside,nearpix,radius*np.pi/180.0)
        near_ra, near_dec = hp.pix2ang(nside, pxls, lonlat=True)
        dist = haversine(np.full_like(near_ra,rai), np.full_like(near_dec,deci), near_ra, near_dec)
        out_ra = np.append(out_ra,  near_ra)
        out_dec = np.append(out_dec, near_dec)
        out_amp = np.append(out_amp, computeGauss(ampi, dist, var))
        #print('Using nearest '+str(len(near_ra))+' pixels to source.')
    print('makeGauss output length: '+str(out_amp.size))
    return out_ra, out_dec, out_amp

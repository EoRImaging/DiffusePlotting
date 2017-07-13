from gaussfuncs import makeGauss
import healpy as hp
import numpy as np
from plothealpix_map import mapping
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from derp import *
import os, errno, re
import cartopy.crs as ccrs

# Hard-coded quantities; to be fixed
fname = './catalogs/1130788624_source_array.sav'
variances = [10000.0, 30000.0]
nside = 1024
obsID = fname[11:21] # TODO: replace with smarter regexp

def crdir(directory):
    try: 
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def numpize(fname):
    sourceList = buildSources(fname)
    print('Sources loaded')
    flag = 0
    for source in sourceList:
        #print('Reading source ' + str(source.ID) + ' with RA ' + str(source.RA) + ', dec ' + str(source.dec) + ', flux ' + str(source.I)) 
        if flag == 0:
            dataArray = source.getSourceData(['RA','dec','I'])
            flag = 1
        else:
            newArray = source.getSourceData(['RA','dec','I'])
            dprint('Appending array of shape '+str(np.shape(newArray))+' to dataArray of shape '+str(np.shape(dataArray)))
            dataArray = np.append(dataArray, newArray,axis=0)

    print('RA, dec, I extracted')
    raout = dataArray[:,0]
    decout = dataArray[:,1]
    fluxout = dataArray[:,2]
    np.savez_compressed(fname[:-3]+'npz',ra=raout,dec=decout,flux=fluxout)
    return raout, decout, fluxout

def gaussify(sra, sdec, sflux, var, nside, obsID):
    # Gaussed data output file name
    gfname='./gaussdata/'+obsID+'_gauss_nside'+str(nside)+'_var'+str(var)+'.npz'

    #if not(os.path.isfile(gfname)):
    #if os.path.isfile(fname):
    #    npzfile = np.load('catrun1.npz')
    #    sra = npzfile['ra']
    #    sdec = npzfile['dec']
    #    sflux = npzfile['flux']
    #else:
    #    sra, sdec, sflux = numpize(fname)
    
    # Normalize to 1.0 to avoid getting negative numbers from log10
    # TODO: Implement +=1, multiplicative normalization switch, compare
    print('Min flux: '+str(min(sflux)))
    print('Max flux: '+str(max(sflux)))
    normConst = 1.0/min(sflux)
    print('Normalizing flux to 1.0 by multiplying by '+str(normConst))
    sflux *= normConst
    bounds = [min(sra)-5.0, max(sra)+5.0, min(sdec)-5.0, max(sdec)+5.0]
    ra, dec, dist, amp = makeGauss(nside,sra,sdec,np.log10(sflux),var,bounds=bounds)
    print('Saving ' + gfname)
    np.savez_compressed(gfname,ra=ra,dec=dec,amp=amp,dist=dist)
    return ra, dec, amp, dist
    #    gaussplot(ra,dec,flux,100.0,var,nside,obsID)
    ### END def gaussify

def gaussplot(ra,dec,amp,cap,var,nside,obsID):
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    plt.tricontourf(ra,dec,np.clip(amp,None,cap),100,cmap='Greys',transform=proj)
    plt.colorbar()
    plt.ylabel('Declination ($^\circ$)')
    plt.xlabel('Right ascension ($^\circ$)')
    plt.title('obsID: '+obsID+'\n var = '+str(var)+' cap = '+str(cap)+' nside='+str(nside))
    ax.set_xticks([45,60,75,90,105,120])
    ax.set_yticks([-20,-30,-40,-50,-60,-70])
    gl = ax.gridlines(crs=proj)
    gl.xlabels_bottom=True
    gl.xlines=True
    #gl.xlocator = mticker.LinearLocator(6)
    gl.ylabels_left=True
    gl.ylines=True
    #gl.ylocator
    #gl.xformatter = LONGITUDE_FORMATTER
    #plt.scatter(srcs['ra'],srcs['dec'],s=srcs['flux'],alpha=0.3)
#plt.savefig('./Plots/OG_run1.png')
    plotname = './gaussplots/run1_prelog'+str(nside)+'_var'+str(var)+'_cap' + str(cap) + '_projected.png'
    print('Saving '+plotname)
    plt.savefig(plotname)
    plt.close()
#   

def makeGaussPlot(fname, nside, var, cap,saveData=True):
    # Create output directories if they don't already exist
    crdir('./gaussdata')
    crdir('./gaussplots')
    
    # Extract obsID from fname using newbie regex
    # DEBUGGING NOTE: Any instance of '113' elsewhere in the file name followed
    # by 7 digits will cause an error!
    searchObj = re.search( r'113\d{7}', fname)
    obsID = searchObj.group()
    print('obsID: ' + obsID)
    # Gauss'd file name
    gfname='./gaussdata/'+obsID+'_gauss_nside'+str(nside)+'_var'+str(var)+'.npz'
    
    # If Gauss'd output file already exists, just plot it
    if os.path.isfile(gfname):
        print('Loading existing Gaussed healpix data')
        npdict = np.load(gfname)
        ra = npdict['ra']
        dec = npdict['dec']
        amp = npdict['amp']
        print('Plotting Gaussians')
        gaussplot(ra, dec, amp, cap, var, nside, obsID)
        print('makeGaussPlot done')
        return
    else: # If Gauss'd data doesn't exist:
        # If ra, dec, flux (I) data is already saved in np format, use it
        if os.path.isfile(fname[:-3]+'npz'):
            print('Loading existing .npz file of RA, dec, flux (I)')
            npdict=np.load(fname[:-3]+'npz')
            sra = npdict['ra']
            sdec = npdict['dec']
            sflux = npdict['flux']
            print('Computing Gaussians')
            ra, dec, amp, dist = gaussify(sra, sdec, sflux, var, nside, obsID)
            print('Plotting Gaussians')
            gaussplot(ra,dec,amp,cap,var,nside,obsID)
            print('makeGaussPlot done')
            return
        else:
            print('Extracting ra, dec, flux (I) from IDL .sav file')
            sra, sdec, sflux = numpize(fname)
            print('Computing Gaussians')
            ra, dec, amp, dist = gaussify(sra, sdec, sflux, var, nside, obsID)
            print('Plotting Gaussians')
            gaussplot(ra, dec, amp, cap, var, nside, obsID)
            print('makeGaussPlot done')
            return
        

#gfname='./catalogs/'+obsID+'_gauss_nside'+str(nside)+'_var'+str(var)+'.npz'

#if not(os.path.isfile(gfname)):
#    if os.path.isfile(fname):
#        npzfile = np.load('catrun1.npz')
#        sra = npzfile['ra']
#        sdec = npzfile['dec']
#        sflux = npzfile['flux']
#    else:
#        sra, sdec, sflux = numpize(fname)
    

    # TODO: Implement +=1, multiplicative normalization switch, evaluatei
#    print('Min flux: '+str(min(sflux)))
#    print('Max flux: '+str(max(sflux)))
#    normConst = 1.0/min(sflux)
#    print('Normalizing to 1.0 by multiplying by '+str(normConst))
#    sflux *= normConst
#    for var in variances:#, 3000.0, 10000.0, 300.0, 100000.0, 300000.0]:
#        ra, dec, dist, amp = makeGauss(nside,sra,sdec,np.log10(sflux),var,bounds=[15,80,-70,-20])
#        print('Pre-logged Gauss grid generated for variance '+str(var))
#        np.savez_compressed(gfname,ra=ra,dec=dec,dist=dist,amp=amp)
#        gaussplot(ra,dec,flux,100.0,var,nside,obsID)
#
#ra[np.where(ra > 180)] -= 360

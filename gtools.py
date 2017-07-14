from gaussfuncs import makeGauss
import healpy as hp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from derp import *
import os, errno, re
import cartopy.crs as ccrs

def crdir(directory):
    """
crdir - Create directory. Creates a directory if it doesn't already exist.
If the directory already exists, it handles the error intelligently.
If the directory cannot be created for any other reason, an error is raised.
Input: string containing path of directory to be created. 
Returns nothing.   
    """
    try: 
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
### END def crdir


def numpize(fname,saveData):
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
    if saveData:
        np.savez_compressed(fname[:-3]+'npz',ra=raout,dec=decout,flux=fluxout)
    return raout, decout, fluxout

def gaussify(sra, sdec, sflux, var, nside, obsID, saveData):
    # Gaussed data output file name
    gfname='./gaussdata/'+obsID+'_gauss_nside'+str(nside)+'_var'+str(var)+'.npz'

    # Normalize to 1.0 to avoid getting negative numbers from log10
    # TODO: Implement +=1, multiplicative normalization switch, compare
    print('Min flux: '+str(min(sflux)))
    print('Max flux: '+str(max(sflux)))
    normConst = 1.0/min(sflux)
    print('Normalizing flux to 1.0 by multiplying by '+str(normConst))
    sflux *= normConst
    
    ra, dec, amp = makeGauss(nside,sra,sdec,np.log10(sflux),var)
    #print('[gaussify] Values returned from makeGauss: ')
    #print('RA: ' + str(ra))
    #print('dec: ' + str(dec))
    #print('amp: ' + str(amp))
    if saveData:
        print('Saving ' + gfname)
        np.savez_compressed(gfname,ra=ra,dec=dec,amp=amp)

    return ra, dec, amp
    ### END def gaussify

def gaussplot(ra,dec,amp,cap,var,nside,obsID,savePlot,colorbarLevels):
    print('[gaussplot] Brightest single point has amplitude ' + str(max(amp)))
    # Set figure size to be nice and big
    mpl.rcParams['figure.figsize'] = [24.0, 18.0]
    
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    #plt.tricontourf(ra,dec,np.clip(amp,None,cap),colorbarLevels,cmap='Greys',transform=proj)
    #print('Parameters passed into gaussplot:')
    #print('RA: ' + str(ra))
    #print('dec: ' + str(dec))
    #print('amp: ' + str(amp))
    plt.scatter(ra,dec,c=np.clip(amp,None,cap),s=0.0001,cmap='Greys',alpha=0.2,transform=proj)
    plt.colorbar()
    plt.gca().invert_xaxis()
    # Label the plot all pretty-like
    #plt.ylabel('Declination ($^\circ$)')
    #plt.xlabel('Right ascension ($^\circ$)')
    # Awful hacky workaround for axis labels because Cartopy sucks
    ax.set_xlabel('Right ascension ($^\circ$)')
    ax.set_ylabel('Declination ($^\circ$)')
    plt.title('obsID: '+obsID+'\n var = '+str(var)+'$^\circ$ cap = '+str(cap)+' nside='+str(nside))
    #ax.set_xticks([45,60,75,90,105,120])
    #ax.set_yticks([-20,-30,-40,-50,-60,-70])
    gl = ax.gridlines(crs=proj)
    gl.xlabels_bottom=True
    gl.xlines=True
    #gl.xlocator = mticker.LinearLocator(6)
    gl.ylabels_left=True
    gl.ylines=True
    #gl.ylocator
    #gl.xformatter = LONGITUDE_FORMATTER
    
    # Plot location of Pictor A galaxy
    #plt.plot((5.0+19.0/60.0+49.7/3600.0)*360.0/24.0,-(45.0+46.0/60.0+44.0/3600.0),'ro',ms=1.0)

    if savePlot:
        plotname = './gaussplots/'+obsID+'_nside'+str(nside)+'_var'+str(var)+'_cap' + str(cap) + '_projected.png'
        print('Saving '+plotname)
        plt.savefig(plotname)
        print('Saved.')
        plt.close()
        return
    else:
        print('ACK NOT IMPLEMENTED DO NOT DO THAT')
        return
        #return plt, ax, gl
#   

def makeFromScratch(fname, nside, var, cap, obsID, saveData, savePlot, 
                    colorbarLevels):
    print('Extracting ra, dec, flux (I) from IDL .sav file')
    sra, sdec, sflux = numpize(fname, saveData)
    print('Computing Gaussians')
    ra, dec, amp = gaussify(sra, sdec, sflux, var, nside, obsID, saveData)
    print('Plotting Gaussians')
    gaussplot(ra, dec, amp, cap, var, nside, obsID, savePlot, colorbarLevels)
    print('makeGaussPlot done')
    return
    
def handleFuckery(fname, nside, var, cap, obsID, saveData, loadData, savePlot, 
                  colorbarLevels):
    print('')
    print('makeGaussPlot started using following parameters:\n')
    print('obsID (GPS seconds): ' + obsID)
    print('nside - HealPIX grid resolution: ' + str(nside))
    print('var (degrees) - Variance for each Gaussian: ' + str(var))
    print('cap (arb. units) - Upper limit on plotted values: ' + str(cap))
    if colorbarLevels > 998:
        raise ValueError('Cannot have more than 998 color levels in output contour plot.')
        #colorbarLevels = 998
    elif colorbarLevels < 0:
        raise ValueError('Cannot have negative number of colorbar levels.')
        #colorbarLevels = 10
    else:
        print('Output colorbar levels: ' + str(colorbarLevels))
    if saveData:
        print('saveData is True; numpy ra/dec/flux arrays for sources and their Gaussian convolution will be saved under ./gaussdata')
    else:
        print('saveData is False; no intermediate data will be saved to disk.')
    if loadData:
        print('loadData is True; existing ra/dec/flux numpy arrays of Gaussian-convolved or catalog data will be loaded if available.')
    else:
        print('loadData is False; the plot will be created from scratch starting with the .sav file.')
    if savePlot:
        print('savePlot is True; the plot will be saved under ./gaussplots')
    else:
        raise NotImplementedError('savePlot is False; this option is not implemented.')
        #savePlot = True
    print('')


def makeGaussPlot(fname, nside, var, cap, saveData=True, loadData=True, 
                  savePlot=True, colorbarLevels=50):
    """Make truncated-Gaussian plot of EoR catalog .sav data.
    
    Arguments:
    fname -- Path to IDL .sav file of source catalog data to plot
    nside -- Resolution of HealPIX grid to convolve & plot upon
    var -- Variance of Gaussian in decimal degrees
    cap -- Upper limit on plotted values, arbitrary units
    
    Keyword arguments: TODO
    """
    # Extract obsID from fname using newbie regex
    # DEBUGGING NOTE: Any instance of '113' elsewhere in the file name followed
    # by 7 digits will cause an error!
    searchObj = re.search( r'113\d{7}', fname)
    obsID = searchObj.group()
    handleFuckery(fname, nside, var, cap, obsID, saveData, loadData, savePlot, colorbarLevels)
    #print('')
    #print('makeGaussPlot started using following parameters:\n')
    #print('obsID (GPS seconds): ' + obsID)
    #print('nside - HealPIX grid resolution: ' + str(nside))
    #print('var (degrees) - Variance for each Gaussian: ' + str(var))
    #print('cap (arb. units) - Upper limit on plotted values: ' + str(cap))
    #if colorbarLevels > 998:
    #    print('Cannot have more than 998 levels in output contour plot. Reducing from ' + str(colorbarLevels) + ' to 998.')
    #    colorbarLevels = 998
    #elif colorbarLevels < 0:
    #    print('Cannot have negative number of colorbar levels. Defaulting to 10.')
    #    colorbarLevels = 10
    #else:
    #    print('Output colorbar levels: ' + str(colorbarLevels))
    #if saveData:
    #    print('saveData is True; numpy ra/dec/flux arrays for sources and their Gaussian convolution will be saved under ./gaussdata')
    #else:
    #    print('saveData is False; no intermediate data will be saved to disk.')
    #if loadData:
    #    print('loadData is True; existing ra/dec/flux numpy arrays of Gaussian-convolved or catalog data will be loaded if available.')
    #else:
    #    print('loadData is False; the plot will be created from scratch starting with the .sav file.')
    #if savePlot:
    #    print('savePlot is True; the plot will be saved under ./gaussplots')
    #else:
    #    print('savePlot is False; this option is not implemented, so savePlot is being set to True.')
    #    savePlot = True
    #print('')
    
    # Create output directories if they don't already exist
    if saveData:
        crdir('./gaussdata')
    if savePlot:
        crdir('./gaussplots')
    
    # Make plots from scratch if loadData is False
    if not(loadData):
        makeFromScratch(fname, nside, var, cap, obsID, saveData, savePlot, colorbarLevels)
        return
    
    
    # Gauss'd file name
    gfname='./gaussdata/'+obsID+'_gauss_nside'+str(nside)+'_var'+str(var)+'.npz'
    
    # If Gauss'd output file already exists, just plot it
    if os.path.isfile(gfname):
        print('Loading existing Gaussian-convolved healpix data')
        npdict = np.load(gfname)
        ra = npdict['ra']
        dec = npdict['dec']
        amp = npdict['amp']
        print('Plotting Gaussians')
        print('Plot data length: ' + str(len(amp)))
        gaussplot(ra, dec, amp, cap, var, nside, obsID, savePlot, colorbarLevels)
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
            ra, dec, amp = gaussify(sra, sdec, sflux, var, nside, obsID, saveData)
            print('Plotting Gaussians')
            gaussplot(ra,dec,amp,cap,var,nside,obsID,savePlot,colorbarLevels)
            print('makeGaussPlot done')
            return
        else: # If np format ra, dec, flux don't exist, extract from sav file
            makeFromScratch(fname, nside, var, cap, obsID, saveData, savePlot, colorbarLevels)
            return

# A cool thing Bryna did that I may need to reference sometime
#ra[np.where(ra > 180)] -= 360

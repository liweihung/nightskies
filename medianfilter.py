#-----------------------------------------------------------------------------#
#medianfilter.py
#
#NPS Night Skies Program
#
#Last updated: 2017/05/01
#
#This script uses multiprocessing to apply median filter to each image. The 
#filter is a circle with 1 degree diameter. This filter size was selected to 
#ensure most (or all) point sources are effectively filtered out. Here, MaxIM DL
#is needed to convert the fits images to tiff images that are compatible with 
#ArcGIS to make mosaics. 
#
#Input: 
#   (1) Calibrated image data
#   (2) Plate scale from image header
#
#Output:
#   (1) Median filtered images in median_ib###.tif format
#
#History:
#	Dan Duriscoe -- Created in java script as "med1.js"; used PixInsight
#	Li-Wei Hung -- Rewrote in python; replaced PixInsight by python functions
#
#-----------------------------------------------------------------------------#

from astropy.io import fits
from glob import glob, iglob
from multiprocessing import Pool
from scipy.ndimage.filters import median_filter
from win32com.client import Dispatch

import itertools
import numpy as n
import os

# Local Source
import filepath  

#-----------------------------------------------------------------------------#    
def FilterImage(arg):
    '''
    Apply the median filter using the given mask size and save the images in 
    .tif format.
    '''
    fn, mask = arg
    T = Dispatch('Maxim.Document')

    m = int(fn[-7:-4])
    if m in range(0,50,5): 
        print 'filtering images %i/45'%m
    
    outFits = 'tiff/median_%s.fit'%fn[-9:-4] #temporary file
    outTiff = 'tiff/median_%s.tif'%fn[-9:-4] #output file
    D = fits.open(fn)[0]
    D.data = median_filter(D.data, footprint=mask)
    D.writeto(fn[:-9]+outFits, clobber=True)
    T.OpenFile(fn[:-9]+outFits)
    T.SaveFile(fn[:-9]+outTiff,5,False,1,0)
    T.Close
    os.remove(fn[:-9]+outFits)
    



def filter(dnight, sets, filter):
    '''
    This module creats a mask and calls the FilterImage module to apply median 
    filter to the calibrated images through multiprocessing.
    '''
        
    #filter paths
    F = {'V':'', 'B':'B/'}
    
    #set the mask radius to be ~0.5 degree
    calsetp = filepath.calibdata+dnight+'/S_0%s/%s' %(sets[0][0],F[filter])
    for fn in iglob(calsetp+'ib???.fit'):
        H = fits.open(fn)[0].header
        if 'CDELT1' in H.keys(): 
            plate_scale = abs(H['CDELT1']) #[deg/pix] X-axis plate scale
            r = n.floor(1./plate_scale/2)  #[pix] radius of the filter mask
            break
    
    #generate the mask
    X, Y = n.meshgrid(n.arange(2*r+1), n.arange(2*r+1))
    R = n.sqrt((X-r)**2+(Y-r)**2)
    mask = n.zeros_like(R)
    mask[n.where(R<=r)] = 1
    
    #multiprocessing the filtering process
    p = Pool()
    
    #loop through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata+dnight+'/S_0%s/%s' %(s[0],F[filter])
        arg = itertools.izip(glob(calsetp+'ib???.fit'), itertools.repeat(mask))
        p.map(FilterImage, arg)
        
    p.close()
    p.join()
    
    #close MaxIm_DL application
    #os.system('taskkill /f /im MaxIm_DL.exe')
    

    
if __name__ == "__main__":
    #import time
    #t1 = time.time()
    #filter('FCNA160803', ['1st',], 'V')
    #t2 = time.time()
    #print 'Total time: %.1f min' %((t2-t1)/60)
    pass

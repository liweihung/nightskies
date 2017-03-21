#-----------------------------------------------------------------------------#
#medianfilter.py
#
#NPS Night Skies Program
#
#Last updated: 2016/12/19
#
#This script applies median filter to each image. The filter is a circle with 1 
#degree diameter. This filter size was selected to ensure most (or all) point 
#sources are effectively filtered out. Here, MaxIM DL is need to convert the 
#fits images to tiff images that are compatible for ArcGIS to make mosaic. 
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
from scipy.ndimage.filters import median_filter
from win32com.client import Dispatch

import numpy as n
import os

# Local Source
import filepath  

#-----------------------------------------------------------------------------#    
def filter(dnight, sets, filter):
    '''
    This module applies the median filter to the calibrated images and saves them.
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
    
    T = Dispatch('Maxim.Document')

    #loop through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata+dnight+'/S_0%s/%s' %(s[0],F[filter])
        
        # loop through each file in the set
        for fn in iglob(calsetp+'ib???.fit'):
            outFits = 'tiff/median_%s.fit'%fn[-9:-4] #temporary file
            outTiff = 'tiff/median_%s.tif'%fn[-9:-4] #output file
            D = fits.open(fn)[0]
            D.data = median_filter(D.data, footprint=mask)
            D.writeto(calsetp+outFits, clobber=True)
            T.OpenFile(calsetp+outFits)
            T.SaveFile(calsetp+outTiff,5,False,1,0)
            T.Close
            os.remove(calsetp+outFits)
            
            m = int(fn[-7:-4])
            if m in range(0,50,5): 
                print 'filtered images %i/45'%m

    
    #close MaxIm_DL application
    os.system('taskkill /f /im MaxIm_DL.exe')
    
if __name__ == "__main__":
    #filter('FCNA160803', ['1st',], 'B')
    pass

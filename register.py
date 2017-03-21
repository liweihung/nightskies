#-----------------------------------------------------------------------------#
#register.py
#
#NPS Night Skies Program
#
#Last updated: 2016/11/10
#
#This script uses the PinPoiot.Plate object (through ACP) to compare the 
#position of the stars in the images to the published Tycho2 catalog position. 
#If the plate object cannot solve the entire image, the input image will be 
#cropped so the plate object will only try to solve the center 200x200 pix. 
#If the images still can't be solved, it will be skipped. The solved X and Y 
#coordinates are stored in the header under 'CRVAL1' and 'CRVAL2'.
#
#Note: In order to use the ACP objects, the Python must be a 32-bit version. 
#
#Input: 
#   (1) Calibrated images
#   (2) Tycho2 catalog
#
#Output:
#   (1) Original calibrated images with updated header
#   (2) Returns a list of files solved with the cropped images and a list of 
#       files that are failed to be solved
#
#History:
#	Dan Duriscoe -- Created in visual basic as "solve_images_v4b.vbs"
#	Li-Wei Hung -- Cleaned and translated to python
#
#-----------------------------------------------------------------------------#

from astropy.io import fits
from datetime import datetime as Dtime
from glob import glob, iglob
from win32com.client import Dispatch

import pdb
import numpy as n

# Local Source
import filepath     

#-----------------------------------------------------------------------------#

def matchstars(dnight, sets, filter):
    p = Dispatch('PinPoint.Plate')
    p.Catalog = 4  # Tycho2 catalog
    p.CatalogPath = filepath.catalog
    
    cropped_fn = []
    failed_fn = []
    l_dir = len(filepath.calibdata+dnight)+1

    #looping through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata + dnight + '/S_0' + s[0] + '/'
        print 'Registering images in', dnight + '/S_0' + s[0] + '...'
        
        #both V and B bands
        if filter == 'V':
            file = glob(calsetp+'ib???.fit')
        elif filter == 'B':
            file = glob(calsetp+'B/ib???.fit')
        
        ti = Dtime.now()
        for fn in file:
            fn_orig = fn
            m = int(fn_orig[-7:-4])
            
            if m in range(0,50,5): 
                print 'Solved images %i/45'%m

            # Masking the area near the horizon in image 0-15
            if m < 16: 
                f = fits.open(fn,uint=False)[0]
                f.data[630:] = 0.
                fn = fn[:-4]+'c.fit'
                f.writeto(fn, clobber=True)
            
            p.attachFits(fn)
            p.ArcsecPerPixelHoriz = 96
            p.ArcsecPerPixelVert = 96
            p.SigmaAboveMean = 3
            p.Minimumbrightness = 2500
            p.FindImageStars()
                        
            p.RightAscension = p.TargetRightAscension
            p.Declination = p.TargetDeclination
            
            if filter == 'B':
                p.colorband = 1 # B-band
            elif filter == 'V':
                p.colorband = 2 # V-band
                
            p.CatalogMaximumMagnitude = 7.
            p.CatalogMinimumMagnitude = 2.
            p.CatalogExpansion = 0.1
            p.Maxsolvetime = 60
            p.FindCatalogStars()
            
            try: 
                p.Solve()
                p.UpdateFITS()
            except:
                # trying to just solve the cropped (200x200 pix) image
                f = fits.open(fn,uint=False)[0]
                l = len(f.data)/2
                f.data = f.data[l-100:l+100,l-100:l+100]
                fn = fn[:-4]+'s.fit'
                f.writeto(fn, clobber=True)
                
                p.DetachFITS()
                p.attachFits(fn)
                p.ArcsecPerPixelHoriz = 96
                p.ArcsecPerPixelVert = 96
                p.SigmaAboveMean = 2
                p.Minimumbrightness = 2000
                p.FindImageStars()
                
                p.RightAscension = p.TargetRightAscension
                p.Declination = p.TargetDeclination
            
                p.CatalogMaximumMagnitude = 9.
                p.CatalogMinimumMagnitude = 2.
                p.CatalogExpansion = 0.3
                p.Maxsolvetime = 40
                p.FindCatalogStars()
                
                try:
                    p.Solve()
                    p.UpdateFITS()
                    cropped_fn.append(fn_orig)
                except:
                    failed_fn.append(fn_orig)
                    
            p.DetachFITS()
            
            #save the updated fits header to the first row horizon images
            if (m<16) & (p.solved==True): 
                f = fits.open(fn_orig,uint=False)
                data = f[0].data.copy()
                f.close()
                header = fits.open(fn)[0].header
                fits.writeto(fn_orig, data, header=header, clobber=True)

        tf = Dtime.now()
        print 'Total time spent on solving plates:', tf-ti
    return(cropped_fn, failed_fn)
    
    
    
    
    
    
    
    
    
    
    
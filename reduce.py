#-----------------------------------------------------------------------------#
#reduce.py
#
#NPS Night Skies Program
#
#Last updated: 2016/11/15
#
#This script performs basic image reduction on the image data collected by the 
#NPS Night Skies Program. The script corrects images for:
#	(1) Bias 
#	(2) Dark 
#	(3) Flat
#   (4) Linearity response of the detector
#
#Note: Camera temperature must be matched to linearity curve
#
#Input: Either from the raw data folder or the info passed by process.py
#   (1) Raw data -- dark, bias, and science frames
#   (2) Master images -- flat
#   (3) Linearity curves
#
#Output:
#   (1) Calibrated science images in fits and tiff formats
#   (2) Measured bias drift in txt file and png image
#
#History:
#	Dan Duriscoe -- Created in 2011 in visual basic as "calibrate images.vbs"
#	Li-Wei Hung -- Cleaned and translated to python
#
#-----------------------------------------------------------------------------#

from astropy.io import fits
from glob import glob, iglob
from scipy.misc import imsave 
from win32com.client import Dispatch

import pdb
import matplotlib.pyplot as plt
import numpy as n
import os
import shutil

# Local Source
import filepath     

#-----------------------------------------------------------------------------#

def reducev(dnight, sets, flatname, curve):
    '''
    This module is for calibrating the V-band data.
    '''
    #read in the linearity curve (ADU, multiplying factor)
    xp, fp = n.loadtxt(filepath.lincurve+curve+'.txt', unpack=True)
    
    #read in the flat
    flat = fits.open(filepath.flats+flatname,unit=False)[0].data
    
    T = Dispatch('Maxim.Document')
    
    #looping through all the sets in that night
    for s in sets:
        rawsetp = filepath.rawdata + dnight + '/' + s + '/'
        calsetp = filepath.calibdata + dnight + '/S_0' + s[0] + '/'
        print dnight + '/' + s + '...'
        if os.path.isdir(calsetp):
            print 'Replacing old calibrated files...'
        else:
            os.makedirs(calsetp)
            os.makedirs(calsetp+'tiff/')
    
        #correct the darks for bias and linearity response; crop the biases
        dark = []
        bias = []
        for i in range(1,6):
            darkraw = fits.open(rawsetp+'dark%i.fit'%i,uint=False)[0]
            biasraw = fits.open(rawsetp+'mbias%i.fit'%i,uint=False)[0]
            biasrawd = biasraw.data
            darkp = (darkraw.data - biasrawd).clip(0) #replace negatives with 0
            darki = darkp * n.interp(darkp,xp,fp) #correct linearity response
            biascrop = biasrawd[486:536,486:536]
            dark.append(darki)
            bias.append(biasrawd)
            fits.writeto(calsetp+'thermal%i.fit'%i, darki, overwrite=True)
            fits.writeto(rawsetp+'biasc%i.fit'%i, biascrop, overwrite=True,
                         header=biasraw.header)
        
        #average combine to generate the master thermal and bias
        corthermal = n.average(dark,axis=0)
        combias = n.average(bias,axis=0)
        fits.writeto(calsetp+'corthermal.fit', corthermal, overwrite=True)
        fits.writeto(calsetp+'combias.fit', combias, overwrite=True)
        
        #measure the bias drift for each frames
        nb = len(glob(rawsetp+'biasc*.fit'))
        biasc = n.empty([nb,50,50])
        Temp = n.empty([nb,1])                    # CCD temperature [C]
        for i in range(1,nb+1):
            biasci = fits.open(rawsetp+'biasc%i.fit'%i,uint=False)[0]
            biasc[i-1] = biasci.data
            Temp[i-1] = biasci.header['CCD-TEMP']
        baseline = n.average(biasc[:6])
        biasdrift_full = n.average(biasc,axis=(1,2)) - baseline
        biasdrift = biasdrift_full[5:]
        n.savetxt(filepath.calibdata+dnight+'/biasdrift_%s.txt'%s[0],
                  biasdrift,fmt='%5.3f',header='delta_bias[ADU]')
        fig = plt.figure('bias')
        plt.plot(n.arange(len(biasc)), n.zeros(len(biasc)), 'k--')
        plt.plot(n.arange(len(biasc)), biasdrift_full, 'o')
        plt.ylim(-5,5)
        plt.title('Bias Drift Compared to the Average of the First 5 Files')
        plt.xlabel('Bias File number')
        plt.ylabel('Delta_Bias [ADU] (bias - %i)'%baseline)
        plt.savefig(filepath.calibdata+dnight+'/biasdrift_%s.png'%s[0])   
        plt.close('bias')
        
        #calibrate the science images
        file = n.hstack((rawsetp+'zenith1.fit',
                         glob(rawsetp+'ib???.fit'),
                         rawsetp+'zenith2.fit'))
            
        for i in range(len(file)):
            f = fits.open(file[i],uint=False)[0]  # science image 
            f.data -= combias+biasdrift[i]        # subtract drift-corrected bias
            f.data *= n.interp(f.data,xp,fp)      # correct for linearity response
            f.data -= corthermal                  # subtract dark
            f.data /= flat                        # divide by flat
            f.header['IMAGETYP'] = 'CALIB_M'
            f.writeto(calsetp+file[i][len(rawsetp):], overwrite=True)
            T.OpenFile(calsetp+file[i][len(rawsetp):])
            T.SaveFile(calsetp+'tiff/'+file[i][len(rawsetp):-4]+'.tif',5,False,1,0)
            T.Close            #imsave(calsetp+'tiff/'+file[i][len(rawsetp):-4]+'.tif', f.data)
        
        for f in iglob(filepath.tiff+'*.tfw'):
            shutil.copy2(f,calsetp+'tiff/')

    #close MaxIm_DL application
    os.system('taskkill /f /im MaxIm_DL.exe')
                    

def reduceb(dnight, sets, flatname, curve):
    '''
    This module is for calibrating the B-band data. Some of the computation is 
    dependent from the output from the reducev module.
    '''
    #read in the linearity curve (ADU, multiplying factor)
    xp, fp = n.loadtxt(filepath.lincurve+curve+'.txt', unpack=True)
    
    #read in the flat
    flat = fits.open(filepath.flats+flatname,unit=False)[0].data
    
    T = Dispatch('Maxim.Document')

    #looping through all the sets in that night
    for s in sets:
        rawsetp = filepath.rawdata + dnight + '/' + s + '/'
        calsetp = filepath.calibdata + dnight + '/S_0' + s[0] + '/B/'
        print dnight + '/' + s + '...'
        if os.path.isdir(calsetp):
            print 'Replacing old calibrated files...'
        else:
            os.makedirs(calsetp)
            os.makedirs(calsetp+'tiff/')
    
        #read in the thermal, bias, and bias drift from the V band directory
        corthermal = fits.open(calsetp[:-2]+'corthermal.fit')[0].data * 1.5
        combias = fits.open(calsetp[:-2]+'combias.fit')[0].data
        biasdrift = n.loadtxt(filepath.calibdata+dnight+'/biasdrift_%s.txt'%s[0], unpack=True)[1:-1]
        
        #calibrate the science images
        file = glob(rawsetp+'ib???b.fit')
        
        for i in range(len(file)):
            f = fits.open(file[i],uint=False)[0]  # science image
            f.data -= combias+biasdrift[i]        # subtract drift-corrected bias
            f.data *= n.interp(f.data,xp,fp)      # correct for linearity response
            f.data -= corthermal                  # subtract dark
            f.data /= flat                        # divide by flat
            f.header['IMAGETYP'] = 'CALIB_M'
            f.writeto(calsetp+file[i][len(rawsetp):-5]+'.fit', overwrite=True)
            T.OpenFile(calsetp+file[i][len(rawsetp):-5]+'.fit')
            T.SaveFile(calsetp+'tiff/'+file[i][len(rawsetp):-5]+'.tif',5,False,1,0)
            T.Close
            #imsave(calsetp+'tiff/'+file[i][len(rawsetp):-5]+'.tif', f.data)
        
        for f in iglob(filepath.tiff+'*.tfw'):
            shutil.copy2(f,calsetp+'tiff/') 
            
    #close MaxIm_DL application
    os.system('taskkill /f /im MaxIm_DL.exe')
                        
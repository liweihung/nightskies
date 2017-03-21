#-----------------------------------------------------------------------------#
#extinction.py
#
#NPS Night Skies Program
#
#Last updated: 2016/12/02
#
#This script finds the best-fit extinction coefficient and the instrumental
#zeropoint by:
#	(1) identifying the standards stars in the images 
#	(2) measuring their background-subtracted flux in [DN/s]
#	(3) computing their elevations to get airmass
#   (4) comparing the measured flux to their absolute magnitude
#   (5) Given the airmass and M-m for each star, fit for the extinction 
#       coefficient and the instrumental zeropoint
#
#Note: In order to use the ACP objects, the Python must be a 32-bit version.
#
#Input: 
#   (1) Calibrated images
#   (2) hipparcos_standards.txt
#   (3) plot_img number (optional) -- if given, the script will display the 
#       bestfit standard stars contours overlaid on the image data
#
#Output:
#   (1) extinction_stars_%s.txt -- list of the standard stars used for fitting
#   (2) extinction_fit_%s.png -- graphical display of the fitting result
#   (3) extinction_fit.txt -- best-fit extinction coefficient and zeropoint
#
#History:
#	Dan Duriscoe -- Created in visual basic as "extinction v4.vbs"
#	Li-Wei Hung -- Cleaned, improved, and translated to python
#
#-----------------------------------------------------------------------------#

from astropy.io import fits
from glob import glob, iglob
from scipy.optimize import curve_fit
from win32com.client import Dispatch

import pdb
import matplotlib.pyplot as plt
import numpy as n


# Local Source
from gaussian import Gaussian_2d
import filepath
import pointing

#-----------------------------------------------------------------------------#

def plot_fit_result(data, x, y, popt_list):
    '''
    This module provides the visual presentation of the bestfit standard stars
    contours overlaid on the image data.
    '''
    # plot the image data
    plt.close()
    fig = plt.figure(1, figsize=(14, 5))
    ax1 = plt.subplot(121)
    im = ax1.imshow(data, interpolation='nearest', vmin=0, vmax=3000)
    fig.colorbar(im)
    plt.title('real image data')

    # add up all the contours of the standard stars in the best-fit image
    ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
    data_fitted = n.zeros_like(data)
    for i in range(len(popt_list)):
        data_fitted += Gaussian_2d((x, y), *popt_list[i]).reshape(1024, 1024)

    # display the contours on both subplots
    ax1.contour(x, y, data_fitted, 8, colors='w')
    ax2.contour(x, y, data_fitted, 8, colors='k')
    fig.colorbar(im)
    plt.title('best-fit contours over the standard stars')
    plt.show(block=False) 


def extinction(dnight, sets, filter, plot_img=0):
    '''
    This module computes the extinction coefficient and the instrumental zero
    point. It returns the number of stars used for the fit and the location of
    the file containing the best-fit extinction coefficient, zeropoint, and
    their uncertainties. 
    '''
    star = Dispatch('NOVAS.Star')
    site = Dispatch('NOVAS.Site')
    util = Dispatch('ACP.Util')
    p = Dispatch('PinPoint.Plate')
    zeropoint_dnight = []
    
    #read in the standard star catalog
    hips = n.loadtxt(filepath.standards+'hipparcos_standards.txt',dtype=object)
    starn = hips[:,0]                                           #star names
    ras, decs, v_mag, bv = n.array(hips[:,1:],dtype=n.float).T  #star properties
    Mag = {'V':v_mag, 'B':v_mag+bv}                    # absolute mag in V and B
    
    #define image xy coordinates
    x = n.arange(0, 1024)
    y = n.arange(0, 1024)
    x, y = n.meshgrid(x, y)
    
    #parameters specific to datasets with different filters
    F = n.array([('/',),('/B/',)],dtype=[('path','S3'),])
    F = F.view(n.recarray)
    k = {'V':0, 'B':1}
    
    #loop through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata+dnight+'/S_0%s%s' %(s[0],F.path[k[filter]])
        bestfit = []
        xscale = []
        yscale = []

        #read in the header to set the site object's parameter
        H = fits.open(calsetp+'ib001.fit')[0].header
        site.longitude = H['LONGITUD']
        site.latitude = H['LATITUDE']
        site.height = 0
        exp = H['exptime'] #[s]
                
        # loop through each file in the set
        for fn in iglob(calsetp+'ib???.fit'):
            H = fits.open(fn)[0].header
            
            #proceed only if the plate is solved
            try:
                if H['PLTSOLVD']: pass
                else: continue
            except KeyError:
                continue

            #find the standard stars within the 24 X 24 deg image
            p.attachFits(fn)
            img_dec = abs(decs-p.Declination) < 12
            img_ra = abs(ras-p.RightAscension)<(12/(15*n.cos(n.deg2rad(decs))))
            w1 = n.where(img_dec & img_ra)   # stars 
            
            #skip the image w/o standard stars
            if len(w1[0])==0:  
                p.DetachFITS()
                continue
            
            #Get the XY pixel coordinates of the given RA/Dec locations
            def SkyToXY(ra, dec):
                p.SkyToXy(ra, dec)
                return p.ScratchX, p.ScratchY
            px1, py1 = n.array(map(SkyToXY, ras[w1], decs[w1])).T
            p.DetachFITS()            
            
            #find the standard stars within 490 pixels of the image center
            w2 = n.where(n.sqrt((px1-512)**2 + (py1-512)**2) < 490)
            w3 = w1[0][w2]
            px, py = px1[w2], py1[w2]
            hip, ra, dec, M = starn[w3], ras[w3], decs[w3], Mag[filter][w3] 

            data = fits.open(fn)[0].data
            popt_plot_list = []
            
            #info needed for calculating the altitude of the stars later
            JD = H['JD']                               #Julian Date
            TJD = util.Julian_TJD(JD)                  #Terrestrial Julian Date
            LAST = pointing.get_last(JD,H['LONGITUD']) #sidereal time [hr]
            ct = util.Newct(H['LATITUDE'],LAST)
            
            #fit 2D Gaussians to standard stars in the image
            for i in range(len(px)):
                #set the aperture radii
                r = ((x-px[i])**2+(y-py[i])**2)**0.5
                w = n.where(r<3)              #source aperture radiu = 3 pix
                b = n.where((r>4) & (r<8))    #background aperture 4-8 pix ring
                
                #subtract background from the fitted data
                bg = n.median(data[b])
                f = data[w].ravel() - bg
                
                #fit
                guess = (px[i], py[i], 0.6, 50000)  #(x,y,std,brightness)
                try:
                    popt = curve_fit(Gaussian_2d, (x[w],y[w]), f, p0=guess)[0]
                except RuntimeError:
                    continue
                
                #calculate the elevation of the star
                star.RightAscension = ra[i]
                star.Declination = dec[i]
                StarTopo = star.GetTopocentricPosition(TJD, site, False)
                ct.RightAscension = StarTopo.RightAscension
                ct.Declination = StarTopo.Declination
                elev = ct.Elevation       #elevation[deg]
                
                #set the acceptance threshold and record the measurement
                delta_position = n.sum(((popt-guess)**2)[0:2])   #position diff
                signal = popt[3]/bg       #brightness over the background level
                sigma = popt[2]           #sigma of the gaussian

                if sigma<2 and signal>25 and delta_position<1:
                    t = [fn[-7:],hip[i],M[i],elev]
                    t.extend(popt[:3])
                    t.append(popt[3]/H['exptime'])
                    bestfit.append(t)
                    popt_plot_list.append(popt)
            
            #plot the image overlaid with the bestfit contours
            if int(fn[-7:-4]) == plot_img:
                plot_fit_result(data, x, y, popt_plot_list)
                
            #reading in the solved plate scale is x and y image plane
            xscale.append(abs(H['CDELT1']))
            yscale.append(abs(H['CDELT2']))
            
        
        #save the list of stars that will be used for calculating the zeropoint
        stars = n.array((bestfit),dtype=object)
        fmt = ['%7s','%8s','%7.2f','%9.2f','%7.1f','%6.1f','%5.2f','%7.f']
        H = 'File    Star   Magnitude Elevation   X      Y   sigma flux[DN/s]'
        fileout = filepath.calibdata+dnight+'/extinction_stars_%s_%s.txt'\
                  %(filter,s[0])
        n.savetxt(fileout,stars,fmt=fmt,header=H)
        
        #fit for the zeropoint and extinction coefficient
        M = n.float64(stars[:,2])            #V_mag, absolute
        elev = n.float64(stars[:,3])         #elevation[deg]
        flux = n.float64(stars[:,7])         #flux, background subtracted [DN]
        airmass = 1/n.sin(n.deg2rad(elev))   #airmass
        m = -2.5*n.log10(flux)               #v_mag, apparent
        
        p, cov = n.polyfit(airmass, M-m, 1, cov=True)
        c, z = p                             #bestfit coefficient and zeropoint
        c_err, z_err = n.sqrt(cov.diagonal())#uncertainties
        
        sx = n.mean(xscale) * 60             #x plate scale ['/pix]
        sy = n.mean(yscale) * 60             #y plate scale ['/pix]
        sa = n.mean(xscale+yscale) * 60      #average plate scale ['/pix]
        
        fit_entry = [int(s[0]), len(stars), z, z_err, c, c_err, sx, sy, sa, exp]
        zeropoint_dnight.append(fit_entry)
                
        #plot the zeropoint and extinction coefficient fitting result
        x = n.arange(8)
        fig = plt.figure('zeropoint')
        plt.plot(airmass, M-m, 'o', label='Hipparcos standard stars')
        plt.plot(x,c*x+z,'-',lw=2,label='Best fit: %.2fx+%.3f' %(c,z))
        plt.errorbar(0,z,z_err,fmt='o',label='zeropoint: %.3f+-%.3f'%(z,z_err))
        plt.legend(loc=0, numpoints=1)
        plt.xlabel('Airmass',fontsize=14)
        plt.ylabel('M-m',fontsize=14)
        plt.title('Zeropoint and Extinction Coefficient',fontsize=14)
        imgout = filepath.calibdata+dnight+'/extinction_fit_%s_%s.png' \
                 %(filter,s[0])
        plt.savefig(imgout)
        plt.close('zeropoint')
    
    #save the bestfit zeropoint and extinction coefficient     
    fileout = filepath.calibdata+dnight+'/extinction_fit_%s.txt' %filter
    fmt = ['%4i', '%9i', '%13.3f', '%10.3f', '%12.3f', '%11.3f', '%11.3f', 
           '%7.3f', '%11.3f', '%13.1f']
    H1 = "set num_star_used zeropoint zeropoint_err extinction extinction_err "
    H2 = "x_scale y_scale avg_scale['/pix], exptime[s]"
    n.savetxt(fileout,n.array((zeropoint_dnight)),fmt=fmt,header=H1+H2)
    
    return len(stars), fileout






#-----------------------------------------------------------------------------#
#pointing.py
#
#NPS Night Skies Program
#
#Last updated: 2016/11/15
#
#This script calculates the actual pointed azimuth (AZ) and altitude (ALT) using
#the solved RA and Dec values from the image headers. If the images are not 
#solved, the RA, Dec, AZ, and ALT values are interpolated.
#	(1) Read in the solved RA and Dec values from the image header.
#	(2) Update the coordinates to the observed date.
#	(3) Translate to the azimuth and altitude given the LAST and the longitude.
#   (4) Write the output to file
#   (5) Insert the interpolated AZ and ALT in the output file 
#   (6) Update the headers with the interpolated RA and Dec if the images are 
#       not solved.
#
#Note: In order to use the ACP objects, the Python must be a 32-bit version. 
#
#Input: 
#   (1) Calibrated images
#
#Output:
#   (1) pointerr_%s.txt
#
#History:
#	Dan Duriscoe -- Created in visual basic as "calc_pointing_error_v4.vbs"
#	Li-Wei Hung -- Cleaned, improved, and translated to Python
#   Davyd Betchkal -- Plotted the pointing error by image number
#
#-----------------------------------------------------------------------------#

from astropy.coordinates import SkyCoord
from astropy.io import fits
from glob import glob, iglob
from scipy.interpolate import UnivariateSpline
from win32com.client import Dispatch

import pdb
import numpy as n
import matplotlib.pyplot as plt

# Local Source
import filepath  

#-----------------------------------------------------------------------------#
def get_last(JD, longitude):

    '''
    This function calculates the local apparent sidereal time given 
    the Julian Date (JD) and the longitude [deg] of the observing site. 
    This calculation is based on the information from http://aa.usno.navy.mil/faq/docs/GAST.php
    The maximum error resulting from the use of these formulas 
    for sidereal time over the period 2000-2100 is 0.432 seconds.

    Parameters
    ----------
    JD: float, the julian day of observation.
    longitude: float, longitude of the observing site (~~which coordinate system??~~)

	Returns
	-------
	float

    '''

    D = JD - 2451545.0                                   #number of days from 2000 January 1, 12h UT
    GMST = 18.697374558 + 24.06570982441908*D            #Greenwich mean sidereal time [hr]
    
    Omega = n.deg2rad(125.04 - 0.052954*D)               #longitude of the ascending node of the Moon [rad]
    L = n.deg2rad(280.47 + 0.98565*D)                    #mean Longitude of the Sun [rad]
    Del_Phi = -0.000319*n.sin(Omega)-0.000024*n.sin(2*L) #nutation in longitude
    eps = n.deg2rad(23.4393 - 0.0000004*D)               #obliquity [rad]
    eqeq = Del_Phi*n.cos(eps)                            #equation of the equinoxes
    
    GAST = GMST + eqeq                                   #Greenwich apparent sidereal time [hr]
    LAST = (GAST+longitude/360*24)%24                     #local apparent sidereal time [hr]
    
    return LAST

    
def interp_coord(filenames, solved_outputs):
    '''
    Interpolate the True_AZ and True_ALT for images that are not solved and 
    update the RA and DEC with the interpolated values in the header.
    '''
    util = Dispatch('ACP.Util')
    solved, Input_AZ, Input_ALT, True_AZ, True_ALT = solved_outputs
    fi = n.array([int(filenames[i][-7:-4]) for i in range(len(filenames))])
    
    w = [0,15,30,40,45] # the number of last image in every elevation row
    
    for i in range(len(w)-1):
        wf = (fi>w[i]) & (fi<=w[i+1])
        
        if not any(filenames[wf]): continue
        if i==0: #using the second row of image in elevation for interpolation
            wi = (solved>w[i+1]) & (solved<=w[i+2])
            k = min(3, sum(wi)-1)
            A = UnivariateSpline(solved[wi]-15, True_AZ[wi], k=1)
            E = UnivariateSpline(solved[wi]-15, True_ALT[wi]-25, k=k)
        
        else:
            wi = (solved>w[i]) & (solved<=w[i+1])
            k = min(3, sum(wi)-1)
            A = UnivariateSpline(solved[wi], True_AZ[wi], k=k)
            E = UnivariateSpline(solved[wi], True_ALT[wi], k=k)
            
        for fn in filenames[wf]:
            #insert the interpolated Obs_AZ and Obs_ALT 
            f = fits.open(fn, mode='update')
            H = f[0].header
            j = int(fn[-7:-4])
            entry = [j,H['AZ'],H['ALT'],float(A(j)),float(E(j))]
            solved_outputs = n.insert(solved_outputs,j-1,entry,axis=1)
    
            #update the RA and DEC in the header with the interpolated values
            LAST = get_last(H['JD'],H['LONGITUD']) #local apparent sidereal time 
        
            ct = util.Newct(H['LATITUDE'],LAST)
            ct.Azimuth = float(A(j))
            ct.Elevation = float(E(j))
            c = SkyCoord(ct.RightAscension, ct.Declination, unit=('hour','deg'))
        
            f[0].header['RA'] = c.ra.to_string(unit='hour',sep=' ',precision=2)
            f[0].header['DEC'] = c.dec.to_string(unit='deg',sep=' ',precision=1)
            f.flush()
            f.close()
            
    return solved_outputs.T

    
def pointing_err(dnight, sets):
    '''
    This module is calculating the pointing error of each image.
    '''

    star = Dispatch('NOVAS.Star')
    site = Dispatch('NOVAS.Site')
    util = Dispatch('ACP.Util')
    p = Dispatch('PinPoint.Plate')
    
    #looping through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata + dnight + '/S_0' + s[0] + '/'
        
        #read in the header to set the site object's parameter
        H = fits.open(calsetp+'ib001.fit',unit=False)[0].header
        site.longitude = H['LONGITUD']
        site.latitude = H['LATITUDE']
        site.height = 0
        
        #calculate the temperture-pressure correction for refraction
        temp = (H['AMTEMP_F']-32)/1.8 + 273                 #temperature [K]
        pres = (1-(0.0065*H['ELEVATIO']/288.15))**5.3       #pressure [atm]
        tpco = pres*(283/temp)                              #correction
        
        #refraction at 7.5 altitude
        refraction = tpco*(1/n.tan(n.deg2rad(7.5+7.31/11.9)))/60 
        
        #just for V band
        solved=[]; notsolved=[]
        True_AZ=[]; True_ALT=[]; Input_AZ=[]; Input_ALT=[]        
        for fn in iglob(calsetp+'ib???.fit'):
            H = fits.open(fn)[0].header
            
            #calculating the pointing error only if the plate is solved
            if 'PLTSOLVD' not in H or H['PLTSOLVD']==False: 
                notsolved.append(fn)
                continue
            
            solved.append(int(fn[-7:-4]))
            p.attachFits(fn)
            star.RightAscension = p.RightAscension
            star.Declination = p.Declination
            
            JD = H['JD']                #Julian Date
            TJD = util.Julian_TJD(JD)   #Terrestrial Julian Date
            
            #updated star's coordinates at the observed date/time and site
            StarTopo = star.GetTopocentricPosition(TJD, site, False)
            
            #local apparent sidereal time [hr]
            LAST = get_last(JD, H['LONGITUD']) 
            
            #new CoordinateTransform object
            ct = util.Newct(H['LATITUDE'],LAST)
            ct.RightAscension = StarTopo.RightAscension
            ct.Declination = StarTopo.Declination
            
            Input_AZ.append(H['AZ'])
            Input_ALT.append(H['ALT'])
            True_AZ.append(ct.Azimuth)

            #correct for atmospheric refraction on images 1-15
            if int(fn[-7:-4]) < 16: 
                True_ALT.append(ct.Elevation + refraction)
            else:
                True_ALT.append(ct.Elevation)
            
            p.DetachFITS()
                          
        
        #interpolate the True_AZ for True_ALT for images that are not solved
        pterr = n.array([solved,Input_AZ,Input_ALT,True_AZ,True_ALT])
        pterr = interp_coord(n.array(notsolved), pterr)

        # calculate errors
        pErr = pterr.T
        pErr[3][n.where((pErr[1]==0)& (pErr[3]>180))] -= 360
        azmErr = (pErr[1] - pErr[3])*n.cos(n.deg2rad(pErr[4]))
        altErr = pErr[2] - pErr[4]
        totErr = n.sqrt(n.power(azmErr,2) + n.power(altErr,2))

        #create a pointing error plot
        errorPlot = plt.figure(figsize=(20,10))
    	ax = errorPlot.add_subplot(111)
    	plt.suptitle("Pointing Error by Image Number", fontsize=25, verticalalignment='top')
    	plt.title("Data Set " + s[0], fontsize=20)
    	plt.plot(pErr[0], azmErr, linestyle="-.", marker="o", markerfacecolor='None', 
    		markersize=4, color = "darkorange", alpha=0.7, label="Azimuth Error")
    	plt.plot(pErr[0], altErr, linestyle="--", marker="o", 
    		markersize=4, color = "darkgreen", alpha=0.7, label="Altitude Error")
    	plt.plot(pErr[0], totErr, linestyle="-", linewidth=2, marker="o", 
    		markersize=6, color = "black", alpha=1, label="Total Error")
    	plt.axhline(0, color="black", linestyle="-", alpha=0.5, zorder=-10)
    	plt.ylim(-3, 3)
    	plt.ylabel("Error in Degrees", fontsize=20, labelpad = 10)
    	plt.xlabel("Image Number", fontsize=20, labelpad = 15)
    	plt.xticks(n.arange(0, 50, 5))
    	plt.legend(loc='upper left', markerscale=1.8, fontsize=18, framealpha=0.3)
    	ax.tick_params(axis='both', which='major', labelsize=15)
    	plt.text(0.5, -2.8, "Average Total Error:   " + '{:.3f}'.format(totErr.mean()) + u'\N{DEGREE SIGN}', fontsize=18)
    	errorPlot.savefig(filepath.calibdata+dnight+'/pointerr_%s.png' %s[0])

        #saving the output file        
        outfile = filepath.calibdata+dnight+'/pointerr_%s.txt' %s[0]
        nformat = ['%4.f','%8.f','%8.1f','%8.2f','%8.2f']
        H = 'file Input_AZ Input_ALT Obs_AZ Obs_ALT' #column names
        n.savetxt(outfile,pterr,fmt=nformat,header=H)

        


if __name__ == "__main__":
    # pass
    print("Hi, running from the console.")
    pointing_err('FCNA160803', ['1st',])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


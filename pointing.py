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
#	(3) Translate to the azimuth and altitude given the LAST and the logitude.
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
#	Li-Wei Hung -- Cleaned, improved, and translated to python
#
#-----------------------------------------------------------------------------#

from astropy.coordinates import SkyCoord
from astropy.io import fits
from glob import glob, iglob
from scipy.interpolate import interp1d
from win32com.client import Dispatch

import pdb
import numpy as n

# Local Source
import filepath  

#-----------------------------------------------------------------------------#
def get_last(JD, logitude):
    '''
    The function calculates the local apparent sidereal time given the Julian 
    Date (JD) and the logitude [deg] of the observing site. This calculation is
    based on the information on http://aa.usno.navy.mil/faq/docs/GAST.php
    The maximum error resulting from the use of these formulas for sidereal 
    time over the period 2000-2100 is 0.432 seconds.
    '''
    D = JD - 2451545.0                                   #number of days from 2000 January 1, 12h UT
    GMST = 18.697374558 + 24.06570982441908*D            #Greenwich mean sidereal time [hr]
    
    Omega = n.deg2rad(125.04 - 0.052954*D)               #longitude of the ascending node of the Moon [rad]
    L = n.deg2rad(280.47 + 0.98565*D)                    #mean Longitude of the Sun [rad]
    Del_Phi = -0.000319*n.sin(Omega)-0.000024*n.sin(2*L) #nutation in longitude
    eps = n.deg2rad(23.4393 - 0.0000004*D)               #obliquity [rad]
    eqeq = Del_Phi*n.cos(eps)                            #equation of the equinoxes
    
    GAST = GMST + eqeq                                   #Greenwich apparent sidereal time [hr]
    LAST = (GAST+logitude/360*24)%24                     #local apparent sidereal time [hr]
    
    return LAST

    
def interp_coord(filenames, solved_outputs):
    '''
    Interpolate the True_AZ and True_ALT for images that are not solved and 
    update the RA and DEC with the interpolated values in the header.
    '''
    util = Dispatch('ACP.Util')
    solved, Input_AZ, Input_ALT, True_AZ, True_ALT = solved_outputs
    f_az = interp1d(solved, True_AZ-Input_AZ, kind='cubic')
    f_al = interp1d(solved, True_ALT-Input_ALT, kind='cubic')

    for fn in filenames:
        #insert the interpolated Obs_AZ and Obs_ALT 
        f = fits.open(fn, mode='update')
        H = f[0].header
        i = int(fn[-7:-4])
        j = i-1
        entry = [i,H['AZ'],H['ALT'],H['AZ']+f_az(i),H['ALT']+f_al(i)]
        solved_outputs = n.insert(solved_outputs,i-1,entry,axis=1)
    
        #update the RA and DEC in the header with the interpolated values
        LAST = get_last(H['JD'], H['LONGITUD']) #local apparent sidereal time 
        
        ct = util.Newct(H['LATITUDE'],LAST)
        ct.Azimuth = H['AZ']+f_az(i)
        ct.Elevation = H['ALT']+f_al(i)
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
            try:
                if H['PLTSOLVD']: pass
                else: notsolved.append(fn); continue
            except KeyError:
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
            
        #interpretation needs to have boundary values from first and last images
        #if the boundary images are not solved, assgin them the default pointing
        for i in [1, 45]:
            if i not in solved:
                fn = calsetp[:-1]+'\\ib%03i.fit' %i
                H = fits.open(fn)[0].header
                solved.insert(i-1, i)
                Input_AZ.insert(i-1, H['AZ'])
                Input_ALT.insert(i-1, H['ALT'])
                True_AZ.insert(i-1, H['AZ'])
                True_ALT.insert(i-1, H['ALT'])
                notsolved.remove(fn)
              
        
        #interpolate the True_AZ for True_ALT for images that are not solved
        pterr = n.array([solved,Input_AZ,Input_ALT,True_AZ,True_ALT])
        pterr = interp_coord(notsolved, pterr)

        #saving the output file        
        outfile = filepath.calibdata+dnight+'/pointerr_%s.txt' %s[0]
        nformat = ['%4.f','%8.f','%8.1f','%8.2f','%8.2f']
        H = 'file Input_AZ Input_ALT Obs_AZ Obs_ALT' #column names
        n.savetxt(outfile,pterr,fmt=nformat,header=H)

        


if __name__ == "__main__":
    pass
    pointing_err('FCNA160803', ['1st',])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


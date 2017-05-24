#-----------------------------------------------------------------------------#
#cooordinates.py
#
#NPS Night Skies Program
#
#Last updated: 2016/12/16
#
#This script calculates ecliptic and galactic coordinates and rotation angles of
#the images. The output will be used in producing natural sky model for a given 
#location, date, and time.
#
#Input: 
#   (1) Calibrated and solved image headers
#   (2) Obs_AZ and Obs_ALT from pointerr_%s.txt
#
#Output:
#   (1) coordinates_%s.txt with all units in degree
#
#History:
#	Dan Duriscoe -- Created in visual basic "compute image coordinates v4.vbs"
#	Li-Wei Hung -- Cleaned, improved, and translated to python
#
#-----------------------------------------------------------------------------#

from astropy import units as u
from astropy.coordinates import get_sun, SkyCoord
from astropy.io import fits
from astropy.time import Time
from glob import glob, iglob
from win32com.client import Dispatch

import numpy as n
import pdb

# Local Source
import filepath
import pointing

#-----------------------------------------------------------------------------#
def bearing_angle(lat1, lon1, lat2, lon2):
    '''
    Calculate the bearing angle of (lat2, lon2) with respect to (lat1, lon1). 
    The bearing angle ranges from 0 to 360 degrees, with zero at due north 
    and increasing clockwise. Both the inputs and outputs are in degrees. 
    '''
    lat1, lon1, lat2, lon2 = n.deg2rad([lat1, lon1, lat2, lon2])
    x = n.cos(lat1)*n.sin(lat2) - n.sin(lat1)*n.cos(lat2)*n.cos(lon1-lon2)
    y = n.sin(lon1-lon2)*n.cos(lat2)
    bearing = n.rad2deg(n.arctan(y/x))
    if x < 0: bearing += 180
    return -bearing % 360.

    
def galactic_ecliptic_coords(dnight, sets):
    '''
    This module computes the galactic and ecliptic coordinates needed for 
    building the natural sky model. 
    '''
    util = Dispatch('ACP.Util')
    star = Dispatch('NOVAS.Star')
    site = Dispatch('NOVAS.Site')
    site.height = 0
    
    ecliptic_pole = [66.56, 18.]            #N pole Dec [deg] and RA [hr]
    galactic_pole = [27.1283, 167.1405]     #N pole latitude and longitude [deg]
    
    #loop through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata+dnight+'/S_0%s/' %s[0]
        outlist = []
        
        #read in the header to set the site object's parameter
        H = fits.open(calsetp+'ib001.fit')[0].header
        site.longitude = H['LONGITUD']
        site.latitude = H['LATITUDE']
        
        #read in the registered images coordinates
        fnum, Obs_AZ, Obs_ALT = n.loadtxt(filepath.calibdata+dnight+
        '/pointerr_%s.txt' %s[0], usecols=(0,3,4)).T

        # loop through each file in the set
        for fn in iglob(calsetp+'ib???.fit'):
            H = fits.open(fn)[0].header
            w = n.where(fnum==int(fn[-7:-4]))
            
            #new CoordinateTransform object
            JD = H['JD']                               #Julian Date
            TJD = util.Julian_TJD(JD)                  #Terrestrial Julian Date
            
            LAST = pointing.get_last(JD, H['LONGITUD']) 
            ct = util.Newct(H['LATITUDE'],LAST)

            #------------- Calculate the galactic coordinates
            ra, dec = H['RA'], H['DEC']

            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), distance=100*u.kpc)
            g = c.galactic
            galactic_l = (-g.l.value+180)%360-180
            galactic_b = g.b.value
            
            if ('PLTSOLVD' in H.keys()) and H['PLTSOLVD']:   #if plate is solved
                b = bearing_angle(c.dec.degree, -c.ra.degree, *galactic_pole)
                galactic_angle = -((b+H['PA'])%360)
            else: 
                star.RightAscension = galactic_pole[1]/15  #ecliptic N pole [hr]
                star.Declination = galactic_pole[0]       #ecliptic N pole [deg]
                StarTopo = star.GetTopocentricPosition(TJD, site, False)
                ct.RightAscension = StarTopo.RightAscension*15
                ct.Declination = StarTopo.Declination 
                b_in = [Obs_ALT[w][0], Obs_AZ[w][0], ct.Elevation, -ct.Azimuth]
                galactic_angle = -bearing_angle(*b_in)

            #------------- Calculate the ecliptic coordinates
            star.RightAscension = ecliptic_pole[1]         #ecliptic N pole [hr]
            star.Declination = ecliptic_pole[0]           #ecliptic N pole [deg]
            StarTopo = star.GetTopocentricPosition(TJD, site, False)

            ct.RightAscension = StarTopo.RightAscension
            ct.Declination = StarTopo.Declination

            ecliptic_angle = bearing_angle(Obs_ALT[w][0], Obs_AZ[w][0], 
            ct.Elevation, ct.Azimuth)

            csun = get_sun(Time(JD, format='jd'))
            ecliptic_l = -(c.heliocentrictrueecliptic.lon.degree-csun.ra.degree)
            ecliptic_l = (ecliptic_l+180)%360-180
            ecliptic_b = c.heliocentrictrueecliptic.lat.degree
            
            outlist.append([int(fn[-7:-4]), galactic_angle, galactic_l, 
            galactic_b , ecliptic_angle, ecliptic_l, ecliptic_b]) #[deg]
            
        #save the coordinates
        fmt = ['%5i','%8.2f','%8.2f','%8.2f','%8.2f','%8.2f','%8.2f']
        H = 'File  Gal_ang   Gal_l    Gal_b   Ecl_ang   Ecl_l    Ecl_b   [deg]'
        fileout = filepath.calibdata+dnight+'/coordinates_%s.txt'%s[0]
        n.savetxt(fileout,n.array(outlist),fmt=fmt,header=H)


if __name__ == "__main__":
    #galactic_ecliptic_coords('FCNA160803', ['1st',])
    pass

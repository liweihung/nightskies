#-----------------------------------------------------------------------------#
#naturalsky.py
#
#NPS Night Skies Program
#
#Last updated: 2018/01/16
#
#This script...
#	(1) 
#	(2)  
#	(3) 
#   (4) 
#
#Note: 
#
#
#Input: 
#   (1) 
#   (2) 
#   (3) 
#
#Output:
#   (1) 
#   (2) 
#
#History:
#	Dan Duriscoe -- Created as "natskyv4.py"
#	Li-Wei Hung -- Made major modification
#
#-----------------------------------------------------------------------------#

#from astropy.io import fits
#from datetime import datetime as Dtime
#from glob import glob, iglob
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIL import Image
#from scipy.misc import imsave 
from scipy import interpolate
from skimage.transform import downscale_local_mean
#from win32com.client import Dispatch

import pdb
import matplotlib.pyplot as plt
import numpy as n
import os
import scipy
#import shutil

# Local Source
import filepath  

#-----------------------------------------------------------------------------#
#conversion
def nl_to_mag(x):
    """
    Converting brightness from nL to magnitude according to Dan's script
    Note: this is inconsistent to mag_to_nL()
    """
    m = 26.3308 - 2.5 * n.log10(x)
    return m
    
def mag_to_nl_dan(m):
    """
    Converting brightness from magnitude to nL according to Dan's script
    Note: this is inconsistent to nL_to_mag()
    """
    x = 34.08*n.exp(20.7233 - 0.92104 * m)
    return x
    
def mag_to_nl_liwei(m):
    """
    Converting brightness from magnitude to nL according to Dan's script
    Note: this is inconsistent to nL_to_mag()
    """
    x = 10**((26.3308-m)/2.5)
    return x

#galactic model, zodiacal model, and median_filtered image
def get_panoramic_raster(dnight, set, band, raster):
    """    
    This function reads in a raster file and converts it into a Python array. 
    Only area above the horizon is preserved to enforce consistency. 

    Arguments:
    dnight -- data night; i.e. 'FCNA160803'
    set -- data set; i.e. 1 
    band -- filter; either 'V' or 'B'
    raster -- input raster file; either 'gal', 'zod', or 'median'

    Returns:
    A -- a 2D Python array of shape (1800,7200)   
    """
    filter = {'V':"",'B':"/B"}
    path = {'gal':"/gal/galtopmags",
            'zod':"/zod/zodtopmags",
            'median':filter[band]+"/median/skybrightmags"}
            
    import arcpy
    
    file = filepath.griddata+dnight+"/S_0"+str(set)+path[raster]
    arcpy_raster = arcpy.sa.Raster(file)  
    A = arcpy.RasterToNumPyArray(arcpy_raster, "#", "#", "#", -9999)[:1800,:]
    
    return A


#-----------------------------------------------------------------------------#
#------------------        Sky Brightness Models           -------------------#
#-----------------------------------------------------------------------------#

class Model(object):

    def __init__(self, *args, **kwargs):
        self.parameters = {}
        self.parameter_list = self.parameters.keys()
        self.dnight = kwargs.get('dnight', 'FCNA160803') #data night
        self.set = kwargs.get('set', 1)                  #data set
        self.filter = kwargs.get('filter', 'V')          #filter used
        self.elevation = kwargs.get('elevation', 0.)     #site elevation [km]
        self.extinction = kwargs.get('extinction', 0.3)  #extinction coefficient
        self.pixscale = kwargs.get('pixscale', 0.05)     #pixscale [deg/pix]
        self.za_min = kwargs.get('za_min', 0.)           #min zenith angle [deg]
        self.za_max = kwargs.get('za_max', 90.)          #max zenith angle [deg]
        self.get_1d_za()                        #self.za, 1D zenith angles [deg]
        self.compute_airmass()                  #self.airmass, 1D airmass [deg]
        
    def get_1d_za(self,):
        """
        This function generates a one-dimensional zenith angles values from min 
        to max with the resolution of pixscale. The output values are centered 
        in between the sampling point boundaries.
        """
        za_range = n.arange(self.za_min,self.za_max,self.pixscale)
        self.za = za_range + self.pixscale/2
        
    def compute_airmass(self,):
        """    
        This function computes the airmass at the given zenith angles according 
        to Pickering, K. A. (2002). This model has taken the atmospheric 
        refraction into account. See the airmass summary on Wikipedia.
        """
        h = 90.-self.za                                # apparent altitude [deg]
        airmass = 1./n.sin(n.deg2rad(h+244./(165+47*h**1.1)))
        self.airmass = airmass[:,n.newaxis]
        
    def get_parameters(self, parameter_list=None):
        if parameter_list is None:
            parameter_list = self.parameter_list
        return [self.parameters[k] for k in parameter_list]

    def set_parameters(self, p, parameter_list=None):
        if parameter_list is None:
            parameter_list = self.parameter_list
        assert len(p) == len(parameter_list)
        for k, v in zip(parameter_list, p):
            self.parameters[k] = v

    def get_input_model(self,):
        """
        gets the input model with default parameters. 
        """
        raise NotImplementedError

    def compute_observed_model(self,):
        """
        computes 1D or 2D array of brightness model as a function of zenith angle. 
        """
        raise NotImplementedError

    def show_input_model(self,):
        """
        show the input model with default parameters
        """
        raise NotImplementedError
        
    def show_observed_model(self,):
        """
        show the observed model image with current parameters
        """
        raise NotImplementedError


#atmospheric diffused light model
_ADLModelBase = Model
class ADL(_ADLModelBase):

    def __init__(self, *args, **kwargs):
        # call base class constructor
        _ADLModelBase.__init__(self, *args, **kwargs)
    
        self.parameters.update({'a':1.20})
        self.parameter_list = self.parameters.keys()
        self.get_input_model()

    def get_input_model(self,):
        """    
        This function gets the atmospheric diffused light [nL] profile as a 
        function of zenith angle [deg]. The ADL is read and interpolated from 
        Dan Duriscoe's raster file called 'adl_05'. See Duriscoe 2013 and Kwon 
        et al. 2004. 
        """
    
        ADL_file = filepath.rasters+'ADL.txt'
    
        if not os.path.exists(ADL_file):
            import arcpy
            ADLraster = arcpy.sa.Raster(filepath.rasters+"adl_05")  
            gridArray = arcpy.RasterToNumPyArray(ADLraster,"#","#","#",-9999)
            za = n.arange(0,95,0.05) + 0.05/2
            za_adl = n.array([za,gridArray[:,0]]).T
            H = 'Zenith Angle [deg], ADL[nL]'
            n.savetxt(ADL_file,za_adl,fmt=['%10.2f','%15.2f'],header=H)
 
        ZA_pts, ADL_pts = n.loadtxt(ADL_file).T  #zenith angles [deg], ADL [nL]
        f = interpolate.interp1d(ZA_pts, ADL_pts, fill_value='extrapolate')
        self.input_model = f(self.za)       #ADL[nL] at the given zenith angles
        
    def compute_observed_model(self,):
        """
        This function computes the atmospheric diffused light [nL] scaled to the
        factor a at the given zenith angles. 
        """
        return self.parameters['a'] * self.input_model
        
    def show_input_model(self,):
        """
        show the input model with default parameters
        """
        fig, ax = plt.subplots()
        ax.plot(self.za, self.input_model)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlabel('Zenith Angle [Degree]',fontsize=14)
        ax.set_ylabel('ADL Brightness [nL]',fontsize=14)
        fig.canvas.set_window_title("ADL_input_model")
        plt.show(block=False)
        
    def show_observed_model(self,):
        """
        show the observed model image with current parameters
        """
        observed_model = self.compute_observed_model()
        fig, ax = plt.subplots()
        ax.plot(self.za, observed_model)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlabel('Zenith Angle [Degree]',fontsize=14)
        ax.set_ylabel('ADL Brightness [nL]',fontsize=14)
        fig.canvas.set_window_title("ADL_observed_model")
        plt.show(block=False)
        

#airglow model
_AirglowModelBase = Model
class Airglow(_AirglowModelBase):

    def __init__(self, *args, **kwargs):
        # call base class constructor
        _AirglowModelBase.__init__(self, *args, **kwargs)
    
        self.parameters.update({'a':20., #zenith airglow 
                                'h':90., #height of the emitting layer above sea
                                'e':0.6})#airglow extinction factor
        self.parameter_list = self.parameters.keys()
        
    def compute_airglow_brightness(self,a,h):
        """
        This function computes the airglow brightness according to the van Rhijn equation (Leinert et al. 1998).
        
        Arguments:
        self.za -- zenith angles [deg]
        self.elevation -- site elevation [km]
        a -- zenith airglow brightness [any linear brightness unit]
        h -- height of the emitting layer above sea level [km]
        
        Returns:
        flux -- Airglow values with the same unit as the input argument a
        """
        R = 6378+self.elevation #Earth's equatorial radius + site elevation [km]
        H = h-self.elevation   #Height of emitting layer above the observer [km]
        flux = a / n.sqrt(1.-(R*n.sin(n.deg2rad(self.za))/(R+H))**2)  #airglow
        return flux[:,n.newaxis]         

    def get_input_model(self,):
        """
        This function compute the airglow brightness [nL] at the given zenith 
        angles [deg] with default parameters.
        """
        return self.compute_airglow_brightness(20,90) #[nL]
        
    def compute_observed_model(self,):
        """
        This function computes the airglow brightness [nL] with the current
        model parameters. 
        """
        a, h, e = self.parameters['a'],self.parameters['h'],self.parameters['e']
        airglow_nl = self.compute_airglow_brightness(a,h)           #[nL]
        airglow_mag = nl_to_mag(airglow_nl)                         #[mag]
        extinction_total = e*self.extinction*self.airmass           #[mag]
        airglow_obs_mag = extinction_total + airglow_mag            #[mag]
        airglow_obs_nl = mag_to_nl_liwei(airglow_obs_mag)           #[nL]
        return airglow_obs_nl                                       #[nL]

    def show_input_model(self,):
        """
        show the input model with default parameters
        """
        fig, ax = plt.subplots()
        ax.plot(self.za, self.get_input_model())
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlabel('Zenith Angle [Degree]',fontsize=14)
        ax.set_ylabel('Ariglow Brightness [nL]',fontsize=14)
        fig.canvas.set_window_title("Airglow_input_model")
        plt.show(block=False)
        
    def show_observed_model(self,):
        """
        show the observed model with the current parameters
        """
        fig, ax = plt.subplots()
        ax.plot(self.za, self.compute_observed_model())
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlabel('Zenith Angle [Degree]',fontsize=14)
        ax.set_ylabel('Ariglow Brightness [nL]',fontsize=14)
        fig.canvas.set_window_title("Airglow_observed_model")
        plt.show(block=False)
                
        
#zodiacal light model
_ZodiacalModelBase = Model
class Zodiacal(_ZodiacalModelBase):

    def __init__(self, *args, **kwargs):
        # call base class constructor
        _ZodiacalModelBase.__init__(self, *args, **kwargs)
    
        self.parameters.update({'e':0.6,}) #Zodiacal light extinction factor 
        self.parameter_list = self.parameters.keys()
        self.get_input_model()
    
    def get_input_model(self,):
        """
        This function reads in the zodiacal light model [mag] specific for this
        set of data. 
        """
        d, s, f = self.dnight, self.set, self.filter
        self.input_model = get_panoramic_raster(d, s, f, 'zod')
        
    def compute_observed_model(self,):
        """
        This function computes the observed brightness model of zodiacal light 
        [mag] with the current parameters. 
        """
        extinction_total = self.parameters['e']*self.extinction*self.airmass
        return self.input_model + extinction_total
        
    def image_template(self, image, title):
        fig = plt.figure(figsize=(12,3.4))
        im = plt.imshow(image, cmap='gist_heat_r', extent=(-180,180,0,90), 
                        vmax=25, vmin=21)
        plt.xticks(n.arange(-180, 181, 30))
        plt.yticks(n.arange(0, 91, 15))
        plt.xlabel('Azimuth (degree)')
        plt.ylabel('Altitude (degree)')
        
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="1.5%", pad="2%")
        cbar = plt.colorbar(im, cax=cax, ticks=n.arange(21,26,1))
        cbar.ax.invert_yaxis()
        cbar.ax.set_yticklabels(['<21','22','23','24','>25'])
        cbar.set_label(r'mag / arcsec$^2$')
        
        fig.canvas.set_window_title(title)

        plt.tight_layout()
        plt.show(block=False)
        
    def show_input_model(self,):
        """
        show the input model with default parameters
        """
        small_img = downscale_local_mean(self.input_model,(10,10))
        self.image_template(small_img, "Zodiacal_light_input_model")

    def show_observed_model(self,):
        """
        show the observed model with the current parameters
        """
        small_img = downscale_local_mean(self.compute_observed_model(),(10,10))
        self.image_template(small_img, "Zodiacal_light_obsersved_model")


#galactic model
_GalacticModelBase = Model 
class Galactic(_GalacticModelBase):

    def __init__(self, *args, **kwargs):
        # call base class constructor
        _GalacticModelBase.__init__(self, *args, **kwargs)
    
        self.parameters.update({'e':0.9,}) #Galactic light extinction factor 
        self.parameter_list = self.parameters.keys()
        self.get_input_model()
        
    def get_input_model(self,):
        """
        This function reads in the galactic light model [mag] specific for this
        set of data. 
        """
        d, s, f = self.dnight, self.set, self.filter
        self.input_model = get_panoramic_raster(d, s, f, 'gal')
        
    def compute_observed_model(self,):
        """
        This function computes the observed brightness model of galactic light 
        [mag] with the current parameters. 
        """
        extinction_total = self.parameters['e']*self.extinction*self.airmass
        return self.input_model + extinction_total
        
    def image_template(self, image, title):
        fig = plt.figure(figsize=(12,3.4))
        im = plt.imshow(image, cmap='Blues', extent=(-180,180,0,90), 
                        vmax=25, vmin=20) #or gist_gray_r
        plt.xticks(n.arange(-180, 181, 30))
        plt.yticks(n.arange(0, 91, 15))
        plt.xlabel('Azimuth (degree)')
        plt.ylabel('Altitude (degree)')
        
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="1.5%", pad="2%")
        cbar = plt.colorbar(im, cax=cax,  ticks=n.arange(20,26,1))
        cbar.ax.invert_yaxis()
        cbar.ax.set_yticklabels(['<20','21','22','23','24','>25'])
        cbar.set_label(r'mag / arcsec$^2$')
        
        fig.canvas.set_window_title(title)

        plt.tight_layout()
        plt.show(block=False)
        
    def show_input_model(self,):
        """
        show the input model with default parameters
        """
        small_img = downscale_local_mean(self.input_model,(10,10))
        self.image_template(small_img, "Galactic_light_input_model")

    def show_observed_model(self,):
        """
        show the observed model with the current parameters
        """
        small_img = downscale_local_mean(self.compute_observed_model(),(10,10))
        self.image_template(small_img, "Galactic_light_obsersved_model")


#mask        
_MaskBase = Model
class Mask(_MaskBase):

    def __init__(self, *args, **kwargs):
        # call base class constructor
        _MaskBase.__init__(self, *args, **kwargs)
    
        self.get_input_model()
        
    def get_input_model(self,):
        """
        This reads in the mask to be used in the fitting process. This mask only
        extends from the zenith to the horizon. Sky = 1; terrain = 0.
        """
        im = Image.open(filepath.griddata+self.dnight+'/mask.tif')
        mask_tif = n.array(im)[:1800]  #only selecting pixels above the horizon
        sky = n.where(mask_tif==65535) #sky pixels
    
        mask = n.zeros_like(mask_tif)  #initialize mask to be opaque
        mask[sky] = 1
    
        self.input_model = mask
        
    def show_input_model(self,):
        """
        show the input mask
        """
        small_img = downscale_local_mean(self.input_model,(10,10))

        fig = plt.figure(figsize=(12,3.4))
        im = plt.imshow(small_img, cmap='binary_r', extent=(-180,180,0,90), 
                        vmax=1, vmin=0)
        plt.xticks(n.arange(-180, 181, 30))
        plt.yticks(n.arange(0, 91, 15))
        plt.xlabel('Azimuth (degree)')
        plt.ylabel('Altitude (degree)')
        
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="1.5%", pad="2%")
        cbar = plt.colorbar(im, cax=cax,  ticks=[0,1])
        cbar.ax.set_yticklabels(['terrain','sky'])
        
        fig.canvas.set_window_title("Terrain_mask")

        plt.tight_layout()
        plt.show(block=False)
############ Start coding from here 

#--------------------------------------------------------------------------- 

A = Airglow()
A.show_input_model()
A.show_observed_model()

D = ADL()
D.show_input_model()
D.show_observed_model()

Z = Zodiacal()
Z.show_input_model()
Z.show_observed_model()

G = Galactic()
G.show_input_model()
G.show_observed_model()

M = Mask()
M.show_input_model()

#mask = get_mask('FCNA160803/mask.tif')
#small_img = downscale_local_mean(mask,(10,10))
#plt.imshow(small_img)
#plt.show(block=False)

#-----------------------------------------------------------------------------#



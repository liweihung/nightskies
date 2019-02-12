#------------------------------------------------------------------------------#
#naturalsky.py
#
#NPS Night Skies Program
#
#Last updated: 2018/07/17
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
#	Li-Wei Hung -- Made major modification and used object-oriented approach
#
#------------------------------------------------------------------------------#

from astropy.io import fits
#from datetime import datetime as Dtime
#from glob import glob, iglob
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIL import Image
#from scipy.misc import imsave 
from scipy import interpolate
from skimage.transform import downscale_local_mean
#from win32com.client import Dispatch

import pdb
import matplotlib
import matplotlib.pyplot as plt
import numpy as n
import os
import scipy
#import shutil
import sys

# Local Source
import filepath  
import colormaps

#------------------------------------------------------------------------------#
#-------------------           Various Functions            -------------------#
#------------------------------------------------------------------------------#

#unit conversion
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
    Converting brightness from magnitude to nL according to nl_to_mag
    """
    x = 10**((26.3308-m)/2.5)
    return x

#read in galactic model, zodiacal model, and median_filtered images
def get_panoramic_raster(dnight, set, band, raster, k=25):
    """    
    This function reads in a raster file and converts it into a Python array. 
    Only area above the horizon is preserved to enforce consistency. This 
    function does not work on full resolution mosaics due limited memory space.

    Arguments:
    dnight -- data night; i.e. 'FCNA160803'
    set -- data set; i.e. 1 
    band -- filter; either 'V' or 'B'
    raster -- input raster; either 'gal', 'zod', or 'median' but not 'fullres'
    k -- downscale factor for making the image smaller

    Returns:
    A -- a 2D Python array of shape [1800,7200]/k   
    """
    filter = {'V':"",'B':"/B"}
    path = {'gal':"/gal/galtopmags",
            'zod':"/zod/zodtopmags",
            'median':filter[band]+"/median/skybrightmags"}
            
    import arcpy
    file = filepath.griddata+dnight+"/S_0"+str(set)+path[raster]
    arcpy_raster = arcpy.sa.Raster(file)  
    A = arcpy.RasterToNumPyArray(arcpy_raster, "#", "#", "#", -9999)[:1800,:7200]
    A_small = downscale_local_mean(A,(k,k))
    return A_small

#read in downscaled models and data from the fits images
def get_downscaled_image(dnight, set, band, image):
    """    
    This function reads in a downscaled fits image. Only area above the horizon 
    is preserved.
    
    Arguments:
    dnight -- data night; i.e. 'FCNA160803'
    set -- data set; i.e. 1 
    band -- filter; either 'V' or 'B'
    image -- input image; either 'gal', 'zod', or 'median' but not 'fullres'

    Returns:
    A_small -- a 2D Python array of shape [72,288]
    """
    filter = {'V':"",'B':"/B"};
    f = {'V':'', 'B':'b'}
    path = {'gal':"/galtopmags%s.fits" %set,
            'zod':"/zodtopmags%s.fits" %set,
            'median':filter[band]+"/skybrightmags%s%s.fits"%(f[band],set)}
            
    file = filepath.griddata+dnight+path[image]
    A_small = fits.open(file,unit=False)[0].data
    return A_small
    
    
#------------------------------------------------------------------------------#
#-------------------    Sky Brightness Model Components     -------------------#
#------------------------------------------------------------------------------#

#base model 
class Model(object):

    def __init__(self, dnight, set, filter, **kwargs):
        self.parameters = {}
        self.parameter_list = self.parameters.keys()
        self.dnight = dnight                           #data night
        self.set = set                                 #data set
        self.filter = filter                           #filter used
        self.pixscale = kwargs.get('pixscale', 0.05)   #pixscale [deg/pix]
        self.downscale = kwargs.get('downscale', 25)   #downscale factor
        self.za_min = kwargs.get('za_min', 0.)         #min zenith angle [deg]
        self.za_max = kwargs.get('za_max', 90.)        #max zenith angle [deg]
        self.mask = kwargs.get('mask', n.array([0,]))  #terrain mask 
        self.get_extinction_coefficient()              #extinction        
        self.get_1d_za()                               #zenith angles 1D [deg]
        self.compute_airmass()                         #airmass 1D [deg]
        

    def get_extinction_coefficient(self,):
        """
        This function reads in the extinction coefficient associated with the 
        data set. 
        """
        F = {'V':'/','B':'/B/'}
        extinctionfile = filepath.calibdata + self.dnight + \
                         '/extinction_fit_%s.txt' %self.filter
        self.extinction = abs(n.loadtxt(extinctionfile, ndmin=2)[self.set-1,4])
        
        
    def get_1d_za(self,):
        """
        This function generates a one-dimensional zenith angles values from min 
        to max with the resolution of pixscale. The output values are centered 
        in between the sampling point boundaries.
        """
        za = n.arange(self.za_min,self.za_max,self.pixscale*self.downscale)
        self.za = za + self.pixscale*self.downscale/2
        
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
        
    def image_template(self, image, title, mask=False, cmapname='NPS_mag', 
                       min=14, max=24):
        plt.rcParams['image.cmap'] = cmapname
        matplotlib.cm.get_cmap().set_bad(color='black') #mask color
        if (self.mask>0).any() and (mask==True):
            image = n.ma.MaskedArray(image,self.mask)
        fig = plt.figure(figsize=(12,3.4))
        im = plt.imshow(image, extent=(-180,180,0,90), vmin=min, vmax=max) 
        plt.xticks(n.arange(-180, 181, 30))
        plt.yticks(n.arange(0, 91, 15))
        plt.xlabel('Azimuth (degree)')
        plt.ylabel('Altitude (degree)')

        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="1.5%", pad="2%")
        cbar = plt.colorbar(im, cax=cax,  ticks=n.arange(min,max+1))
        cbar.ax.invert_yaxis()
        
        labels = map(str, n.arange(min,max+1))
        labels[0] = '<'+labels[0]
        labels[-1] = '>'+labels[-1]
        
        cbar.ax.set_yticklabels(labels)
        cbar.set_label(r'mag / arcsec$^2$')
        
        fig.canvas.set_window_title(title)

        plt.tight_layout()
        plt.show(block=False)    

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
        
        self.get_site_elevation()
        
    def get_site_elevation(self,):
        """
        This function reads in the elevation of the observing site from the first image in the data set.
        """
        F = {'V':'/','B':'/B/'}
        headerfile = filepath.calibdata + self.dnight + '/S_0' + \
                     str(self.set) + F[self.filter] + 'ib001.fit'
        self.elevation = fits.open(headerfile)[0].header['ELEVATIO']/1000. #[km]
        
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
        self.input_model = f(self.za)[:,n.newaxis] #ADL[nL] at the zenith angles
        
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
        self.input_model = get_downscaled_image(d, s, f, 'gal')
        
    def compute_observed_model(self,unit='nl'):
        """
        This function computes the observed brightness model of galactic light 
        [mag] with the current parameters. 
        """
        extinction_total = self.parameters['e']*self.extinction*self.airmass
        if unit=='mag':
            return self.input_model + extinction_total
        else:
            return mag_to_nl_liwei(self.input_model + extinction_total)
            
    def show_input_model(self,):
        """
        show the input model with default parameters
        """
        self.image_template(self.input_model, "Galactic_light_input_model")

    def show_observed_model(self,):
        """
        show the observed model with the current parameters
        """
        img = self.compute_observed_model(unit='mag')
        self.image_template(img, "Galactic_light_obsersved_model", mask=True)

                            
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
        self.input_model = get_downscaled_image(d, s, f, 'zod')
        
    def compute_observed_model(self, unit='nl'):
        """
        This function computes the observed brightness model of zodiacal light 
        [mag] with the current parameters. 
        """
        extinction_total = self.parameters['e']*self.extinction*self.airmass
        if unit=='mag':
            return self.input_model + extinction_total
        else:
            return mag_to_nl_liwei(self.input_model + extinction_total)
        
    def show_input_model(self,):
        """
        show the input model with default parameters
        """
        self.image_template(self.input_model, "Zodiacal_light_input_model")

    def show_observed_model(self,):
        """
        show the observed model with the current parameters
        """
        img = self.compute_observed_model(unit='mag')
        self.image_template(img,"Zodiacal_light_obsersved_model", mask=True)
        

#------------------------------------------------------------------------------#
#-------------------              Terrain Mask              -------------------#
#------------------------------------------------------------------------------#

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
        extends from the zenith to the horizon. Sky = 0; terrain = 1.
        """
        im = Image.open(filepath.griddata+self.dnight+'/mask.tif')
        mask_tif = n.array(im)[:1800]  #only selecting pixels above the horizon
        
        mask = downscale_local_mean(mask_tif,(self.downscale,self.downscale))
        terrain = n.where(mask!=65535) #terrain pixels
        mask[:,:] = 0                  #initialize mask to sky (zero)
        mask[terrain] = 1              #set terrain pixels to one
        self.input_model = mask
        
    def show_input_model(self,):
        """
        show the input mask
        """
        self.image_template(self.input_model, "Terrain_mask", cmapname='binary',
                            min=0, max=1)

        
#------------------------------------------------------------------------------#
#-------------------             Combine Models             -------------------#
#------------------------------------------------------------------------------# 

_AggregateModelBase = Model
class AggregateModel(_AggregateModelBase):
    """
    This class is for combining the input models 
    - "fix_param" will keep parameters fixed during the fitting process. It
      takes the form {model:[fixed params]}
    """

    def __init__(self, model_list, *args, **kwargs):
        # call base class constructor
        _AggregateModelBase.__init__(self, *args, **kwargs)

        self.model_list = model_list
        self.fix_param = kwargs.get('fix_param',{}) #{model:[fixed params]}

    def get_parameters(self, parameter_list=None):
        """
        Get the values of floating parameters from all the input models. 
        """
        assert parameter_list is None

        all_parms = []
        for model in self.model_list:
            model_params = model.parameters.copy()
            if model in self.fix_param:
                [model_params.pop(k) for k in self.fix_param[model]]

            floating_params = model_params.values()
            all_parms.extend(floating_params)
        return all_parms

    def nfree(self,):
        """
        Returns the number of free parameter.
        """
        n_free = len(self.get_parameters())
        return n_free

    def set_parameters(self, p, parameter_list=None):
        assert parameter_list is None

        # assert that the number of parameters is correct
        assert len(p) == self.nfree(),"the number of the parameter given does \
                                      not match the number of free parameters"

        # set the free parameters to the new values 
        i = 0
        for model in self.model_list:
            model_params = model.parameters.copy()
            if model in self.fix_param:
                [model_params.pop(k) for k in self.fix_param[model]]
            n_parms = len(model_params)
            parms = p[i:i+n_parms]
            model.set_parameters(parms, parameter_list=model_params.keys())
            i += n_parms
    
    def compute_observed_model(self, unit='nl'):
        """
        This function computes the combined brightness from all the input models  
        """
        im = self.model_list[0].compute_observed_model()
        for model in self.model_list[1:]:
            im += model.compute_observed_model()
        if unit=='mag':
            return nl_to_mag(im)
        else:
            return im
            
    def show_observed_model(self,):
        """
        show the observed model with the current parameters
        """
        img = self.compute_observed_model(unit='mag')
        self.image_template(img, "Natural_sky_model", mask=True)
        
        
#------------------------------------------------------------------------------#
#-----------  Calculate Difference between the Model and the Data  ------------#
#------------------------------------------------------------------------------# 

class ModelComparator(object):

    def __init__(self, model, *args, **kwargs):
        """image format: {lam: data}, where data[0] = image and data[1] = sigma
           sed format: array([lam, flux, sigma])"""

        self.model = model
        self.image = kwargs.get('image', n.array(()))
        
    def compute_difference(self, limitedkey=[], return_2D=False):
        """
        compute the difference in flux between the data and the model. 
        limitedkey is a list of parameters that are constrained to be positive
        """

        #----limit the range of the values of parameters in limitedkey
        badp = 0
        for modeli in self.model.model_list:
            if any([modeli.parameters[key]<0. for key in limitedkey if key in \
                    modeli.parameters]):
                badp = 1

        #----Compute the image error
        if badp: 
            artificial_light_nl = 1000*n.ones_like(self.image)
        else:
            #assume the uncertainty is the same across the pixels
            data_nl = mag_to_nl_liwei(self.image)
            natural_sky_nl = self.model.compute_observed_model() 
            artificial_light_nl = data_nl - natural_sky_nl
            assert self.model.mask.any()     #assert the mask is applied

        if return_2D:
            return artificial_light_nl
        else:
            return n.ravel(artificial_light_nl)

    def chi2(self,):
        imgerr = self.compute_difference()
        chi2img = n.sum(imgerr**2)
        nfree = self.model.nfree()
        chi2aimg = chi2img/(len(imgerr)-nfree)
        print "Chi^2: %i" %chi2img
        print "Reduced Chi^2: %.2f" %chi2aimg

    def showimg(self, p=None): 
        if p != None:
            self.model.set_parameters(p)
            
        natural_sky = self.model.compute_observed_model(unit='mag')
        
        fig = plt.figure(self.model.dnight+' Set %i'%self.model.set, 
                         figsize=(11.7,9), dpi=80)
        fig.clear()
        plt.rcParams['image.cmap'] = 'NPS_mag'
        fsize = 12
     
        #------------ (a)Data: observed
        ax = fig.add_subplot(311)
        plt.imshow(self.image, extent=(-180,180,0,90), vmin=14, vmax=24)
        ax.tick_params(axis='both', labelsize=fsize)
        ax.set_ylabel('Altitude (degree)', fontsize=fsize)

        #------------ (b)Model: natural sky model
        ax = fig.add_subplot(312, sharex=ax)
        plt.imshow(natural_sky, extent=(-180,180,0,90), vmin=14, vmax=24)
        ax.tick_params(axis='both', labelsize=fsize)
        ax.set_ylabel('Altitude (degree)', fontsize=fsize)

        #------------ (c)Residual: light pollution only
        ax = fig.add_subplot(313, sharex=ax)
        artificial_light = nl_to_mag(self.compute_difference(return_2D=True))
        print n.sum(artificial_light**2)
        im = plt.imshow(artificial_light, extent=(-180,180,0,90), vmin=14, vmax=24)
        ax.tick_params(axis='both', labelsize=fsize)
        ax.set_ylabel('Altitude (degree)', fontsize=fsize)
        ax.set_xlabel('Azimuth (degree)', fontsize=fsize)

        plt.tight_layout(rect=[0.02,0,0.95,1], h_pad=0.1)
        cbar_ax = fig.add_axes([0.93, 0.072, 0.01, 0.907])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label(r'mag / arcsec$^2$',fontsize=fsize)
        plt.show(block=False)

    '''
    def showparams(self, fitoutput):
        p = fitoutput[0]
        ps = fitoutput[0].astype(str)
        if type(fitoutput[1]) == type(None):
            print "The covariance matrix is 'None'."
            for i in range(len(p)):
                ps[i] = format(p[i],'.3g')
        else: 
            p_err = n.sqrt(n.diagonal(fitoutput[1]))
            for i in range(len(p)):
                ps[i] = format(p[i],'.3g') + ' +- ' + format(p_err[i],'.3g')
        self.model.set_parameters(p)
        self.model.set_parameters(ps)
        unit = {'a0':'[um]','amax':'[um]','amin':'[um]','b':'[]',
                'diskcen':'[pix]','distance':'[pc]','fg':'[]','fi':'[]',
                'fo':'[]','fp':'[]','fs':'[]','Fstar':'[erg/s/cm^2]',
                'gamma':'[]','inc':'[deg]','k':'[]','L':'[ergs/sec]',
                'n0':'[#/um^2(/um)]','p':'[]','PA':'[deg]','radius':'[AU]',
                'ri':'[AU]','ro':'[AU]','starpos':'[pix]','T':'[K]'}
        unit = defaultdict(lambda:'[?]', unit)
        for m in self.model.model_list:
            print 'model', m.__class__
            for key in m.parameters.keys():
                print key,':',m.parameters[key],unit[key]
        self.model.set_parameters(p)
    '''

       
#------------------------------------------------------------------------------#
dnight = 'FCNA160803' #data night
set = 1               #data set
filter = 'V'          #filter used

Pa = [dnight, set, filter]
Pk = {'pixscale':0.05, #unit?
      'downscale':25, 
      'za_min':0., 
      'za_max':90.}

K = Mask(*Pa, **Pk)
#print(K.parameters)
#K.show_input_model()

Pk['mask'] = K.input_model

A = Airglow(*Pa, **Pk)
#A.parameters['a'] = 30.
print('Airglow: ', A.parameters)
#A.show_input_model()
#A.show_observed_model()

D = ADL(*Pa, **Pk)
print('ADL: ', D.parameters)
#D.show_input_model()
#D.show_observed_model()

G = Galactic(*Pa, **Pk)
#G.parameters['e'] = 1.
print('Galactic: ', G.parameters)
G.show_input_model()
G.show_observed_model()

Z = Zodiacal(*Pa, **Pk)
#Z.parameters['e'] = 0.1
print('Zodiacal: ', Z.parameters)
Z.show_input_model()
Z.show_observed_model()

M = AggregateModel([G,Z,A,D],*Pa,**Pk)
#M.show_observed_model()


S = get_downscaled_image(K.dnight, K.set, K.filter, 'median')

C = ModelComparator(M, image=S)
C.showimg()
#K.image_template(S, "Median Filtered Data")
'''
D_nl = mag_to_nl_liwei(S) - M.compute_observed_model()
D_mag = nl_to_mag(D_nl)
D_mag[n.isnan(D_mag)] = 30
M.image_template(D_mag, "Difference")
'''

#small_mask = downscale_local_mean(M.input_model,(25,25))
#small_img = downscale_local_mean(S-M.compute_observed_model(),(25,25))
#masked_img = n.ma.masked_array(small_img, mask=small_mask)
#plt.imshow(S,interpolation='nearest')
#plt.imshow(masked_img, cmap='binary')
#cbar = plt.colorbar()
#plt.show(block=False)

#-----------------------------------------------------------------------------#



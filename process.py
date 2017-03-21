#-----------------------------------------------------------------------------#
#process.py
#
#NPS Night Skies Program
#
#Last updated: 2017/03/12
#
#This script uses multiprocessing to reduce, calibrate, and process the images 
#collected by the NPS Night Skies Program. This script does the following:
#	(1) Reduction --- bias, dark, flat, and linearity response correction
#	(2) Registration --- image pointing determination
#	(3) Photometric Calibration --- fitting for the zeropoint and extinction
#	(4) Star Removal --- through median filtering
#	(5) Mosaicing --- producing panoramic images of the galactic model, zodiacal
#       model, median-filtered data, and full-resolution data
#
#
#Note:  
#   Importing arcpy will alter some default function name spaces (ie, time and
#   datetime). Thus, be aware of where arcpy is imported and the potential 
#   namespace conflicts. Here, arcpy is imported through importing a local 
#   mosaic source file. 
#
#
#Input:
#   (1) filelist.txt - input list of file folders to be processed and the 
#       associated calibration files. 
#
#
#Output:
#   (1) Calibdata folder (calibrated images, processlog.txt, and other outputs) 
#   (2) Griddata folder (4 photometrically calibrated panoramic images: galactic 
#       model, zodiacal model, median-filtered data, and full-resolution data)
#
#
#History:
#	Dan Duriscoe -- Created (named firstbatchv4vb.py originally)
#	Li-Wei Hung -- Cleaned, modified, and changed some processing methods 
#
#-----------------------------------------------------------------------------#
print ' '
print '--------------------------------------------------------------'
print ' '
print '        NPS NIGHT SKIES PROGRAM RAW IMAGE PROCESSING'
print ' '
print '--------------------------------------------------------------'

print ' '
#-----------------------------------------------------------------------------#

#Standard libraries
import matplotlib.pyplot as plt
#import multiprocessing
import numpy as n
import os
import pdb
import time
#from astropy.io import fits
from datetime import datetime as Dtime

#-----------------------------------------------------------------------------#

##########################  Definitions  ######################################
def loghistory(string):
    '''To log the processing history and print status on the screen'''
    if type(string) is list:
        for s in string:
            history.write(s+'\n')
            print s
    else:
        history.write(string+'\n')
        print string    


def log_inputs(p):
    '''starting the history file and log the input files used'''
    print 'Recording the reduction process in: '
    print filepath.calibdata+dnight+'/processlog.txt'
    print ' '
    loghistory('%s processed on %s by %s' %(p[0],str(Dtime.now())[:19],p[6]))
    loghistory('')
    loghistory('Input files:')
    loghistory('Reducing V-band data? ' + p[1])
    loghistory('Reducing B-band data? ' + p[2])
    loghistory('vflat = ' + p[3])
    loghistory('bflat = ' + p[4])
    loghistory('curve = ' + p[5])
    loghistory('')
    loghistory('____________________ Processing history: _____________________')
    loghistory('')
   

def update_progressbar(x,y,t):
    '''update the progress bar plot'''
    if t < 1: texty = '%.1f' %t
    else: texty = '%d' %round(t)
    if t < 6: c = 'k'
    else: c = 'w'
    barax.pcolor([x,x+1], [y+5,y+6], [[t],], cmap='Greens', vmin=0, vmax=10)
    barax.text(x+0.5, y+5.5, texty, color=c, 
               horizontalalignment='center',
               verticalalignment='center', 
               size='medium')
    plt.pause(0.05)  #draw the new data and run the GUI's event loop
    Z[y+5,x] = t   #recording the time in a master array


def timing_function(some_function):
    '''Outputs the time a function takes to execute.'''
    def wrapper(*args, **kwargs):
        t1 = time.time()
        j = function_number[some_function.__name__]
        barax.pcolor([j,j+1], [i+5,i+6], [[4],], cmap='gray', vmin=0, vmax=5)
        plt.pause(0.05)
        some_function(*args, **kwargs)
        t2 = time.time()
        dt = round((t2-t1)/60,1) #minute
        update_progressbar(j,i,dt)
        loghistory('Time it took to run this function: %s min \n' %str(dt))
    return wrapper
    
@timing_function    
def reduce_images(*args):
    '''Basic image reduction'''
    loghistory('--------------------------------------------------------------')
    if 'V' in args[2]:
        #V-band Basic image reduction
        loghistory('Calibrating V-band images ('+R.__name__+'.reducev)')
        R.reducev(args[0],args[1],args[2]['V'],args[3])
        loghistory('Finished calibrating images in V-band successfully.')
        loghistory('')
    if 'B' in args[2]:
        #B-band Basic image reduction
        loghistory('Calibrating B-band images ('+R.__name__+'.reduceb)')
        R.reduceb(args[0],args[1],args[2]['B'],args[3])
        loghistory('Finished calibrating images in B-band successfully.')
        loghistory('')
    

@timing_function  
def register_coord(*args):
    '''Registering coordinates of each image by matching up star positions'''
    loghistory('--------------------------------------------------------------')
    loghistory('Registering the image coordinates (%s.matchstars)' %RE.__name__)
    loghistory('')
    for filter in args[2]:
        loghistory('%s band' %filter)
        arg = list(args)
        arg[2] = filter
        cropped_fn, failed_fn = RE.matchstars(*arg)
        a = 'Images are registered using full frames unless noted otherwise.'
        b = 'These images are registered using the central 200x200 pixels only:'
        loghistory([a,'',b])
        loghistory(cropped_fn)
        loghistory('')
        loghistory('The following images are failed to be registered:')
        loghistory(failed_fn)
        loghistory('')


@timing_function
def calculate_pointing_err(*args):
    '''Calculating the pointing error'''
    loghistory('--------------------------------------------------------------')
    loghistory('Calculating the pointing error ( %s.pointing_err)' \
                %pointing.__name__)
    loghistory('')
    pointing.pointing_err(*args)
    loghistory('Output file pointerr.txt')
    loghistory('')

    
@timing_function
def fit_zeropoint(*args, **kwargs):
    '''Fitting for the extinction coefficient and zeropoint'''
    loghistory('--------------------------------------------------------------')
    loghistory('Fitting for extinction and zeropoint (%s.extinction)' 
               %EX.__name__)
    loghistory('')
    for filter in args[2]:
        loghistory('%s band' %filter)
        arg = list(args)
        arg[2] = filter
        nstar, bestfit_file = EX.extinction(*arg, **kwargs)
        bestfitp = list(n.loadtxt(bestfit_file, dtype=str, comments='', 
                        delimiter='!'))
        loghistory('The best fit parameters are: ')
        loghistory(bestfitp)
        loghistory('')


@timing_function
def compute_coord(*args):
    '''Computig the galactic and ecliptic coordinates of the images'''
    loghistory('--------------------------------------------------------------')
    a = 'Computing galactic and ecliptic coordinates'
    b = '(%s.galactic_ecliptic_coords)' %(CO.__name__)
    loghistory(a+b)
    loghistory('')
    CO.galactic_ecliptic_coords(*args)
    loghistory('Output file coordinate_#.txt')
    loghistory('')


@timing_function
def apply_filter(*args):
    '''Apply ~1 degree (in diameter) median filter to the images'''
    loghistory('--------------------------------------------------------------')
    loghistory('Applying the median filter to images (%s.filter)' %MF.__name__)
    loghistory('')
    for filter in args[2]:
        loghistory('%s band' %filter)
        arg = list(args)
        arg[2] = filter
        MF.filter(*arg)
        loghistory('Output median filtered images tiff/median_*.tif')
        loghistory('')


@timing_function
def mosaic_galactic(*args):
    '''Creates the mosaic of the galactic model'''
    loghistory('--------------------------------------------------------------')
    loghistory('Creating the mosaic of the galactic model (%s.mosaic)' 
               %(galactic.__name__))
    loghistory('')
    galactic.mosaic(*args)
    loghistory('Output galactic mosaic files Griddata/dnight/galtopmags%s.lyr')
    loghistory('')


@timing_function
def mosaic_zodiacal(*args):
    '''Creates the mosaic of the zodiacal model'''
    loghistory('--------------------------------------------------------------')
    loghistory('Creating the mosaic of the zodiacal model (%s.mosaic)' 
               %(zodiacal.__name__))
    loghistory('')
    zodiacal.mosaic(*args)
    loghistory('Output zodiacal mosaic files Griddata/dnight/zodtopmags%s.lyr')
    loghistory('')


@timing_function  
def mosaic_median(*args):
    '''Creates the mosaic from the median-filtered data'''
    loghistory('--------------------------------------------------------------')
    loghistory('Creating the mosaic from median-filtered images (%s.mosaic)' 
               %medianmosaic.__name__)
    loghistory('')
    for filter in args[2]:
        loghistory('%s-band' %filter)
        arg = list(args)
        arg[2] = filter
        medianmosaic.mosaic(*arg)
        a = 'Griddata/dnight/skybrightmags%s.lyr'
        loghistory('Output median-filtered data mosaic files in' + a)
        loghistory('')


@timing_function
def mosaic_fullres(*args):
    '''Creates the mosaic from the full-resolution data'''
    loghistory('--------------------------------------------------------------')
    loghistory('Creating the mosaic from full-resolution images (%s.mosaic)' \
               %fullmosaic.__name__)
    loghistory('')
    for filter in args[2]:
        loghistory('%s-band' %filter)
        arg = list(args)
        arg[2] = filter
        fullmosaic.mosaic(*arg)
        a = 'Griddata/dnight/skytopomags%s.lyr'
        loghistory('Output full-resolution mosaic layer files' + a)
        loghistory('')
    

function_number = {'reduce_images':4, 
                   'register_coord':5, 
                   'calculate_pointing_err':6,
                   'fit_zeropoint':7,
                   'compute_coord':8,
                   'apply_filter':9,
                   'mosaic_galactic':10,
                   'mosaic_zodiacal':11,
                   'mosaic_median':12,
                   'mosaic_fullres':13}
    
#------------ Read in the processing list and initialize ---------------------#
#Local sources
import filepath
import progressbars

#Read in the processing dataset list and the calibration file names 
filelist = n.loadtxt(filepath.processlist+'filelist.txt', dtype=str, ndmin=2)
Dataset, V_band, B_band, Flat_V, Flat_B, Curve, Processor = filelist.T

#Check the calibration files exist    
for i in range(len(filelist)):
    if V_band == 'Yes':
        open(filepath.flats+Flat_V[i])
    if B_band == 'Yes':
        open(filepath.flats+Flat_B[i])
    open(filepath.lincurve+Curve[i]+'.txt')

#Determine the number of data sets collected in each night 
img_sets = set(['1st','2nd','3rd','4th','5th','6th','7th','8th'])
dnight_sets = {}
nsets = []
for dnight in Dataset:
    n_path = filepath.rawdata + dnight + '/'
    dnight_sets[dnight] = []
    for f in os.listdir(n_path):
        if os.path.isdir(n_path+f) & (f in img_sets):
            dnight_sets[dnight].append(f)
    nsets.append(len(dnight_sets[dnight]))
    
    #Make calibration folders
    if not os.path.exists(filepath.calibdata+dnight):
        os.makedirs(filepath.calibdata+dnight)


#Plot the progress bar template
barfig, barax = progressbars.bar(Dataset, nsets)
if all(Processor == 'L_Hung'):
    barfig.canvas.manager.window.move(2755,0)  
else:
    print 'You have 5 seconds to adjust the position of the progress bar window'
    plt.pause(5) #users have 5 seconds to adjust the figure position

#Progress bar array (to be filled with processing time)
Z = n.empty((5+len(filelist),14))*n.nan



#------------ Main data processing code --------------------------------------#

#Import local source code; takes about 12 seconds
print 'Importing local source code... \n'
import coordinates as CO
import extinction as EX
import filepath
import fullmosaic
import galactic
import medianfilter as MF
import medianmosaic
import pointing
import progressbars
import reduce as R
import register as RE
import zodiacal

#Looping through multiple data nights
for i in range(len(filelist)):
    history = open(filepath.calibdata+Dataset[i]+'/processlog.txt', 'w')
    sets = dnight_sets[Dataset[i]]
    
    Filter = []; Filterset = {}
    if V_band[i] == 'Yes': 
        Filter.append('V')
        Filterset['V'] = Flat_V[i]
    if B_band[i] == 'Yes': 
        Filter.append('B')
        Filterset['B'] = Flat_B[i]

    P0 = [Dataset[i],sets,Filterset,Curve[i]]
    P1 = [Dataset[i],sets,Filter] 
    P2 = [Dataset[i],sets]
    
    log_inputs(filelist[i])
    reduce_images(*P0)             #reduce images     
    register_coord(*P1)            #register image positions
    calculate_pointing_err(*P2)    #calculate pointing error
    fit_zeropoint(*P1)             #fit for zeropoint
    compute_coord(*P2)             #compute galactic & ecliptic coordinates
    apply_filter(*P1)              #apply median filter 
    mosaic_galactic(*P2)           #make galactic mosaic
    mosaic_zodiacal(*P2)           #make zodiacal mosaic
    mosaic_median(*P1)             #make median mosaic
    mosaic_fullres(*P1)            #make full-resolution mosaic

    #save the timing records for running the script
    n.savetxt(filepath.calibdata+Dataset[i]+'/processtime.txt', Z, fmt='%4.1f')
    barfig.savefig(filepath.calibdata+Dataset[i]+'/processtime.png')
  
#Todo: Multiprocessing
history.close()
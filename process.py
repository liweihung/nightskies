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


##########################  Definitions  ######################################
def logm(m,string):
    '''To log the processing history and print status on the screen'''
    if type(string) is list:
        for s in string:
            m.append(s+'\n')
            print s
    else:
        m.append(string+'\n')
        print string    


def loghistory(message):
    '''To log the processing history and print status on the screen'''
    [history.write(line) for line in message]
    

def log_inputs(p):
    '''starting the history file and log the input files used'''
    print 'Recording the reduction process in: '
    print filepath.calibdata+dnight+'/processlog.txt'
    print ' '
    m = []
    logm(m,'%s processed on %s by %s' %(p[0],str(Dtime.now())[:19],p[6]))
    logm(m,'')
    logm(m,'Input files:')
    logm(m,'Reducing V-band data? ' + p[1])
    logm(m,'Reducing B-band data? ' + p[2])
    logm(m,'vflat = ' + p[3])
    logm(m,'bflat = ' + p[4])
    logm(m,'curve = ' + p[5])
    logm(m,'')
    logm(m,'____________________ Processing history: _____________________')
    logm(m,'')
    loghistory(m)
   

def update_progressbar(x,y,t=0):
    '''update the progress bar plot'''
    if t == 0:          #gray out for the in-progress status 
        barax.pcolor([x+4,x+5],[y+5,y+6],[[4],],cmap='gray',vmin=0,vmax=5)
    else:               #update for the completed status
        t/=60           #[min]
        if t < 1: texty = '%.1f' %t
        else: texty = '%d' %round(t)
        if t < 6: c = 'k'
        else: c = 'w'
        Z[y+5,x] = t    #record the time [min] in a master array
        barax.pcolor([x+4,x+5],[y+5,y+6],[[t],],cmap='Greens',vmin=0,vmax=10)
        barax.text(x+4.5, y+5.5, texty, color=c, horizontalalignment='center',
                   verticalalignment='center', size='medium')
    plt.pause(0.05)     #draw the new data and run the GUI's event loop
    

def reduce_images(*args):
    '''Basic image reduction'''
    update_progressbar(0,i)
    t1 = time.time()
    print 'Reducing images... '
    import reduce as R
    if 'V' in args[2]:
        R.reducev(args[0],args[1],args[2]['V'],args[3])
    if 'B' in args[2]:
        R.reduceb(args[0],args[1],args[2]['B'],args[3])
    t2 = time.time()
    update_progressbar(0,i,t2-t1)


def register_coord(*args):
    '''Registering coordinates of each image by matching up star positions'''
    update_progressbar(1,i)
    t1 = time.time()
    m = ['Images are normally registered using full frames.','',]
    import register as RE
    for filter in args[2]:
        cropped_fn, failed_fn = RE.matchstars(args[0],args[1],filter)
        logm(m,'These %s images were registered using central 200 pix:'%filter)
        logm(m,cropped_fn)
        logm(m,'')
        logm(m,'These %s images failed to be registered:'%filter)
        logm(m,failed_fn)
        logm(m,'')
    loghistory(m)
    t2 = time.time()
    update_progressbar(1,i,t2-t1)


def pointing_error(*args):
    '''Calculating the pointing error'''
    t1 = time.time()
    import pointing
    pointing.pointing_err(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)


def fit_zeropoint(*args):
    '''Fitting for the extinction coefficient and zeropoint'''
    t1 = time.time()
    m = []
    import extinction as EX
    for filter in args[2]:
        arg = list(args[:-1])
        arg[2] = filter
        nstar, bestfit_file = EX.extinction(*arg)
        bestfitp = list(n.loadtxt(bestfit_file, dtype=str, comments='', 
                        delimiter='!'))
        logm(m,'%s-band best fit zeropoint and extinction: ' %filter)
        logm(m, bestfitp)
        logm(m, '\n')        
    t2 = time.time()
    args[-1].put(t2-t1)
    args[-1].put(m)


def apply_filter(*args):
    '''Apply ~1 degree (in diameter) median filter to the images'''
    t1 = time.time()
    import medianfilter
    for filter in args[2]:
        print 'Applying median filter to %s-band images' %filter
        medianfilter.filter(args[0],args[1],filter)
    t2 = time.time()
    args[-1].put(t2-t1)


def compute_coord(*args):
    '''Computig the galactic and ecliptic coordinates of the images'''
    t1 = time.time()
    import coordinates as CO
    CO.galactic_ecliptic_coords(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)


def mosaic_galactic(*args):
    '''Creates the mosaic of the galactic model'''
    t1 = time.time()
    import galactic
    print 'Creating the mosaic of the galactic model'
    galactic.mosaic(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)


def mosaic_zodiacal(*args):
    '''Creates the mosaic of the zodiacal model'''
    t1 = time.time()
    import zodiacal
    print 'Creating the mosaic of the zodiacal model' 
    zodiacal.mosaic(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)


def mosaic_full(*args):
    '''Creates the mosaic from the full-resolution data'''
    t1 = time.time()
    import fullmosaic
    for filter in args[2]:
        print 'Creating the mosaic from full-resolution %s-band images' %filter
        fullmosaic.mosaic(args[0],args[1],filter)
    t2 = time.time()
    args[-1].put(t2-t1)

        
def mosaic_median(*args):
    '''Creates the mosaic from the median-filtered data'''
    t1 = time.time()
    import medianmosaic
    for filter in args[2]:
        print 'Creating the mosaic from median-filtered %s-band images' %filter
        medianmosaic.mosaic(args[0],args[1],filter)
    t2 = time.time()
    args[-1].put(t2-t1)
    

import time
import numpy as n

if __name__ == '__main__':
    t1 = time.time()
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
    import os
    import pdb
    import warnings
    
    from datetime import datetime as Dtime
    from multiprocessing import Process, Queue
    
    warnings.filterwarnings("ignore",".*GUI is implemented.*")
        
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
    if all(Processor == 'L_Hung2'):
        barfig.canvas.manager.window.move(2755,0)  
    else:
        print 'You have 5 seconds to adjust the position of the progress bar window'
        plt.pause(5) #users have 5 seconds to adjust the figure position
    
    #Progress bar array (to be filled with processing time)
    Z = n.empty((5+len(filelist),14))*n.nan
    
    
    
    #------------ Main data processing code --------------------------------------#
    import filepath
    import progressbars
        
    #Looping through multiple data nights
    for i in range(len(filelist)):
        history = open(filepath.calibdata+Dataset[i]+'/processlog.txt', 'w')
        log_inputs(filelist[i])
        
        Filter = []; Filterset = {}
        if V_band[i] == 'Yes': 
            Filter.append('V')
            Filterset['V'] = Flat_V[i]
        if B_band[i] == 'Yes': 
            Filter.append('B')
            Filterset['B'] = Flat_B[i]
        
        sets = dnight_sets[Dataset[i]]
        K0 = (Dataset[i],sets,Filterset,Curve[i])
        K1 = (Dataset[i],sets,Filter) 
        K2 = (Dataset[i],sets)  

        q2=Queue(); Q2=(q2,); p2=Process(target=pointing_error,args=K2+Q2)
        q3=Queue(); Q3=(q3,); p3=Process(target=fit_zeropoint,args=K1+Q3)
        q4=Queue(); Q4=(q4,); p4=Process(target=apply_filter,args=K1+Q4)
        q5=Queue(); Q5=(q5,); p5=Process(target=compute_coord,args=K2+Q5)  
        q6=Queue(); Q6=(q6,); p6=Process(target=mosaic_galactic,args=K2+Q6)
        q7=Queue(); Q7=(q7,); p7=Process(target=mosaic_zodiacal,args=K2+Q7)
        q8=Queue(); Q8=(q8,); p8=Process(target=mosaic_full,args=K1+Q8)
        q9=Queue(); Q9=(q9,); p9=Process(target=mosaic_median,args=K1+Q9)
        
        reduce_images(*K0)                            #image reduction   
        register_coord(*K1)                           #pointing 
        p2.start(); update_progressbar(2,i)           #pointing error
        p3.start(); update_progressbar(3,i)           #zeropoint & extinction
        p4.start(); update_progressbar(4,i)           #median filter
        p2.join() ; update_progressbar(2,i,q2.get())
        p5.start(); update_progressbar(5,i)           #galactic & ecliptic coord
        p5.join() ; update_progressbar(5,i,q5.get())
        p6.start(); update_progressbar(6,i)           #galactic mosaic
        p7.start(); update_progressbar(7,i)           #zodiacal mosaic
        p3.join() ; update_progressbar(3,i,q3.get())
        p8.start(); update_progressbar(8,i)           #full mosaic
        p4.join() ; update_progressbar(4,i,q4.get())
        p9.start(); update_progressbar(9,i)           #median mosaic
        p6.join() ; update_progressbar(6,i,q6.get())
        p7.join() ; update_progressbar(7,i,q7.get())
        p8.join() ; update_progressbar(8,i,q8.get())
        p9.join() ; update_progressbar(9,i,q9.get())
        
        #log the processing history
        q_all = [q2,q3,q4,q5,q6,q7,q8,q9]
        for q in q_all:
            while not q.empty():
                loghistory(q.get())
        
    
        #save the timing records for running the script
        n.savetxt(filepath.calibdata+Dataset[i]+'/processtime.txt', Z, fmt='%4.1f')
        barfig.savefig(filepath.calibdata+Dataset[i]+'/processtime.png')
        
    t2 = time.time()
    print 'Total processing time: %.1f min' %((t2-t1)/60)
    loghistory(['Total processing time: %.1f min' %((t2-t1)/60),])
    history.close()
    
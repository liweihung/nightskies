#-----------------------------------------------------------------------------#
#fullmosaic.py
#
#NPS Night Skies Program
#
#Last updated: 2017/03/08
#
#This script makes the whole sky mosaic from the full resolution images 
#according to the location of observed sky. The temporary files generated during
#the processing stage are stored in filepath.rasters/scratch_fullres folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) full resolution tiff files in the filepath.calibdata
#   (2) pointerr_%s.txt
#   (3) extinction_fit_%s.txt
#   (4) raster files in the filepath.rasters folder
#
#Output:
#   (1) layer files skytopomags%s.lyr for full-resolution mosaic
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#
#-----------------------------------------------------------------------------#
from glob import glob, iglob
from scipy.misc import imread

import arcpy
import pdb
import numpy as n
import os
import shutil

# Local Source
import filepath  

#-----------------------------------------------------------------------------#
if not os.path.exists(filepath.rasters+'scratch_fullres/'):
    os.makedirs(filepath.rasters+'scratch_fullres/')
    
geogcs = "GEOGCS['GCS_Sphere_EMEP',\
          DATUM['D_Sphere_EMEP',SPHEROID['Sphere_EMEP',6370000.0,0.0]],\
          PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
          
#set arcpy environment variables part 1/2
arcpy.env.rasterStatistics = "NONE"
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = "NONE"

# define source control points
source_pnt = "'0 0';'0 296039.8';'0 590759.1';'0 884157.9';'0 1176236';\
'0 1466994';'0 -296039.8';'0 -590759.1';'0 -884157.9';'0 -1176236';\
'0 -1466994';'-296039.8 0';'-590759.1 0';'-884157.9 0';'-1176236 0';\
'-1466994 0';'296039.8 0';'590759.1 0';'884157.9 0';'1176236 0';'1466994 0';\
'1241985 1241985';'-1241985 -1241985';'-1241985 1241985';'1241985 -1241985';\
'1445714 1445714';'-1445714 1445714';'-1445714 -1445714';'1445714 -1445714';\
'1037322 1037322';'-1037322 1037322';'-1037322 -1037322';'1037322 -1037322';\
'417730 417730';'-417730 417730';'-417730 -417730';'417730 -417730'"

# define target control points
target_pnt = "'0 0';'0 296708';'0 593400';'0 890100';'0 1186800';'0 1483500';\
'0 -296700';'0 -593400';'0 -890100';'0 -1186800';'0 -1483500';'-296700 0';\
'-593400 0';'-890100 0';'-1186800 0';'-1483500 0';'296700 0';'593400 0';\
'890100 0';'1186800 0';'1483500 0';'1258791 1258791';'-1258791 -1258791';\
'-1258791 1258791';'1258791 -1258791';'1468590 1468590';'-1468590 1468590';\
'-1468590 -1468590';'1468590 -1468590';'1048993 1048993';'-1048993 1048993';\
'-1048993 -1048993';'1048993 -1048993';'419597 419597';'-419597 419597';\
'-419597 -419597';'419597 -419597'"
          
#-----------------------------------------------------------------------------#
def clip_envelope(AZ, ALT, i):
    if i < 15:
        bond = [AZ[i]-13,-6,AZ[i]+13,ALT[i]+12.9] 
    elif i < 30: 
        bond = [AZ[i]-13,ALT[i]-12.7,AZ[i]+13,ALT[i]+12.6] 
    elif i < 40:
        bond = [AZ[i]-18.6,ALT[i]-12.7,AZ[i]+18.6,ALT[i]+12.6] 
    elif i < 45:
        bond = [AZ[i]-39.6,ALT[i]-12.7,AZ[i]+39.6,ALT[i]+12.7] 
    return ' '.join(str(i) for i in bond) #'xmin ymin xmax ymax' 
    
    
def tc(lon,lat):
    '''Returns the topocentric coordinate setting'''
    return "PROJCS['gnomonic',%s,PROJECTION['Gnomonic'],\
    PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],\
    PARAMETER['Longitude_Of_Center',%s],PARAMETER['Latitude_Of_Center',%s],\
    UNIT['Meter',1.0]]"%(geogcs,str(lon),str(lat))
    

def mosaic(dnight, sets, filter):
    '''
    This module creates the mosaic of full-resolution images for each data set.
    '''
    #set arcpy environment variables part 2/2
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.workspace = filepath.rasters+'scratch_fullres/'
    arcpy.env.scratchWorkspace = filepath.rasters+'scratch_fullres'

    #filter paths
    F = {'V':'', 'B':'B/'}
    f = {'V':'', 'B':'b'}
    
    for s in sets:
        #file paths
        calsetp = filepath.calibdata+dnight+'/S_0%s/%s' %(s[0],F[filter])
        gridsetp = filepath.griddata+dnight+'/S_0%s/%sfullres/' %(s[0],F[filter])
        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp)
        os.makedirs(gridsetp)
                
        #read in the registered images coordinates
        file = filepath.calibdata+dnight+'/pointerr_%s.txt' %s[0]
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        
        #read in the best-fit zeropoint and plate scale
        file = filepath.calibdata+dnight+'/extinction_fit_%s.txt' %filter
        zeropoint, platescale, exptime = n.loadtxt(file, usecols=(2,8,9), unpack=True, ndmin=2)
        
        #loop through each file in the set
        for w in range(len(Obs_AZ)+1):

            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360
            
            if v in range(0,50,5): print 'Generating fullres image %i/45'%v
            
            arcpy.CopyRaster_management(calsetp+'/tiff/ib%03d.tif' %(w+1), 'ib%03d.tif' %v,"DEFAULTS","","","","","16_BIT_UNSIGNED")
            
            #re-define projection to topocentric coordinates
            arcpy.DefineProjection_management("ib%03d.tif" %v,tc(Obs_AZ[w],Obs_ALT[w]))
            
            #warp image to remove barrel distortion image
            arcpy.Warp_management("ib%03d.tif"%v, source_pnt, target_pnt, 'ibw%03d.tif'%v, "POLYORDER3", "BILINEAR")

            #reproject into GCS
            arcpy.ProjectRaster_management('ibw%03d.tif' %v, 'fwib%03d.tif' %v, geogcs, "BILINEAR", "0.0261")
                                       
            #clip to image boundary
            rectangle = clip_envelope(Obs_AZ, Obs_ALT, w)
            arcpy.Clip_management("fwib%03d.tif"%v, rectangle, "fcib%03d"%v)
            
        #mosaic raster list must start with an image with max pixel value > 256
        v=1; mstart=1
        while v < (len(Obs_AZ)+1):
            im = imread(filepath.rasters+'scratch_fullres/ib%03d.tif' %v)
            if n.max(im) > 255:
                mstart = v
                break
            v+=1
                        
        #mosaic raster list
        R1 = ';'.join(['fcib%03d' %i for i in range(mstart,47)])
        R2 = ';'.join(['fcib%03d' %i for i in range(1,mstart)])
        R = R1+';'+R2
        
        #mosaic to topocentric coordinate image; save in Griddata\
        print "Mosaicking into all sky full-resolution image"
        arcpy.MosaicToNewRaster_management(R, gridsetp, 'skytopo', geogcs, 
                                        "32_BIT_FLOAT", "0.0261", "1", "BLEND", 
                                        "FIRST")
        
        #convert to magnitudes per square arc second
        print "Converting the mosaic to mag per squard arcsec"
        psa = 2.5*n.log10((platescale[int(s[0])-1]*60)**2) # platescale adjustment
        
        # break the arcpy calculations into a few steps
        stm1 = arcpy.sa.raster(gridsetp + skytopo)
        stm2 = arcpy.sa.log10(stm1)
        stm3 = 2.5 * stm2 / esptime[0]
        
        skytopomags = zeropoint[int(s[0])-1] + psa - stm3
        
        #save mags mosaic to disk
        skytopomags.save(gridsetp+'skytopomags')
    
        print "Creating layer files for full-resolution mosaic"
        layerfile = filepath.griddata+dnight+'/skytopomags%s%s.lyr' %(f[filter],s[0])
        arcpy.MakeRasterLayer_management(gridsetp+'skytopomags', dnight+'_%s_fullres%s'%(s[0],f[filter]))
        arcpy.SaveToLayerFile_management(dnight+'_%s_fullres%s'%(s[0],f[filter]), layerfile, "RELATIVE")
    
        #Set layer symbology to magnitudes layer
        symbologyLayer = filepath.rasters+'magnitudes.lyr'
        arcpy.ApplySymbologyFromLayer_management(layerfile, symbologyLayer)
        lyrFile = arcpy.mapping.Layer(layerfile)
        lyrFile.replaceDataSource(gridsetp,'RASTER_WORKSPACE','skytopomags',
                                  'FALSE')
        lyrFile.save()

    
if __name__ == "__main__":
    #mosaic('FCNA160803', ['1st',],'B')
    pass

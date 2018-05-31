#-----------------------------------------------------------------------------#
#galactic.py
#
#NPS Night Skies Program
#
#Last updated: 2016/12/23
#
#This script makes the whole sky mosaic of the galactic model according to the 
#time and location of the observed sky. The temporary files generated during the 
#processing stage are stored in filepath.rasters/scratch_galactic folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) raster files in the filepath.rasters folder
#   (2) coordinates_%s.txt
#   (3) pointerr_%s.txt
#
#Output:
#   (1) layer files galtopmags%s.lyr for galactic mosaic
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#
#-----------------------------------------------------------------------------#

from glob import glob, iglob

import arcpy
#import pdb
import numpy as n
import os
import shutil

# Local Source
import filepath  

#-----------------------------------------------------------------------------#
if not os.path.exists(filepath.rasters+'scratch_galactic/'):
    os.makedirs(filepath.rasters+'scratch_galactic/')

#set arcpy environment variables part 1/2
arcpy.env.rasterStatistics = "STATISTICS 2 2 (-999)"
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = "NONE"

#input rasters for galactic model
galraster  = arcpy.sa.Raster(filepath.rasters+'galmagsnew')
galraster1 = arcpy.sa.Raster(filepath.rasters+'galmagsnew180')
galrastern = arcpy.sa.Raster(filepath.rasters+'galnewnorth')
galrasters = arcpy.sa.Raster(filepath.rasters+'galnewsouth')

geogcs = "GEOGCS['GCS_Sphere_EMEP',\
          DATUM['D_Sphere_EMEP',SPHEROID['Sphere_EMEP',6370000.0,0.0]],\
          PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
    
#-----------------------------------------------------------------------------#
def gal_envelope(lon,lat):
    if abs(lat)<71:
        #expension angle
        if abs(lat)<55: ang = 22/n.cos(n.deg2rad(lat))
        else: ang = 89
        lon = (lon+90) % 180 - 90
        bond = [lon-ang, lat-19, lon+ang, lat+19] 
    else:
        bond = [-2223549.5, -2223549.5, 2223549.5, 2223549.5]  
    return ' '.join(str(i) for i in bond) #'xmin ymin xmax ymax' 

    
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
    
    
def get_galgn(lon,lat):
    rectangle = gal_envelope(lon,lat)
    if abs(lat)<71:
        if abs(lon)<90: arcpy.Clip_management(galraster,rectangle,"galclip.tif")
        else: arcpy.Clip_management(galraster1, rectangle, "galclip.tif")
        lon = (lon+90) % 180 - 90
        p = [tc(lon,lat),"BILINEAR", "6000"]
        arcpy.ProjectRaster_management("galclip.tif", "galgn.tif", *p)
    else: 
        p = [tc(lon,lat),"BILINEAR", "6000"]
        if lat>=71:
            arcpy.ProjectRaster_management(galrastern, "galtemp.tif", *p)
        else:
            arcpy.ProjectRaster_management(galrasters, "galtemp.tif", *p)
        arcpy.Clip_management('galtemp.tif', rectangle, 'galgn.tif')
        

def mosaic(dnight, sets):
    '''
    This module creates the mosaic of the galactic model for each data set.
    '''
    #set arcpy environment variables part 2/2
    arcpy.env.workspace = filepath.rasters+'scratch_galactic/'
    arcpy.env.scratchWorkspace = filepath.rasters+'scratch_galactic'

    for s in sets:
        #file paths
        calsetp = filepath.calibdata+dnight+'/S_0%s/' %s[0]
        gridsetp = filepath.griddata+dnight+'/S_0%s/gal/' %s[0]
        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp)
        os.makedirs(gridsetp)
        
        #read in the galactic coordinates from coordinates_%s.txt
        file = filepath.calibdata+dnight+'/coordinates_%s.txt'%s[0]
        Gal_ang, Gal_l, Gal_b = n.loadtxt(file,usecols=(1,2,3)).T
        
        #read in the registered images coordinates
        file = filepath.calibdata+dnight+'/pointerr_%s.txt' %s[0]
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        
        #loop through each file in the set
        for w in range(len(Obs_AZ)+1):
            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360
                
            get_galgn(Gal_l[w], Gal_b[w])
            if v in range(0,50,5): print 'Generating galactic image %i/45'%v
        
            #rotate by galctic angle
            arcpy.Rotate_management('galgn.tif', 
                                    'rotaterasterg.tif', 
                                    str(Gal_ang[w]), 
                                    "0 0",
                                    "BILINEAR")
                                    
            #re-define projection to topocentric coordinates
            arcpy.DefineProjection_management('rotaterasterg.tif',
                                            tc(Obs_AZ[w],Obs_ALT[w]))
                                            
            #reproject into GCS
            arcpy.ProjectRaster_management('rotaterasterg.tif', 
                                        'gal%02d.tif'%v, 
                                        geogcs,
                                        "BILINEAR",
                                        "0.05")
                                        
            #clip to image boundary
            rectangle = clip_envelope(Obs_AZ, Obs_ALT, w)
            arcpy.Clip_management("gal%02d.tif"%v, rectangle, "gali%02d"%v)
    
        #Mosaic to topocentric coordinate model; save in Griddata\
        print "Mosaicking into all sky galactic model"
        R = ';'.join(['gali%02d' %i for i in range(1,47)])
        arcpy.MosaicToNewRaster_management(R, gridsetp, 'galtopmags', geogcs, 
                                        "32_BIT_FLOAT", "0.05", "1", "BLEND", 
                                        "FIRST")
    
        print "Creating layer files for galactic mosaic"
        layerfile = filepath.griddata+dnight+'/galtopmags%s.lyr' %s[0]
        arcpy.MakeRasterLayer_management(gridsetp+'galtopmags', 'galtoplyr')
        arcpy.SaveToLayerFile_management('galtoplyr', layerfile, "ABSOLUTE")
    
        #Set layer symbology to magnitudes layer
        symbologyLayer = filepath.rasters+'magnitudes.lyr'
        arcpy.ApplySymbologyFromLayer_management(layerfile, symbologyLayer)
        lyrFile = arcpy.mapping.Layer(layerfile)
        lyrFile.replaceDataSource(gridsetp,'RASTER_WORKSPACE','galtopmags',
                                  'FALSE')
        lyrFile.save()
    
if __name__ == "__main__":
    #mosaic('FCNA160803', ['1st',])
    pass









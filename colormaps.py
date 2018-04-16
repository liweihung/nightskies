#-----------------------------------------------------------------------------#
#colormaps.py
#
#NPS Night Skies Program
#
#Last updated: 2018/04/16
#
#This script generates a customized color map for displaying NPS night sky 
#panoramic data. This color map replicated the existing ArcGIS color map by 
#using the same RGB values. The range of the color map should be set from 14 to 
#24 magnitudes. 
#
#Input: 
#   (1) colormap_magnitudeslyr.txt -- the RGB values from magnitudes.lyr
#
#Output:
#   (1) registered plt color map 'NPS_mag'
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
import matplotlib
import matplotlib.pyplot as plt
import numpy as n

# Local Source
import filepath  

#-----------------------------------------------------------------------------#

#RGB values originally from magnitudes.lyr
colormap_file = filepath.rasters+'colormap_magnitudeslyr.txt'
mag_start, mag_end, R, G, B = n.loadtxt(colormap_file).T    #RGB in 0-255 scale 

#color positions in 0-1 scale
mag_range = mag_start[-1]-mag_start[0]
mag_percent = (mag_start-mag_start[0])/mag_range

#RGB values in 0-1 scale
R /= 255.
G /= 255.
B /= 255.

#color map lists 
red = []; green = []; blue = []
for i in range(len(mag_start)):
    red.append((mag_percent[i],R[i],R[i]))
    green.append((mag_percent[i],G[i],G[i]))
    blue.append((mag_percent[i],B[i],B[i]))

#declare color map setting   
cdict = {'red':red,'green':green,'blue':blue}

#register the color map
plt.register_cmap(name='NPS_mag', data=cdict)

if __name__ == '__main__':
    plt.rcParams['image.cmap'] = 'NPS_mag'
    img = n.arange(10000).reshape(100,100)
    plt.imshow(img, interpolation='nearest')#, cmap=nps_mag
    cbar = plt.colorbar()
    plt.show(block=False)
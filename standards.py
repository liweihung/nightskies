#-----------------------------------------------------------------------------#
#standards.py
#
#NPS Night Skies Program
#
#Last updated: 2016/11/22
#
#This script puts the numbers in the hip3.xlsx to a plain txt file. hip3.txt 
#contains 371 standard stars from the Hipparcos Catalog. This list serves as a 
#standard star library for doing zeropoint and extinction correction. This 
#script should only need to be run once; once the text file is generated, this 
#script does not need to be used again. The starndard stars are ploted; the 
#three output plots are identical except for the colors used. 
#
#Input: 
#   (1) hip3.xlsx
#
#Output:
#   (1) hipparcos_standards.txt
#   (2) standards_b.png
#   (3) standards_c.png
#   (4) standards_w.png
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#

from xlrd import open_workbook
import matplotlib.pyplot as plt
import numpy as n
# Local Source
import filepath  


#Read in the first sheet (process_list) from filelist.xlsx
F = open_workbook(filepath.standards+'hip3.xlsx').sheets()[0]
star = F.col_values(0) #From the Hipparcos Catalog (ESA 1997)
RA = F.col_values(1)
Dec = F.col_values(2)
Mag = F.col_values(3) #Magnitude
B_V = F.col_values(4) #B-V magnitude

standards = n.array((star,RA,Dec,Mag,B_V), dtype = object).T
F = ['%7s','%10.6f','%11.6f','%6.3f','%7.3f']
H = ' Star     RA[h]     Dec[deg]  V_mag    B-V'
n.savetxt(filepath.standards+'hipparcos_standards.txt',
          standards,fmt=F,header=H)


def plot_standards():
    '''common setting for the plots'''
    cb = plt.colorbar()
    cb.ax.invert_yaxis()
    cb.set_label('V magnitude')
    plt.xlim(0,24)
    plt.ylim(-90,90)
    plt.title('Standard star library used in calculating extinction and ZP')
    plt.xlabel('RA (hr)', fontsize=14)
    plt.ylabel('Dec (degree)', fontsize=14)

#colored figure
fig1, ax1 = plt.subplots()
plt.scatter(RA, Dec, s=50, c=Mag)
plot_standards()
plt.savefig(filepath.standards+'standards_c.png',dpi=200)      

#black and white figure with white background
fig2, ax2 = plt.subplots()
plt.scatter(RA, Dec, s=50,cmap='gist_gray_r', c=Mag)
plot_standards()
plt.savefig(filepath.standards+'standards_w.png',dpi=200)      

#black and white figure with black background
fig3, ax3 = plt.subplots()
plt.scatter(RA, Dec, s=50,cmap='gist_gray_r', c=Mag, edgecolor='gray')
ax3.patch.set_facecolor('black')
plot_standards()
plt.savefig(filepath.standards+'standards_b.png',dpi=200)      

plt.show(block=False)

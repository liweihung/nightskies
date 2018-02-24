# NSNSD Night Skies CCD Camera Data Reduction Pipeline
![false color CCD image](https://github.com/liweihung/nightskies/blob/master/static/FalseColor_example.png)

<!-- MarkdownTOC autolink=true depth=3 bracket=round -->

- [Purpose](#purpose)
- [Required Software and Install Procedures](#required-software)
- [Preparing Data For Processing](#preparing-data)
- [Processing Flow Chart](#processing-flow-chart)
- [Module Documentation](#module-documentation)
  - [1. Reduction](#1-reduction)
  - [2. Registration](#2-registration)
  - [3. Pointing Error](#3-pointing-error)
  - [4. Zeropoint and Extinction](#4-zeropoint-and-extinction)
  - [5. Median Filter](#5-median-filter)
  - [6. Galactic and Zodiacal Coordinates](#6-galactic-and-zodiacal-coordinates)
  - [7. Galactic Mosaic](#7-galactic-mosaic)
  - [8. Zodiacal Mosaic](#8-galactic-and-zodiacal-coordinates)
  - [9. Full-resolution Mosaic](#9-galactic-and-zodiacal-coordinates)
  - [10. Median-filtered Mosaic](#10-galactic-and-zodiacal-coordinates)

- [Installation](#installation)
  - [Dependencies](#dependencies)
  - [Public domain](#public-domain)

<!-- /MarkdownTOC -->


----------------------------------


## Purpose

## Required Software and Install Procedures

## Preparing Data For Processing

## Processing Flow Chart

![processing flow chart](https://github.com/liweihung/nightskies/blob/master/static/FlowChart.png)

## Module Documentation

### 1. Reduction

#### Purpose: 
This script performs basic image reduction, including corrections for bias, dark, flat, and linearity response of the detector.

#### Source code
`process.py` > `Reduce images()` > `reduce.py` > `reducev()` and `reduceb()`  

#### Methods
In the beginning of data collection, for each data set, 5 dark (**D**) images showing the thermal noise and 5 bias (**B**) images showing the read noise were taken alternately for calibration purposes. Each dark image is then processed here as following to obtain the calibrated dark (**D<sub>c</sub>**) image:  

![D_c equation](https://github.com/liweihung/nightskies/blob/master/static/D_c.png)

where **L** is the linearity curve for correcting the detector response. Then, the script creates the master dark image **D<sub>m</sub>** through averaging the 5 calibrated darks and the master bias image **B<sub>m</sub>** through averaging the 5 biases. While collecting science images, an additional 50 pixels by 50 pixels small bias image was taken immediately after each science image. We use these small bias images to track the bias drift over the course of the observation. We compute the bias drift **B<sub>d</sub>** by subtracting the average value of the central 50 pixels by 50 pixels of the master bias from the average pixel value of each one of the small bias images. This measured bias drift information is saved in baisdrift.txt and biasdrift.png in the calibdata folder. To obtain scientific, calibrated images **S<sub>c</sub>**, we use the following equation:

![S_c equation](https://github.com/liweihung/nightskies/blob/master/static/S_c.png)
 
where **S** is the raw science image and **F** is the flat image taken and processed in the lab. All of the terms in the above equation are 2D image arrays except for **B<sub>d</sub>** and **L** which are single-value scale factors. This script outputs calibrated science images in both fits and tiff formats. 


### 2. Registration

#### Purpose: 
This script will register the pointing of the images to the sky coordinates by matching the position of standard stars captured in the image to the position from the Tycho2 catalog. This script uses the PinPoiot.Plate object (through ACP Observatory Software) for the matching process. 

#### Source code: 
`process.py` > `register_coord()` > `register.py` > `matchstars()`

#### Methods: 
The script first reads in the reduced images. If the image contains the horizon, the script will mask the part of the image below ~2 degrees in elevation. Then, the script uses the PinPoiot.Plate object (through ACP Observatory Software) to compare the position of the standard stars in the images to the Tycho2 catalog. If the image cannot be solved, it will try solving the image with only the central 200x200 pixels. If the images still can't be solved, it will be skipped. The script will return a list of files solved with the cropped images and a list of files that are failed to be solved. If an image is solved successfully, the image header will be updated with the solved coordinate information.


### 3. Pointing Error

#### Purpose: 
This script calculates the pointed azimuth and altitude using the solved RA and Dec values from the image headers. If the images are not solved, the RA, Dec, AZ, and ALT values are interpolated. The output from this script will be used later for making the mosaic. 

#### Source code: 
`process.py` > `pointing_error()` > `pointing.py` > `pointing_err()`

#### Methods: 
It reads in the solved RA and Dec values from the image header, updates the coordinates to the observed date, translates to the azimuth and altitude given the LAST and the longitude, interpolates RA, Dec, AZ, and ALT values if the images are not solved, updates the headers, and records these values in a text file.


### 4. Zeropoint & Extinction

#### Purpose: 
This script finds the best-fit extinction coefficient and the instrumental zeropoint. These values are needed for the photometric calibration later. 

#### Source code: 
`process.py` > `fit_zeropoint` > `extinction.py` > `extinction()`

#### Methods: 
This script finds the best-fit extinction coefficient and the instrumental zeropoint by identifying the standards stars in the images using PinPoiot.Plate object, measuring their background-subtracted flux in [DN/s] through aperture photometry with source aperture radius = 3 pixels and ring-shaped background with radius extending from 4 to 8 pixels, computing their airmass through their elevation, comparing the measured flux to their absolute magnitude, and finally fitting for the extinction given the airmass and M-m for all of the standard stars. The best-fit was found through numpy.polyfit() with the uncertainties of the best-fit parameters estimated using the covariance matric. This script outputs a list of the standard stars used for fitting, a graphical display of the fitting result, and the best-fit extinction coefficient and zeropoint.


### 5. Median Filter

#### Purpose: 
This script applies median filter to each image to remove the stars. The filtered images will show the sky background brightness.

#### Source code: 
`process.py` > `apply_filter()` > `medianfilter.py` > `filter()`

#### Methods: 
We first read in the plate scale from the header of a solved image. Then we set the filter mask radius to be ~ 0.5 degree, which translates to a radius of 18 pixels for ML3 and ML4 CCD cameras. This filter size was selected to ensure most (or all) point sources are effectively filtered out, and the median value in that filtered window is a good representation for the background brightness. This script uses multiprocessing to speed up the filtering process for all the images. Here, MaxIM DL is used to convert the filtered fits images to tiff images so that they are compatible with ArcGIS to make mosaics in the future.


### 6. Galactic & Zodiacal Coordinates

#### Purpose: 
This script calculates ecliptic and galactic coordinates and rotation angles of the images. These output coordinates will be used in producing natural sky model for the given location, date, and time specific to the data set.

#### Source code: 
`process.py` > `compute_coord()` > `coordinates.py` > `galactice_ecliptic_coords()`

#### Methods: 
The script reads in the longitude and latitude from the image header. It also reads in the registered image coordinates in azimuth and altitude from the input file. Then, it uses ACP to calculate the corresponding galactic and zodiacal coordinates. The output coordinates are saved in a text file for the future use of making the mosaic images of galactic and zodiacal model.


### 7. Galactic Mosaic

#### Purpose: 
This script makes the whole sky mosaic of the galactic model according to the time and location of the observed sky.

#### Source code: 
`process.py` > `mosaic_galactic()` > `galactic.py` > `mosaic()`

#### Methods: 
This script reads in some premade raster templates from the raster folder.  Then it reads in the galactic coordinates and the pointing error from the input files. We use arcpy (ArcGIS) to manipulate the images with rotation, projection, clipping, and mosaic. The output raster and layer files are stored in the Griddata folder.


### 8. Zodiacal Mosaic

#### Purpose: 
This script makes the whole sky mosaic of the zodiacal model according to the time and location of the observed sky.

#### Source code: 
`process.py` > `mosaic_zodiacal()` > `zodiacal.py` > `mosaic()`

#### Methods: 
This script reads in some premade raster templates from the raster folder.  Then it reads in the zodiacal coordinates and the pointing error from the input files. We use arcpy (ArcGIS) to manipulate the images with rotation, projection, clipping, and mosaic. The output raster and layer files are stored in the Griddata folder.


### 9. Full-resolution Mosaic

#### Purpose: 
This script makes the whole sky mosaic from the full resolution images according to the location of observed sky.

#### Source code: 
`process.py` > `mosaic_full()` > `fullmosaic.py` > `mosaic()`

#### Methods: 
This script reads in the full resolution tiff files and the pointing error file. It also uses some premade raster templates from the raster folder. We use arcpy (ArcGIS) to manipulate the images with projection, removal of distortion, clipping, and mosaic. Note that the mosaic raster list must start with an image with maximum pixel value being greater than 256 to avoid _“no data”_ in the final mosaic. We then read in the zeropoint from the input file and use it to convert the mosaic from the raw unit in data number (**DN**) to calibrated unit in magnitudes per square arcsecond using the following equation: 

![M_c equation](https://github.com/liweihung/nightskies/blob/master/static/M_c.png)

where **M<sub>c</sub>** is the mosaic in calibrated unit of magnitudes per square arc second, **Z** is the instrumental zeropoint magnitude, Mr  is the mosaic in raw unit of **DN**, **t** is the exposure time, and **P** is the plate scale in arcsecond per pixel. The photometrically calibrated raster and layer files are stored in the Griddata folder.


### 10. Median-filtered Mosaic

#### Purpose: 
This script makes the whole sky mosaic from the median filtered images according to the location of observed sky.

#### Source code: 
`process.py` > `mosaic_median()` > `medianmosaic.py` > `mosaic()`

#### Methods: 
Same as the method section under [_“full-resolution mosaic”_](### Full-resolution Mosaic)  with the exception of using median filtered images as the input. 


## Installation

### Dependencies

### Public domain

This project is in the worldwide [public domain](LICENSE.md). As stated in [CONTRIBUTING](CONTRIBUTING.md):

> This project is in the public domain within the United States,
> and copyright and related rights in the work worldwide are waived through the
> [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
>
> All contributions to this project will be released under the CC0 dedication.
> By submitting a pull request, you are agreeing to comply with this waiver of copyright interest.

I. DESCRIPTION

Precise registration and orthorectification of remote sensing images are the basic 
processes for quantitative remote sensing applications, especially for multi-temporal 
image analysis. We developed an automated precise registration and orthorectification 
package (AROP) for Landsat and Landsat-like data processing. The Landsat and 
Landsat-like satellite images acquired from different sensors at different spatial 
resolutions and projections can be re-projected, co-registered, and orthorectified to the 
same projection, geographic extent, and spatial resolution using a common base image; 
this allows us to perform pixel-by-pixel image analysis directly. This package has been 
tested on the Landsat MSS/TM/ETM+, TERRA ASTER, CBERS CCD and IRS-P6 
AWiFS data. The detailed description of algorithm can be cited and found in:

Gao., F., Masek, J. G., Wolfe, R. F., An automated registration and orthorectification 
package for Landsat and Landsat-like data processing, Journal of Applied Remote 
Sensing, Vol.3, 033515, doi: 10.1117/1.3104620, 2009. 


II. INSTALLATION (tested under Linux bash)

gzip -d public_AROP_v2.2.6.tar.gz
tar -xvf public_AROP_v2.2.6.tar
cd public_AROP_v2.2.6
source env.sh
make
make install

You will see an executable program named "ortho" under sub-directory 
"public_AROP_v2.2.6/bin" if programs were compiled successfully.


III. USAGE

Usage: ortho [-r|-o|-b|-v] <parameter_file>
       -r  do registration only 
       -o  do orthorectification only (assume warp image has been registrated)
       -b  do both registration and othorectification (recommend for all L1G data)
       -v  verify image matching

Input parameter file defines the detailed information for base, warp and output image. 
The base image and the matching band from warp images will be used for the area-
based tie point searching base on the tie point control file. Please check appendix A.1 
and A.2 published in the Journal of Applied Remote Sensing by Gao et al. (2009) for 
descriptions of input parameter file and tie point control file. The updated input 
parameter file is provided in the end of this page.

IV. TIPS

1. Follow examples under "tests" subdirectory to familiar with different AROP options 
in handling different data sources. 

2. The AROP package assumes the warp image covers the COMPLETE swath (not 
subset) when running orthorectification options ("-b" and "-o"). The AROP uses scene 
boundaries (and pointing angle if not zero) to determine the satellite nadir track 
function. 

3. If the input file is in GeoTIFF format and the input parameters (e.g., NSAMPLE, 
NLINE, PIXEL_SIZE, or UPPER_LEFT_CORNER) are different from what stored 
in the GeoTIFF file, the AROP package will use the values retrieved from the 
GeoTIFF file. 

4. The MATCH_BAND for warp and out images should use a single band with similar 
spectral bandwidth to the BASE_LANDSAT. We recommend selecting NIR or SWIR 
band (e.g., TM/ETM+ band 4 or 5) for a better control points searching process. 

5. The AROP package allows input warp image, GeoCover Landsat and SRTM DEM 
data in different spatial resolutions and projection. Fine resolution DEM data will 
produce a more accurate result. Use 1 arcsec complete (gap-filled) SRTM DEM data if 
it is available.

6. The AROP package accepts inputs in per band per file format. Outputs include 
image files in binary format each associated an ENVI file header. You can open output 
images in ENVI software and link to the base image for visual overlay checking.

7. You may run "-v" option first and see if how the images match each other and then 
decide if further co-registration process is needed. 


IV. TESTING

Landsat and DEM data for the following tests (under sub-directory tests/) can be 
downloaded from the Global Land Cover Facility (GLCF) at the university of 
Maryland through http://glcfapp.glcf.umd.edu:8080/esdi/index.jsp
 
1. ETM+ testing

This test runs orthorectification for a Landsat ETM+ scene which only requires one 
step precise registration by shifting upper left coordinates.
 
search ETM+ scene with path = 17 and row = 34
* download GeoCover ETM+ data: p017r034_7t20000610_z17_nn50.tif.gz 
(6/10/2000, 
ETM+)
* download 1 arcsec SRTM elevation data: SRTM_u01_p017r034.tif.gz
* download USGS L1G Landsat ETM+ data: L7?017034_03420001101_B?0.L1G.gz 
(11/1/2000, ETM+) 

unzip files and then run 
../source/ortho -b ETM_example.inp
You can open orthorectified outputs directly using ENVI.

2. TM testing

This test runs orthorectification for a Landsat TM scene which involves re-
orthorectification process by using first-order polynomial transformation and
orthorectification after geo-location checking fails.

search TM scene with path = 15 and row = 33
* download GeoCover TM data: p015r33_5t870516_nn5.tif  (05/16/1987, TM)
* download 1 arcsec SRTM elevation data: SRTM_u01_p015r033.tif
* download USGS L1G Landsat TM data: LT5015033009022410.I?.gz (08/12/1990, 
TM)
(this TM image comes with -11.722845 rotation)

unzip files and then run
../source/ortho -b TM_example.inp

3. MSS testing

This test performs orthorectification for a Landsat MSS scene using GeoCover MSS 
image as a base scene. 

search MSS scene with path 16 and row 33 (in WRS-1)
* download GeoCover MSS data: p016r33_1m19721011_03.tif.gz (10/11/1972, MSS
* download MSS L1G data: LM2016033008119290.I?.gz (7/11/1981, MSS)
search SRTM DEM data with path 15 and row 33 (in WRS-2)
* download 1 arcsec SRTM DEM data: SRTM_u01_p015r033.tif

unzip files and then run
../source/ortho -b MSS_example1.inp 
(use GeoCover MSS data as base)
 
OR

use GeoCover ETM+ data as base and 
dowload GeoCover ETM+ data: p015r033_7t20011005_z18_nn40.tif (10/05/2001, 
ETM+) 

../source/ortho -b MSS_example2.inp

4. Orthorectification only testing

This test orthorectify image without registration processing
(not recommend for the image that has not been precisely registrated)

search ETM+ scene with path = 17 and row = 34
* download 1 arcsec SRTM elevation data: SRTM_u01_p017r034.tif.gz
* download USGS L1G Landsat ETM+ data: L7?017034_03420001101_B?0.L1G.gz 
(11/01/2000) 

unzip files and then run 
../source/ortho -o ETM_ortho_only.inp

5. Precise registration only testing

This test co-register two orthorectified (or unorthorectified) images

GeoCover WRS1 MSS  p16r33  (10/11/1972, MSS)
GeoCover WRS2 ETM+ p15r33  (10/05/2001, ETM+)
Register ETM+ to MSS and aggregte output from 28.5m to 57m resolution
../source/ortho -r ETM_MSS_coreg.inp

Register MSS to ETM+ and resample from 57m to 28.5m resolution
../source/ortho -r MSS_ETM_coreg.inp


6. Rotation, reprojection and orthorectification testing

This test shows a combined series operations on the warp image that is saved in a 
rotated SOM projection.

search TM scene with path = 15 and row = 33
* download GeoCover TM data: p015r33_5t870516_nn5.tif (05/16/1987, TM)
* download 1 arcsec SRTM elevation data: SRTM_u01_p015r033.tif
* download USGS L1G Landsat TM data: LT5015033009628910.B* (10/15/1996, 
TM)
 
unzip files and then run
../source/ortho -b TM_example_reproject.inp


7. ASTER data testing

This test (ASTER_example.inp) includes the necessary parameters for 
* map rotation (from swath orientation to North)
       WARP_ORIENTATION_ANGLE = -9.049154
* map reprojection (in different UTM zones)
       WARP_PROJECTION_CODE = 1
       WARP_UTM_ZONE = 1
       WARP_UNIT = 2
* ASTER sensor pointing angle (central line is not nadir)
       WARP_SATELLITE_POINTINGANGLE = 5.674000
* upper left coordinate in degrees
       WARP_UPPER_LEFT_CORNER_DEGREE = -78.509630, 37.575253
(note that ASTER L1B HDF file need to be converted to binary format
first and be saved in one file per band. All VNIR and SWIR bands are saved
in 15m resolution)


8. Mixed data type (ASTER)

This example (ASTER_mixed_example.inp) demonstrates the process on the mixed 
data type. In the example, we have ASTER VNIR and SWIR bands (1-9) stored in 1-
byte data type and TIR bands (10-14) stored in 2-byte data type. For a 
ortho process dealing with the mixed data type, we need to define
        WARP_BAND_DATA_TYPE = 1 1 1 1 1 1 1 1 1 2 2 2 2 2
in the input parameter file. The output band will be stored in the same data type as 
input band.


9. DEM data with different projection

This example (DEM_reprojection.inp) shows how to process DEM data in different 
projection (different UTM zone here). The input file needs define 
	DEM_PROJECTION_CODE = 1
	DEM_UTM_ZONE = 16 
For input with different projection, the projection parameters need to be defined
	DEM_PROJECTION_PARAM = (15 gctpc parameters)

10. Image matching verification

User may run program with option -v for image matching verification. This example 
(verification_example.inp) will not create actual output image files. Instead program 
will verify images matching from four quadrants using automatically selected tie points 
and display them on screen.


Appendix

Complete List of Parameters 

PARAMETER_FILE
###################
# define base image 
###################
# define input file type, use GEOTIFF or BINARY
BASE_FILE_TYPE = 
# number of samples (columns) or use value from GEOTIFF file (-1)
BASE_NSAMPLE = 
# number of lines (rows) or use value from GEOTIFF file (-1)
BASE_NLINE = 
# pixel size in meters
BASE_PIXEL_SIZE = 
# upper left coordinate in meters or use value from GEOTIFF file (-1)
BASE_UPPER_LEFT_CORNER = 
# base image for "r" and "b" options (one band)
BASE_LANDSAT = 
# base image must be in UTM projection for this version
# use positive UTM for North and negative for South
BASE_UTM_ZONE = 
# define datum for input base image
BASE_DATUM = 
# supported satellites: Landsat 1-5 and Landsat 7
BASE_SATELLITE =
####################
# define warp images
####################
# define input file type, use GEOTIFF or BINARY
WARP_FILE_TYPE = 
# number of samples (columns)
WARP_NSAMPLE = 
# number of lines (rows)
WARP_NLINE = 
# pixel size in meters
WARP_PIXEL_SIZE = 
# upper left coordinate can be in meters or degrees (only need one)
# define upper left coordinate in meters
# WARP_UPPER_LEFT_CORNER =
# define upper left coordinate in degrees
#WARP_UPPER_LEFT_CORNER_DEGREE = 
# supported satellites: Landsat 1-5, Landsat 7, TERRA, CBERS1, CBERS2, AWIFS
WARP_SATELLITE = 
# sensor pointing angle in degrees (for ASTER)
WARP_SATELLITE_POINTINGANGLE = 
# map orientation angle in degrees for warp image
WARP_ORIENTATION_ANGLE = 
# number of input bands 
WARP_NBANDS = 
# list each band filename separated by comma or space
WARP_LANDSAT_BAND = 
# define data type for each band in: 1=8-bit; 2=16-bit
WARP_BAND_DATA_TYPE = 
# define matching warp band to BASE_LANDSAT
WARP_BASE_MATCH_BAND = 
# define projection for warp image if different from base image
# don't need these information if warp image projection is same as base 
# define projection in GCTPC format (see GCTPC documents for details)
# 0=GEO; 1=UTM (default: UTM)
WARP_PROJECTION_CODE = 
# use positive UTM for North and negative for South
WARP_UTM_ZONE = 
# 15 GCTPC projection parameters (default: 0.0; not required for UTM)
# WARP_PROJECTION_PARAM = 
# 0=radians; 1=US feet; 2=meters; 3=seconds of arc; 4=degree of arc
# 5=international feet (default: meters)
WARP_UNIT = 
# 0=Clarke 1866; 8=GRS 1980; 12=WGS 84 (default: WGS 84)
WARP_DATUM = 
######################
# define output images
######################
# define output pixel resolution
OUT_PIXEL_SIZE = 
# define resampling approach (NN, BI, CC, AGG) 
# NN for nearest neighbor, 
# BI for bilinear interpolation
# CC for cubic convolution 
# AGG for pixel aggregation
RESAMPLE_METHOD =
# define image extent for the output (BASE, WARP, DEF) 
# BASE uses base map extent
# WARP uses warp map extent
# DEF takes user defined map extent 
OUT_EXTENT = 
# define background or fill value (default is 0)
BK_VALUE = 
# if DEF is defined, the following two lines need to be defined (in meters)
# OUT_UPPER_LEFT_CORNER =	   
# OUT_LOWER_RIGHT_CORNER =	   
# define corresponding output files for each band separated by comma
OUT_LANDSAT_BAND = 
# define one corresponding output matching band for geolocation verification
# define one band matching to BASE_LANDSAT  
OUT_BASE_MATCH_BAND = 
# the maximum degree of polynomial transformation (0, 1, 2)
# note that 2nd degree is not recommended unless have to 
OUT_BASE_POLY_ORDER = 
# ancillary inputs for orthorectification process
# define terrain elevation file (must be in GeoTIFF format)
# the SRTM DEM data in GeoTIFF format can be downloaded from the UMD GLCF  
INPUT_DEM_FILE = 
# define projection for DEM data if it's different from base image
# projection information is not needed if projection for DEM data is same as base 
image 
# define projection in GCTPC format (see GCTPC documents for details) (default: 
UTM)
# DEM_PROJECTION_CODE =
# use positive UTM for North and negative for South
# DEM_UTM_ZONE =
# 15 GCTPC projection parameters (default: 0.0)
# DEM_PROJECTION_PARAM =
# 0=radians; 1=US feet; 2=meters; 3=seconds of arc; 4=degree of arc
# 5=international feet (default: meters)
# DEM_UNIT = 
# 0=Clarke 1866; 8=GRS 1980; 12=WGS 84 (default: WGS 84)
# DEM_DATUM =
# tie point searching control parameters file (defined separately)
CP_PARAMETERS_FILE = lndortho.cps_par.ini
END


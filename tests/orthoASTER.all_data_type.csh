#!/bin/csh -f

# convert ASTER L1B HDF file to binary file and do orthorectification 

# set command 
set sds2bin = "/home/fgao/Tools/LDOPE_Linux/MODAT/bin/sds2bin"
set enlarge_sds = "/home/fgao/Tools/LDOPE_Linux/MODAT/bin/enlarge_sds"
set ortho = "/home/fgao/Tools/ASTER_ortho_v2.2.2/lndortho/source/ortho"
set ncdump = "/usr/local/bin/ncdump"
set listgeo = "/usr/local/bin/listgeo"

# parse input parameters
if $#argv != 2 then
    echo "Usage: orthoASTER.csh <ASTER_L1B_hdf_file> <BASE_Ortho_Landsat_file>"
else
    set aster = $argv[1]
    set landsat = $argv[2]
endif

$ncdump -h $aster > tmp.hdr
set nrows = `grep "ImageLine_VNIR_Swath" tmp.hdr | awk 'NR==1 {printf("%d", $3)}'`
set ncols = `grep "ImagePixel_VNIR_Swath" tmp.hdr | awk 'NR==1 {printf("%d", $3)}'`

set angle = `grep "MAPORIENTATIONANGLE" -C 1 tmp.hdr | awk 'NR==4 {printf("%f", $4)}'`
set ul = `grep -C 1 "UPPERLEFT" tmp.hdr | awk 'NR==4 {print $4 $5}' | sed 's/(//' | sed 's/)\\n",//' | sed 's/,/ /'`
set uly = `echo $ul | awk '{printf("%f", $1)}'`
set ulx = `echo $ul | awk '{printf("%f", $2)}'`
set warp_zone = `grep -C 1 "UTMZONECODE1" tmp.hdr | awk 'NR==4 {printf("%d", $4)}'`
set pangle = `grep "POINTINGANGLE\\n" -C 1 tmp.hdr | awk 'NR==5 {printf("%f", $4)}'`
rm -f tmp.hdr

# create binary file from ASTER L1B HDF file in 15m resolution
set out_stem = `echo $aster | sed 's/.hdf//'`
$sds2bin -sds="ImageData1" -of=$out_stem.b1 $aster
$sds2bin -sds="ImageData2" -of=$out_stem.b2 $aster
$sds2bin -sds="ImageData3N" -of=$out_stem.b3 $aster
# need enlarge SWIR bands from 30m resolution to 15m for one ortho processing
set tempout="ImageData_x2.hdf"
set swir = "ImageData4,ImageData5,ImageData6,ImageData7,ImageData8,ImageData9"
$enlarge_sds -sds=$swir -sf=2 -of=$tempout $aster
$sds2bin -sds="ImageData4_2x" -of=$out_stem.b4 $tempout
$sds2bin -sds="ImageData5_2x" -of=$out_stem.b5 $tempout
$sds2bin -sds="ImageData6_2x" -of=$out_stem.b6 $tempout
$sds2bin -sds="ImageData7_2x" -of=$out_stem.b7 $tempout
$sds2bin -sds="ImageData8_2x" -of=$out_stem.b8 $tempout
$sds2bin -sds="ImageData9_2x" -of=$out_stem.b9 $tempout
rm -f $tempout
# need enlarge TIR bands from 90m resolution to 15m for one ortho processing
set tempout="ImageData_x6.hdf"
set tir = "ImageData10,ImageData11,ImageData12,ImageData13,ImageData14"
$enlarge_sds -sds=$tir -sf=6 -of=$tempout $aster
$sds2bin -sds="ImageData10_6x" -of=$out_stem.b10 $tempout
$sds2bin -sds="ImageData11_6x" -of=$out_stem.b11 $tempout
$sds2bin -sds="ImageData12_6x" -of=$out_stem.b12 $tempout
$sds2bin -sds="ImageData13_6x" -of=$out_stem.b13 $tempout
$sds2bin -sds="ImageData14_6x" -of=$out_stem.b14 $tempout
rm -f $tempout

# create ortho input card
set tempout="$out_stem.ortho.inp"

set res=`$listgeo $landsat | awk 'NR==9 {printf("%f", $1)}'`
set base_zone=`$listgeo $landsat | grep "UTM zone" | awk 'NR==1 {printf("%d",$9)}'`

echo "PARAMETER_FILE" > $tempout
echo "# LANDSAT base image" >> $tempout
echo "BASE_FILE_TYPE = GEOTIFF" >> $tempout
echo "BASE_SATELLITE = Landsat7" >> $tempout
echo "BASE_NSAMPLE = -1" >> $tempout
echo "BASE_NLINE = -1" >> $tempout
echo "BASE_PIXEL_SIZE = $res" >> $tempout
echo "BASE_UPPER_LEFT_CORNER = -1, -1" >> $tempout
echo "BASE_LANDSAT = $landsat" >> $tempout
echo "UTM_ZONE = $base_zone[1]" >> $tempout
echo "BASE_SATELLITE = Landsat7" >> $tempout
echo "" >> $tempout

echo "# Landsat warp images" >> $tempout
echo "WARP_FILE_TYPE = BINARY" >> $tempout
echo "WARP_NSAMPLE = $ncols" >> $tempout
echo "WARP_NLINE = $nrows" >> $tempout
echo "WARP_PIXEL_SIZE = 15.0" >> $tempout
echo "WARP_UPPER_LEFT_CORNER_DEGREE = $ulx, $uly" >> $tempout
echo "# Landsat 1-5, Landsat 7, TERRA, CBERS2 #" >> $tempout
echo "WARP_SATELLITE = TERRA" >> $tempout
echo "WARP_SATELLITE_POINTINGANGLE = $pangle" >> $tempout
echo "WARP_ORIENTATION_ANGLE = -$angle" >> $tempout
echo "WARP_PROJECTION_CODE = 1" >> $tempout
echo "WARP_UTM_ZONE = $warp_zone" >> $tempout
echo "WARP_UNIT = 2" >> $tempout
echo "WARP_NBANDS = 14" >> $tempout
echo "WARP_LANDSAT_BAND = $out_stem.b1 $out_stem.b2 $out_stem.b3 $out_stem.b4 $out_stem.b5 $out_stem.b6 $out_stem.b7 $out_stem.b8 $out_stem.b9 $out_stem.b10 $out_stem.b11 $out_stem.b12 $out_stem.b13 $out_stem.b14" >> $tempout
echo "WARP_BAND_DATA_TYPE = 1 1 1 1 1 1 1 1 1 2 2 2 2 2" >> $tempout
echo "WARP_BASE_MATCH_BAND = $out_stem.b4" >> $tempout

echo "" >> $tempout
echo "# Landsat orthorectied output images" >> $tempout
echo "OUT_PIXEL_SIZE = $res" >> $tempout
echo "# NN-nearest neighbor; BI-bilinear interpolation; CC-cubic convolution; AGG-aggregation #" >> $tempout
echo "RESAMPLE_METHOD = AGG" >> $tempout
echo "# BASE-use base map extent; WARP-use warp map extent; DEF-user defined extent #" >> $tempout
echo "OUT_EXTENT = WARP" >> $tempout
echo "OUT_LANDSAT_BAND = Ortho.$out_stem.b1 Ortho.$out_stem.b2 Ortho.$out_stem.b3 Ortho.$out_stem.b4 Ortho.$out_stem.b5 Ortho.$out_stem.b6 Ortho.$out_stem.b7 Ortho.$out_stem.b8 Ortho.$out_stem.b9 Ortho.$out_stem.b10 Ortho.$out_stem.b11 Ortho.$out_stem.b12 Ortho.$out_stem.b13 Ortho.$out_stem.b14" >> $tempout
echo "OUT_BASE_MATCH_BAND = Ortho.$out_stem.b4" >> $tempout

set dem_file=`ls SRTM*.tif`

listgeo $dem_file > tmp.hdr
set dem_zone = `grep "UTM zone" tmp.hdr | sed "s/N)//" | awk '{printf("%s", $9)}'`
rm -f tmp.hdr

echo "" >> $tempout
echo "# ancillary input for orthorectification process" >> $tempout
echo "INPUT_DEM_FILE = $dem_file" >> $tempout
echo "DEM_PROJECTION_CODE = 1" >> $tempout
echo "DEM_UTM_ZONE = $dem_zone" >> $tempout
echo "CP_PARAMETERS_FILE = aster_ortho.cps_par.ini" >> $tempout
echo "" >> $tempout
echo "END" >> $tempout

$ortho -b $tempout
foreach band(1 2 3 4 5 6 7 8 9 10 11 12 13 14)
rm -f $out_stem.b$band
end

# write a ENVI meta to include all ortho bands for easier loading
set out_meta = Ortho.$out_stem.meta
set samples = `grep "samples" Ortho.$out_stem.b1.hdr | awk '{printf("%d", $3)}'`
set lines = `grep "lines" Ortho.$out_stem.b1.hdr | awk '{printf("%d", $3)}'`
echo "ENVI META FILE" > $out_meta
foreach band(1 2 3 4 5 6 7 8 9 10 11 12 13 14)
echo "FILE : Ortho.$out_stem.b$band" >> $out_meta
echo "Bands: 1" >> $out_meta
echo "Dims : 1-$samples,1-$lines" >> $out_meta
echo "" >> $out_meta
end 




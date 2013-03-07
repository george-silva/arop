/**   
 * ! Description
 *   
 *   This is the main program of the automated registration and orthorectification package
 *   (AROP) for Landsat and Landsat-like data processing. It has been tested on the Landsat
 *   MSS/TM/ETM+, TERRA ASTER, CBERS CCD and IRS-P6 AWiFS data. The algorithm details
 *   are described in the following paper as well as readme file in this package  
 *
 *   Gao., F., Masek, J. G., Wolfe, R. F., 2008, "An automated registration and orthorectification
 *   package (AROP) for Landsat and Landsat-like data processing," Journal of Applied Remote
 *   Sensing, submitted. 
 *
 *   The AROP has been tested in Linux system.
 *
 * ! Input
 *   - Geocover Landsat data as base image
 *   - warp mid-resolution (Landsat) satellite images
 *   - SRTM DEM data ( Landsat WRS2 scene in geotiff format from UMD GLCF)
 *   - Landsat WRS2 coordinates file (a text file included in the package)
 *   - extent of output orthorectified image (default use input base image)
 *   - output spatial resolution
 *   - output resampling method
 *
 * ! Ouput
 *   - registered and/or orthorectified image that matches GeoCover base image
 *
 * ! Usage
 *
 *   lndortho [-r|-o|-b|-v] input_parameter_file
 *   -r registration only
 *   -o orthorectification only
 *   -b combined registration and orthorectification operations
 *   -v image matching verification
 * 
 * ! Credits
 *
 *   Feng Gao (Feng.Gao@nasa.gov), Developer
 *   Jeff Masek (Jeffrey.G.Masek@nasa.gov)
 *   Robert Wolfe (Robert.E.Wolfe@nasa.gov)
 *
 * ! Revision History
 *
 *   Revision 2.2.6  8/22/2011  Feng Gao
 *   - included image matching checking option "-v" (no image output) 
 *   - added parameter for HJ-1 satellite (WARP_SATELLITE = HY1)
 *
 *   Revision 2.2.5  1/27/2010  Feng Gao
 *   - included datum input for base image, so this package will revoke reprojection if datum 
 *     from base, warp or dem are different
 *   - adjusted initial tie point to the nearby point with the maximum stdev to avoid selcting 
 *     control chip on a "flat" homogeneous area 
 *   - A tip for choosing tie points over a homogeneous area: increase size of the control chip 
 *   - fixed a bug in using warp file with absolute directory (temporary file errors)
 *   - adjusted CP_SEED_WIN according to the image size to avoid too few tie points
 *
 *   Revision 2.2.4  8/17/2009  Feng Gao
 *   - fixed potential leaks/bugs in reading 1-byte or 2-byte data type 
 *   - adjusted reduction scale automatically for coarse image if coarse image is too small
 *     for pyramid registration which may cause too few tie points for registration.  
 *
 *   Revision 2.2.3  2/27/2009  Feng Gao
 *   - excluded saturated (and fill value) pixels in BI, CC and AGG resmapling
 *
 *   Revision 2.2.2  10/10/2008 Feng Gao
 *   - allow 1 byte or 2 bytes data type as input data format (the matching base and warp band
 *     should be still 1 byte data type). This is revised for the ASTER data processing which
 *     includes both 1 byte (VNIR and SWIR, band 1-9) and 2 bytes (TIR, band 10-14) data type.
 *   - improved result verification by using absolute error instead of relative error for the
 *     tie points (more direct and reliable)
 *   - added capability in accepting SRTM DEM data in different projection from the BASE image 
 *
 *   Revision 2.2.1  6/03/2008  Feng Gao
 *   - removed Landsat WRS2 definition; all nadir track will be detected from the warp image
 *     This requires the warp image in the original swath (SUBSET of the warp image will NOT work)  
 *   - removed Landsat specific parameters from input file (better for Landsat-like data processing)
 *     LANDSAT_WRS2_FILE = 
 *     PATH =
 *     ROW =
 *   - Increased the processing limit of image size 
 * 
 *   Revision 2.2.0  6/02/2008  Feng Gao
 *   - implemented combined single resampling strategy for warp image that needs rotation 
 *     and/or reprojection for both "r" and "b" options. 
 *   - fixed a bug in computing the nadir track of the warp image in the base image coordinate 
 *      
 *   Revision 2.1.0  5/09/2008  Feng Gao
 *   - fixed a bug in computing nadir track of warp image in output space
 *     (previous versions assume the nadir track of warp image in output space always passes 
 *      the center of output image which may not true if output image extent is not same as
 *      WARP image. In this version, we convert central pixel of warp image to output space 
 *      and compute nadir track in output space using central location of warp image)
 *   - write DATUM (WGS-84 in SRTM DEM) information in ENVI header
 *   - implemented upper left corner input in degrees (for ASTER processing)  
 *     WARP_UPPER_LEFT_CORNER_DEGREE = $lon, $lat
 *   - implemented sensor pointing angle for ASTER data in extracting nadir track function 
 *     WARP_SATELLITE_POINTINGANGLE = 
 *   - added satellite information for BASE image such allowing base and warp images from 
 *     different satellite platforms (different satellite altitude)
 *     BASE_SATELLITE =
 *   - used actual center pixel of scene (both base and warp) to compute earth radius 
 *     instead of predefined WRS2 coordinates which is not appropriate for non-WRS2 scenes 
 *   - always computed nadir track of output/warp image detected from warp image 
 *     (previous versione used WRS2 coordinates if supplied. This is good for Landsat WRS2 
 *     scenes but may not good for other sensors) 
 *   - revised tie points distribution checking from base image to output image
 *     (better explanation on the ortho and registration results)
 * 
 *   Revision 2.0.0  4/15/2008  Feng Gao
 *   - clean code and document for public release
 *
 *   Revision 2.0.0 (Beta1) 12/21/2007  Feng Gao
 *   - fixed a bug while defined area is much smaller than input area (ouput subset operation)
 *   - normalized coordinate system for inputs and output to the UL of the UL pixel
 *     so output uses upper left coordinate of a pixel  
 *     (GeoCover data used center of pixel while matching SRTM DEM used UL of pixel)
 *
 *   Revision 2.0.0 (Beta0) 11/30/2007  Feng Gao
 *   - restructured code for easier maintainence
 *   - added option for 
 *       -r do registration only
 *       -o do orthorectification only
 *       -b do both (same as version 1.9.0 and before)
 *
 *   Revision 1.9.0  11/23/2007  Feng Gao
 *   - accepts warp image that needs rotation by defining
 *     WARP_ORIENTATION_ANGLE = 
 *   - allows applying to different satellite type by defining 
 *     WARP_SATELLITE = 
 *   - allows using different resolution image as base (such as using GeoCover TM/ETM+ as base for MSS)
 *   - allows defining output image in a different spatial resolution 
 *     so base, warp and output images can be in different spatial resolutions
 *   - added aggregation option for output for fine input warp images 
 *
 *   Revision 1.8.3  11/09/2007  Feng Gao
 *   - adjusted maximum searching distance (MAX_SHIFT) after registration on coarse images 
 *   - forced to run preliminary registration if MAX_SHIFT > 100
 *
 *   Revision 1.8.2  11/06/2007  Feng Gao
 *   - fixed memory leak problem for BI and CC resampling options 
 *
 *   Revision 1.8.1  10/29/2007  Feng Gao
 *   - fixed junk line output when resampled data is outside the input imagery
 *
 *   Revision 1.8.0  7/31/1007   Feng Gao
 *   - added capability to define output extent from input using following key words
 *     OUT_EXTENT = DEF
 *     OUT_UPPER_LEFT_CORNER = 
 *     OUT_LOWER_RIGHT_CORNER = 
 *
 *   Revision 1.7.0  11/12/2006  Feng Gao
 *   - added capability to do preliminary registration based on the aggregated coarse resolution images
 *     This solves problem of precise registration on images with large geolocation error such as some MSS images
 *   - added definition of control point searching for preliminary registration in lndortho.cps_par.ini  
 *
 *   Revision 1.6.1  10/26/2006  Feng Gao
 *   - if no enough control points can be found in the first attemp, quit program and ask to revise
 *     initial parameter file to increase number of control point seeds and/or increase searching 
 *     range
 *   - strict test of matching points with revised prompt 
 * 
 *   Revision 1.6    10/5/2006   Feng Gao
 *   - removed thresholds definition from program to parameter file "lndortho.cps_par.ini"
 *     and thus allows users to process specific scenes that do not work with default values 
 *     by revising parameter file without changing and re-compiling program. 
 *   - relax geolocation checking (exclude checking on the regions which has less than 5 gcps)
 *
 *   Revision 1.5    4/11/2006   Feng Gao
 *   - added second order polynomial precise registration for MSS scene processing (default use first order)
 *   - added analysis of nadir path track from image directly if scene coordinate file (WRS-1) is not included 
 *   - fixed bugs in computing correlation coefficient & redo processing loop
 *
 *   Revision 1.4    2/16/2006   Feng Gao
 *   - combined precise registration process (lndregortho) and orthorectification process (lndortho)
 *   - revised code to accept L1G landsat (binary or GeoTIFF) input for public release 
 *
 *   Revision 1.3    11/23/2005   Feng Gao
 *   - implemented sub-pixel accuracy matches using neighbor pixel cross correlation values
 *     (used approach in EOSAT IDPS Theoretical Supplement by Robert Wolfe)
 *
 *   Revision 1.2    11/18/2005   Feng Gao
 *   - added recheck function after orthorectification to see if output is really match with base image
 *   - added first order geolocation rectification if output doesn't match with base image 
 *   - do once resampling for both first order geolocation rectification and orthorectification if necessary
 *   - include thermal band (brightness temperature) in output orthorectified HDF file
 *
 *   Revision 1.1    11/1/2005    Feng Gao
 *   - added three sampling methods (NN, BI and CC)
 *   - added variable output resolutions 
 *   - changed processing to deal with output from forward sampling to backward sampling, so it can process
 *     different output sampling methods and variable output spatial resolution
 *
 *   Revision 1.0     09/2005     Feng Gao
 *   Original version  
 *   - accepts different resolution DEM data from SRTM
 *   - can produce output in defined extent
 *   - output has same resolution as input Landsat data 
 */

#include "lndortho.h"

int main(int argc, char *argv[])
{
  int  i, ret=SUCCESS;

  /* all input parameters */
  IN_PARAMETERS *inPars;

  /* SRTM dem data covering whole area */
  DEM   *sdem;

  /* define actual input and output space */
  BASE_LANDSAT  *base;
  WARP_LANDSAT  *warp;
  OUT_LANDSAT   *out;
  
  /* temporary warp image that has been rotated and reprojected */ 
  /* contains only one matching band */
  WARP_LANDSAT  *temp_warp, *rotate_warp, *proj_warp;

  /* define working space */
  BASE_LANDSAT  *wbase;
  WARP_LANDSAT  *wwarp;
  OUT_LANDSAT   *wout;

  if(argc!=3 || (strcmp(argv[1], "-r")!=0 && strcmp(argv[1], "-o")!=0 && strcmp(argv[1], "-b")!=0 && strcmp(argv[1], "-v")!=0)) {
    usage(argv);
    exit(1);
  }

  if(strcmp(argv[1], "-r") == 0)
    fprintf(stderr, "\nStarting Automated Registration for %s\n", argv[2]);
  else if(strcmp(argv[1], "-o") == 0)
    fprintf(stderr, "\nStarting Orthorectification for %s\n", argv[2]);
  else if(strcmp(argv[1], "-b") == 0)
    fprintf(stderr, "\nStarting Combined Automated Registration and Orthorectification for %s\n", argv[2]);
  else if(strcmp(argv[1], "-v") == 0)
    fprintf(stderr, "\nStarting Image Matching Verification for %s\n", argv[2]);

  /* allocate early variables */
  if(!(inPars = malloc(sizeof(IN_PARAMETERS)))) exit(1);
  if(!(base = malloc(sizeof(BASE_LANDSAT)))) exit(1);
  if(!(warp  = malloc(sizeof(WARP_LANDSAT)))) exit(1);
  if(!(out = malloc(sizeof(OUT_LANDSAT)))) exit(1);
  if(!(temp_warp  = malloc(sizeof(WARP_LANDSAT)))) exit(1);
  if(!(rotate_warp  = malloc(sizeof(WARP_LANDSAT)))) exit(1);
  if(!(proj_warp  = malloc(sizeof(WARP_LANDSAT)))) exit(1);
  if(!(wbase = malloc(sizeof(BASE_LANDSAT)))) exit(1);
  if(!(wwarp  = malloc(sizeof(WARP_LANDSAT)))) exit(1);
  if(!(wout = malloc(sizeof(OUT_LANDSAT)))) exit(1);
  if(!(sdem = malloc(sizeof(DEM)))) exit(1);


  /*** READ AND PREPARE DATA ***/  
  /* get input parameters from command line */
  if(getInParameter(inPars, argc, argv)==FAILURE) {
    fprintf(stderr, "\nRetrieve input parameters error!\n");
    exit(1);
  }

  printf("\n- retrieving parameters from inputs");

  /* retrieve metadata from SRTM DEM data */
  if(inPars->option == ORTHO_ONLY || inPars->option == DO_BOTH) {
    strcpy(sdem->fileName, inPars->dem_fileName);
    if(getSRTMMetaInfo(sdem, inPars)==FAILURE) {
      fprintf(stderr, "Retrieve %s SRTM GeoTiff file error!\n", sdem->fileName);
      exit(1);
    }
    printf("\n\tDEM file:\n\t\t%s", sdem->fileName);
    printf("\n\t\tnrows = %d\tncols = %d", sdem->nrows, sdem->ncols);
    printf("\n\t\tulx = %8.1f\tuly = %8.1f\n\t\tres = %5.1f", sdem->ulx, sdem->uly, sdem->res);
    if(sdem->proj.sys == UTM)
      printf("\n\t\tprojection = UTM\n\t\tutm_zone = %d", sdem->proj.utm_zone);
    else {
      printf("\n\t\tprojection = %ld", sdem->proj.sys);
      printf("\n\t\tparameters[15] = ");
      for(i=0; i<15; i++)
	printf("%f ", sdem->proj.param[i]);
    }
  } /* endif ORTHO and BOTH */

  /* retrieve Landsat metadata from input parameters */
  if(getLandsatMetaInfo(inPars, base, warp, out)==FAILURE) {
    fprintf(stderr, "Retrieve input parameters error!\n");
    exit(1);
  } 

  if(inPars->option == REG_ONLY || inPars->option == DO_BOTH || inPars->option == VERIFY_ONLY) {
    /* print out base file info */
    printf("\n\tInput BASE file:\n\t\t%s",base->lnd.fileName);
    printf("\n\t\tnrows = %d;\tncols = %d", base->lnd.nrows, base->lnd.ncols); 
    printf("\n\t\tulx = %8.1f;\tuly = %8.1f", base->lnd.ulx, base->lnd.uly);
    printf("\n\t\tcenx =%8.3f;\tceny =%8.3f (degrees)", base->lnd.centerLonLat[0], base->lnd.centerLonLat[1]);
    printf("\n\t\tres = %5.1f;\tfillV = %d", base->lnd.res, base->lnd.fillValue);
    printf("\n\t\tutm_zone = %d; \tdatum = %ld", base->utm_zone, base->datum);
  } else {
    /* copy warp to base for ORTHO_ONLY option if base is not defined 
       (useless just keeps program going) */
    if(strlen(base->lnd.fileName)==0) {
      copyLandsat(&(base->lnd), &(warp->lnd));
      strcpy(base->lnd.fileName, warp->fileName[0]);
      strcpy(warp->lnd.fileName, warp->fileName[0]);
    }
  }

  /* print out warp file info */
  printf("\n\tInput WARP file:");
  for(i=0; i<warp->nbands; i++) 
    printf("\n\t\t%s", warp->fileName[i]);
  printf("\n\t\tnrows = %d;\tncols = %d", warp->lnd.nrows, warp->lnd.ncols); 
  printf("\n\t\tulx = %8.1f;\tuly = %8.1f", warp->lnd.ulx, warp->lnd.uly);
    printf("\n\t\tcenx =%8.3f;\tceny =%8.3f (degrees)", warp->lnd.centerLonLat[0], warp->lnd.centerLonLat[1]);
  printf("\n\t\tres = %5.1f;\tfillV = %d;", warp->lnd.res, warp->lnd.fillValue);
  printf("\n\t\tsatellite = %s", inPars->warp_satellite);
  printf("\n\t\torientation_angle = %f", warp->orientationAngle);
  printf("\n\t\tsensor_pointing_angle = %f", warp->pointingAngle);
  if(warp->proj.sys == UTM)
    printf("\n\t\tprojection = UTM\n\t\tutm_zone = %d\tdatum = %ld", warp->proj.utm_zone, warp->proj.datum);
  else {
    printf("\n\t\tprojection = %ld", warp->proj.sys);
    printf("\n\t\tparameters[15] = ");
    for(i=0; i<15; i++)
      printf("%f ", warp->proj.param[i]);
  }

  if(inPars->option != VERIFY_ONLY) {
    printf("\n\t\tresample_method = %s", warp->resampleMethod);
    printf("\n\t\thigher order polynomial transformation = %d", out->lnd.porder);

    /* print out output file info */
    printf("\n\tOutput orthorectified OUT file:");
    for(i=0; i<warp->nbands; i++)
      printf("\n\t\t%s", out->fileName[i]);
    printf("\n\t\tnrows = %d;\tncols = %d", out->lnd.nrows, out->lnd.ncols); 
    printf("\n\t\tulx = %8.1f;\tuly = %8.1f", out->lnd.ulx, out->lnd.uly);
    printf("\n\t\tres = %5.1f;\tfillV = %d", out->lnd.res, out->lnd.fillValue);
  }

  /* create temporay warp image to store rotated and reprojected matching band */
  if(createTempWarp(temp_warp, warp)==FAILURE) {
    fprintf(stderr, "creating temporary warp image error\n");
    exit(1);
  }
 
  /*** DO PRE_PROCESSING ON ROTATION IF REQUIRED ***/
  if(inPars->rotation_flag) {
    printf("\n\n- rotating warp images");
    /* do rotation for warp image if necessary */
    if(rotateWarp(temp_warp)==FAILURE) {
      fprintf(stderr, "Rotate warp image error!\n");
      exit(1);
    }
    printf("\n\t\tnrows = %d;\tncols = %d", temp_warp->lnd.nrows, temp_warp->lnd.ncols); 
    printf("\n\t\tulx = %8.1f;\tuly = %8.1f", temp_warp->lnd.ulx, temp_warp->lnd.uly);
    /* if extent is defined same as warp, then extent must be updated */
    if(strcasecmp(inPars->extentSource, "WARP")==0)
      updateOutCoor(&(temp_warp->lnd), &(out->lnd));
    /* save rotated warp image for later combined resampling */
    if(createTempWarp(rotate_warp, temp_warp)==FAILURE) {
      fprintf(stderr, "creating temporary warp image error\n");
      exit(1);
    }
  }

  /*** DO PRE_PROCESSING ON REPROJECTION IF REQUIRED ***/
  if(inPars->projection_flag) {
    printf("\n\n- converting projection for warp images");
    /* do reprojection for warp image if necessary */
    if(reproject(temp_warp, base)==FAILURE) {
      fprintf(stderr, "Reproject warp image error!\n");
      exit(1);
    }
    printf("\n\t\tnrows = %d;\tncols = %d", temp_warp->lnd.nrows, temp_warp->lnd.ncols); 
    printf("\n\t\tulx = %8.1f;\tuly = %8.1f", temp_warp->lnd.ulx, temp_warp->lnd.uly);
    /* if extent is defined same as warp, then extent must be updated */
    if(strcasecmp(inPars->extentSource, "WARP")==0)
      updateOutCoor(&(temp_warp->lnd), &(out->lnd));
    /* save reprojected warp image for later combined resampling */
    if(createTempWarp(proj_warp, temp_warp)==FAILURE) {
      fprintf(stderr, "creating temporary warp image error\n");
      exit(1);
    }
  }

  /* if projection for DEM data is different from BASE image - Feng (10/08) */
  if(inPars->projection_dem_flag && inPars->option != REG_ONLY) {
    printf("\n\n- converting projection for DEM data to BASE image\n");
    /* do reprojection for warp image if necessary */
    if(reprojectDEM(sdem, base)==FAILURE) {
      fprintf(stderr, "Reproject DEM image error!\n");
      exit(1);
    }
    printf("\t\tnrows = %d;\tncols = %d", sdem->nrows, sdem->ncols); 
    printf("\n\t\tulx = %8.1f;\tuly = %8.1f", sdem->ulx, sdem->uly);
  }
  
  
  /*** CONVERT ACTUAL INPUT/OUTPUT SPACE TO WORKING SPACE ***/
  printf("\n\n- converting input space to working space");
  /* so base, warp and output images all have one same spatial resolution 
     for control point searching and image matching validation  */
  if(convertSpace(base, temp_warp, out, wbase, wwarp, wout, inPars)==FAILURE) {
    fprintf(stderr, "Convert input/output space to working space error!\n");
    exit(1);
  }

  /*** PREPARE DATA IN WORKING SPACE ***/
  if(inPars->option == ORTHO_ONLY || inPars->option == DO_BOTH) {
    printf("\n\n- computing parameters for orthorectification");

    /* use central point to compute Earth radius based on WGS84 ellipsoid and geoid */
    computeRadius(&(wbase->lnd));
    computeRadius(&(wwarp->lnd));
    computeRadius(&(wout->lnd));
    printf("\n\tEarth Radius = %9.1f m\n\tSatellite Altitude = %8.1f m", 
	   wbase->lnd.radius, wbase->lnd.altitude);
    

    /* extract the nadir track from the warp image */
    if(extractNadirPath(&(wwarp->lnd))==FAILURE) {
      fprintf(stderr, "Extract base image nadir view track from image error!\n");
      exit(1);
    }

    /* compute the warp image nadir track projected in output image */  
    computeNadirViewfromWarp(&(wout->lnd), wwarp);

    /* compute the warp image nadir track projected in base image */
    computeNadirViewfromWarp(&(wbase->lnd), wwarp);
    
    /* need nadir track function for computing distance from a given pixel to nadir point */
    printf("\n\tWarp Image Nadir Track in Output: y = %f * x + %f", wout->lnd.lpara[0], wout->lnd.lpara[1]);
    printf("\n\tWarp Image Nadir Track in Base  : y = %f * x + %f", wbase->lnd.lpara[0], wbase->lnd.lpara[1]);
  } /* endif ORTHO and BOTH */

  /* initialize transform polynomial function */
  initializePoly(&(wbase->lnd), &(wout->lnd));
  
  if(inPars->option == REG_ONLY || inPars->option == DO_BOTH) {

    /* force preliminary registration if MAX_SHIFT > 100  to avoid 
       1) too much computation
       2) too many false matches and thus ends up without enough control points (11/9/07) */
    if(CP_KEYS.MAX_SHIFT > 100)
      CP_KEYS.PRELIMINARY_REGISTRATION = 1;
    
    /*** DO PRELIMINARY REGISTRATION IF NECESSARY (TRY SHIFT ONLY) ***/
    if(CP_KEYS.PRELIMINARY_REGISTRATION == 1) {
      /* perform preliminary registration based on coarse resolution images */
      printf("\n\n- doing automated preliminary registration based on coarse resolution images\n");
      ret = preliminaryRegistration(wbase, wwarp);
      if(ret == FAILURE) {
	fprintf(stderr, "\nError in preliminary registration based on coarse resolution images!\n");
	exit(1);
      } 
 
      /* update output extent if it is defined same as WARP image */
      if(strcasecmp(inPars->extentSource, "WARP")==0)
	updateOutCoor(&(wwarp->lnd), &(out->lnd));
      
      /* reset maximum shift for fine resolution image to 3 times coarse pixel resolution */
      /* it should be large enough for precise registration */
      CP_KEYS.MAX_SHIFT = 3 * CP_KEYS.COARSE_SCALE;
    } 
  }

  if(inPars->option == DO_BOTH) {
    /*** DO PRECISE REGISTRATION (FIRST TRY SHIFT ONLY) ***/
    printf("\n\n- doing automated precise registration between fine BASE and WARP image\n");
    ret = preciseRegistration(wbase, wwarp, sdem);
    if(ret == FAILURE) {
      fprintf(stderr, "\nPrecise Registartion Error!\n");
      exit(1);
    }  
 
    /* update output extent if it is defined same as WARP image */
    if(strcasecmp(inPars->extentSource, "WARP")==0)
      updateOutCoor(&(wwarp->lnd), &(out->lnd));

    /*** TRY INITIAL ORTHORECTIFICATION IN WORKING SPACE ***/
    printf("\n\n- start initial orthorectification processing\n");
    ret = initialOrthoCheck(wbase, wwarp, wout, sdem);
  } /* endif BOTH */


  /*** DO FINAL ORTHORECTIFICATION AND PRECISE REGISTRATION IN REAL INPUT/OUTPUT SPACE ***/
  if(inPars->option == ORTHO_ONLY || inPars->option == DO_BOTH) {
    if(ret == SUCCESS) {
      printf("\n\n- doing final orthorectification processing\n");
      ret = finalOrtho(wbase, warp, temp_warp, wwarp, wout, out, sdem, inPars);
    }
    else {
      fprintf(stderr, "\nOrthorectification try failed!\n");
      exit(1);
    }
  }
  else if(inPars->option == REG_ONLY) {
    printf("\n\n- doing precise registration and resampling\n");
    ret = preciseRegistrationWithSampling(wbase, warp, temp_warp, wwarp, out, wout, inPars);
  }
  else if(inPars->option == VERIFY_ONLY) {
    printf("\n\n- searching tie points for image matching verification ...\n");
    ret = matchingVerification(wbase, wwarp, wout);    
  }

  if(inPars->option != VERIFY_ONLY) {
    /*** WRITE AND PRINT OUT ENVI HAEDER FILES FOR EACH BAND ***/
    printf("\n\n- writing ENVI header file for each band");
    ret = writeENVIheader(temp_warp, out);
  }

  /*** CLEAN UP MEMORY AND FILE SPACE ***/
  if(inPars->option == ORTHO_ONLY || inPars->option == DO_BOTH) {
    if(strcasecmp(sdem->fileType, "GEOTIFF")==0)
      XTIFFClose(sdem->fp_tiff);
    else
      fclose(sdem->fp);
    free(sdem->buf);
  }
 
  if(inPars->option == REG_ONLY || inPars->option == DO_BOTH || inPars->option == VERIFY_ONLY) {
    close_LANDSAT_file(&(base->lnd)); 
    close_LANDSAT_file(&(wbase->lnd)); 
    close_LANDSAT_file(&(warp->lnd));
    close_LANDSAT_file(&(temp_warp->lnd));
    close_LANDSAT_file(&(wwarp->lnd));
  }

#ifndef DEBUG
  /* remove temporary working file */ 
  remove(wbase->lnd.fileName);
  remove(wwarp->lnd.fileName);
  remove(wout->lnd.fileName);
  /* remove temporary rotation files */
  if(inPars->rotation_flag)
    remove(temp_warp->tempRotateFileName[0]);
  /* remove temporary reprojection files */
  if(inPars->projection_flag)
    remove(temp_warp->tempProjectFileName[0]);
  if(inPars->projection_dem_flag)
    remove(sdem->fileName);
#endif  

  if(inPars->option == DO_BOTH) {
    if(!wout->checking_passed  && wout->num_redo<= CP_KEYS.MAX_NUM_ITER-1) {
      for(i=0; i<wout->nbands; i++)
	fclose(wout->fp[i]);
    }
    else 
      close_LANDSAT_file(&(wout->lnd));
  }
   
  if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0)
    for(i=0; i<warp->nbands; i++) 
      XTIFFClose(warp->fp_tiff[i]);
  else 
    for(i=0; i<warp->nbands; i++) 
      fclose(warp->fp[i]);

  if(strcasecmp(temp_warp->lnd.fileType, "GEOTIFF")==0) 
    XTIFFClose(temp_warp->fp_tiff[0]);
  else
    fclose(temp_warp->fp[0]);

  if(strcasecmp(wwarp->lnd.fileType, "GEOTIFF")==0) 
    XTIFFClose(wwarp->fp_tiff[0]);
  else
    fclose(wwarp->fp[0]);
  
  free(base);
  free(warp);
  free(temp_warp);
  free(out);
  free(wbase);
  free(wwarp);
  free(wout);
  free(sdem);

  printf("\n\n- finished successfully!\n\n");
  return SUCCESS;

}


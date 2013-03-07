/**
 * ! Description
 *   other utility subroutines for Landsat orthorectification program 
 *
 * ! Credits
 *   Feng Gao (Feng.Gao@nasa.gov), Developer
 *   Jeff Masek (Jeffrey.G.Masek@nasa.gov)
 *   Robert Wolfe (Robert.E.Wolfe@nasa.gov)
 *
 * ! Revision 2.2.6  8/22/2011 
 * ! Revision 2.2.5  1/6/2010 
 * ! Revision 2.2.3  2/27/2009
 * ! Revision 2.2    5/22/2008
 */

#include "lndortho.h"

/* get input parameters from command line */
int getInParameter(IN_PARAMETERS *inPars, int argc, char *argv[])
{
  int   i;
  char  buffer[MAX_STRLEN] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *seperator = "= ,";
  FILE  *in;

  /* determine input option */
  if(strcasecmp(argv[1], "-r")==0)
    inPars->option = REG_ONLY;
  else if (strcasecmp(argv[1], "-o")==0)
    inPars->option = ORTHO_ONLY;
  else if(strcasecmp(argv[1], "-b")==0) 
    inPars->option = DO_BOTH;
  else if(strcasecmp(argv[1], "-v")==0) 
    inPars->option = VERIFY_ONLY;
  else {
    fprintf(stderr, "\nWrong option %s", argv[1]);
    usage(argv);
    return FAILURE;
  }
  
  /* read from input parameter file */
  if((in=fopen(argv[2],"r"))==NULL) {
    fprintf(stderr, "Can't open input %s\n", argv[2]);
    return FAILURE;
  }

  /* set default values */
  inPars->porder = 1;
  inPars->bkValue = 0;
  inPars->satValue = 255;
  inPars->satValue2 = 32768; 
  inPars->data_type = 1; 
  inPars->datum = 12;  /* base default in WGS84 */
  inPars->warp_angle = 0;
  inPars->warp_satellite_pangle = 0;
  inPars->warp_sys = -1;   
  inPars->warp_utm_zone = 999;
  inPars->warp_unit = 2;   /* default in meter */
  inPars->warp_datum = 12; /* warp default in WGS84 */
  for(i=0; i<15; i++) inPars->warp_param[i] = 0.0;
  inPars->warp_ulx_degree = -9999;
  inPars->warp_uly_degree = -9999;
  /* default data type in 1 byte */
  for(i=0; i<MAX_NBANDS; i++) inPars->warp_nbyte[i] = 1;
  inPars->dem_sys = -1;  
  inPars->dem_utm_zone = 999;
  inPars->dem_unit = 2;   /* default in meter */
  inPars->dem_datum = 12; /* output default in WGS84 */
  for(i=0; i<15; i++) inPars->dem_param[i] = 0.0;
  label = malloc(MAX_STRLEN);

  /* process line by line */
  while(fgets(buffer, MAX_STRLEN, in) != NULL) {

    /* get string token */
    tokenptr = strtok(buffer, seperator);
    label=tokenptr;
 
    /* skip comment line */
    if(strcmp(label,"#") == 0) continue;

    while(tokenptr != NULL) {
      
      tokenptr = strtok(NULL, seperator);
      
      /* get input from base file */
      if(strcasecmp(label, "BASE_FILE_TYPE") == 0)
	sscanf(tokenptr, "%s", inPars->base_fileType);
      else if(strcasecmp(label, "BASE_NSAMPLE") == 0)
	inPars->base_ncols = atoi(tokenptr);
      else if(strcasecmp(label, "BASE_NLINE") == 0)
	inPars->base_nrows = atoi(tokenptr);
      else if(strcasecmp(label, "BASE_PIXEL_SIZE") == 0)
	inPars->base_res = atof(tokenptr);
      else if(strcasecmp(label, "BASE_UPPER_LEFT_CORNER") == 0) {
	inPars->base_ulx = atof(tokenptr);
	tokenptr = strtok(NULL, seperator);	
	inPars->base_uly = atof(tokenptr);
      }
      else if(strcasecmp(label, "BASE_LANDSAT") == 0) 
	sscanf(tokenptr, "%s", inPars->base_fileName);
      else if(strcasecmp(label, "UTM_ZONE") == 0 || strcasecmp(label, "BASE_UTM_ZONE") == 0)
	inPars->utm_zone = atoi(tokenptr);
      else if(strcasecmp(label, "BASE_DATUM") == 0) {
	inPars->datum = atoi(tokenptr);
	/*printf("**** %s = %ld ****\n", label, inPars->datum);*/
      }
      else if(strcasecmp(label, "BASE_SATELLITE")==0)
	sscanf(tokenptr, "%s", inPars->base_satellite);
      
      /* get input from warp image */
      else if(strcasecmp(label, "WARP_FILE_TYPE") == 0 || strcasecmp(label, "WRAP_FILE_TYPE") == 0)
	sscanf(tokenptr, "%s", inPars->warp_fileType);
      else if(strcasecmp(label, "WARP_NSAMPLE") == 0 || strcasecmp(label, "WRAP_NSAMPLE") == 0)
	inPars->warp_ncols = atoi(tokenptr);
      else if(strcasecmp(label, "WARP_NLINE") == 0 || strcasecmp(label, "WRAP_NLINE") == 0)
	inPars->warp_nrows = atoi(tokenptr);
      else if(strcasecmp(label, "WARP_PIXEL_SIZE") == 0 || strcasecmp(label, "WRAP_PIXEL_SIZE") == 0)
	inPars->warp_res = atof(tokenptr);
      else if(strcasecmp(label, "WARP_ORIENTATION_ANGLE") == 0 || strcasecmp(label, "WRAP_ORIENTATION_ANGLE") == 0)
	inPars->warp_angle = atof(tokenptr);
      else if(strcasecmp(label, "WARP_UPPER_LEFT_CORNER") == 0 || strcasecmp(label, "WRAP_UPPER_LEFT_CORNER") == 0) {
	inPars->warp_ulx = atof(tokenptr);
	tokenptr = strtok(NULL, seperator);	
	inPars->warp_uly = atof(tokenptr);
      }
      else if(strcasecmp(label, "WARP_UPPER_LEFT_CORNER_DEGREE") == 0 || strcasecmp(label, "WRAP_UPPER_LEFT_CORNER_DEGREE") == 0) {
	inPars->warp_ulx_degree = atof(tokenptr);
	tokenptr = strtok(NULL, seperator);	
	inPars->warp_uly_degree = atof(tokenptr);
      }
      else if(strcasecmp(label, "WARP_NBANDS")==0 || strcasecmp(label, "WRAP_NBANDS")==0)
	inPars->warp_nbands = atoi(tokenptr);
      else if(strcasecmp(label, "WARP_LANDSAT_BAND") == 0 || strcasecmp(label, "WRAP_LANDSAT_BAND") == 0)
	for(i=0; i<inPars->warp_nbands; i++) {
	  sscanf(tokenptr, "%s", inPars->warp_fileName[i]);
	  tokenptr = strtok(NULL, seperator);
	}
      else if(strcasecmp(label, "WARP_BAND_DATA_TYPE") == 0 || strcasecmp(label, "WRAP_BAND_DATA_TYPE") == 0)
	for(i=0; i<inPars->warp_nbands; i++) {
	  sscanf(tokenptr, "%d", &(inPars->warp_nbyte[i]));
	  tokenptr = strtok(NULL, seperator);
	}
      else if(strcasecmp(label, "WARP_BASE_MATCH_BAND")==0 || strcasecmp(label, "WRAP_BASE_MATCH_BAND")==0)
	sscanf(tokenptr, "%s", inPars->warp_match_fileName);
      else if(strcasecmp(label, "WARP_SATELLITE")==0 || strcasecmp(label, "WRAP_SATELLITE")==0)
	sscanf(tokenptr, "%s", inPars->warp_satellite);
      else if(strcasecmp(label, "WARP_SATELLITE_POINTINGANGLE")==0 || strcasecmp(label, "WRAP_SATELLITE_POINTINGANGLE")==0)
	inPars->warp_satellite_pangle= atof(tokenptr);
      else if(strcasecmp(label, "WARP_PROJECTION_CODE")==0 || strcasecmp(label, "WRAP_PROJECTION_CODE")==0)
	inPars->warp_sys = atoi(tokenptr);
      else if(strcasecmp(label, "WARP_UTM_ZONE")==0 || strcasecmp(label, "WRAP_UTM_ZONE")==0)
	inPars->warp_utm_zone = atoi(tokenptr);     
      else if(strcasecmp(label, "WARP_PROJECTION_PARAM")==0 || strcasecmp(label, "WRAP_PROJECTION_PARAM")==0)
	for(i=0; i<15; i++) {
	  inPars->warp_param[i] = atof(tokenptr);
	  tokenptr = strtok(NULL, seperator);
	}
      else if(strcasecmp(label, "WARP_UNIT")==0 || strcasecmp(label, "WRAP_UNIT")==0)
	inPars->warp_unit = atoi(tokenptr);    
      else if(strcasecmp(label, "WARP_DATUM")==0 || strcasecmp(label, "WRAP_DATUM")==0)
	inPars->warp_datum = atoi(tokenptr);    

      else if(strcasecmp(label, "INPUT_DEM_FILE") == 0)
	sscanf(tokenptr, "%s", inPars->dem_fileName);
      else if(strcasecmp(label, "DEM_PROJECTION_CODE")==0)
	inPars->dem_sys = atoi(tokenptr);
      else if(strcasecmp(label, "DEM_UTM_ZONE")==0)
	inPars->dem_utm_zone = atoi(tokenptr);     
      else if(strcasecmp(label, "DEM_PROJECTION_PARAM")==0)
	for(i=0; i<15; i++) {
	  inPars->dem_param[i] = atof(tokenptr);
	  tokenptr = strtok(NULL, seperator);
	}
      else if(strcasecmp(label, "DEM_UNIT")==0)
	inPars->dem_unit = atoi(tokenptr);    
      else if(strcasecmp(label, "DEM_DATUM")==0)
	inPars->dem_datum = atoi(tokenptr);    

      else if(strcasecmp(label, "OUT_PIXEL_SIZE") == 0)
	inPars->out_res = atof(tokenptr);
      else if(strcasecmp(label, "RESAMPLE_METHOD") == 0)
	sscanf(tokenptr, "%s", inPars->resampleMethod);
      else if(strcasecmp(label, "BK_VALUE") == 0)
	inPars->bkValue = atoi(tokenptr);
      else if(strcasecmp(label, "OUT_EXTENT") == 0)
	sscanf(tokenptr, "%s", inPars->extentSource);
      else if(strcasecmp(label, "OUT_UPPER_LEFT_CORNER") == 0) {
	inPars->out_ulx = atof(tokenptr);
	tokenptr = strtok(NULL, seperator);	
	inPars->out_uly = atof(tokenptr);
      }
      else if(strcasecmp(label, "OUT_LOWER_RIGHT_CORNER") == 0) {
	inPars->out_lrx = atof(tokenptr);
	tokenptr = strtok(NULL, seperator);	
	inPars->out_lry = atof(tokenptr);
      }
      else if(strcasecmp(label, "OUT_LANDSAT_BAND") == 0) {
	for(i=0; i<inPars->warp_nbands; i++) {
	  sscanf(tokenptr, "%s", inPars->out_fileName[i]);
	  tokenptr = strtok(NULL, seperator);
	}
      }
      else if(strcasecmp(label, "OUT_BASE_MATCH_BAND")==0)
	sscanf(tokenptr, "%s", inPars->out_match_fileName);
      else if(strcasecmp(label, "OUT_BASE_POLY_ORDER")==0)
	inPars->porder = atoi(tokenptr);      
      else if(strcasecmp(label, "CP_PARAMETERS_FILE")==0)
	sscanf(tokenptr, "%s", inPars->cp_parameters_file);
      
      /* in case label (key words) is no the first word in a line */
      label = tokenptr;

    }  /* while token */
  } /* while line */
 
  fclose(in);

  /* use base image projection if warp image projection is not defined */
  if(inPars->warp_sys == -1) 
    inPars->warp_sys = UTM;    /* use UTM */
  if(inPars->warp_utm_zone == 999)
    inPars->warp_utm_zone = inPars->utm_zone;

   /* use base image projection if DEM image projection is not defined */
  if(inPars->dem_sys == -1)
    inPars->dem_sys = UTM;    /* use UTM */
  if(inPars->dem_utm_zone == 999)
    inPars->dem_utm_zone = inPars->utm_zone;

  inPars->rotation_flag = 0;
  if(fabs(inPars->warp_angle) > 0.00001)
    inPars->rotation_flag = 1;
  
  inPars->projection_flag = 0;
  if(!(inPars->warp_sys == UTM && inPars->warp_utm_zone == inPars->utm_zone && inPars->warp_datum == inPars->datum))
    inPars->projection_flag = 1;

  inPars->projection_dem_flag = 0;
  if(!(inPars->dem_sys == UTM && inPars->dem_utm_zone == inPars->utm_zone && inPars->dem_datum == inPars->datum))
    inPars->projection_dem_flag = 1;

  CP_KEYS.CHIP_SIZE = 11;          
  CP_KEYS.CP_SEED_WIN = 50;         
  CP_KEYS.MAX_SHIFT = 50;         
  CP_KEYS.MAX_SIZE = 2*(CP_KEYS.MAX_SHIFT+CP_KEYS.CHIP_SIZE);         
  CP_KEYS.MIN_STDEV = 1;
  CP_KEYS.MAX_NUM_HIGH_CORR = 3;   
  CP_KEYS.MIN_ACCEPTABLE_NCP = 10;  
  CP_KEYS.ACCEPTABLE_CORR = 0.60;
  CP_KEYS.MAX_AVE_ERROR = 0.5;   
  CP_KEYS.MAX_NUM_ITER = 1;
  /* preliminary registration settings for coarse resolution images */
  CP_KEYS.PRELIMINARY_REGISTRATION = 1;
  CP_KEYS.COARSE_SCALE = 10;
  CP_KEYS.COARSE_CP_SEED_WIN = 5;
  CP_KEYS.COARSE_MAX_SHIFT = 100; /* 10*100 pixels in fine resolution */
  CP_KEYS.MAX_RMSE = 0.75;
  CP_KEYS.DIFF_RMSE = 0.01;

  if((in=fopen(inPars->cp_parameters_file, "r"))==NULL) {
    printf("\nControl point parameter file %s doesn't exist.", inPars->cp_parameters_file);
    printf("\nDefault values will be used for this run\n");
  }
  else {
    /* process line by line */
    while(fgets(buffer, MAX_STRLEN, in) != NULL) {
      
      /* get string token */
      tokenptr = strtok(buffer, seperator);
      label=tokenptr;
 
      /* skip comment line */
      if(strcmp(label,"#") == 0) continue;
      while(tokenptr != NULL) {
	tokenptr = strtok(NULL, seperator);
	if(strcasecmp(label, "CHIP_SIZE") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.CHIP_SIZE));
	else if(strcasecmp(label, "CP_SEED_WIN") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.CP_SEED_WIN));
	else if(strcasecmp(label, "MAX_SHIFT") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.MAX_SHIFT));
	else if(strcasecmp(label, "MAX_NUM_HIGH_CORR") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.MAX_NUM_HIGH_CORR));
 	else if(strcasecmp(label, "MIN_ACCEPTABLE_NCP") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.MIN_ACCEPTABLE_NCP));
 	else if(strcasecmp(label, "MAX_NUM_ITER") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.MAX_NUM_ITER));
	else if(strcasecmp(label, "ACCEPTABLE_CORR") == 0)
	  sscanf(tokenptr, "%f", &(CP_KEYS.ACCEPTABLE_CORR));
	else if(strcasecmp(label, "MAX_AVE_ERROR") == 0)
	  sscanf(tokenptr, "%f", &(CP_KEYS.MAX_AVE_ERROR));
 	else if(strcasecmp(label, "PRELIMINARY_REGISTRATION") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.PRELIMINARY_REGISTRATION));
 	else if(strcasecmp(label, "COARSE_SCALE") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.COARSE_SCALE));
 	else if(strcasecmp(label, "COARSE_MAX_SHIFT") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.COARSE_MAX_SHIFT));
 	else if(strcasecmp(label, "COARSE_CP_SEED_WIN") == 0)
	  sscanf(tokenptr, "%d", &(CP_KEYS.COARSE_CP_SEED_WIN));
	else if(strcasecmp(label, "MAX_ACCEPTABLE_RMSE") == 0)
	  sscanf(tokenptr, "%f", &(CP_KEYS.MAX_RMSE));
	
	/* maximum array size for control point searching */
	CP_KEYS.MAX_SIZE  = 2*(CP_KEYS.MAX_SHIFT+CP_KEYS.CHIP_SIZE);         
	/* in case label (key words) is no the first word in a line */
	label = tokenptr;
      }  /* while token */
    } /* while line */
    fclose(in);
  }


#ifdef DEBUG
  printf("PRELIMINARY_REGISTRATION = %d\n", CP_KEYS.PRELIMINARY_REGISTRATION);
  printf("COARSE_SCALE = %d\n", CP_KEYS.COARSE_SCALE);
  printf("COARSE_MAX_SHIFT = %d\n", CP_KEYS.COARSE_MAX_SHIFT);
  printf("COARSE_CP_SEED_WIN = %d\n", CP_KEYS.COARSE_CP_SEED_WIN);
  printf("CHIP_SIZE  = %d\n", CP_KEYS.CHIP_SIZE);
  printf("CP_SEED_WIN = %d\n", CP_KEYS.CP_SEED_WIN);
  printf("MAX_SHIFT = %d\n", CP_KEYS.MAX_SHIFT);
  printf("MAX_NUM_HIGH_CORR = %d\n", CP_KEYS.MAX_NUM_HIGH_CORR); 
  printf("ACCEPTABLE_CORR = %4.2f\n", CP_KEYS.ACCEPTABLE_CORR); 
  printf("MIN_ACCEPTABLE_NCP = %d\n", CP_KEYS.MIN_ACCEPTABLE_NCP); 
  printf("MAX_AVE_ERROR = %4.2f\n", CP_KEYS.MAX_AVE_ERROR);
  printf("MAX_NUM_ITER = %d\n", CP_KEYS.MAX_NUM_ITER);
#endif

  return SUCCESS;
}


/* get DEM metadata from SRTM data (in GeoTiff format) */ 
int getSRTMMetaInfo(DEM *sdem, IN_PARAMETERS *inPars) {

  int i;
  uint16 count, coor_sys;
  double *tiePoint, *pixelScale;
  GTIF *gtif;

  strcpy(sdem->fileType, "GEOTIFF");

  if((sdem->fp_tiff = XTIFFOpen(sdem->fileName, "r"))==NULL) {
    fprintf(stderr, "Can't open file %s\n", sdem->fileName);
    return FAILURE;
  }

  /* get SRTM metadata from tiff file */
  if(TIFFGetField(sdem->fp_tiff, TIFFTAG_IMAGEWIDTH, &(sdem->ncols))==0) {
    printf("Retrieve SRTM file %s error\n", sdem->fileName);
    return FAILURE;
  }
  if(TIFFGetField(sdem->fp_tiff, TIFFTAG_IMAGELENGTH, &(sdem->nrows))==0) {
    printf("Retrieve SRTM file %s error\n", sdem->fileName);
    return FAILURE;
  }
  count=6;
  if(TIFFGetField(sdem->fp_tiff, TIFFTAG_GEOTIEPOINTS, &count, &tiePoint)==0) {
    printf("Retrieve SRTM file %s error\n", sdem->fileName);
    return FAILURE;
  }
  count=3;
  if(TIFFGetField(sdem->fp_tiff, TIFFTAG_GEOPIXELSCALE, &count, &pixelScale)==0) {
    printf("Retrieve SRTM file %s error\n", sdem->fileName);
    return FAILURE;
  }

  sdem->res = pixelScale[0];

  /* GeoKey 1025 (GTRasterTypeGeoKey) dictates whether the reference
     coordinate is the UL (*RasterPixelIsArea*, code 1) or center
     (*RasterPixelIsPoint*, code 2) of the UL pixel. If this key is missing,
     the default (as defined by the specification) is to be
     *RasterPixelIsArea*, which is the UL of the UL pixel. */
  gtif = GTIFNew(sdem->fp_tiff);
  if (GTIFKeyGet(gtif, GTRasterTypeGeoKey, &coor_sys, 0, 1) != 1) {
    printf("Coordinate system is not defined in %s\n", sdem->fileName);
    printf("assume used UL of the UL pixel\n");
  }
  if (coor_sys == RasterPixelIsPoint){
    sdem->ulx = tiePoint[3] - 0.5 * sdem->res;
    sdem->uly = tiePoint[4] + 0.5 * sdem->res;
  }
  else {  /* default use RasterPixelIsArea */
    sdem->ulx = tiePoint[3];
    sdem->uly = tiePoint[4];
  }
  GTIFFree(gtif);

  if(inPars->dem_sys != -1) {
    /* if defined new projection system */
    sdem->proj.sys = inPars->dem_sys;
    sdem->proj.utm_zone = inPars->dem_utm_zone;
    sdem->proj.unit = inPars->dem_unit;
    sdem->proj.datum = inPars->dem_datum;
    for(i=0; i<15; i++)
      sdem->proj.param[i] = inPars->dem_param[i];
  }
  else {
    sdem->proj.sys = UTM;
    if(inPars->dem_utm_zone == 999)
      sdem->proj.utm_zone = inPars->utm_zone;
    else
      sdem->proj.utm_zone = inPars->dem_utm_zone;
  }

  alloc_1dim_contig((void **) (&sdem->buf), sdem->ncols, sizeof(int16));

  return SUCCESS;
}



/* get inputs from input parameters */
int getLandsatMetaInfo(IN_PARAMETERS *inPars, BASE_LANDSAT *base, WARP_LANDSAT *warp, OUT_LANDSAT *out) 
{
  int i, iband;
  int nrows, ncols, min_rc;
  double *tiePoint, *pixelScale;
  uint16 count, coor_sys;
  GTIF *gtif;

  /* variables for GCTPC */
  double incoor[2];
  double outcoor[2];
  long   insys = 0;
  long   inzone;
  double inparm[15];
  long   inunit;
  long   indatum;
  long   ipr = 0;
  long   jpr = 999;
  long   outsys;
  long   outzone;
  double outparm[15];
  long   outunit;
  long   outdatum;
  long   flg;
  char   file27[MAX_STRLEN],file83[MAX_STRLEN],efile[MAX_STRLEN],file1[MAX_STRLEN];

  /*** prepare parameters for BASE image ***/
  strcpy(base->lnd.fileName, inPars->base_fileName);
  strcpy(base->lnd.fileType, inPars->base_fileType);
  base->lnd.nrows = inPars->base_nrows;
  base->lnd.ncols = inPars->base_ncols;
  base->lnd.ulx = inPars->base_ulx;
  base->lnd.uly = inPars->base_uly;
  base->lnd.res = inPars->base_res;
  base->lnd.fillValue = inPars->bkValue;
  base->lnd.satValue = inPars->satValue;
  base->lnd.nbyte = inPars->data_type;
  base->utm_zone = inPars->utm_zone;
  base->datum = inPars->datum;

  /* open file for read */
  if(strcasecmp(base->lnd.fileType, "GEOTIFF")==0) {
    if((base->lnd.fp_tiff = XTIFFOpen(base->lnd.fileName, "r"))==NULL) {
      fprintf(stderr, "Can't open base GEOTIFF file %s\n", base->lnd.fileName);
      return FAILURE;
    } 
    /* get BASE Landsat metadata from tiff file */
    if(TIFFGetField(base->lnd.fp_tiff, TIFFTAG_IMAGEWIDTH, &ncols)==0) {
      printf("Retrieve BASE Landsat file %s error\n", base->lnd.fileName);
      return FAILURE;
    }
    if(ncols != base->lnd.ncols) {
      /* printf("\n\tWARNING: Number of samples in %s does not match your input.", base->lnd.fileName);
	 printf("\n\t         The setting in GEOTIFF file will be used instead.");*/
      base->lnd.ncols = ncols;
    }
    if(TIFFGetField(base->lnd.fp_tiff, TIFFTAG_IMAGELENGTH, &nrows)==0) {
      printf("Retrieve BASE Landsat file %s error\n", base->lnd.fileName);
      return FAILURE;
    }
    if(nrows != base->lnd.nrows) {
      /* printf("\n\tWARNING: Number of rows in %s does not match your input.", base->lnd.fileName);
	 printf("\n\t         The setting in GEOTIFF file will be used instead.");*/
      base->lnd.nrows = nrows;
    }
    count=6;
    if(TIFFGetField(base->lnd.fp_tiff, TIFFTAG_GEOTIEPOINTS, &count, &tiePoint)==0) {
      printf("Retrieve BASE Landsat file %s error\n", base->lnd.fileName);
      return FAILURE;
    }
    if(fabs(tiePoint[3]-base->lnd.ulx)>1.0||fabs(tiePoint[4]-base->lnd.uly)>1.0) {
      /* printf("\n\tWARNING: The upper left coordinate in %s does not match your input.", base->lnd.fileName);
	 printf("\n\t         The setting in GEOTIFF file will be used instead.");*/
      base->lnd.ulx = tiePoint[3];
      base->lnd.uly = tiePoint[4];
    }
    count=3;
    if(TIFFGetField(base->lnd.fp_tiff, TIFFTAG_GEOPIXELSCALE, &count, &pixelScale)==0) {
      printf("Retrieve BASE Landsat file %s error\n", base->lnd.fileName);
      return FAILURE;
    }
    if(fabs(pixelScale[0]-base->lnd.res)>0.01) {
      /* printf("\n\tWARNING: The pixel resolution in %s does not match your input.", base->lnd.fileName);
	 printf("\n\t         The setting in GEOTIFF file will be used instead.");*/
      base->lnd.res = pixelScale[0];
    }    

    gtif = GTIFNew(base->lnd.fp_tiff);
    if (GTIFKeyGet(gtif, GTRasterTypeGeoKey, &coor_sys, 0, 1) != 1) {
      printf("Coordinate system is not defined in %s\n", base->lnd.fileName);
      printf("assume used UL of the UL pixel\n");
    }
    if (coor_sys == RasterPixelIsPoint){
      base->lnd.ulx = tiePoint[3] - 0.5 * base->lnd.res;
      base->lnd.uly = tiePoint[4] + 0.5 * base->lnd.res;
    }
    GTIFFree(gtif);

  }
  else {
    if((base->lnd.fp=fopen(base->lnd.fileName, "rb"))==NULL) {
      /* it's ok without base image for ortho only option */
      if(inPars->option != ORTHO_ONLY) 
	{
	  fprintf(stderr, "Can't open base binary file %s\n", base->lnd.fileName);
	  return FAILURE;
	}
    }
  }
  
  /*** prepare parameters for WARP images ***/
  warp->nbands = inPars->warp_nbands;
  strcpy(warp->lnd.fileType, inPars->warp_fileType);
  warp->lnd.nrows = inPars->warp_nrows;
  warp->lnd.ncols = inPars->warp_ncols;
  for(i=0; i< warp->nbands; i++) {
    warp->nbyte[i] = inPars->warp_nbyte[i];
    /* use maximum value as saturation value */
    if(warp->nbyte[i] == 1)
      warp->satValue[i] = inPars->satValue;
    if(warp->nbyte[i] == 2)
      warp->satValue[i] = inPars->satValue2; 
    out->nbyte[i] = inPars->warp_nbyte[i];
  }
  
  /* convert degree to meters if necessary */
  if(inPars->warp_ulx_degree > -180 && inPars->warp_ulx_degree < 180 &&
     inPars->warp_uly_degree > -90 && inPars->warp_uly_degree < 90) {
    insys = GEO; 
    inzone = 62;
    for(i=0; i<15; i++) inparm[i] = 0;
    inunit = 4;
    indatum = 0;

    outsys = inPars->warp_sys;
    outzone = inPars->warp_utm_zone;
    for(i=0; i<15; i++) outparm[i] = inPars->warp_param[i];
    outunit = inPars->warp_unit;
    outdatum = inPars->warp_datum;

    /* convert upperleft coordinate from degree to meters in WARP image projection */
    incoor[0] = inPars->warp_ulx_degree;
    incoor[1] = inPars->warp_uly_degree;
    gctp(incoor,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,file1,
	 outcoor,&outsys,&outzone,outparm,&outunit,&outdatum,file27,file83,&flg);
    inPars->warp_ulx = outcoor[0];
    inPars->warp_uly = outcoor[1];
  }

  /* if warp UL is not trustable, use base instead */
  if(fabs(inPars->warp_ulx+1)<0.00001) {
    warp->lnd.ulx = base->lnd.ulx;
    warp->lnd.uly = base->lnd.uly;
  }
  else {
    warp->lnd.ulx = inPars->warp_ulx;
    warp->lnd.uly = inPars->warp_uly;
  }
  warp->lnd.res = inPars->warp_res;
  if(inPars->warp_sys != -1) {
    /* if defined new projection system */
    warp->proj.sys = inPars->warp_sys;
    warp->proj.utm_zone = inPars->warp_utm_zone;
    warp->proj.unit = inPars->warp_unit;
    warp->proj.datum = inPars->warp_datum;
    for(i=0; i<15; i++)
      warp->proj.param[i] = inPars->warp_param[i];
  }
  else {
    warp->proj.sys = UTM;
    if(inPars->warp_utm_zone == 999)
      warp->proj.utm_zone = inPars->utm_zone;
    else
      warp->proj.utm_zone = inPars->warp_utm_zone;
  }

  warp->lnd.fillValue = inPars->bkValue;
  warp->lnd.satValue = inPars->satValue;
  warp->lnd.nbyte = inPars->data_type;
  warp->orientationAngle = inPars->warp_angle; 
  warp->pointingAngle = inPars->warp_satellite_pangle; 

  strcpy(warp->resampleMethod, inPars->resampleMethod);
  
  /* select the right satellite major axis */
  getSatMajorAxis(inPars->warp_satellite, &(warp->lnd));
  getSatMajorAxis(inPars->base_satellite, &(base->lnd));
  
  /* open all input warp files */
  for(iband=0; iband<warp->nbands; iband++) {
    strcpy(warp->fileName[iband], inPars->warp_fileName[iband]);
    if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0) {
      if((warp->fp_tiff[iband] = XTIFFOpen(warp->fileName[iband], "r"))==NULL) {
	fprintf(stderr, "Can't open warp GEOTIFF file %s\n", warp->fileName[iband]);
	return FAILURE;
      } 
    }
    else {
      if((warp->fp[iband]=fopen(warp->fileName[iband], "rb"))==NULL) {
	fprintf(stderr, "Can't open warp binary file %s\n", warp->fileName[iband]);
	return FAILURE;
      }
    }
  }  
  /* open warp matching file */
  strcpy(warp->lnd.fileName, inPars->warp_match_fileName);
  if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0) {
    if((warp->lnd.fp_tiff = XTIFFOpen(warp->lnd.fileName, "r"))==NULL) {
      fprintf(stderr, "Can't open matching GEOTIFF file %s\n", warp->lnd.fileName);
      return FAILURE;
    } 
    /* get WARP Landsat metadata from tiff file */
    if(TIFFGetField(warp->lnd.fp_tiff, TIFFTAG_IMAGEWIDTH, &ncols)==0) {
      fprintf(stderr, "Retrieve WARP Landsat file %s error\n", warp->lnd.fileName);
      return FAILURE;
    }
    if(ncols != warp->lnd.ncols) {
      /*printf("\n\tWARNING: Number of samples in %s does not match your input.", warp->fileName[iband]);
	printf("\n\t         The setting in GEOTIFF file will be used instead.");*/
      warp->lnd.ncols = ncols;
    }
    if(TIFFGetField(warp->lnd.fp_tiff, TIFFTAG_IMAGELENGTH, &nrows)==0) {
      printf("Retrieve WARP Landsat file %s error\n", warp->lnd.fileName);
      return FAILURE;
    }
    if(nrows != warp->lnd.nrows) {
      /*printf("\n\tWARNING: Number of rows in %s does not match your input.", warp->fileName[iband]);
	printf("\n\t         The setting in GEOTIFF file will be used instead.");*/
      warp->lnd.nrows = nrows;
    }
    count=6;
    if(TIFFGetField(warp->lnd.fp_tiff, TIFFTAG_GEOTIEPOINTS, &count, &tiePoint)==0) {
      printf("Retrieve WARP Landsat file %s error\n", warp->lnd.fileName);
      return FAILURE;
    }
    if(fabs(tiePoint[3]-warp->lnd.ulx)>1.0||fabs(tiePoint[4]-warp->lnd.uly)>1.0) {
      /*printf("\n\tWARNING: The upper left coordinate in %s does not match your input.", warp->fileName[iband]);
	printf("\n\t         The setting in GEOTIFF file will be used instead.");*/
      warp->lnd.ulx = tiePoint[3];
      warp->lnd.uly = tiePoint[4];
    }
    count=3;
    if(TIFFGetField(warp->lnd.fp_tiff, TIFFTAG_GEOPIXELSCALE, &count, &pixelScale)==0) {
      printf("Retrieve WARP Landsat file %s error\n", warp->lnd.fileName);
      return FAILURE;
    }
    if(fabs(pixelScale[0]-warp->lnd.res)>0.01) {
      /*printf("\n\tWARNING: The pixel resolution in %s does not match your input.", warp->fileName[iband]);
	printf("\n\t         The setting in GEOTIFF file will be used instead.");*/
      warp->lnd.res = pixelScale[0];
    }    

    gtif = GTIFNew(warp->lnd.fp_tiff);
    if (GTIFKeyGet(gtif, GTRasterTypeGeoKey, &coor_sys, 0, 1) != 1) {
      printf("Coordinate system is not defined in %s\n", warp->lnd.fileName);
      printf("assume used UL of the UL pixel\n");
    }
    if (coor_sys == RasterPixelIsPoint){
      warp->lnd.ulx = tiePoint[3] - 0.5 * warp->lnd.res;
      warp->lnd.uly = tiePoint[4] + 0.5 * warp->lnd.res;
    }
    GTIFFree(gtif);
  }
  else {
    if((warp->lnd.fp = fopen(warp->lnd.fileName, "rb"))==NULL) {
      /* it's ok without base image for ortho only option */
      if(inPars->option != ORTHO_ONLY) {
	fprintf(stderr, "Can't open matching binary file %s\n", warp->lnd.fileName);
	return FAILURE;
      }
    }
  }
  
  /*** Prepare parameters for orthorectified output file ***/
  if(inPars->option!=VERIFY_ONLY) {
    out->nbands = warp->nbands;
    out->utm_zone = inPars->utm_zone;
    out->lnd.res = inPars->out_res;
    out->lnd.porder = inPars->porder;
    out->lnd.fillValue = inPars->bkValue;
    out->lnd.satValue = inPars->satValue;
    out->lnd.nbyte = inPars->data_type;
    for(iband=0; iband<out->nbands; iband++) {
      strcpy(out->fileName[iband], inPars->out_fileName[iband]);
      if((out->fp[iband]=fopen(out->fileName[iband],"w+b"))==NULL) {
	fprintf(stderr, "Can't open output file %s\n", out->fileName[iband]);
	return FAILURE;
      }
    }

    if(strcasecmp(inPars->extentSource, "BASE")==0) {
      out->lnd.ulx = base->lnd.ulx;
      out->lnd.uly = base->lnd.uly;
      out->lnd.nrows = (int)(base->lnd.nrows * base->lnd.res / out->lnd.res + 0.5);
      out->lnd.ncols = (int)(base->lnd.ncols * base->lnd.res / out->lnd.res + 0.5);
    }
    else if(strcasecmp(inPars->extentSource, "WARP")==0 || strcasecmp(inPars->extentSource, "WRAP")==0) {
      out->lnd.ulx = warp->lnd.ulx;
      out->lnd.uly = warp->lnd.uly;
      out->lnd.nrows = (int)(warp->lnd.nrows * warp->lnd.res / out->lnd.res + 0.5);
      out->lnd.ncols = (int)(warp->lnd.ncols * warp->lnd.res / out->lnd.res + 0.5);
    }
    else if(strcasecmp(inPars->extentSource, "DEF")==0) {
      out->lnd.ulx = inPars->out_ulx;
      out->lnd.uly = inPars->out_uly;
      out->lnd.nrows = (int)((inPars->out_uly - inPars->out_lry) / out->lnd.res + 0.5);
      out->lnd.ncols = (int)((inPars->out_lrx - inPars->out_ulx) / out->lnd.res + 0.5);
    } 
    else {
      fprintf(stderr, "Wrong OUT_EXTENT option\n");
      return FAILURE;
    }
    /* save output matching file */
    strcpy(out->lnd.fileName, inPars->out_match_fileName);
    strcpy(out->lnd.fileType, "BINARY");

    /* copy satellite semi-major axis from warp to out */
    out->lnd.satMajorAxis = warp->lnd.satMajorAxis;
  }

  /* compute lat/lon for central pixels for BASE, WARP and OUT */
  outsys = GEO; 
  outzone = 62;
  for(i=0; i<15; i++) outparm[i] = 0;
  outunit = 4;
  outdatum = 0;

  /* for BASE */
  insys = UTM;
  inzone = inPars->utm_zone;
  for(i=0; i<15; i++) inparm[i] = 0;
  inunit = 2;
  indatum = 0;

  incoor[0] = base->lnd.ulx + base->lnd.ncols/2.0*base->lnd.res;
  incoor[1] = base->lnd.uly - base->lnd.nrows/2.0*base->lnd.res;
  gctp(incoor,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,file1,
       outcoor,&outsys,&outzone,outparm,&outunit,&outdatum,file27,file83,&flg);
  base->lnd.centerLonLat[0] = outcoor[0];
  base->lnd.centerLonLat[1] = outcoor[1];

  /*for OUT (projection should be same as BASE) */
  incoor[0] = out->lnd.ulx + out->lnd.ncols/2.0*out->lnd.res;
  incoor[1] = out->lnd.uly - out->lnd.nrows/2.0*out->lnd.res;
  gctp(incoor,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,file1,
       outcoor,&outsys,&outzone,outparm,&outunit,&outdatum,file27,file83,&flg);
  out->lnd.centerLonLat[0] = outcoor[0];
  out->lnd.centerLonLat[1] = outcoor[1];

  /* for WARP */
  insys = inPars->warp_sys;
  inzone = inPars->warp_utm_zone;
  for(i=0; i<15; i++) inparm[i] = inPars->warp_param[i];
  inunit = inPars->warp_unit;
  indatum = inPars->warp_datum;
  incoor[0] = warp->lnd.ulx + warp->lnd.ncols/2.0*warp->lnd.res;
  incoor[1] = warp->lnd.uly - warp->lnd.nrows/2.0*warp->lnd.res;
  gctp(incoor,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,file1,
       outcoor,&outsys,&outzone,outparm,&outunit,&outdatum,file27,file83,&flg);
  warp->lnd.centerLonLat[0] = outcoor[0];
  warp->lnd.centerLonLat[1] = outcoor[1];

  /* adjust number of tie points (seed distance) according to base and warp image size */
  /* at least 100 tie points for the smallest edge */
  min_rc = warp->lnd.nrows;
  if(warp->lnd.ncols < min_rc) min_rc = warp->lnd.ncols;
  if(base->lnd.nrows < min_rc) min_rc = base->lnd.nrows;
  if(base->lnd.ncols < min_rc) min_rc = base->lnd.ncols;
  CP_KEYS.CP_SEED_WIN = (int) (min_rc / 100.0);

  return SUCCESS;
}


/* convert coordinates in different projections 
 * flg = FOWARD_PROJ or BACKWARD_PROJ 
 */
int convertCoor(BASE_LANDSAT *base, double *bx, double *by, PROJECTION *proj, double *wx, double *wy, int ptype)
{
  int i, ret;
  /* variables for GCTPC */
  double incoor[2];
  double outcoor[2];
  long   insys = 0;
  long   inzone;
  double inparm[15];
  long   inunit;
  long   indatum;
  long   ipr = 0;
  long   jpr = 999;
  long   outsys;
  long   outzone;
  double outparm[15];
  long   outunit;
  long   outdatum;
  long   flg;
  char   file27[MAX_STRLEN],file83[MAX_STRLEN],efile[MAX_STRLEN],file1[MAX_STRLEN];

  /* assign input projection parameters */
  insys = proj->sys;
  inzone = proj->utm_zone;
  for(i=0; i<15; i++) inparm[i] = proj->param[i];
  inunit = proj->unit;
  indatum = proj->datum;

  /* assign output (base image) projection parameters */
  outsys = UTM;    /* in UTM: 1 */
  outzone = base->utm_zone;
  for(i=0; i<15; i++) outparm[i] = 0;
  outunit = 2;   /* in meter: 2 */
  outdatum = base->datum; 

  /* set values to change the datum and radius */
  indatum=0;
  if ((inparm[0] != 0) && (inparm[0] != 6370997))
    if ((insys != 1) && (insys != 2))
      indatum = -1;
  
  outdatum=0;
  if ((outparm[0] != 0) && (outparm[0] != 6370997))
    if ((outsys != 1) && (outsys != 2))
      outdatum = -1;
  
  if(ptype == FORWARD_PROJ) {
    incoor[0] = *wx;
    incoor[1] = *wy;

    gctp(incoor,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,file1,
	 outcoor,&outsys,&outzone,outparm,&outunit,&outdatum,
	 file27,file83,&flg);

    if(flg==0) {
      *bx = outcoor[0];
      *by = outcoor[1];
      ret = SUCCESS;
    }
    else
      ret = FAILURE;
  }
  else if(ptype == BACKWARD_PROJ) {
    incoor[0] = *bx;
    incoor[1] = *by;

    gctp(incoor,&outsys,&outzone,outparm,&outunit,&outdatum,&ipr,efile,&jpr,file1,
	 outcoor,&insys,&inzone,inparm,&inunit,&indatum,
	 file27,file83,&flg);

    if(flg==0) {
      *wx = outcoor[0];
      *wy = outcoor[1];
      ret = SUCCESS;
    }
    else
      ret = FAILURE;
  }
  else ret = FAILURE;
  return ret;
}


/* find the right satellite major axis */
void getSatMajorAxis(char *satellite, LANDSAT *lnd)
{
  if(strcasecmp(satellite, "Landsat1")==0)
    lnd->satMajorAxis = Landsat1;
  else if(strcasecmp(satellite, "Landsat2")==0)
    lnd->satMajorAxis = Landsat2; 
  else if(strcasecmp(satellite, "Landsat3")==0)
    lnd->satMajorAxis = Landsat3; 
  else if(strcasecmp(satellite, "Landsat4")==0)
    lnd->satMajorAxis = Landsat4; 
  else if(strcasecmp(satellite, "Landsat5")==0)
    lnd->satMajorAxis = Landsat5; 
  else if(strcasecmp(satellite, "Landsat7")==0)
    lnd->satMajorAxis = Landsat7; 
  else if(strcasecmp(satellite, "CBERS1")==0)
    lnd->satMajorAxis = CBERS1; 
  else if(strcasecmp(satellite, "CBERS2")==0)
    lnd->satMajorAxis = CBERS2; 
  else if(strcasecmp(satellite, "TERRA")==0)
    lnd->satMajorAxis = TERRA;
  else if(strcasecmp(satellite, "AWIFS")==0)
    lnd->satMajorAxis = AWIFS;
  else if(strcasecmp(satellite, "HJ1")==0)
    lnd->satMajorAxis = HJ1;
  else if(strcasecmp(satellite, "MySensor")==0)
    lnd->satMajorAxis = MySensor;
  else {
    printf("\nDidn't find any matching satellite for %s, use Landsat7 instead %d\n", 
	   satellite, Landsat7);
    lnd->satMajorAxis = Landsat7;
  }
}



/* make coarse resolution images */
int makeCoarse(LANDSAT *fine, LANDSAT *coarse)
{
  coarse->nrows = fine->nrows / CP_KEYS.COARSE_SCALE;
  coarse->ncols = fine->ncols / CP_KEYS.COARSE_SCALE;
  strcpy(coarse->fileType, "BINARY");
  coarse->ulx = fine->ulx;
  coarse->uly = fine->uly;
  coarse->res = fine->res * CP_KEYS.COARSE_SCALE;
  coarse->fillValue = fine->fillValue;
  coarse->satValue = fine->satValue;
  coarse->nbyte = fine->nbyte; 

  aggregateImage(fine, coarse);

  if((coarse->fp=fopen(coarse->fileName, "rb"))==NULL) {
    fprintf(stderr, "Error in opening temporary coarse resolution image for reading\n");
    return FAILURE;
  }  

  return SUCCESS;
}


/**
   CONVERT ACTUAL INPUT/OUTPUT SPACE TO WORKING SPACE 
   so base, warp and output images all have one same spatial resolution 
   (=max(base, warp)) 
   for control point searching and image matching validation 
**/
int convertSpace(BASE_LANDSAT *base, WARP_LANDSAT *warp, OUT_LANDSAT *out, 
		 BASE_LANDSAT *wbase, WARP_LANDSAT *wwarp, OUT_LANDSAT *wout, 
		 IN_PARAMETERS *inPars)
{
  int  choice, nrows, old_scale;
  char command[MAX_STRLEN];
  enum RESOLUTION {NO_CHANGE=0, USE_BASE=1, USE_WARP=2};

  /* if resolution are close enough (<10%), then we don't need extra process */
  if(fabs(warp->lnd.res-base->lnd.res)/base->lnd.res < 0.10) {
    choice = NO_CHANGE;
    nrows = base->lnd.nrows;
  }
  else if(warp->lnd.res > base->lnd.res) {
    choice = USE_WARP;
    /* adjust scale for coarse resolution output (always based on the finer res. input) */
    CP_KEYS.COARSE_SCALE = (int) (CP_KEYS.COARSE_SCALE*base->lnd.res/warp->lnd.res); 
    nrows = warp->lnd.nrows;
  }
  else {
    choice = USE_BASE;
    /* adjust scale for coarse resolution output (always based on the finer res. input) */
    CP_KEYS.COARSE_SCALE = (int) (CP_KEYS.COARSE_SCALE*warp->lnd.res/base->lnd.res); 
    nrows = base->lnd.nrows;
  }

  /* adjust reduction scale for pyramid registration if coarse image size is too small */
  if(nrows/MIN_COARSE_NROWS < CP_KEYS.COARSE_SCALE && 
     nrows/MIN_COARSE_NROWS > 1) {
    old_scale = CP_KEYS.COARSE_SCALE;
    CP_KEYS.COARSE_SCALE = nrows/MIN_COARSE_NROWS;
    CP_KEYS.COARSE_MAX_SHIFT = CP_KEYS.COARSE_MAX_SHIFT * old_scale / CP_KEYS.COARSE_SCALE;
    printf("\n\tCoarse image is too small, adjusted COARSE_SCALE = %d", CP_KEYS.COARSE_SCALE); 
    printf("\n\tAdjusted COARSE_MAX_SHIFT = %d", CP_KEYS.COARSE_MAX_SHIFT);
  }

  /* create working file for base, warp and output */
  /* for base file */
  copyLandsat(&(wbase->lnd), &(base->lnd));
  strcpy(wbase->lnd.fileName, "temp_working_base.dat");
  wbase->utm_zone = base->utm_zone;
  wbase->datum = base->datum;

  /* for warp file (need only one band for working file) */
  copyLandsat(&(wwarp->lnd), &(warp->lnd));
  strcpy(wwarp->lnd.fileName, "temp_working_warp.dat");
  wwarp->nbands = 1;
  wwarp->nbyte[0] = 1;
  strcpy(wwarp->fileName[0], wwarp->lnd.fileName);
  /* always use NN in working space operation since output resolution same as input */
  strcpy(wwarp->resampleMethod, "NN");  
  wwarp->pointingAngle = warp->pointingAngle;

  /* for out file (need only one band for output validation) */
  copyLandsat(&(wout->lnd), &(out->lnd));
  wout->nbands = 1;
  wout->nbyte[0] = 1;
  strcpy(wout->fileName[0], "temp_working_out.dat");
  strcpy(wout->lnd.fileName, wout->fileName[0]);
  if((wout->fp[0] = fopen(wout->fileName[0], "w+b"))==NULL) {
    fprintf(stderr, "\nError in opening working space output file");
    return FAILURE;
  }    

  /* open wbase or wwarp depending on the choice */
  if(choice == USE_BASE || choice == NO_CHANGE) {
    /* copy base and open wbase */
    sprintf(command, "cp -f %s %s", base->lnd.fileName, wbase->lnd.fileName);
    system(command);
    if(strcasecmp(wbase->lnd.fileType, "GEOTIFF")==0) {
      if((wbase->lnd.fp_tiff = XTIFFOpen(wbase->lnd.fileName, "r"))==NULL) {
	fprintf(stderr, "Can't open base GEOTIFF file %s\n", wbase->lnd.fileName);
	return FAILURE;
      } 
    }
    else 
      if((wbase->lnd.fp=fopen(wbase->lnd.fileName, "rb"))==NULL) {
	fprintf(stderr, "Can't open base binary file %s\n", wbase->lnd.fileName);
	return FAILURE;
      }
  }

  if(choice == USE_WARP || choice == NO_CHANGE) {
    /* copy warp and open wwarp */
    sprintf(command, "cp -f %s %s", warp->lnd.fileName, wwarp->lnd.fileName);
    system(command);
    if(strcasecmp(wwarp->lnd.fileType, "GEOTIFF")==0) {
      if((wwarp->fp_tiff[0] = XTIFFOpen(wwarp->fileName[0], "r"))==NULL) {
	fprintf(stderr, "Can't open warp GEOTIFF file %s\n", wwarp->fileName[0]);
	return FAILURE;
      }     
      if((wwarp->lnd.fp_tiff = XTIFFOpen(wwarp->lnd.fileName, "r"))==NULL) {
	fprintf(stderr, "Can't open matching GEOTIFF file %s\n", wwarp->lnd.fileName);
	return FAILURE;
      } 
    }
    else {
      if((wwarp->fp[0]=fopen(wwarp->fileName[0], "rb"))==NULL) {
	fprintf(stderr, "Can't open warp binary file %s\n", wwarp->fileName[0]);
	return FAILURE;
      }
      if((wwarp->lnd.fp = fopen(wwarp->lnd.fileName, "rb"))==NULL) {
	fprintf(stderr, "Can't open matching binary file %s\n", wwarp->lnd.fileName);
	return FAILURE;
      }
    }
  }
 
  switch (choice) {

  case NO_CHANGE: 
    /* set output resolution same as base for easy geolocation validation in working space */
    wout->lnd.res = base->lnd.res;
    break;
    
  case USE_WARP:
    printf("\n\taggregating base image with resolution from %6.2f to %6.2f", base->lnd.res, warp->lnd.res); 
    /* deal with base image */
    wbase->lnd.res = warp->lnd.res;
    wbase->lnd.nrows = base->lnd.nrows * base->lnd.res / wbase->lnd.res;
    wbase->lnd.ncols = base->lnd.ncols * base->lnd.res / wbase->lnd.res;
    strcpy(wbase->lnd.fileType, "BINARY");
    aggregateImage(&(base->lnd), &(wbase->lnd));
    /* reopen newly created wbase file */
    if((wbase->lnd.fp = fopen(wbase->lnd.fileName, "rb"))==NULL) {
      fprintf(stderr, "\nError in opening working base file");
      return FAILURE;
    }
    /* reset working output resolution use warp image */
    wout->lnd.res = warp->lnd.res;
    break;

  case USE_BASE:
    printf("\n\taggregating warp image with resolution from %6.2f to %6.2f", warp->lnd.res, base->lnd.res); 
    /* deal with warp image */
    wwarp->nbands = 1;
    wwarp->lnd.res = base->lnd.res;
    wwarp->lnd.nrows = warp->lnd.nrows * warp->lnd.res / wwarp->lnd.res;
    wwarp->lnd.ncols = warp->lnd.ncols * warp->lnd.res / wwarp->lnd.res;
    strcpy(wwarp->lnd.fileType, "BINARY");
    aggregateImage(&(warp->lnd), &(wwarp->lnd));
    /* reopen newly created wwarp file */
    if((wwarp->fp[0] = fopen(wwarp->fileName[0], "rb"))==NULL) {
      fprintf(stderr, "\nError in opening working warp file");
      return FAILURE;
    }
    if((wwarp->lnd.fp = fopen(wwarp->lnd.fileName, "rb"))==NULL) {
      fprintf(stderr, "\nError in opening working warp file");
      return FAILURE;
    }
    /* use base resolution for output working file */
    wout->lnd.res = base->lnd.res;
    break;
  default:;
  }  

  if(strcasecmp(inPars->extentSource, "DEF")==0) {
    /* use warp image extent for working output image so the whole image are covered */
    wout->lnd.ulx = wwarp->lnd.ulx;
    wout->lnd.uly = wwarp->lnd.uly;
    /* recompute nrows and ncols for working output image */
    wout->lnd.nrows = wwarp->lnd.nrows *  wwarp->lnd.res / wout->lnd.res;
    wout->lnd.ncols = wwarp->lnd.ncols *  wwarp->lnd.res / wout->lnd.res;
  }  
  else {
    wout->lnd.nrows = out->lnd.nrows *  out->lnd.res / wout->lnd.res;
    wout->lnd.ncols = out->lnd.ncols *  out->lnd.res / wout->lnd.res;
  }

  return SUCCESS;
}


/** 
    aggregate image from fine resolution to coarse resolution 
    sum(data_i*area_i)/sum(area_i) 
**/
int aggregateImage(LANDSAT *fine, LANDSAT *coarse)
{
  uint8 *fbuf8;
  int16 *fbuf, cbuf;
  int   irow, icol, m, n;
  float **data, **count;
  float x, y, dx, dy, scale;

  scale = coarse->res / fine->res;

  alloc_1dim_contig ((void **) (&fbuf), fine->ncols, sizeof(int16));
  alloc_1dim_contig ((void **) (&fbuf8), fine->ncols, sizeof(uint8));
  alloc_2dim_contig ((void ***) (&count), coarse->nrows+1, coarse->ncols+1, sizeof(float));
  alloc_2dim_contig ((void ***) (&data), coarse->nrows+1, coarse->ncols+1, sizeof(float));

  if((coarse->fp=fopen(coarse->fileName, "wb"))==NULL) {
    fprintf(stderr, "Error in opening aggregate image for writing\n");
    return FAILURE;
  }  

  for(irow=0; irow<coarse->nrows+1; irow++)
    for(icol=0; icol<coarse->ncols+1; icol++) {
      count[irow][icol] = 0;
      data[irow][icol] = 0.0;
    }

  for(irow=0; irow<fine->nrows; irow++) {
     if(strcasecmp(fine->fileType, "GEOTIFF")==0) {
      /* read data array from GeoTiff file */
      if(fine->nbyte == 1) {
	if (!TIFFReadScanline(fine->fp_tiff, fbuf8, irow, 0)) {
	  fprintf(stderr, "Read line %d error\n", irow);
	  return FAILURE;
	} 
	for(icol=0; icol<fine->ncols; icol++)
	  fbuf[icol] = fbuf8[icol];
      }
      else
	if (!TIFFReadScanline(fine->fp_tiff, fbuf, irow, 0)) {
	  fprintf(stderr, "Read line %d error\n", irow);
	  return FAILURE;
	} 
    }
    else {
      /* read data array from binary file */
      if(fine->nbyte == 1) {
	fread(fbuf8, fine->nbyte, fine->ncols, fine->fp);
	for(icol=0; icol<fine->ncols; icol++)
	  fbuf[icol] = fbuf8[icol];
      }
      else	
	fread(fbuf, fine->nbyte, fine->ncols, fine->fp);
    }

    for(icol=0; icol<fine->ncols; icol++) {
      x = icol/scale;
      y = irow/scale;
      m = (int) x;
      n = (int) y;

      /* exclude fill value */ 
      if(fbuf[icol] == fine->fillValue) {
	count[n][m] = -1;  
	continue;
      }
      /* exclude saturation value */
      if(fbuf[icol] == fine->satValue) {
	count[n][m] = -2;
	continue;
      }

      /* discard it if cell box contains fill or saturation value */
      if(count[n][m]>-0.000001) {
	/* overlap with upper left pixel */
	dx = (m+1)*coarse->res - x;
	dy = (n+1)*coarse->res - y;
	if(dx > fine->res) dx = fine->res;
	if(dx < 0) dx = 0;
	if(dy > fine->res) dy = fine->res;
	if(dy < 0) dy = 0;
	data[n][m] += dx*dy*fbuf[icol];
	count[n][m] += dx*dy;
	/* overlap with upper right pixel */
	dx = x + fine->res - (m+1)*coarse->res;
	dy = (n+1)*coarse->res - y;
	if(dx > fine->res) dx = fine->res;
	if(dx < 0) dx = 0;
	if(dy > fine->res) dy = fine->res;
	if(dy < 0) dy = 0;
	data[n][m] += dx*dy*fbuf[icol];
	count[n][m] += dx*dy;
	/* overlap with lower left pixel */
	dx = (m+1)*coarse->res - x;
	dy = y + fine->res - (n+1)*coarse->res;
	if(dx > fine->res) dx = fine->res;
	if(dx < 0) dx = 0;
	if(dy > fine->res) dy = fine->res;
	if(dy < 0) dy = 0;
	data[n][m] += dx*dy*fbuf[icol];
	count[n][m] += dx*dy;
	/* overlap with lower right pixel */
	dx = x + fine->res - (m+1)*coarse->res;
	dy = y + fine->res - (n+1)*coarse->res;
	if(dx > fine->res) dx = fine->res;
	if(dx < 0) dx = 0;
	if(dy > fine->res) dy = fine->res;
	if(dy < 0) dy = 0;
	data[n][m] += dx*dy*fbuf[icol];
	count[n][m] += dx*dy;
      }  
    }
  }

  for(irow=0; irow<coarse->nrows; irow++)
    for(icol=0; icol<coarse->ncols; icol++) {
      if(count[irow][icol] > 0.000001) 
	cbuf = (int) (data[irow][icol]/count[irow][icol] + 0.5);
      else if(fabs(count[irow][icol]+2) < 0.000001) 
	cbuf = fine->satValue;
      else 
	cbuf = fine->fillValue;

      fwrite(&cbuf, fine->nbyte, 1, coarse->fp);
    }
  
  free(fbuf);
  free_2dim_contig((void **) data);
  free_2dim_contig((void **) count);

  fclose(coarse->fp);  
  return SUCCESS;
}



/* rotate warp image and replace old one */
int rotateWarp(WARP_LANDSAT *warp)
{
  uint8 *fbuf, **data, tmp;
  int   i, j, k, m, n, irow, icol;
  int   start_row, start_col, end_row, end_col;
  int   new_nrows, new_ncols;
  double x, y, x1, y1, xc[4], yc[4];
  double minx, miny, maxx, maxy, new_ulx, new_uly;
  char temp_name[MAX_STRLEN];
  FILE *temp_out;
 
  alloc_1dim_contig ((void **) (&fbuf), warp->lnd.ncols, sizeof(uint8));
  alloc_2dim_contig ((void ***) (&data), warp->lnd.nrows, warp->lnd.ncols, sizeof(uint8));
 
  /* compute new coordinates for four corners after rotation */
  x = 0;
  y = 0;
  rotation(x, y, &x1, &y1, -1.0*warp->orientationAngle);
  xc[0] = x1; yc[0] = y1;
  
  x = warp->lnd.ncols;
  y = warp->lnd.nrows;
  rotation(x, y, &x1, &y1, -1.0*warp->orientationAngle);
  xc[1] = x1; yc[1] = y1;
  
  x = 0;
  y = warp->lnd.nrows;
  rotation(x, y, &x1, &y1, -1.0*warp->orientationAngle);
  xc[2] = x1; yc[2] = y1;

  x = warp->lnd.ncols;
  y = 0;
  rotation(x, y, &x1, &y1, -1.0*warp->orientationAngle);
  xc[3] = x1; yc[3] = y1;
  
  minx = maxx = xc[0];
  miny = maxy = yc[0];
  for(i=1; i<4; i++) {
    if(xc[i] < minx) minx = xc[i];
    if(xc[i] > maxx) maxx = xc[i];
    if(yc[i] < miny) miny = yc[i];
    if(yc[i] > maxy) maxy = yc[i];
  }    
  
  /* compute new image size */
  start_row = (int) miny;
  end_row = (int) maxy;
  start_col = (int) minx;
  end_col = (int) maxx;
  new_ncols = end_col - start_col;
  new_nrows = end_row - start_row;

  /* compute new upper left coordinate */
  new_ulx = warp->lnd.ulx + start_col * warp->lnd.res; 
  new_uly = warp->lnd.uly - start_row * warp->lnd.res;
  
  /* create rotated warp image for each band */
  for(k=0; k<warp->nbands; k++) {

    printf("\n\tfor %s", warp->fileName[k]);
    /* create temporary files for warp image */
    sprintf(temp_name, "temp_rotated_warp.b%d.dat", k+1);
    if((temp_out=fopen(temp_name, "wb"))==NULL) {
      fprintf(stderr, "Error in writing temporary rotation file for warp image\n");
      return FAILURE;
    }

    /* load whole scene to memory for one band */
    for(irow=0; irow<warp->lnd.nrows; irow++) {
      if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0) {
	/* read data array from GeoTiff file */
	if (!TIFFReadScanline(warp->fp_tiff[k], fbuf, irow, 0)) {
	  fprintf(stderr, "Read line %d error\n", irow);
	  return FAILURE;
	} 
      }
      else 
	/* read data array from binary file */
	fread(fbuf, sizeof(uint8), warp->lnd.ncols, warp->fp[k]);
      for(icol=0; icol<warp->lnd.ncols; icol++)
	data[irow][icol] = fbuf[icol];
    }
    
    for(i=start_row; i<end_row; i++)
      for(j=start_col; j<end_col; j++) {
	x = j;
	y = i;

	/* find location in input image */
	rotation(x, y, &x1, &y1, warp->orientationAngle);
	m = (int)(y1+0.5);
	n = (int)(x1+0.5);

	if(m>0 && m<warp->lnd.nrows && n>0 && n<warp->lnd.ncols)
	  tmp = data[m][n];
	else
	  tmp = warp->lnd.fillValue;	
	fwrite(&tmp, sizeof(uint8), 1, temp_out);
      }

    /* replace warp with the rotated warp map */
    fclose(temp_out);
    if(strcasecmp(warp->fileName[k], warp->lnd.fileName)==0)
      strcpy(warp->lnd.fileName, temp_name);
    strcpy(warp->fileName[k], temp_name);
    strcpy(warp->tempRotateFileName[k], temp_name);
    /* close and reopen the rotated images */
    if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0) 
      XTIFFClose(warp->fp_tiff[k]);
    else
      fclose(warp->fp[k]);
    if((warp->fp[k] = fopen(warp->fileName[k], "rb"))==NULL) {
      fprintf(stderr, "Error in opening rotated warp image\n");
      return FAILURE;
    }
  }
  
  /* close and reopen warp matching band */
  if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0) 
    XTIFFClose(warp->lnd.fp_tiff);
  else
    fclose(warp->lnd.fp);
  strcpy(warp->lnd.fileType, "BINARY");
  if((warp->lnd.fp = fopen(warp->lnd.fileName, "rb"))==NULL) {
    fprintf(stderr, "Can't open matching binary file %s\n", warp->lnd.fileName);
    return FAILURE;
  }
  
  warp->lnd.ulx = new_ulx;
  warp->lnd.uly = new_uly;
  warp->lnd.nrows = new_nrows;
  warp->lnd.ncols = new_ncols;
  
  free(fbuf);  
  free_2dim_contig((void **) data);
  return SUCCESS;
}



/* compute new coordinate after rotation: x = icol and y = irow */
void rotation(double x, double y, double *x1, double *y1, double angle)
{
  double sin_cita, cos_cita; 
  sin_cita = sin(angle * PI/180.0);
  cos_cita = cos(angle * PI/180.0);
  *x1 = x * cos_cita - y * sin_cita;
  *y1 = x * sin_cita + y * cos_cita; 
}



/**
 * reproject warp image to the projection used for the base image 
 * keep same warp image spatial resolution  
 * and replace warp files with new ones
 */
int reproject(WARP_LANDSAT *warp, BASE_LANDSAT *base)
{
  char   temp_name[MAX_STRLEN];
  int    i, j, k, percent_done, ret;
  int    new_ncols, new_nrows;
  short int irow, icol;
  double MinX, MaxX, MinY, MaxY;
  double wx, wy, bx, by;    
  uint8  *fbuf, **data, tmp;
  FILE   *temp_out;
 
  alloc_1dim_contig ((void **) (&fbuf), warp->lnd.ncols, sizeof(uint8));
  alloc_2dim_contig ((void ***) (&data), warp->lnd.nrows, warp->lnd.ncols, sizeof(uint8));

  /* Get the coordinate of four coners first */
  MinX=1.0e12;
  MaxX=-1.0e12;
  MinY=1.0e12;
  MaxY=-1.0e12;
  
  printf("\n\tcomputing output image boundary\n");
  for(i=0; i<warp->lnd.nrows+100; i+=100)
    for(j=0; j<warp->lnd.ncols+100; j+=100) {

      if(i>warp->lnd.nrows) irow = warp->lnd.nrows;
      else irow = i;
      if(j>warp->lnd.ncols) icol = warp->lnd.ncols; 
      else icol = j;

      wx = warp->lnd.ulx + icol * warp->lnd.res;
      wy = warp->lnd.uly - irow * warp->lnd.res;
      
      if(convertCoor(base, &bx, &by, &(warp->proj), &wx, &wy, FORWARD_PROJ)==SUCCESS) {
	if(bx>MaxX) MaxX=bx;
        if(bx<MinX) MinX=bx;
        if(by>MaxY) MaxY=by;
        if(by<MinY) MinY=by;
      }
    }    	  
  
  new_ncols = (int) ((MaxX-MinX)/warp->lnd.res+0.5);
  new_nrows = (int) ((MaxY-MinY)/warp->lnd.res+0.5);

  /* create reprojected warp image for each band */
  for(k=0; k<warp->nbands; k++) {
    fprintf(stderr, "\tonly for matching band --- ");
    /* create temporary files for warp image */
    sprintf(temp_name, "temp_reprojected_warp.b%d.dat", k+1);
    if((temp_out=fopen(temp_name, "wb"))==NULL) {
      fprintf(stderr, "Error in writing temporary projection file for warp image\n");
      return FAILURE;
    }

    /* load whole scene to memory for one band */
    for(irow=0; irow<warp->lnd.nrows; irow++) {
      if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0) {
	/* read data array from GeoTiff file */
	if (!TIFFReadScanline(warp->fp_tiff[k], fbuf, irow, 0)) {
	  fprintf(stderr, "Read line %d error\n", irow);
	  return FAILURE;
	} 
      }
      else 
	/* read data array from binary file */
	fread(fbuf, sizeof(uint8), warp->lnd.ncols, warp->fp[k]);
      for(icol=0; icol<warp->lnd.ncols; icol++)
	data[irow][icol] = fbuf[icol];
    }

    for(i=0; i<new_nrows; i++) {

      if(new_nrows % (new_nrows/100)) {
	percent_done = round((i+1)*100.0/new_nrows);
	fprintf(stderr, "%3d\b\b\b", percent_done);
      }
      if(i == new_nrows-1) fprintf(stderr, "100");
      
      for(j=0; j<new_ncols; j++) {
   
	/* use nearest neighbor resampling approach, convert coordinate to the center of pixel  */
	bx = MinX + (j+0.5) * warp->lnd.res;
	by = MaxY - (i+0.5) * warp->lnd.res;

	ret = convertCoor(base, &bx, &by, &(warp->proj), &wx, &wy, BACKWARD_PROJ);

	icol = (int)((wx - warp->lnd.ulx)/warp->lnd.res);
	irow = (int)((warp->lnd.uly - wy)/warp->lnd.res); 

	/* check if this pixel is in the range of input image */
	if(icol>=0 && icol<warp->lnd.ncols && irow>=0 && irow<warp->lnd.nrows && ret==SUCCESS) 
	  tmp = data[irow][icol];
	else 
	  tmp = warp->lnd.fillValue;	

	fwrite(&tmp, sizeof(uint8), 1, temp_out);  
      }
    } /* end of one image processing */

    /* replace warp image with the reprojected warp map */
    fclose(temp_out);
    if(strcasecmp(warp->fileName[k], warp->lnd.fileName)==0)
      strcpy(warp->lnd.fileName, temp_name);
    strcpy(warp->fileName[k], temp_name);
    strcpy(warp->tempProjectFileName[k], temp_name);
    /* close and reopen the reprojected images */
    if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0) 
      XTIFFClose(warp->fp_tiff[k]);
    else
      fclose(warp->fp[k]);
    if((warp->fp[k] = fopen(warp->fileName[k], "rb"))==NULL) {
      fprintf(stderr, "Error in opening reprojected warp image\n");
      return FAILURE;
    }

  } /* end of bands processing */
  
  /* close and reopen warp matching band */
  if(strcasecmp(warp->lnd.fileType, "GEOTIFF")==0) 
    XTIFFClose(warp->lnd.fp_tiff);
  else
    fclose(warp->lnd.fp);
  strcpy(warp->lnd.fileType, "BINARY");
  if((warp->lnd.fp = fopen(warp->lnd.fileName, "rb"))==NULL) {
    fprintf(stderr, "Can't open matching binary file %s\n", warp->lnd.fileName);
    return FAILURE;
  }
  
  /* update upper left coordinates and image size */
  warp->lnd.ulx = MinX;
  warp->lnd.uly = MaxY;
  warp->lnd.nrows = new_nrows;
  warp->lnd.ncols = new_ncols;
  /* reset projection info */
  /*warp->sys = 1;
    warp->utm_zone = base->utm_zone;*/
  
  free(fbuf);  
  free_2dim_contig((void **) data);

  return SUCCESS;
  
}



/**
 * reproject SRTM DEM to the projection used for the base image 
 * keep same DEM spatial resolution  
 * and replace DEM files with new ones
 */
int reprojectDEM(DEM *sdem, BASE_LANDSAT *base)
{
  char   temp_name[MAX_STRLEN];
  int    i, j, ret;
  int    new_ncols, new_nrows;
  short int irow, icol;
  double MinX, MaxX, MinY, MaxY;
  double wx, wy, bx, by;    
  int16  *fbuf, **data, tmp;
  FILE   *temp_out;
 
  alloc_1dim_contig ((void **) (&fbuf), sdem->ncols, sizeof(int16));
  alloc_2dim_contig ((void ***) (&data), sdem->nrows, sdem->ncols, sizeof(int16));

  /* Get the coordinate of four coners first */
  MinX = 1.0e12;
  MaxX = -1.0e12;
  MinY = 1.0e12;
  MaxY = -1.0e12;
  
  for(i=0; i<sdem->nrows+100; i+=100)
    for(j=0; j<sdem->ncols+100; j+=100) {

      if(i>sdem->nrows) irow = sdem->nrows;
      else irow = i;
      if(j>sdem->ncols) icol = sdem->ncols; 
      else icol = j;

      wx = sdem->ulx + icol * sdem->res;
      wy = sdem->uly - irow * sdem->res;
      
      if(convertCoor(base, &bx, &by, &(sdem->proj), &wx, &wy, FORWARD_PROJ)==SUCCESS) {
	if(bx>MaxX) MaxX=bx;
        if(bx<MinX) MinX=bx;
        if(by>MaxY) MaxY=by;
        if(by<MinY) MinY=by;
      } 
    }    	  
  
  new_ncols = (int) ((MaxX-MinX)/sdem->res+0.5);
  new_nrows = (int) ((MaxY-MinY)/sdem->res+0.5);

  /* create temporary files for DEM */
  sprintf(temp_name, "temp_reprojected_dem.dat");
  if((temp_out=fopen(temp_name, "wb"))==NULL) {
    fprintf(stderr, "Error in writing temporary projection file for DEM\n");
    return FAILURE;
  }

  /* load whole scene to memory for one band */
  for(irow=0; irow<sdem->nrows; irow++) {
    /* read data array from GeoTiff file */
    if (!TIFFReadScanline(sdem->fp_tiff, fbuf, irow, 0)) {
      fprintf(stderr, "Read line %d error\n", irow);
      return FAILURE;
    } 
    for(icol=0; icol<sdem->ncols; icol++)
      data[irow][icol] = fbuf[icol];
  }
  
  for(i=0; i<new_nrows; i++) {
    
    for(j=0; j<new_ncols; j++) {
   
      /* use nearest neighbor resampling approach, convert coordinate to the center of pixel  */
      bx = MinX + (j+0.5) * sdem->res;
      by = MaxY - (i+0.5) * sdem->res;
      
      ret = convertCoor(base, &bx, &by, &(sdem->proj), &wx, &wy, BACKWARD_PROJ);
      
      icol = (int)((wx - sdem->ulx)/sdem->res);
      irow = (int)((sdem->uly - wy)/sdem->res); 

      /* check if this pixel is in the range of input image */
      if(icol>=0 && icol<sdem->ncols && irow>=0 && irow<sdem->nrows && ret==SUCCESS) 
	tmp = data[irow][icol];
      else 
	tmp = DEM_FILL;	

      fwrite(&tmp, sizeof(int16), 1, temp_out);  
    }
  }

  /* replace warp image with the reprojected warp map */
  fclose(temp_out);
  XTIFFClose(sdem->fp_tiff);

  strcpy(sdem->fileName, temp_name);
  strcpy(sdem->fileType, "BINARY");
  if((sdem->fp = fopen(sdem->fileName, "rb"))==NULL) {
    fprintf(stderr, "Error in opening reprojected DEM image\n");
    return FAILURE;
  }
    
  /* update upper left coordinates and image size */
  sdem->ulx = MinX;
  sdem->uly = MaxY;
  sdem->nrows = new_nrows;
  sdem->ncols = new_ncols;
  /* reset projection info */
  sdem->proj.sys = UTM;
  sdem->proj.utm_zone = base->utm_zone;
  sdem->proj.datum = base->datum;
    
  free(fbuf);  
  free_2dim_contig((void **) data);

  free(sdem->buf);
  alloc_1dim_contig((void **) (&sdem->buf), sdem->ncols, sizeof(int16));

  return SUCCESS;
  
}



/* copy necessary Landsat data structure */ 
void copyLandsat(LANDSAT *dest, LANDSAT *src)
{
  int i;

  strcpy(dest->fileName, src->fileName);
  strcpy(dest->fileType, src->fileType);
  dest->ulx = src->ulx;
  dest->uly = src->uly;
  for(i=0; i<2; i++)
    dest->centerLonLat[i] = src->centerLonLat[i];
  dest->nrows = src->nrows;
  dest->ncols = src->ncols;
  dest->res = src->res;
  dest->nbyte = src->nbyte;
  dest->fillValue = src->fillValue;
  dest->satValue = src->satValue;
  for(i=0; i<2; i++)
    dest->lpara[i] = src->lpara[i];
  dest->radius = src->radius;
  dest->satMajorAxis = src->satMajorAxis;
  dest->altitude = src->altitude;
  dest->porder = src->porder;
  for(i=0; i<6; i++) {
    dest->a[i] = src->a[i];
    dest->b[i] = src->b[i];
  }
}


/** Copy orignal warp image to temporary warp image for rotation and reprojection
 *  It only contains matching band for ortho testing and verification 
 */
int createTempWarp(WARP_LANDSAT  *twarp, WARP_LANDSAT  *warp)
{
  int i;

  twarp->nbands = 1;
  copyLandsat(&(twarp->lnd), &(warp->lnd));  
  strcpy(twarp->fileName[0], twarp->lnd.fileName);
  strcpy(twarp->resampleMethod, "NN"); 
  twarp->orientationAngle = warp->orientationAngle;
  twarp->pointingAngle = warp->pointingAngle;
  twarp->proj.sys = warp->proj.sys;
  twarp->proj.utm_zone = warp->proj.utm_zone;
  for (i=0; i<15; i++)
    twarp->proj.param[i] = warp->proj.param[i];
  twarp->proj.unit = warp->proj.unit;
  twarp->proj.datum = warp->proj.datum;
  for(i=0; i<warp->nbands; i++) twarp->nbyte[i] = warp->nbyte[i];

  if(strcasecmp(twarp->lnd.fileType, "GEOTIFF")==0) {
    if((twarp->fp_tiff[0] = XTIFFOpen(twarp->fileName[0], "r"))==NULL) {
      fprintf(stderr, "Can't open warp GEOTIFF file %s\n", twarp->fileName[0]);
      return FAILURE;
    }     
    if((twarp->lnd.fp_tiff = XTIFFOpen(twarp->lnd.fileName, "r"))==NULL) {
      fprintf(stderr, "Can't open matching GEOTIFF file %s\n", twarp->lnd.fileName);
      return FAILURE;
    } 
  }
  else {
    if((twarp->fp[0]=fopen(twarp->fileName[0], "rb"))==NULL) {
      fprintf(stderr, "Can't open warp binary file %s\n", twarp->fileName[0]);
      return FAILURE;
    }
    if((twarp->lnd.fp = fopen(twarp->lnd.fileName, "rb"))==NULL) {
      fprintf(stderr, "Can't open matching binary file %s\n", twarp->lnd.fileName);
      return FAILURE;
    }
  }

  return SUCCESS;  
}


/* update output metadata from source */
void updateOutCoor(LANDSAT *src, LANDSAT *dest) {
  dest->ulx = src->ulx;
  dest->uly = src->uly;
  dest->nrows = (int)(src->nrows * src->res / dest->res + 0.5);
  dest->ncols = (int)(src->ncols * src->res / dest->res + 0.5);
}


/*** WRITE AND PRINT OUT ENVI HAEDER FILES FOR EACH BAND ***/
int writeENVIheader(WARP_LANDSAT *warp, OUT_LANDSAT *out)
{
  char fname[MAX_STRLEN];
  int  iband;
  FILE *fp;

  for(iband=0; iband<out->nbands; iband++) {
    sprintf(fname, "%s.hdr", out->fileName[iband]);
    if((fp=fopen(fname, "w"))==NULL) {
      fprintf(stderr, "can't write header for output file %s\n", out->fileName[iband]);
      return FAILURE;
    }
    fprintf(fp, "ENVI\n");
    fprintf(fp, "description = {\n");
    fprintf(fp, "orthorectified Landsat image from %s\n", warp->fileName[iband]);
    fprintf(fp, "resampling method: %s}\n", warp->resampleMethod);
    fprintf(fp, "samples = %d\n", out->lnd.ncols);
    fprintf(fp, "lines = %d\n", out->lnd.nrows);
    fprintf(fp, "bands = 1\n");
    fprintf(fp, "file type = ENVI Standard\n");
    fprintf(fp, "data type = %d\n", out->nbyte[iband]);
    fprintf(fp, "interleave = bsq\n");
    fprintf(fp, "sensor type = Landsat\n");
    fprintf(fp, "byte order = 0\n");
    fprintf(fp, "map info = {UTM, 1.000, 1.000, %f, %f, %f, %f, ", out->lnd.ulx, out->lnd.uly, out->lnd.res, out->lnd.res); 
    if(out->utm_zone > 0)
      fprintf(fp, "%d, North, WGS-84, units=Meters}\n", abs(out->utm_zone));
    else
      fprintf(fp, "%d, South, WGS-84, units=Meters}\n", abs(out->utm_zone));   
    fclose(fp);
  }
  return SUCCESS;
}

/* display usage information */
void usage(char *argv[]) 
{
  fprintf(stderr, "\nUsage: %s [-r|-o|-b|-v] <parameter_file>", argv[0]);
  fprintf(stderr, "\n       -r  do registration only (assume warp image has been orthorectified)");
  fprintf(stderr, "\n       -o  do orthorectification only (assume warp image has been registrated)");
  fprintf(stderr, "\n       -b  do both registration and othorectification (recommend for L1G data)");
  fprintf(stderr, "\n       -v  do image matching verification for two input images\n\n");
}


/* close and remove coarse resolution file */
void clean_coarse(LANDSAT *sr)
{
  close_LANDSAT_file(sr);
#ifndef DEBUG
  remove (sr->fileName);
#endif
}

/* close file defined in LANDSAT structure */
void close_LANDSAT_file(LANDSAT *sr)
{
  if(strcasecmp(sr->fileType,"GEOTIFF")==0) 
    XTIFFClose(sr->fp_tiff);
  else
    fclose(sr->fp);
}


/* 1D array allocation */
void alloc_1dim_contig (void **ptr, int d1, int elsize)
{
   void *p = NULL;

   p = calloc (d1, elsize);
   if (!p) {
      printf ("Memory allocation error in alloc_1dim_contig\n");
      exit(1);
   }
   *ptr = p;
   return;
}


/* 2D array allocation */
void alloc_2dim_contig (void ***ptr, int d1, int d2, int elsize)
{
   void *p = NULL;
   void **pp = NULL;
   int i = 0;

   /* alloc array for data */
   alloc_1dim_contig ((void **) (&p), d1 * d2, elsize);

   /* alloc array for pointers */
   alloc_1dim_contig ((void **) (&pp), d1, sizeof (void *));

   /* Set the pointers to indicate the appropriate elements of the data array. */
   for (i = 0; i < d1; i++) {
      pp[i] = (char *) p + (i * d2 * elsize);
   }

   *ptr = pp;

   return;
}


/* 3D array allocation */
void alloc_3dim_contig (void ****ptr, int d1, int d2, int d3, int elsize)
{
   void *p = NULL;
   void **pp = NULL;
   void ***ppp = NULL;
   int i = 0;

   /* allocate the data array */
   alloc_1dim_contig ((void **) (&p), d1 * d2 * d3, elsize);

   /* alloc the double pointers */
   alloc_1dim_contig ((void **) (&pp), d1 * d2, sizeof (void *));

   /* and again for the first dimensional pointers */
   alloc_1dim_contig ((void **) (&ppp), d1, sizeof (void **));

   /* first set the d1 pointers */
   for (i = 0; i < d1; i++) {
      ppp[i] = pp + (i * d2);
   }

   /* next set all of the d2 pointers */
   for (i = 0; i < d1 * d2; i++) {
      pp[i] = (char *) p + (i * d3 * elsize);
   }

   *ptr = ppp;

   return;
}


void free_2dim_contig (void **a)
{
  free (a[0]);
  free (a);
  return;
} 


void free_3dim_contig (void ***a)
{
   free (a[0][0]);
   free (a[0]);
   free (a);
   return;
}

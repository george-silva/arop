/**
 * ! Description
 *   geolocation checking and process subroutines for Landsat orthorectification program 
 *
 * ! Credits
 *   Feng Gao (Feng.Gao@nasa.gov), Developer
 *   Jeff Masek (Jeffrey.G.Masek@nasa.gov)
 *   Robert Wolfe (Robert.E.Wolfe@nasa.gov)
 *
 * ! Revision 2.2.5  08/22/2011
 * ! Revision 2.2.5  01/12/2010
 * ! Revision 2.2.3  02/27/2009
 * ! Revision 2.2.2  10/30/2008
 * ! Revision 2.2    05/20/2008
 */

#include "lndortho.h"


/**
 * do preliminary registration based on coarse images and apply shifts to warp image if needs 
 */
int preliminaryRegistration(BASE_LANDSAT  *wbase, WARP_LANDSAT  *wwarp)
{
  int    num_usable, num_total, ret;
  float  **match_p;
  CP_PARAMETERS CP_KEYS_BAK;
  /* coarse resolution base and warp image for preliminary registration */
  LANDSAT cbase, cwarp;

  printf("\tcreating temporary coarse resolution images ...\n");  
  strcpy(cbase.fileName, "temp_working_base_coarse.dat");
  if(makeCoarse(&(wbase->lnd), &cbase) == FAILURE) {
    fprintf(stderr, "Make coarse resolution file %s error!\n", cbase.fileName);
    return FAILURE;
  }

  strcpy(cwarp.fileName, "temp_working_warp_coarse.dat");
  if(makeCoarse(&(wwarp->lnd), &cwarp) == FAILURE) {
    fprintf(stderr, "Make coarse resolution file %s error!\n", cwarp.fileName);
    return FAILURE;
  }

  fprintf(stderr, "\tsearching tie points between two coarse resolution images ...");
  /* save copy for restore */
  CP_KEYS_BAK = CP_KEYS;
  /* change settings for coarse resolution image */
  CP_KEYS.CP_SEED_WIN = CP_KEYS.COARSE_CP_SEED_WIN; 
  CP_KEYS.MAX_SHIFT = CP_KEYS.COARSE_MAX_SHIFT;
  CP_KEYS.MAX_CPS = ((int) (cbase.nrows/(CP_KEYS.CP_SEED_WIN)) + 1) * ((int) (cbase.ncols/(CP_KEYS.CP_SEED_WIN)) + 1);
  CP_KEYS.MAX_SIZE = 2*(CP_KEYS.MAX_SHIFT+CP_KEYS.CHIP_SIZE);  
  alloc_2dim_contig((void ***) (&match_p), CP_KEYS.MAX_CPS, 7, sizeof(float));

#ifdef REG_DEBUG
    printf("\nCoarse Image Info: (max_n_cps=%d)", CP_KEYS.MAX_CPS);
    printf("\nCBASE: nrows=%d ncols=%d res=%f", cbase.nrows, cbase.ncols, cbase.res);
    printf("\nCWARP: nrows=%d ncols=%d res=%f", cwarp.nrows, cwarp.ncols, cwarp.res);
    printf("\nseed_window_size=%d  max_shift=%d", CP_KEYS.CP_SEED_WIN, CP_KEYS.MAX_SHIFT);
    printf("\n\n=== Control Point Candidates Based on Coarse Resolution Images ===\n");
    printf("Base_x,y\tProjected_x,y\tMatched_x,y\tSub_Pixel_Match\t\tCorr.\tNum_High_Corr.\n");
#endif
    /* find control points from coarse images */
    num_total = searchCPs(&cbase, &cwarp, match_p);
    /* filter out outliers */
    num_usable = filterCandidates(&cbase, &cwarp, match_p, num_total);
    if(num_usable >= CP_KEYS.MIN_ACCEPTABLE_NCP) {
      /* if difference is larger than one coarse resolution pixel, then update it */
      /* or in other word, we trust control point accuracy is within one pixel */
      if(fabs(wwarp->lnd.ulx-cwarp.ulx)/wwarp->lnd.res > CP_KEYS.COARSE_SCALE ||
	 fabs(wwarp->lnd.uly-cwarp.uly)/wwarp->lnd.res > CP_KEYS.COARSE_SCALE ) {
	printf("\n\tpreliminary registration performed!!!");
	printf("\n\t\tupper left coordinate has been updated");
	printf("\n\t\tfrom (%7.1f, %7.1f)", wwarp->lnd.ulx, wwarp->lnd.uly);
	printf("\n\t\tto   (%7.1f, %7.1f)", cwarp.ulx, cwarp.uly);
	wwarp->lnd.ulx = cwarp.ulx;
	wwarp->lnd.uly = cwarp.uly;
      }
      else 
	printf("\n\tpreliminary registration is not performed since current coordinates are close enough");
      ret = SUCCESS;
    }
    else {
      printf("\n\tpreliminary registration can not be performed due to lack of good control points");
      ret = FAILURE;
    }

    /* restore initial settings for fine resolution image */
    CP_KEYS = CP_KEYS_BAK;

    /* clean up all work for coarse resolution images */
    free_2dim_contig((void **)match_p);
    clean_coarse(&cbase);
    clean_coarse(&cwarp);

    return ret;
}



/**
 * do precise registration and apply shifts to warp image if needs
 */
int preciseRegistration (BASE_LANDSAT  *wbase, WARP_LANDSAT  *wwarp, DEM *sdem)
{
  int    num_usable, num_total;
  float  **match_p;
  double tx, ty;
#ifdef REG_DEBUG  
  int i;
#endif

  fprintf(stderr, "\tsearching tie points between BASE and WARP image -- ");
#ifdef REG_DEBUG
  printf("\n\t\tfound tie points: ");
  printf("\n\n=== Control Point Candidates ===\n");
  printf("Base_x,y\tProjected_x,y\tMatched_x,y\tSub_Pixel_Match\t\tCorr.\tNum_High_Corr.\n");
#endif

  CP_KEYS.MAX_CPS = ((int) (wbase->lnd.nrows/CP_KEYS.CP_SEED_WIN) + 1) * ((int) (wbase->lnd.ncols/CP_KEYS.CP_SEED_WIN) + 1);
  alloc_2dim_contig((void ***) (&match_p), CP_KEYS.MAX_CPS, 7, sizeof(float));
  num_total = searchCPs(&(wbase->lnd), &(wwarp->lnd), match_p);
  printf("\n\tcomputing new upper left corner for WARP images based on control points");
#ifdef REG_DEBUG  
  printf("max_cps==%d\n",CP_KEYS.MAX_CPS); 
  printf("\n\n=== Compute Terrain Displacement ===\n");
#endif
  
  /* convert to orthorectified coordinates */
  if(computeOrthoShifts(&(wbase->lnd), sdem, match_p, num_total)==FAILURE) {
    fprintf(stderr, "compute orthorectification shiftments error\n");
    return FAILURE;
  }

#ifdef REG_DEBUG  
  printf("\n\n=== Filtered Control Points (Final Selection) ===\n");
#endif

  /* save current uppper left corner */
  tx = wwarp->lnd.ulx;
  ty = wwarp->lnd.uly;
  
  /* filter out outliers */
  num_usable = filterCandidates(&(wbase->lnd), &(wwarp->lnd), match_p, num_total);

#ifdef REG_DEBUG
  for(i=0; i<num_total; i++)
    if(fabs(match_p[i][4]-1)<0.0000001)  
      printf("%6.1f %6.1f\t%6.1f %6.1f\t%6.1f %6.1f\n", match_p[i][0], match_p[i][1], match_p[i][2], match_p[i][3], match_p[i][5], match_p[i][6]);
#endif

  if(num_total <  CP_KEYS.MIN_ACCEPTABLE_NCP) {
    printf("\n- Orthorectification Processing Fails!");
    printf("\n\tPrecise registration can't be applied due to the limited numbers of control points");
    printf("\n\tPlease check initial parameter file");
    printf("\n\tYou may increase MAX_SHIFT to increase searching range to account for the large geolocation error");
    printf("\n\tAND/OR decrease control point seed window size CP_SEED_WIN to increase number of control points"); 
    printf("\n\tAND/OR use preliminary registration for images with large geolocation error\n");
    return FAILURE;
  }  
 
  printf("\n\t\tNum_CP_Candidates: %d\tNum_CP_Selected: %d",num_total, num_usable);
  printf("\n\t\tbase image UL:  (%10.1f, %10.1f)", wbase->lnd.ulx, wbase->lnd.uly);
  printf("\n\t\twarp image UL:  (%10.1f, %10.1f)", tx, ty);
 
  /* add new upleft corner metadata */
  if(num_usable >= CP_KEYS.MIN_ACCEPTABLE_NCP) {
    printf("\n\t\twarp image adjusted UL: (%10.1f, %10.1f)", wwarp->lnd.ulx, wwarp->lnd.uly);   
    printf("\n\tprecise registration is applied to all warp images");
  }
  else {
    wwarp->lnd.ulx = tx;
    wwarp->lnd.uly = ty;
    printf("\n\tsimple bias precise registration can NOT be applied");
    printf("\n\tdue to great distortion or limited usable control points");
  }
  free_2dim_contig((void **)match_p);

  return SUCCESS;
}


/** 
 * do precise registration with sampling (assume both base and warp images have been orthorectified
 * or both not orthorectified
 */
int  preciseRegistrationWithSampling (BASE_LANDSAT  *wbase, WARP_LANDSAT  *warp, WARP_LANDSAT  *temp_warp, WARP_LANDSAT  *wwarp, OUT_LANDSAT *out, OUT_LANDSAT *wout, IN_PARAMETERS *inPars)
{
  char   command[MAX_STRLEN];
  int    i, j, k, iband, tmp, num_usable, num_total;
  int    current_row=0, percent_done, checking_passed;
  float  **match_p;
  double tx, ty;
  double ai, aj, wi, wj, out_scale;
  double bx, by, wx, wy, rx, ry;
  double shift_x, shift_y;
  int16  **outRow, ***inRows;        /* hold on input and output lines */
  PIXEL   *px;
  LANDSAT fine_out, agg_out;         /* define output for aggregation option */  

  /* compute shift from preliminary and precise registration */
  shift_x = temp_warp->lnd.ulx - wwarp->lnd.ulx;
  shift_y = wwarp->lnd.uly - temp_warp->lnd.uly;

  if(!(px   = malloc(sizeof(PIXEL)))) return FAILURE;

  fprintf(stderr, "\tsearching tie points between BASE and WARP image -- ");
#ifdef REG_DEBUG
  printf("\n\t\tfound tie points: ");
  printf("\n\n=== Control Point Candidates ===\n");
  printf("Base_x,y\tProjected_x,y\tMatched_x,y\tSub_Pixel_Match\t\tCorr.\tNum_High_Corr.\n");
#endif

  CP_KEYS.MAX_CPS = ((int) (wbase->lnd.nrows/CP_KEYS.CP_SEED_WIN) + 1) * ((int) (wbase->lnd.ncols/CP_KEYS.CP_SEED_WIN) + 1);
  alloc_2dim_contig((void ***) (&match_p), CP_KEYS.MAX_CPS, 7, sizeof(float));
  num_total = searchCPs(&(wbase->lnd), &(wwarp->lnd), match_p);
  printf("\n\tcomputing new upper left corner for WARP images based on control points");

#ifdef REG_DEBUG  
  printf("\n\n=== Filtered Control Points (Final Selection) ===\n");
#endif

  /* save current uppper left corner */
  tx = wwarp->lnd.ulx;
  ty = wwarp->lnd.uly;
  
  /* filter out outliers */
  num_usable = filterCandidates(&(wbase->lnd), &(wwarp->lnd), match_p, num_total);

#ifdef REG_DEBUG
  for(i=0; i<num_total; i++)
    if(fabs(match_p[i][4]-1)<0.0000001) { 
      printf("%6.1f %6.1f\t%6.1f %6.1f\n", match_p[i][0], match_p[i][1], match_p[i][2], match_p[i][3]);
      break;
    }
#endif

  if(num_total <  CP_KEYS.MIN_ACCEPTABLE_NCP) {
    printf("\n- Orthorectification Processing Fails!");
    printf("\n\tPrecise registration can't be applied due to the limited numbers of control points");
    printf("\n\tPlease check initial parameter file");
    printf("\n\tYou may increase MAX_SHIFT to increase searching range to account for the large geolocation error");
    printf("\n\tAND/OR decrease control point seed window size CP_SEED_WIN to increase number of control points"); 
    printf("\n\tAND/OR use preliminary registration for images with large geolocation error\n");
    return FAILURE;
  }  
 
  printf("\n\t\tNum_CP_Candidates: %d\tNum_CP_Selected: %d",num_total, num_usable);
  printf("\n\t\tbase image UL:  (%10.1f, %10.1f)", wbase->lnd.ulx, wbase->lnd.uly);
  printf("\n\t\twarp image UL:  (%10.1f, %10.1f)", tx, ty);
 
  /* add new upleft corner metadata */
  if(num_usable >= CP_KEYS.MIN_ACCEPTABLE_NCP) {
    printf("\n\t\twarp images may be adjusted with new UL: (%10.1f, %10.1f)", wwarp->lnd.ulx, wwarp->lnd.uly);
    /* don't update UL, I will try 1st order instead */
    wwarp->lnd.ulx = tx;
    wwarp->lnd.uly = ty;;
  }   
  else 
    printf("\n\t\tneed try higher order polynomial transformation");

  /* use all control points */
  for(i=0; i<num_total; i++) match_p[i][4] = 1;
  redoPreciseGeoCorrection(match_p, num_total, &(wout->lnd));
  /* check registered results */
  for(i=0; i<num_total; i++) {
    match_p[i][0] = match_p[i][5];
    match_p[i][1] = match_p[i][6];
  }
  printf("\tverifying tie points");
  checking_passed = checkingTPs(&(wbase->lnd), &(wout->lnd), &(wwarp->lnd), match_p, num_total);
  if(checking_passed)
    printf("\n\t\tmatching test passed");
  else
    printf("\n\t\tmatching test failed, use with CAUTION!!!");

  printf("\n\tapplying precise registration to all warp images\n");
  
  /* allocate memory to save one row of output */
  /* choice the largest nrows since output image maybe finer or coarser than working file */
  tmp = wout->lnd.ncols * wout->lnd.res / warp->lnd.res;
  if(out->lnd.ncols > tmp) tmp = out->lnd.ncols;
  alloc_2dim_contig((void ***) (&outRow), tmp, out->nbands, sizeof(int16));
  /* allocate memory to save a block of warp input */
  alloc_3dim_contig((void ****) (&inRows), BLOCK_NROWS, warp->lnd.ncols, warp->nbands, sizeof(int16));

  px->res = wout->lnd.res;
  px->ulx = wout->lnd.ulx;
  px->uly = wout->lnd.uly;

  if(strcasecmp(warp->resampleMethod, "AGG")==0) {
    /* use NN for fine resolution output and then aggregate to coarse output */
    copyLandsat(&agg_out, &(out->lnd));
    /* use original spatial resolution */
    out->lnd.nrows = out->lnd.nrows * out->lnd.res / warp->lnd.res;
    out->lnd.ncols = out->lnd.ncols * out->lnd.res / warp->lnd.res;
    out->lnd.res = warp->lnd.res; 
  }
  
  fprintf(stderr, "\tregistration processing (in percent) --");

  /* in case output resolution different in working and output space */      
  out_scale = out->lnd.res / wout->lnd.res;

  /* foreach row in the orthorectified image space (output) */
  for(i=0; i<out->lnd.nrows; i++) {
    
    if(out->lnd.nrows % (out->lnd.nrows/100)) {
      percent_done = round((i+1)*100.0/out->lnd.nrows);
      fprintf(stderr, "%3d\b\b\b", percent_done);
    }
    if(i == out->lnd.nrows-1) fprintf(stderr, "100");

    /* convert i to working space (use center of pixel) */
    px->irow = (double) (wout->lnd.uly - (out->lnd.uly - (i+0.5) * out->lnd.res)) / wout->lnd.res;
    wi = px->irow;
    
    /* foreach pixel in output file */
    for(j=0; j<out->lnd.ncols; j++) {
	
      /* initialize output arrays */
      for(iband=0; iband<out->nbands; iband++)
	outRow[j][iband] = out->lnd.fillValue;
	
      /* convert j to working space (use center of pixel) */
      px->icol = (double) ((out->lnd.ulx + (j+0.5) * out->lnd.res) - wout->lnd.ulx) / wout->lnd.res;
      wj = px->icol;
	
      /* convert from output coordinate to base coordinate */
      ai = (wbase->lnd.uly - wout->lnd.uly + wi * wout->lnd.res)/wbase->lnd.res;
      aj = (wout->lnd.ulx - wbase->lnd.ulx + wj * wout->lnd.res)/wbase->lnd.res;

      /* precise registration using polynomial transformation function 
	 convert coordinate from base to warp in working space */
      if(wout->lnd.porder == 1) { 
	px->irow = wout->lnd.b[0]+wout->lnd.b[1]*aj+wout->lnd.b[2]*ai;
	px->icol = wout->lnd.a[0]+wout->lnd.a[1]*aj+wout->lnd.a[2]*ai;
      }
      else {
	px->irow = wout->lnd.b[0]+wout->lnd.b[1]*aj+wout->lnd.b[2]*ai+
	  wout->lnd.b[3]*aj*aj+wout->lnd.b[4]*aj*ai+wout->lnd.b[5]*ai*ai;
	px->icol = wout->lnd.a[0]+wout->lnd.a[1]*aj+wout->lnd.a[2]*ai+
	  wout->lnd.a[3]*aj*aj+wout->lnd.a[4]*aj*ai+wout->lnd.a[5]*ai*ai;
      }

      /* convert to actual coordinates of wwarp in working space  */
      px->ux = wwarp->lnd.ulx + px->icol * wwarp->lnd.res;
      px->uy = wwarp->lnd.uly - px->irow * wwarp->lnd.res;
      
      /* adjust to coordinate before shift registration (to temp_warp) */
      px->ux += shift_x;
      px->uy -= shift_y;

      /* convert coordinate if it is reprojected (to proj_warp)*/
      if(inPars->projection_flag) {
	bx = px->ux;
	by = px->uy;
	convertCoor(wbase, &bx, &by, &(warp->proj), &wx, &wy, BACKWARD_PROJ);
	px->ux = wx;
	px->uy = wy;
	/* convert coordinate if it is rotated (to rotate_warp)*/
	if(inPars->rotation_flag) {
	  wx = wx - warp->lnd.ulx;
	  wy = warp->lnd.uly - wy;
	  rotation(wx, wy, &rx, &ry, warp->orientationAngle);
	  px->ux = warp->lnd.ulx + rx;
	  px->uy = warp->lnd.uly - ry;
	}
      }
      else 
	/* convert coordinate if it is rotated */
	if(inPars->rotation_flag) {
	  wx = px->ux - warp->lnd.ulx;
	  wy = warp->lnd.uly - px->uy;
	  rotation(wx, wy, &rx, &ry, warp->orientationAngle);
	  px->ux = warp->lnd.ulx + rx;
	  px->uy = warp->lnd.uly - ry;
	}

      /* convert working space back to the input space  (use warp) */
      px->urow = (int)((warp->lnd.uly - px->uy)/warp->lnd.res);
      px->ucol = (int)((px->ux - warp->lnd.ulx)/warp->lnd.res);

      /* process first pixel and load data block for next uses for this row */
      if(j == 0) {	  
	/* load data block into memory */
	current_row = px->urow;
	if(loadDataBlock(warp, current_row, inRows)==FAILURE) {
	  fprintf(stderr, "load data block error!\n");
	  return FAILURE;
	}	
      }

      /* do resample and get value */
      resampleData(px, warp, current_row, inRows, outRow[j]); 
	
#ifdef DEBUG
      if(i*out_scale==DEBUG_irow && j*out_scale==DEBUG_icol) {
	printf("\nin real output (%d %d)\n", j, i);
	printf("in working base (%f %f)\n", aj, ai);
	printf("in working warp (%f %f)\n", px->icol, px->irow);
	printf("working warp UL (%8.1f, %8.1f) res=%5.1f (%8.1f, %8.1f)\n", 
	       wwarp->lnd.ulx, wwarp->lnd.uly, wwarp->lnd.res, px->ux, px->uy);
	printf("in real warp UL (%8.1f, %8.1f) res=%5.1f\n", 
	       warp->lnd.ulx, warp->lnd.uly, warp->lnd.res);		
	printf("in real warp (%4d,%4d)\n", px->ucol, px->urow);
	for(k=0; k<out->nbands; k++)
	  printf("band%d = %d\n", k+1, outRow[j][k]);
      }
#endif      
	
    }   /* end of column looping */
      
    /* write one row to output file in output space */
    if(writeOneRow(out, outRow)==FAILURE) {
      fprintf(stderr, "Write line error iline=%d\n!", i);
      return FAILURE;
    }
  } /* end of line looping */
  
  /* close output files */
  for(iband=0; iband<out->nbands; iband++)
    fclose(out->fp[iband]);    

  /* do aggregation for outputs that are coarser than inputs */
  if(strcasecmp(warp->resampleMethod, "AGG")==0) {
    copyLandsat(&fine_out, &(out->lnd));
    for(iband=0; iband<out->nbands; iband++) {
      printf("\n\taggregating output for band%d", iband+1);
      /* copy fine resolution file to a temporary file and open it */
      sprintf(command, "cp %s %s", out->fileName[iband], wout->lnd.fileName);
      system(command);
      strcpy(fine_out.fileName, wout->lnd.fileName);
      if((fine_out.fp=fopen(fine_out.fileName, "rb"))==NULL) {
	fprintf(stderr, "Open output file for aggregation error!\n");
	return FAILURE;
      }
      /* aggregate temporary file and rewrite output file */
      strcpy(agg_out.fileName, out->fileName[iband]);
      fine_out.satValue = warp->satValue[iband];
      fine_out.nbyte = out->nbyte[iband];
      aggregateImage(&fine_out, &agg_out);
      fclose(fine_out.fp);
    }
    /* restore original output metadata */
    copyLandsat(&(out->lnd), &agg_out);
  }
  
  /* also print output file information to screen */
  printf("\n\tprecise registration output filename:");
  for(iband=0; iband<out->nbands; iband++) 
    printf("\n\t\t%s", out->fileName[iband]);
  printf("\n\tnumber of samples = %d", out->lnd.ncols);
  printf("\n\tnumber of lines   = %d", out->lnd.nrows);
  printf("\n\tupper left X = %f", out->lnd.ulx);
  printf("\n\tupper left Y = %f", out->lnd.uly);
  printf("\n\toutput pixel resolution = %f", out->lnd.res);
  printf("\n\tUTM zone number = %d", out->utm_zone);
  printf("\n\tresampling method: %s", warp->resampleMethod);

  free(px);
  free_2dim_contig((void **) outRow);
  free_2dim_contig((void **) match_p);
  free_3dim_contig((void ***) inRows);
  
  return SUCCESS;
  
}



/** 
 * do image matching verification 
 */
int  matchingVerification (BASE_LANDSAT  *wbase, WARP_LANDSAT  *wwarp, OUT_LANDSAT *wout)
{
  int    i, num_usable, num_total, checking_passed;
  float  **match_p;
  PIXEL  *px;

  if(!(px   = malloc(sizeof(PIXEL)))) return FAILURE;

  /*printf("\tsearching tie points between BASE and WARP image -- ");
  printf("\n\t\tfound matching tie points: ");
  printf("\n\n=== Control Point Candidates ===\n");
  printf("Base_x,y\tProjected_x,y\tMatched_x,y\tSub_Pixel_Match\t\tCorr.\tNum_High_Corr.\n");*/

  CP_KEYS.MAX_CPS = ((int) (wbase->lnd.nrows/CP_KEYS.CP_SEED_WIN) + 1) * ((int) (wbase->lnd.ncols/CP_KEYS.CP_SEED_WIN) + 1);
  alloc_2dim_contig((void ***) (&match_p), CP_KEYS.MAX_CPS, 7, sizeof(float));
  num_total = searchCPs(&(wbase->lnd), &(wwarp->lnd), match_p);

  /*printf("\n\n=== Filtered Control Points ===\n");
    printf("Base_x Base_y\tWarp_x Warp_y\n");*/
  /* filter out outliers */
  /* num_usable = filterCandidates(&(wbase->lnd), &(wwarp->lnd), match_p, num_total);*/

  /*for(i=0; i<num_total; i++)
    if(fabs(match_p[i][4]-1)<0.0000001)
      printf("%6.1f %6.1f\t%6.1f %6.1f\n", match_p[i][0], match_p[i][1], match_p[i][2], match_p[i][3]);
  */

  printf("\n- verifying tie points");

  checking_passed = checkingTPs(&(wbase->lnd), &(wwarp->lnd), &(wwarp->lnd), match_p, num_total);

  if(checking_passed)
    printf("\n\t\tmatching test passed");
  else
    printf("\n\t\tmatching test failed, use with CAUTION!!!");

  free(px);
  free_2dim_contig((void **) match_p);
  
  return SUCCESS;
 
}


/**
 * do orthorectification and image match verification (checking)
 */
int initialOrthoCheck(BASE_LANDSAT  *wbase, WARP_LANDSAT  *wwarp, OUT_LANDSAT *wout, DEM *sdem)
{
  int    i, j, k, iband, npars, num_total;
  int    num_redo, checking_passed, cps_passed;
  int    current_row=0, percent_done;
  float  **match_p;
  double ai, aj, wi, wj;
  int16  **outRow, ***inRows;       /* hold on input and output lines */
  PIXEL  *px;

  if(!(px = malloc(sizeof(PIXEL)))) return FAILURE;

  alloc_2dim_contig((void ***) (&match_p), CP_KEYS.MAX_CPS, 7, sizeof(float));
  /* allocate memory to save one row of output */
  alloc_2dim_contig((void ***) (&outRow), wout->lnd.ncols, wout->nbands, sizeof(int16));
  /* allocate memory to save a block of warp input */
  alloc_3dim_contig((void ****) (&inRows), BLOCK_NROWS, wwarp->lnd.ncols, wwarp->nbands, sizeof(int16));

  num_redo = 0;
  cps_passed = 1; 
  checking_passed = 1;
  px->res = wout->lnd.res;
  px->ulx = wout->lnd.ulx;
  px->uly = wout->lnd.uly;
  npars = (wout->lnd.porder+1)*(wout->lnd.porder+2)/2;

  do {    

    fprintf(stderr, "\ttrying orthorectification (in percent) --");
 
   /* foreach row in the orthorectified image space (output) */
    for(i=0; i<wout->lnd.nrows; i++) {

      if(wout->lnd.nrows % (wout->lnd.nrows/100)) {
	percent_done = (int)((i+1)*100.0/wout->lnd.nrows+0.5);
	fprintf(stderr, "%3d\b\b\b", percent_done);
      }
      if(i == wout->lnd.nrows-1) fprintf(stderr, "100");
      /* use center of pixel */
      px->irow = i + 0.5;
      wi = px->irow;

      /* read a corresponding row from SRTM DEM file (resolution can be different) */ 
      if(readHeightLine(px, sdem)==FAILURE) exit(1);

      /* foreach pixel in the orthorectified image space (output) */
      for(j=0; j<wout->lnd.ncols; j++) {

	/* initialize output arrays */
	for(iband=0; iband<wout->nbands; iband++)
	  outRow[j][iband] = wout->lnd.fillValue;
	/* use center of pixel */
	px->icol = j + 0.5;
	wj = px->icol;

	/* find DEM data */
	k=(int)((wout->lnd.ulx + px->icol*wout->lnd.res - sdem->ulx)/sdem->res);
	if(k < 0)             k = 0;              /* if out of left boundary */ 
	if(k >= sdem->ncols)  k = sdem->ncols-1;  /* if out of right boundary */
	sdem->icol = k;

	/* convert output location (i,j) to master image location (ai, aj) */
	ai = (wbase->lnd.uly - wout->lnd.uly + wi * wout->lnd.res)/wbase->lnd.res;
	aj = (wout->lnd.ulx - wbase->lnd.ulx + wj * wout->lnd.res)/wbase->lnd.res;

	/* adjust location USING the first order transformation 
	   between previous orthorectified image and master image (ai, aj) */
	if(wout->lnd.porder == 1) { 
	  px->irow = wout->lnd.b[0]+wout->lnd.b[1]*aj+wout->lnd.b[2]*ai;
	  px->icol = wout->lnd.a[0]+wout->lnd.a[1]*aj+wout->lnd.a[2]*ai;
	}
	else {
	  px->irow = wout->lnd.b[0]+wout->lnd.b[1]*aj+wout->lnd.b[2]*ai+
	    wout->lnd.b[3]*aj*aj+wout->lnd.b[4]*aj*ai+wout->lnd.b[5]*ai*ai;
	  px->icol = wout->lnd.a[0]+wout->lnd.a[1]*aj+wout->lnd.a[2]*ai+
	    wout->lnd.a[3]*aj*aj+wout->lnd.a[4]*aj*ai+wout->lnd.a[5]*ai*ai;
	}

	/* process first pixel and load data block for next uses for this row */
	if(j == 0) {

	  /* find next valid height value */
	  while(sdem->buf[k] == DEM_FILL) k++;
	  px->height = sdem->buf[k];

	  /* compute distance from nadir view */
	  computeOffDis(&(wout->lnd), px);	

	  /* compute displacement caused by elevation */
	  computeDisplacement(&(wout->lnd), px);

	  /* convert center of pixel to un-orthorectified space  */
	  px->ux = wout->lnd.ulx + px->icol*wout->lnd.res + px->dx;
	  px->uy = wout->lnd.uly - px->irow*wout->lnd.res + px->dy;
	  px->urow = (int)((wwarp->lnd.uly - px->uy)/wwarp->lnd.res);
	  px->ucol = (int)((px->ux - wwarp->lnd.ulx)/wwarp->lnd.res);

	  /* load data block into memory */
	  current_row = px->urow;
	  if(loadDataBlock(wwarp, current_row, inRows)==FAILURE) {
	    fprintf(stderr,"load data block error!\n");
	    exit(1);
	  }	

	  /* do resample and get value */
	  resampleData(px, wwarp, current_row, inRows, outRow[j]); 

	} 
	
	else {
	  if(sdem->buf[k] != DEM_FILL || px->height != DEM_FILL) {

	    /* if current elevation data is valid then use it */
	    if(sdem->buf[k] != DEM_FILL)
	      px->height = sdem->buf[k];
	    /* else use previous height value (in px->height) */
	  
	    /* compute distance from nadir view */
	    computeOffDis(&(wout->lnd), px);

	    /* compute displacement caused by elevation */
	    computeDisplacement(&(wout->lnd), px);          

	    /* convert center of pixel to un-orthorectified space  */
	    px->ux = wout->lnd.ulx + px->icol*wout->lnd.res + px->dx;
	    px->uy = wout->lnd.uly - px->irow*wout->lnd.res + px->dy;
	    px->urow = (int)((wwarp->lnd.uly - px->uy)/wwarp->lnd.res);
	    px->ucol = (int)((px->ux - wwarp->lnd.ulx)/wwarp->lnd.res);

	    /* do resample and get value */
	    resampleData(px, wwarp, current_row, inRows, outRow[j]); 
	  }  /* end if for valid elevation */
	
#ifdef DEBUG
	  if(i==DEBUG_irow && j==DEBUG_icol) {
	    printf("\nh=%f ndis=%8.1f ds=%6.1f dx=%6.1f dy=%6.1f\n", 
		   px->height, px->off_nadir_dis, px->ds, px->dx, px->dy);	
	    printf("base_If(%6.1f, %6.1f) out_loc(%8.1f, %8.1f) UL(%8.1f, %8.1f) res=%5.1f\n", 
		   px->icol, px->irow, px->ux, px->uy, wwarp->lnd.ulx, wwarp->lnd.uly, wwarp->lnd.res);	
	    printf("Of(%4d,%4d) Dem(%4d,%4d) If(%4d,%4d) h=%f\n", 
		   j, i, sdem->icol, sdem->irow, px->ucol, px->urow, px->height);
	    for(k=0; k<wwarp->nbands; k++)
	      printf("band%d = %d\n", k+1, outRow[j][k]);
	  }
#endif      
	} /* end of else (column from 1 to ncols) */

      }   /* end of column looping */
      
      /* write one row to output file */
      if(writeOneRow(wout, outRow)==FAILURE) {
	fprintf(stderr, "Write line error iline=%d !\n", i);
	return FAILURE;
      }
    } /* end of line looping */

    /* close all output files */
    for(iband=0; iband<wout->nbands; iband++)
      fclose(wout->fp[iband]);

    /*** VALIDATE GEOLOCATION (THIS TIME TRY HIGHER ORDER POLYNOMIAL) ***/    
    /* check orthorectified output with base image */    
    printf("\n\n- validating geolocation between orthorectified output and base image\n"); 
    
    /* open output matching file for read */
    if((wout->lnd.fp = fopen(wout->lnd.fileName, "rb"))==NULL) {
      fprintf(stderr, "Can't open matching binary file %s\n", wout->lnd.fileName);
      return FAILURE;
    }

    fprintf(stderr, "\tsearching control points between BASE and OUT image --  ");
#ifdef REG_DEBUG
    printf("\n\n=== Control Point Candidates ===\n");
    printf("Base_x,y\tProjected_x,y\tMatched_x,y\tSub_Pixel_Match\t\tCorr.\tNum_High_Corr.\n");
#endif
    /* search control points based on output and base (master) image */
    num_total = searchCPs(&(wbase->lnd), &(wout->lnd), match_p); 

    printf("\n\tchecking if higher order geolocation transformation is needed");
    /* check to see if control points are match */
    /*checking_passed = checkingCPs(&(wbase->lnd), &(wout->lnd), match_p, num_total);*/
    checking_passed = checkingTPs(&(wbase->lnd), &(wout->lnd), &(wwarp->lnd), match_p, num_total);

    /* if not match then we need redo ortho with new higher order transformation */
    if(!checking_passed && num_redo<=CP_KEYS.MAX_NUM_ITER-1) { 
      printf("\n\n- checking failed, redo orthorectification along with higher order polynomial transformation");

      /* use all control points redo precise registration */    
      for(k=0; k<num_total; k++) 
	match_p[k][4] = 1;
      
      /* number of control points must be big enough to make a meaningful higher order transformation */
      if(num_total>npars+4) {  
	redoPreciseGeoCorrection(match_p, num_total, &(wout->lnd));
	num_redo++;
	cps_passed = 1;
#ifdef REG_DEBUG
	printf("\n");
	for(k=0; k<npars; k++)
	  printf("%f ", wout->lnd.a[k]);
	printf("\n");
	for(k=0; k<npars; k++)
	  printf("%f ", wout->lnd.b[k]);
	printf("\n");
#endif
	close_LANDSAT_file(&(wout->lnd));
	for(i=0; i<wout->nbands; i++) {
	  if((wout->fp[i]=fopen(wout->fileName[i], "w+b"))==NULL) {
	    printf("Can't open output file %s for reading\n", wout->fileName[i]);
	    exit(1);
	  }
	}
      }
      else 
	cps_passed = 0;      
    }
    else {
      if(checking_passed)
	printf("\n\tmatching test passed, no need to redo precise registration");
      else
	printf("\n\tmatching test failed, use orthorectified image with CAUTION!!!");
      break;
    }
    /* only iterate process while
       process didn't pass the control points checking &
       number of iteration is less than the limited numbers of redo processing & 
       there are enough number of control points are available */
    /*** END OF GEOLOCATION VALIDATION ***/

  } while (!checking_passed && num_redo<=CP_KEYS.MAX_NUM_ITER && cps_passed);

  wout->num_redo = num_redo;
  wout->checking_passed = checking_passed;
  wout->cps_passed = cps_passed;

  free(px);
  free_2dim_contig((void **) match_p);
  free_2dim_contig((void **) outRow);
  free_3dim_contig((void ***) inRows);

  if(cps_passed || checking_passed) 
    return SUCCESS;
  else 
    return FAILURE;
}


/**
 *  do final orthorectification along with precise registration function 
 *  wbase - base image in working space
 *  warp - original warp image including all bands
 *  rotate_warp - warp image after rotation (matching band only)
 *  proj_warp - warp image after reprojection (matching band only)
 *  temp_warp - warp image after rotation and reprojection, 
 *              input for working space coversion (matching band only)
 *              temp_warp will be same as warp (but one band) if rotattion 
 *              and reprojection are not necessary
 *  wwarp - warp image in working space
 *  wout - output in working space
 *  out - actual output including all bands
 *  sdem - SRTM DEM data
 *
 *           --------warp     (oroginal warp image, include all bands)
 *          |         |    
 *          |         V       if need rotation (for matching band)
 *          |   -rotate_warp
 *          |  |      |    
 *          |  |      V       if need reprojection (for matching band)
 *          |  |   proj_warp
 *          |  |      |
 *          |  |      V
 *           -- --->temp_warp (only matching band)
 *                    |
 *                    V       convert to working space in same resolution
 *                  wwarp     all registration and ortho processing are based on                  
 *
 */

int finalOrtho(BASE_LANDSAT  *wbase, WARP_LANDSAT  *warp, WARP_LANDSAT *temp_warp, WARP_LANDSAT *wwarp, OUT_LANDSAT *wout, OUT_LANDSAT *out, DEM *sdem, IN_PARAMETERS *inPars)
{
  char   command[MAX_STRLEN];
  int    i, j, k, iband, tmp;
  int    current_row=0, percent_done;
  double ai, aj, wi, wj, out_scale;
  int16  **outRow, ***inRows;        /* hold on input and output lines */

  PIXEL   *px;
  LANDSAT fine_out, agg_out;         /* define output for aggregation option */

  /* allocate memory to save one row of output */
  /* choice the largest nrows since output image maybe finer or coarser than working file */
  tmp = wout->lnd.ncols * wout->lnd.res / warp->lnd.res;
  if(out->lnd.ncols > tmp) tmp = out->lnd.ncols;
  alloc_2dim_contig((void ***) (&outRow), tmp, out->nbands, sizeof(int16));
  /* allocate memory to save a block of warp input */
  alloc_3dim_contig((void ****) (&inRows), BLOCK_NROWS, warp->lnd.ncols, warp->nbands, sizeof(int16));

  if(!(px   = malloc(sizeof(PIXEL)))) return FAILURE;
  px->res = wout->lnd.res;
  px->ulx = wout->lnd.ulx;
  px->uly = wout->lnd.uly;

  if(strcasecmp(warp->resampleMethod, "AGG")==0) {
    /* use NN for fine resolution output and then aggregate to coarse output */
    copyLandsat(&agg_out, &(out->lnd));
    /* use original spatial resolution */
    out->lnd.nrows = out->lnd.nrows * out->lnd.res / warp->lnd.res;
    out->lnd.ncols = out->lnd.ncols * out->lnd.res / warp->lnd.res;
    out->lnd.res = warp->lnd.res; 
  }

  fprintf(stderr, "\tfinal ortho processing (in percent) --");

  /* in case output resolution different in working and output space */      
  out_scale = out->lnd.res / wout->lnd.res;

  /* foreach row in the orthorectified image space (output) */
  for(i=0; i<out->lnd.nrows; i++) {

    if(out->lnd.nrows % (out->lnd.nrows/100)) {
      percent_done = (int)((i+1)*100.0/out->lnd.nrows+0.5);
      fprintf(stderr, "%3d\b\b\b", percent_done);
    }
    if(i == out->lnd.nrows-1) fprintf(stderr, "100");

    /* convert center of pixel to working space */
    px->irow = (double) (wout->lnd.uly - (out->lnd.uly - (i+0.5) * out->lnd.res)) / wout->lnd.res;
    wi = px->irow;
      
    /* read a corresponding row from SRTM DEM file (resolution can be different) */ 
    if(readHeightLine(px, sdem)==FAILURE) return FAILURE;
    
    /* foreach pixel in the orthorectified image space (output) */
    for(j=0; j<out->lnd.ncols; j++) {
	
      /* initialize output arrays */
      for(iband=0; iband<out->nbands; iband++)
	outRow[j][iband] = out->lnd.fillValue;
	
      /* convert center of pixel to working space */
      px->icol = (double) ((out->lnd.ulx + (j+0.5) * out->lnd.res) - wout->lnd.ulx) / wout->lnd.res;
      wj = px->icol;
	
      /* find DEM data */
      k=(int)((wout->lnd.ulx + px->icol*wout->lnd.res - sdem->ulx)/sdem->res);
      if(k < 0)             k = 0;              /* if out of left boundary */ 
      if(k >= sdem->ncols)  k = sdem->ncols-1;  /* if out of right boundary */
      sdem->icol = k;
	
      if(inPars->option == DO_BOTH) {

	/* convert output location (i,j) to master image location (ai, aj) 
	   for polynomial adjustment */
	ai = (wbase->lnd.uly - wout->lnd.uly + wi * wout->lnd.res)/wbase->lnd.res;
	aj = (wout->lnd.ulx - wbase->lnd.ulx + wj * wout->lnd.res)/wbase->lnd.res;

	/* precise registration using polynomial transformation function 
	   between previous orthorectified image and master image (ai, aj) 
	   and then convert coordinates back to the correponding location in 
	   output image.
	   in this way, we combines precise registration and orthorectification 
	   in one processing and avoid double resampling processing 
	*/
	if(wout->lnd.porder == 1) { 
	  px->irow = wout->lnd.b[0]+wout->lnd.b[1]*aj+wout->lnd.b[2]*ai;
	  px->icol = wout->lnd.a[0]+wout->lnd.a[1]*aj+wout->lnd.a[2]*ai;
	}
	else {
	  px->irow = wout->lnd.b[0]+wout->lnd.b[1]*aj+wout->lnd.b[2]*ai+
	    wout->lnd.b[3]*aj*aj+wout->lnd.b[4]*aj*ai+wout->lnd.b[5]*ai*ai;
	  px->icol = wout->lnd.a[0]+wout->lnd.a[1]*aj+wout->lnd.a[2]*ai+
	    wout->lnd.a[3]*aj*aj+wout->lnd.a[4]*aj*ai+wout->lnd.a[5]*ai*ai;
	}
      }
      else {
	/* ORTHO only option by assuming warp image has been precise registrated */
	px->irow = wi;
	px->icol = wj;
      }

      /* process first pixel and load data block for next uses for this row */
      if(j == 0) {
	  
	/* find next valid height value */
	while(sdem->buf[k] == DEM_FILL) k++;
	px->height = sdem->buf[k];

	/* trace location backward from output to original warp image */	  
	traceLocBackward(px, wbase, warp, temp_warp, wwarp, wout, inPars);

	/* load data block into memory */
	current_row = px->urow;

	if(loadDataBlock(warp, current_row, inRows)==FAILURE) {
	  fprintf(stderr, "load data block error!\n");
	  return FAILURE;
	}	

	/* do resample and get value */
	resampleData(px, warp, current_row, inRows, outRow[j]); 
      } 
	
      else {
	if(sdem->buf[k] != DEM_FILL || px->height != DEM_FILL) {
	  
	  /* if current elevation data is valid then use it */
	  if(sdem->buf[k] != DEM_FILL)
	    px->height = sdem->buf[k];
	  /* else use previous height value (in px->height) */

	  /* trace location backward from output to original warp image */	  
	  traceLocBackward(px, wbase, warp, temp_warp, wwarp, wout, inPars);

	  /* do resample and get value */
	  resampleData(px, warp, current_row, inRows, outRow[j]);

	} /* end if for valid elevation */
	
#ifdef DEBUG
	if(i*out_scale==DEBUG_irow && j*out_scale==DEBUG_icol) {
	  printf("\nh=%f ndis=%8.1f ds=%6.1f dx=%6.1f dy=%6.1f\n", 
		 px->height, px->off_nadir_dis, px->ds, px->dx, px->dy);
	  printf("base_If(%6.1f, %6.1f) out_loc(%8.1f, %8.1f) UL(%8.1f, %8.1f) res=%5.1f\n", 
		 px->icol, px->irow, px->ux, px->uy, warp->lnd.ulx, warp->lnd.uly, warp->lnd.res);	
	  printf("Of(%4d,%4d) Dem(%4d,%4d) If(%4d,%4d) h=%f\n", 
		 j, i, sdem->icol, sdem->irow, px->ucol, px->urow, px->height);
	  for(k=0; k<out->nbands; k++)
	    printf("band%d = %d\n", k+1, outRow[j][k]);
	}
#endif      
      } /* end of else (column from 1 to ncols) */
	
    }   /* end of column looping */
      
    /* write one row to output file in output space */
    if(writeOneRow(out, outRow)==FAILURE) {
      fprintf(stderr, "Write line error iline=%d\n!", i);
      return FAILURE;
    }
  } /* end of line looping */

    /* close output files */
  for(iband=0; iband<out->nbands; iband++)
    fclose(out->fp[iband]);    

  /* do aggregation for outputs that are coarser than inputs */
  if(strcasecmp(warp->resampleMethod, "AGG")==0) {
    copyLandsat(&fine_out, &(out->lnd));
    for(iband=0; iband<warp->nbands; iband++) {
      printf("\n\t\taggregating output for %s", out->fileName[iband]);
      /* copy fine resolution file to a temporary file and open it */
      sprintf(command, "cp %s %s", out->fileName[iband], wout->lnd.fileName);
      system(command);
      strcpy(fine_out.fileName, wout->lnd.fileName);
      if((fine_out.fp=fopen(fine_out.fileName, "rb"))==NULL) {
	fprintf(stderr, "Open output file for aggregation error!\n");
	return FAILURE;
      }
      /* aggregate temporary file and rewrite output file */
      strcpy(agg_out.fileName, out->fileName[iband]);
      fine_out.satValue = warp->satValue[iband];
      fine_out.nbyte = out->nbyte[iband];
      aggregateImage(&fine_out, &agg_out);
      fclose(fine_out.fp);
    }
    /* restore original output metadata */
    copyLandsat(&(out->lnd), &agg_out);
  }
  
  /* also print output file information to screen */
  printf("\n\torthorectified output filename:");
  for(iband=0; iband<warp->nbands; iband++) 
    printf("\n\t\t%s", out->fileName[iband]);
  printf("\n\tnumber of samples = %d", out->lnd.ncols);
  printf("\n\tnumber of lines   = %d", out->lnd.nrows);
  printf("\n\tupper left X = %f", out->lnd.ulx);
  printf("\n\tupper left Y = %f", out->lnd.uly);
  printf("\n\toutput pixel resolution = %f", out->lnd.res);
  printf("\n\tUTM zone number = %d", out->utm_zone);
  if(inPars->option == DO_BOTH) {
    if(wout->num_redo>0)
      printf("\n\thigher order geolocation transformation was applied with orthorectification process");
    if(!(wout->checking_passed))
      printf("\n\tmatching test between BASE and OUT images did NOT pass");
    if(!(wout->cps_passed) && !(wout->checking_passed))
      printf("\n\tredo orthorectification was not performed due to the limited number of control points");
  }
  printf("\n\tresampling method: %s", warp->resampleMethod);

  free(px);
  free_2dim_contig((void **) outRow);
  free_3dim_contig((void ***) inRows);

  return SUCCESS;
}



/**
 *  trace pixel location backward (after registration) in original warp image
 *  reg_out -> wwarp -> temp_warp -> (proj_warp) -> (rotate_warp) -> warp   
 */
void traceLocBackward(PIXEL *px, BASE_LANDSAT  *wbase, WARP_LANDSAT  *warp, WARP_LANDSAT *temp_warp, WARP_LANDSAT *wwarp, OUT_LANDSAT *wout, IN_PARAMETERS *inPars) 
{
  double bx, by, wx, wy, rx, ry;
  double shift_x, shift_y;

  /* compute shift from preliminary and precise registration */
  shift_x = temp_warp->lnd.ulx - wwarp->lnd.ulx;
  shift_y = wwarp->lnd.uly - temp_warp->lnd.uly;
  /*printf("\n%f %f\n", shift_x, shift_y);*/

  /* compute distance from nadir view */
  computeOffDis(&(wout->lnd), px);	

  /* compute displacement caused by elevation */
  computeDisplacement(&(wout->lnd), px);
	  
  /* convert to un-orthorectified space in working warp image (to wwarp) */
  px->ux = wout->lnd.ulx + px->icol*wout->lnd.res + px->dx;
  px->uy = wout->lnd.uly - px->irow*wout->lnd.res + px->dy;

  /* adjust to coordinate before shift registration (to temp_warp) */
  px->ux += shift_x;
  px->uy -= shift_y;

  /* convert coordinate if it is reprojected (to proj_warp)*/
  if(inPars->projection_flag) {
    bx = px->ux;
    by = px->uy;
    convertCoor(wbase, &bx, &by, &(warp->proj), &wx, &wy, BACKWARD_PROJ);
    px->ux = wx;
    px->uy = wy;
    /* convert coordinate if it is rotated (to rotate_warp)*/
    if(inPars->rotation_flag) {
      wx = wx - warp->lnd.ulx;
      wy = warp->lnd.uly - wy;
      rotation(wx, wy, &rx, &ry, warp->orientationAngle);
      px->ux = warp->lnd.ulx + rx;
      px->uy = warp->lnd.uly - ry;
    }
  }
  else 
    /* convert coordinate if it is rotated */
    if(inPars->rotation_flag) {
      wx = px->ux - warp->lnd.ulx;
      wy = warp->lnd.uly - px->uy;
      rotation(wx, wy, &rx, &ry, warp->orientationAngle);
      px->ux = warp->lnd.ulx + rx;
      px->uy = warp->lnd.uly - ry;
    }
	
  /* convert working space back to the input space  (use warp) */
  px->urow = (int)((warp->lnd.uly - px->uy)/warp->lnd.res);
  px->ucol = (int)((px->ux - warp->lnd.ulx)/warp->lnd.res);

}


/** 
 *  initialize first order polynomial coefficients 
 *  so coordinate from output can be converted to base image with these parameters  
 *  only scale+offset are allowed for first ortho iteration 
 */
void initializePoly(LANDSAT *base, LANDSAT *out)
{
  int i, npars;

  npars = (out->porder+1)*(out->porder+2)/2;
  out->a[0] = (base->ulx - out->ulx) / out->res;   
  out->a[1] = base->res / out->res;
  out->a[2] = 0.0; 
  for(i=3; i<npars; i++) out->a[i] = 0.0;
  out->b[0] = (out->uly - base->uly) / out->res;   
  out->b[1] = 0.0;
  out->b[2] = base->res / out->res;
  for(i=3; i<npars; i++) out->b[i] = 0.0;
}


/**
 * Search for control point condidates based on cross correlation between 
 * base and warp Landsat images
 * return the number of total candidates 
 */
int searchCPs(LANDSAT *base, LANDSAT *warp, float **match_p) 
{
  int i,j,m,n,index,tn,num;
  int irow, icol;
  int VAR_SHIFT;

  uint8 *buf;
  uint8 **bdata, **sdata;
  uint8 **base_image, **warp_image;
  int **r_bdata, **r_sdata;
 
  int32 start[2];
  int32 length[2];
  /*long int offset;*/

  int **control_p; /* index 0 for x and 1 for y */
  int num_total;
  
  int shift_x, shift_y, total_shift;
  int sul_x, sul_y;
  int match_x=0, match_y=0;
  int srow, erow, scol, ecol, valid_flag;
  float maxc_loc[2];

  short int B=0, S=0;
  float sum_base, sum_slave, sum2_base;
  float ave_base, ave_slave, stdev, max_stdev;
  float sum_BS, sum_BB, sum_SS;

  int num_high_correlation[4];  /* 0:>0.5; 1:>0.6; 2:>0.7; 3:>0.8 */
  float correlation, max_correlation;
  float **neighbor_cc;

  alloc_1dim_contig ((void **) (&buf), MAX_NCOLS, sizeof(uint8));
  alloc_2dim_contig((void ***) (&bdata), CP_KEYS.CHIP_SIZE, CP_KEYS.CHIP_SIZE, sizeof(uint8));
  alloc_2dim_contig((void ***) (&sdata), CP_KEYS.MAX_SIZE, CP_KEYS.MAX_SIZE, sizeof(uint8));
  alloc_2dim_contig((void ***) (&r_bdata), CP_KEYS.CHIP_SIZE, CP_KEYS.CHIP_SIZE, sizeof(int));
  alloc_2dim_contig((void ***) (&r_sdata), CP_KEYS.CHIP_SIZE, CP_KEYS.CHIP_SIZE, sizeof(int));
  alloc_2dim_contig((void ***) (&control_p), CP_KEYS.MAX_CPS, 2, sizeof(int));
  alloc_2dim_contig((void ***) (&neighbor_cc), 2*CP_KEYS.MAX_SHIFT, 2*CP_KEYS.MAX_SHIFT, sizeof(float));
  alloc_2dim_contig((void ***) (&base_image), base->nrows, base->ncols, sizeof(uint8));
  alloc_2dim_contig((void ***) (&warp_image), warp->nrows, warp->ncols, sizeof(uint8));

  /* set up control point chips (up-left corner) uniformly in base map */
  index=0;
  for(i=0; i<base->nrows; i+=CP_KEYS.CP_SEED_WIN)
    for(j=0; j<base->ncols; j+=CP_KEYS.CP_SEED_WIN) {
      control_p[index][0] = j + CP_KEYS.CP_SEED_WIN/2;   /* X coordinate */
      control_p[index][1] = i + CP_KEYS.CP_SEED_WIN/2;   /* Y coordinate */ 
      index++;
    }

  num_total = 0;   /* number of matching control points */
  VAR_SHIFT = CP_KEYS.MAX_SHIFT;  /* VAR_SHIFT will be adjusted later to accelerate searching process */
  total_shift = 0.0;

  /* load whole matching image to memory to reduce I/O operations in the searching processing */
  if(strcasecmp(warp->fileType, "BINARY")==0) rewind(warp->fp);
  for(irow=0; irow<warp->nrows; irow++) {
    if(strcasecmp(warp->fileType, "GEOTIFF")==0) {
      /* read data array from GeoTiff file */
      if (!TIFFReadScanline(warp->fp_tiff, buf, irow, 0)) {
	fprintf(stderr, "Read line %d error\n", irow);
	return FAILURE;
      } 
    }
    else 
      /* read data array from binary file */
      fread(buf, sizeof(uint8), warp->ncols, warp->fp);
    for(icol=0; icol<warp->ncols; icol++)
      warp_image[irow][icol] = buf[icol];
  }

  if(strcasecmp(base->fileType, "BINARY")==0) rewind(base->fp);
  for(irow=0; irow<base->nrows; irow++) {
    if(strcasecmp(base->fileType, "GEOTIFF")==0) {
      /* read data array from GeoTiff file */
      if (!TIFFReadScanline(base->fp_tiff, buf, irow, 0)) {
	fprintf(stderr, "Read line %d error\n", irow);
	return FAILURE;
      } 
    }
    else 
      /* read data array from binary file */
      fread(buf, sizeof(uint8), base->ncols, base->fp);
    for(icol=0; icol<base->ncols; icol++)
      base_image[irow][icol] = buf[icol];
  }

  /* find nearby "best" initial tie points with maximum stdev from base image (01/2010) */
  num = 0;

  for(i=0; i<index; i++) {
    srow = control_p[i][1] - CP_KEYS.CP_SEED_WIN/2;
    scol = control_p[i][0] - CP_KEYS.CP_SEED_WIN/2;
    erow = control_p[i][1] + CP_KEYS.CP_SEED_WIN/2;
    ecol = control_p[i][0] + CP_KEYS.CP_SEED_WIN/2;
    length[0] = CP_KEYS.CHIP_SIZE;
    length[1] = CP_KEYS.CHIP_SIZE;
    max_stdev = 0.0;
    for(m=srow; m<erow; m++)
      for(n=scol; n<ecol; n++) {
	sum_base = 0.0;
	sum2_base = 0.0;
	valid_flag = 1;
	for(irow=m; irow<m+length[0]; irow++) 
	  for(icol=n; icol<n+length[1]; icol++) {
	    if(irow<=0||irow>=base->nrows||icol<0||icol>=base->ncols) 
	      {valid_flag = 0; break;}  
	    B = base_image[irow][icol];
	    if(B == base->fillValue || B <= 0 || B == base->satValue) 
	      {valid_flag = 0; break;}
	    sum_base += B; 
	    sum2_base += B*B;
	  }
	/* if input is not a valid number then do next point */
	if(valid_flag == 0)
	  continue;
	else {
	  tn = CP_KEYS.CHIP_SIZE*CP_KEYS.CHIP_SIZE;
	  ave_base = sum_base/(float)(tn);
	  stdev = sqrt(sum2_base/tn-ave_base*ave_base);
	}

	if(stdev > max_stdev) {
	  max_stdev = stdev;
	  control_p[num][1] = m;
	  control_p[num][0] = n;
	}	
      }

    /* max stdev must be valid and significant enough */
    if(max_stdev >= CP_KEYS.MIN_STDEV) {
      /*printf("\n%d %d %d %d",scol, srow, ecol, erow);
	printf("\nnum:%d (%d %d) -> (%d %d) %f", num, control_p[i][0], control_p[i][1], control_p[num][0], control_p[num][1], max_stdev);*/
      num++;
      } 
  }

  /* start to search candidates from warp image for each control point chip */
  for(i=0; i<num; i++) {

    start[0] = control_p[i][1];
    start[1] = control_p[i][0];
    length[0] = CP_KEYS.CHIP_SIZE;
    length[1] = CP_KEYS.CHIP_SIZE;

    /* if loading sub-image is out of range, then do next */
    if(start[0]<0||start[0]>=base->nrows-length[0]||
       start[1]<0||start[1]>=base->ncols-length[1])
      continue;

    /* predict possible matching location in slave image based on their current corridinates*/
    shift_x = (base->ulx - warp->ulx + control_p[i][0]* base->res) / warp->res;
    shift_y = (warp->uly - base->uly + control_p[i][1]* base->res) / warp->res;
    
#ifdef REG_DEBUG     
    /*    printf("%5d: (%4d,%4d)\t", i, control_p[i][0], control_p[i][1]);
	  printf("(%4d,%4d)\n", shift_x, shift_y); */
#endif

    if(shift_y<0||shift_y>=warp->nrows-length[0]||
       shift_x<0||shift_x>=warp->ncols-length[1])
      continue;

    for(irow=start[0]; irow<start[0]+length[0]; irow++) {

      /* if(strcasecmp(base->fileType, "GEOTIFF")==0) {
	if (!TIFFReadScanline(base->fp_tiff, buf, irow, 0)) {
	  fprintf(stderr, "Read line %d error\n", irow);
	  return FAILURE;
	} 
      }
      else {
	offset = (long) irow * base->ncols; 
	fseek(base->fp, offset, 0);
	fread(buf, sizeof(uint8), base->ncols, base->fp);
	}*/

      for(icol=start[1]; icol<start[1]+length[1]; icol++) 
	bdata[irow-start[0]][icol-start[1]] = base_image[irow][icol];
    }    

    /* load slave image around (shift_y, shift_x) with expanded range of +-VAR_SHIFT */
    start[0] = shift_y - VAR_SHIFT;
    start[1] = shift_x - VAR_SHIFT;
    length[0] = 2 * VAR_SHIFT + CP_KEYS.CHIP_SIZE;
    length[1] = 2 * VAR_SHIFT + CP_KEYS.CHIP_SIZE;

    if(start[0]<0||start[0]>=warp->nrows-length[0]||
       start[1]<0||start[1]>=warp->ncols-length[1])
      continue;

    for(irow=start[0]; irow<start[0]+length[0]; irow++) {

      /*if(strcasecmp(warp->fileType, "GEOTIFF")==0) {
	if (!TIFFReadScanline(warp->fp_tiff, buf, irow, 0)) {
	  fprintf(stderr, "Read line %d error\n", irow);
	  return FAILURE;
	} 
      }
      else {
	offset = (long) irow * warp->ncols; 
	fseek(warp->fp, offset, 0);
	fread(buf, sizeof(uint8), warp->ncols, warp->fp);
	}*/

      for(icol=start[1]; icol<start[1]+length[1]; icol++) 
	sdata[irow-start[0]][icol-start[1]] = warp_image[irow][icol];
    }

    /* compute average value of a control point chip in base image */
    sum_base = 0.0;
    valid_flag = 1;
     for(irow=0; irow<CP_KEYS.CHIP_SIZE; irow++)
      for(icol=0; icol<CP_KEYS.CHIP_SIZE; icol++) {
	B = bdata[irow][icol];
	/* discard this point if there is a fill value within control point chip */
	if(B == base->fillValue || B <= 0 || B == base->satValue) {valid_flag=0; break;}
	else sum_base += B;
      }
     /* if input is not a valid number then do next point */
     if(valid_flag ==0)
       continue;
     else {
       tn = CP_KEYS.CHIP_SIZE*CP_KEYS.CHIP_SIZE;
       ave_base = sum_base/(float)(tn);
     }
 
     for(irow=0; irow<CP_KEYS.CHIP_SIZE; irow++)
       for(icol=0; icol<CP_KEYS.CHIP_SIZE; icol++) 
	 r_bdata[irow][icol] = bdata[irow][icol] - ave_base;
     
     for(j=0; j<4; j++) 
       num_high_correlation[j] = 0;
     max_correlation = -1.0;

     /* search matching point in a loaded slave subregion (-VAR_SHIFT, VAR_SHIFT) */
     for(m=-VAR_SHIFT; m<VAR_SHIFT; m++)
       for(n=-VAR_SHIFT; n<VAR_SHIFT; n++) {
	 
	 sul_x = n + VAR_SHIFT;
	 sul_y = m + VAR_SHIFT;

	 /* if out of slave image range, then do next searching */
	 if(shift_x+n<0||shift_y+m<0||shift_x+n>=warp->ncols||shift_y+m>=warp->nrows)
	   continue;

	/* compute average value of a candidate control point chip in slave image */
	sum_slave = 0.0;
	for(irow=0; irow<CP_KEYS.CHIP_SIZE; irow++)
	  for(icol=0; icol<CP_KEYS.CHIP_SIZE; icol++) { 
	    /* discard this point if there is a fill value within control point chip */
	    S = sdata[irow+sul_y][icol+sul_x];
	    if(S == warp->fillValue || S<=0 || S==warp->satValue) break;
	    else sum_slave += S;
	  }
	if(S == warp->fillValue || S<=0 || S==warp->satValue) break;  /* discard search for this position */
	else ave_slave = sum_slave/(CP_KEYS.CHIP_SIZE*CP_KEYS.CHIP_SIZE);

	for(irow=0; irow<CP_KEYS.CHIP_SIZE; irow++)
	  for(icol=0; icol<CP_KEYS.CHIP_SIZE; icol++)
	    r_sdata[irow][icol] = sdata[irow+sul_y][icol+sul_x] - ave_slave;

	/* compute normalized cross correlation */
	sum_BS = 0.0;
	sum_BB = 0.0;
	sum_SS = 0.0;
	for(irow=0; irow<CP_KEYS.CHIP_SIZE; irow++)
	  for(icol=0; icol<CP_KEYS.CHIP_SIZE; icol++) {
	    B = r_bdata[irow][icol];
	    S = r_sdata[irow][icol];
	    sum_BS += B*S;
	    sum_BB += B*B;
	    sum_SS += S*S;
	  }
	if(sum_BB == 0 || sum_SS == 0) break;  /* discard searching if invalid */ 
	correlation = sum_BS / (sqrt(sum_BB) * sqrt(sum_SS));
	/* save correlation to array for later use in sub-pixel matching */ 
	neighbor_cc[sul_y][sul_x] = correlation;

	if(correlation > max_correlation) {
	  match_x = shift_x + n;
	  match_y = shift_y + m;
	  max_correlation = correlation;
	  maxc_loc[0] = sul_x;
	  maxc_loc[1] = sul_y;
	}

	/* count numbers of high correlation */
	if(correlation > 0.50) num_high_correlation[0]++;
	if(correlation > 0.60) num_high_correlation[1]++;
	if(correlation > 0.70) num_high_correlation[2]++;
	if(correlation > 0.80) num_high_correlation[3]++;
 
      }  /* end of searching from slave image */

    /** acceptable control point must:
     *  1) has a high correlation larger than ACCEPTABLE_CORR
     *  2) has a small number of high correlation less than MAX_NUM_HIGH_CORR
     */

    j = (int)(max_correlation*10) - 6;  /* index for num_high_correlation[] */
    if(j<0) j = 0;
    if(j>3) j = 3;    
    
    if(max_correlation >= CP_KEYS.ACCEPTABLE_CORR && num_high_correlation[j] <= CP_KEYS.MAX_NUM_HIGH_CORR) {

      /* save control points to a array */
      match_p[num_total][0] = control_p[i][0];
      match_p[num_total][1] = control_p[i][1];
      /* use sub-pixel searching function */
      seekSubPixelLoc(neighbor_cc, maxc_loc);
      match_p[num_total][2] = match_x + maxc_loc[0];
      match_p[num_total][3] = match_y + maxc_loc[1];
      match_p[num_total][4] = 1;  /* 1=usable; 0=invalid */
      num_total++;

#ifdef REG_DEBUG     
      printf("(%4.0f,%4.0f)\t", match_p[num_total-1][0], match_p[num_total-1][1]);
      printf("(%4d,%4d)\t", shift_x, shift_y); 
      printf("(%4d,%4d)\t", match_x, match_y);
      printf("(%6.1f,%6.1f)\t", match_p[num_total-1][2], match_p[num_total-1][3]);
      printf("\t%5.3f\t", max_correlation);
      printf("%d (>%3.1f)\n", num_high_correlation[j], 0.5+j*0.1);
#endif

      total_shift += sqrt((match_x-shift_x)*(match_x-shift_x)+(match_y-shift_y)*(match_y-shift_y));
	
      fprintf(stderr, "%4d\b\b\b\b", num_total);

      /* estimate offset and adjust VAR_SHIFT to accelerate searching process 
	 after 10 GCPs are found (for reliable estimation) */
      if(num_total == 10) {
	VAR_SHIFT = (int)total_shift/num_total+5;
	if(VAR_SHIFT > CP_KEYS.MAX_SHIFT) VAR_SHIFT = CP_KEYS.MAX_SHIFT; 
#ifdef REG_DEBUG
	printf("Landsat CP searching distance have been adjusted from %d to %d to accelerate process\n", CP_KEYS.MAX_SHIFT, VAR_SHIFT); 
#endif
      }
	
      /*  printf("---base image---\n");
      for(irow=0; irow<CP_KEYS.CHIP_SIZE; irow++) {
	for(icol=0; icol<CP_KEYS.CHIP_SIZE; icol++) 
	  printf("%d ",bdata[irow][icol]); 
	printf("\n");
      }

      printf("---shift slave---\n");
      for(irow=0; irow<CP_KEYS.CHIP_SIZE; irow++) {
	for(icol=0; icol<CP_KEYS.CHIP_SIZE; icol++) 
	  printf("%d ",sdata[irow+VAR_SHIFT][icol+VAR_SHIFT]); 
	printf("\n");
      }

      printf("---match slave---\n");
      for(irow=match_y-shift_y+VAR_SHIFT; irow<CP_KEYS.CHIP_SIZE+match_y-shift_y+VAR_SHIFT; irow++) {
	for(icol=match_x-shift_x+VAR_SHIFT; icol<CP_KEYS.CHIP_SIZE+match_x-shift_x+VAR_SHIFT; icol++) 
	  printf("%d ",sdata[irow][icol]); 
	printf("\n");
	}*/
      
    }    
 
  } /* end of searching of control point candidate */

  fprintf(stderr, "%4d", num_total);
  
  free (buf);
  free_2dim_contig((void **) bdata);
  free_2dim_contig((void **) sdata);
  free_2dim_contig((void **) r_bdata);
  free_2dim_contig((void **) r_sdata);
  free_2dim_contig((void **) control_p);
  free_2dim_contig((void **) neighbor_cc);
  free_2dim_contig((void **) base_image);
  free_2dim_contig((void **) warp_image);

  return num_total;
}



/**
 * compute displacement caused by terrain effects
 * those displacments should be considered while matching orthorectified and 
 * un-orthorectified Landsat images
 */
int computeOrthoShifts(LANDSAT *base, DEM *sdem, float **match_p, int num_total)
{
  int i, k;
  PIXEL *px;
  
  px = malloc(sizeof(PIXEL));
  px->res = base->res;
  px->ulx = base->ulx;
  px->uly = base->uly;
  for(i=0; i<num_total; i++) {
    px->irow = match_p[i][1];
    px->icol = match_p[i][0];
    if(readHeightLine(px, sdem)==FAILURE) return FAILURE;
    k=(int)((base->ulx + px->icol * base->res - sdem->ulx)/sdem->res);
    if(k < 0)             k = 0;              /* if out of left boundary */ 
    if(k >= sdem->ncols)  k = sdem->ncols-1;  /* if out of right boundary */
    sdem->icol = k;
    if(sdem->buf[k] == DEM_FILL) {
      match_p[i][4] = 0;
      continue;
    }
    px->height = sdem->buf[k];
 
    /* compute distance from nadir view */
    computeOffDis(base, px);
    
    /* compute displacement caused by elevation */
    computeDisplacement(base, px);

    match_p[i][5] = px->dx;
    match_p[i][6] = px->dy;
 
    /*  printf("Base(%f %f), Dem(%d %d) height=%f\n",px->icol, px->irow, sdem->icol, sdem->irow, px->height);
      printf("(%6.1f,%6.1f)\t", match_p[i][0], match_p[i][1]);
      printf("(%6.1f,%6.1f)\t", match_p[i][2], match_p[i][3]);
      printf("(%6.1f,%6.1f)\t", match_p[i][5], match_p[i][6]);
      getchar();
    */
    
  }
  free(px);
  return SUCCESS;
}


/**
 * filter control point candidate 
 * remove one control point each time until maximum error is less than 0.5 pixel
 * return the number of final selection
 */
int filterCandidates(LANDSAT *base, LANDSAT *slave, float **match_p, int num_total) 
{
  int index, i, j=0;
  int num_usable;

  float shift_x, shift_y;  
  float sum_shift_x, sum_shift_y;
  float ave_shift_x, ave_shift_y;
  float predict_x, predict_y;
  float dis_x, dis_y;
  float dis, max_dis;


  /* do loop check on control point candidates, filter out outliers */  
  do {
    index = 0;
    sum_shift_x = 0;
    sum_shift_y = 0;
    for(i=0; i<num_total; i++) 
      if(fabs(match_p[i][4]-1)<0.0001) {
	shift_x = match_p[i][0]*base->res + match_p[i][5] - match_p[i][2]*slave->res;
	shift_y = match_p[i][1]*base->res - match_p[i][6] - match_p[i][3]*slave->res;
	sum_shift_x += shift_x;
	sum_shift_y += shift_y;
	index++;
      }
    ave_shift_x = (float)sum_shift_x/index;
    ave_shift_y = (float)sum_shift_y/index;
    num_usable = index;
    max_dis = -1.0;
    for(i=0; i<num_total; i++)
      if(fabs(match_p[i][4]-1)<0.0001) {
	/* compute error for each control point */
	predict_x = match_p[i][0]*base->res + match_p[i][5] - ave_shift_x;
	predict_y = match_p[i][1]*base->res - match_p[i][6] - ave_shift_y;
	dis_x = fabs(predict_x - match_p[i][2]*slave->res);
	dis_y = fabs(predict_y - match_p[i][3]*slave->res);
	dis = sqrt(dis_x*dis_x+dis_y*dis_y);
	if(dis > max_dis) {
	  max_dis = dis;
	  j = i;
	}
      }
    /* if the maximum error larger than half pixel, then discard it */
    if(max_dis > CP_KEYS.MAX_AVE_ERROR * base->res) { 
      match_p[j][4] = 0;
      num_usable--;
    }
  } while (num_usable != index);

  /* compute adjusted up-left corner coordinate */
  slave->ulx = base->ulx + ave_shift_x;
  slave->uly = base->uly - ave_shift_y;

#ifdef REG_DEBUG  
  printf("\nnum_tps=%d\tp_shift_x: %4.1f\tp_shift_y: %4.1f\tmax_err: %4.2f (in meters)\n", 
	 num_usable, ave_shift_x, ave_shift_y, max_dis);
#endif

  return num_usable;
}



/**
 * checks control point candidates and returns if two matches or not
 */
int checkingCPs(LANDSAT *base, LANDSAT *slave, float **match_p, int num_total) 
{
  int index, i, j=0, k, ret;
  int num_usable, num_cp_sec[3], num_cp_bak[3];

  float shift_x, shift_y;  
  float sum_shift_x, sum_shift_y;
  float ave_shift_x, ave_shift_y;
  float predict_x, predict_y;
  float dis_x, dis_y;
  float dis, max_dis;
  double x0, y0, a, b, c;

  /* first divide control points into three sections (left offnadir (0), nadir (1) and right offnadir(2)) */
  /* use line passes central point with nadir direction */ 
  a = slave->lpara[0];
  b = -1.0;
  c = slave->ncols/2.0 - slave->nrows/2.0 * a;
  /* initialize number of control points in each section */
  for(i=0; i<3; i++) num_cp_sec[i] = 0;
  for(i=0; i<num_total; i++) {
    /* check tie point distribution in warp image */
    x0 = (match_p[i][2]+0.5);
    y0 = (match_p[i][3]+0.5);
    /* distance from point (x0, y0) to ax+by+c=0 */
    dis = fabs(a*x0+b*y0+c)/sqrt(a*a+b*b);
    if(dis < slave->ncols/6.0) {
      match_p[i][4] = 1;
      num_cp_sec[1]++;
    }
    else {
      /* decides sign of displacement (right to nadir track use 2; left to nadir track use 0) */
      if(x0 * a + y0 * b + c > 0) {
	match_p[i][4] = 2;
	num_cp_sec[2]++;
      }
      else {
	match_p[i][4] = 0;
	num_cp_sec[0]++;
      }
    }
  }  
  for(i=0; i<3; i++)  num_cp_bak[i] = num_cp_sec[i];

  /* do loop check on control point candidates, filter out outliers */  
  do {
    index = 0;
    sum_shift_x = 0;
    sum_shift_y = 0;
    for(i=0; i<num_total; i++) 
      if(match_p[i][4] > -0.0001 ) {
	shift_x = match_p[i][0]*base->res - match_p[i][2]*slave->res;
	shift_y = match_p[i][1]*base->res - match_p[i][3]*slave->res;
	sum_shift_x += shift_x;
	sum_shift_y += shift_y;
	index++;
      }
    ave_shift_x = (float)sum_shift_x/index;
    ave_shift_y = (float)sum_shift_y/index;
    num_usable = index;
    max_dis = -1.0;
    for(i=0; i<num_total; i++)
      if(match_p[i][4] > -0.0001 ) {
	/* compute error for each control point */
	predict_x = match_p[i][0]*base->res - ave_shift_x;
	predict_y = match_p[i][1]*base->res - ave_shift_y;
	dis_x = fabs(predict_x - match_p[i][2]*slave->res);
	dis_y = fabs(predict_y - match_p[i][3]*slave->res);
	dis = sqrt(dis_x*dis_x+dis_y*dis_y);
	if(dis > max_dis) {
	  max_dis = dis;
	  j = i;
	}
      }
    /* if the maximum error larger than a pixel, then discard it */
    if(max_dis > base->res) {
      k = (int)(match_p[j][4]+0.01);
      num_cp_sec[k]--;
      match_p[j][4] = -1;
      num_usable--;
    }
  } while (num_usable != index);

  /* if half of control points in each section doesn't match, then return false */
  if((num_usable < 0.5*num_total)|| 
     (num_cp_sec[0] < 0.5*num_cp_bak[0]) ||
     (num_cp_sec[1] < 0.5*num_cp_bak[1]) ||  
     (num_cp_sec[2] < 0.5*num_cp_bak[2]))
    ret = 0;
  else
    ret = 1;

  printf("\n\t\tOVERALL   TOTAL_CPS=%4d  NUM_MATCHES=%4d", num_total, num_usable);
  printf("\n\t\tSECTION1(left)  TOTAL_CPS=%4d  NUM_MATCHES=%4d", num_cp_bak[0], num_cp_sec[0]);
  printf("\n\t\tSECTION2(nadir) TOTAL_CPS=%4d  NUM_MATCHES=%4d", num_cp_bak[1], num_cp_sec[1]);
  printf("\n\t\tSECTION3(right) TOTAL_CPS=%4d  NUM_MATCHES=%4d", num_cp_bak[2], num_cp_sec[2]);

  return ret;
}



/**
 * checks tie points in four regions (UL, UR, LR, LL) in valid output image 
 * return success if match testing passes
 */
int checkingTPs(LANDSAT *base, LANDSAT *out, LANDSAT *warp, float **match_p, int num_total) 
{
  int i, k, ret;
  int num_usable, num_cp_sec[4], num_cp_bak[4];
  int x0, y0, cx, cy, ulx, uly, lrx, lry;

  /*  int j=0, index;
      float max_dis;
      float shift_x, shift_y;  
      float sum_shift_x, sum_shift_y;
      float ave_shift_x, ave_shift_y; */
  float predict_x, predict_y;
  float dis_x, dis_y;
  float dis;

  /* compute center pixel of valid region in output space */
  ulx = (int)((warp->ulx - out->ulx) / out->res + 0.5);
  if(ulx < 0) ulx = 0;
  uly = (int)((out->uly - warp->uly) / out->res + 0.5);
  if(uly < 0) uly = 0;
  lrx = (int)((warp->ulx + warp->ncols * warp->res - out->ulx) / out->res + 0.5);
  if(lrx > out->ncols) lrx = out->ncols;
  lry = (int)((out->uly - (warp->uly - warp->nrows * warp->res)) / out->res + 0.5);
  if(lry > out->nrows) lry = out->nrows;
  cx = (ulx + lrx) / 2.0;
  cy = (uly + lry) / 2.0;


  /* initialize number of control points in each section */
  /*   0 | 1
   *   --+--
   *   2 | 3
   *  where + is (cx, cy)
   */
  for(i=0; i<4; i++) num_cp_sec[i] = 0;

  for(i=0; i<num_total; i++) {
    /* check tie point distribution in output image */
    x0 = match_p[i][2];
    y0 = match_p[i][3];

    if(x0 < cx && y0 < cy) {
      match_p[i][4] = 0;
      num_cp_sec[0]++;
    }
    else if(x0 >= cx && y0 < cy) {
      match_p[i][4] = 1;
      num_cp_sec[1]++;
    }     
    else if(x0 < cx && y0 >= cy) {
      match_p[i][4] = 2;
      num_cp_sec[2]++;
    }     
    else if(x0 >= cx && y0 >= cy) {
      match_p[i][4] = 3;
      num_cp_sec[3]++;
    }
  }  

  for(i=0; i<4; i++)  num_cp_bak[i] = num_cp_sec[i];

  /* use absolute coordinate instead of relative - Feng (10/2008) */ 
  num_usable = num_total;
  for(i=0; i<num_total; i++)
    if(match_p[i][4] > -0.0001 ) {
      /* compute error for each control point */
      predict_x = base->ulx + match_p[i][0]*base->res;
      predict_y = base->uly - match_p[i][1]*base->res;
      dis_x = fabs(predict_x - (out->ulx + match_p[i][2]*out->res));
      dis_y = fabs(predict_y - (out->uly - match_p[i][3]*out->res));
      dis = sqrt(dis_x*dis_x+dis_y*dis_y);
      /* if the error larger than a pixel, then discard it */
      if(dis > base->res) {
	k = (int)(match_p[i][4]+0.01);
	num_cp_sec[k]--;
	num_usable--;
      }
    }


  /* do loop check on control point candidates, filter out outliers */  
  /*do {
    index = 0;
    sum_shift_x = 0;
    sum_shift_y = 0;
    for(i=0; i<num_total; i++) 
      if(match_p[i][4] > -0.0001 ) {
	shift_x = match_p[i][0]*base->res - match_p[i][2]*out->res;
	shift_y = match_p[i][1]*base->res - match_p[i][3]*out->res;
	shift_x = match_p[i][0]*base->res - match_p[i][2]*out->res;
	sum_shift_x += shift_x;
	sum_shift_y += shift_y;
	index++;
      }
    ave_shift_x = (float)sum_shift_x/index;
    ave_shift_y = (float)sum_shift_y/index;
    num_usable = index;
    max_dis = -1.0;
    for(i=0; i<num_total; i++)
    if(match_p[i][4] > -0.0001 ) {*/
	/* compute error for each control point */
	/*predict_x = match_p[i][0]*base->res - ave_shift_x;
	predict_y = match_p[i][1]*base->res - ave_shift_y;
	dis_x = fabs(predict_x - match_p[i][2]*out->res);
	dis_y = fabs(predict_y - match_p[i][3]*out->res);
	dis = sqrt(dis_x*dis_x+dis_y*dis_y);
	if(dis > max_dis) {
	  max_dis = dis;
	  j = i;
	}
	}*/
    /* if the maximum error larger than a pixel, then discard it */
    /*if(max_dis > base->res) {
      k = (int)(match_p[j][4]+0.01);
      num_cp_sec[k]--;
      match_p[j][4] = -1;
      num_usable--;
    }
    } while (num_usable != index);*/

  /* if half of control points in each section doesn't match, then return false */
  if((num_usable < 0.6*num_total) && ( 
     (num_cp_sec[0] < 0.5*num_cp_bak[0]) ||
     (num_cp_sec[1] < 0.5*num_cp_bak[1]) ||  
     (num_cp_sec[2] < 0.5*num_cp_bak[2]) ||
     (num_cp_sec[3] < 0.5*num_cp_bak[3])))
    ret = 0;
  else
    ret = 1;

  printf("\n\t\tOverall (num_match / num_total) = (%4d / %4d)", num_usable, num_total);
  printf("\n\t\t(%4d/%4d) | (%4d/%4d)", num_cp_sec[0], num_cp_bak[0], num_cp_sec[1], num_cp_bak[1]);
  printf("\n\t\t------------+------------");
  printf("\n\t\t(%4d/%4d) | (%4d/%4d)", num_cp_sec[2], num_cp_bak[2], num_cp_sec[3], num_cp_bak[3]);

  return ret;
}



/**
 * redo precise registration between orthorectified image and master image (GeoCover)
 *  ! Input
 *    double match_p[][5]: control point pairs
 *    0=base_x, 1=base_y, 2=warp_x, 3=warp_y, 4=valid_flag, 5=predict_warp_x, 6=predict_warp_y
 *    int tnum: total number of control points
 *  ! Output
 *    float a[]: first order polynomial coefficient for X
 *    float b[]: first order polynomial coefficient for Y
 */

void redoPreciseGeoCorrection(float **match_p, int tnum, LANDSAT *out)
{
  int i, num, mi, npars;
  float **cp, a[6], b[6];
  float dis, rmse=0.0, max_dis, pre_rmse;
  float *x, *y;
  float tx, ty, txy, tx2, ty2;
  
  alloc_1dim_contig((void **) (&x),CP_KEYS.MAX_CPS, sizeof(float));
  alloc_1dim_contig((void **) (&y),CP_KEYS.MAX_CPS, sizeof(float));
  alloc_2dim_contig((void ***) (&cp), CP_KEYS.MAX_CPS, 3, sizeof(float));

  npars = (out->porder+1)*(out->porder+2)/2;

  do {
      num = 0;
      for(i=0; i<tnum; i++) {
	if(fabs(match_p[i][4]-1)<0.0001) {
	  cp[num][0]=match_p[i][0];
	  cp[num][1]=match_p[i][1];    
	  cp[num][2]=match_p[i][2];
	  num++;
	}
      }    
      getAffine(cp, num, out->porder, a);

      num=0;
      for(i=0; i<tnum; i++) {
	if(fabs(match_p[i][4]-1)<0.0001) {
	  cp[num][0]=match_p[i][0];
	  cp[num][1]=match_p[i][1];   
	  cp[num][2] = match_p[i][3];
	  num++;
	}
      }
      getAffine(cp, num, out->porder, b);

      pre_rmse = rmse;
      rmse = 0.0;
      max_dis = 0.0;
      mi = 0;
      for(i=0; i<tnum; i++) {
	if(fabs(match_p[i][4]-1)<0.0001) {
	  tx = match_p[i][0];
	  ty = match_p[i][1];
	  tx2 = tx * tx;
	  txy = tx * ty;
	  ty2 = ty * ty;
	  if(out->porder == 1) {
	    x[i]=a[0]+a[1]*tx+a[2]*ty;
	    y[i]=b[0]+b[1]*tx+b[2]*ty;
	  }
	  else {
	    x[i]=a[0]+a[1]*tx+a[2]*ty+a[3]*tx2+a[4]*txy+a[5]*ty2;
	    y[i]=b[0]+b[1]*tx+b[2]*ty+b[3]*tx2+b[4]*txy+b[5]*ty2;
	  }
	  match_p[i][5] = x[i];
	  match_p[i][6] = y[i];
	  dis = sqrt(pow(x[i]-match_p[i][2],2)+pow(y[i]-match_p[i][3],2));
	  if(dis>max_dis) {max_dis = dis; mi = i;}
	  rmse += dis*dis;
	  /*printf("(%7.2f,%7.2f) (%7.2f,%7.2f) %7.2f (%7.2f,%7.2f) %7.4f\n", match_p[i][0], match_p[i][1], match_p[i][2], match_p[i][3], sqrt(pow(match_p[i][0]-match_p[i][2],2)+pow(match_p[i][1]-match_p[i][3],2)), x[i], y[i], dis);*/
	} 
      }

      rmse = sqrt(rmse/num);
      match_p[mi][4] = 0;

#ifdef REG_DEBUG
      printf("\n");
      for(i=0; i<npars; i++) 
	printf("%f ", a[i]);
      printf("\n");
      for(i=0; i<npars; i++) 
	printf("%f ", b[i]);
      printf("\n");
      printf("%d: (%7.2f,%7.2f) (%7.2f,%7.2f) (%7.2f,%7.2f) %7.4f\n", mi, match_p[mi][0], match_p[mi][1], match_p[mi][2], match_p[mi][3], x[mi], y[mi], max_dis);
      printf("RMSE = %f  NUM = %d\n", rmse, num);
#endif

      /* add more limitations on the tie points, continue eliminating "bad" tie points 
	 if one of following conditions is not satisfied:
	 - the rmse MUST be smaller than CP_KEYS.MAX_RMSE (defined in *.ini file)
	 - the maximum error among the tie point MUST less than 2 pixels
	 - the difference of rmse between this loop and previous one MUST small than predefined value (0.01)
	 and of course, we need to maintain the minumum number of tie points for regression
	 Feng (10/2008)
      */
  } while((rmse>CP_KEYS.MAX_RMSE || max_dis>2 || fabs(pre_rmse-rmse)>CP_KEYS.DIFF_RMSE) && num>npars+2); 
  
  for(i=0; i<npars; i++) { out->a[i]=a[i]; out->b[i]=b[i]; }
 
  if(out->porder == 1) {
    printf("\n\ttrying 1st order polynomial transformation");
    printf("\n\t\tNum_CP_Selected=%d, RMSE=%6.3f", num, rmse);
    printf("\n\t\tx' = %9.3f + %6.3f * x + %6.3f * y", a[0], a[1], a[2]);
    printf("\n\t\ty' = %9.3f + %6.3f * x + %6.3f * y\n", b[0], b[1], b[2]);
  }
  else if(out->porder == 2) {
    printf("\n\ttrying 2nd order polynomial transformation"); 
    printf("\n\t\tNum_CP_Selected=%d, RMSE=%6.3f", num, rmse); 
    printf("\n\t\tx' = %9.3f + %6.3f * x + %6.3f * y + %6.3f * x^2 + %6.3f * xy + %6.3f *y^2", 
	   a[0], a[1], a[2], a[3], a[4], a[5]);
    printf("\n\t\ty' = %9.3f + %6.3f * x + %6.3f * y + %6.3f * x^2 + %6.3f * xy + %6.3f *y^2\n", 
	   b[0], b[1], b[2], b[3], b[4], b[5]);
  }

  free(x);
  free(y);
  free_2dim_contig((void **) cp);
}



/**
 * find subpixel location based on the neighbor pixel cross correlation 
 * ! Input
 *     float neighbor_cc[][]: cross correlation matrix 
 *     float maxc_loc[]: location (in the matrix) of maximum correlation
 * ! Ouput
 *     float maxc_loc[]: adjusted subpixel location [-2, 2]    
 */ 
void seekSubPixelLoc(float **neighbor_cc, float *maxc_loc)
{
  int i, j, k, m, n;
  float **cp;
  float c[6], temp;

  alloc_2dim_contig((void ***) (&cp), 25, 3, sizeof(float));
  
  k=0;
  for(i=-2; i<=2; i++) {
    /*printf("\n");*/
    for(j=-2; j<=2; j++) {
      cp[k][0] = j;
      cp[k][1] = i;
      m = i+maxc_loc[1];
      n = j+maxc_loc[0];
      if(m<0||m>=CP_KEYS.MAX_SHIFT*2||n<0||n>=CP_KEYS.MAX_SHIFT*2) {
	maxc_loc[0] = 0.0;
	maxc_loc[1] = 0.0;
	return;
      }
      cp[k][2] = neighbor_cc[m][n];
      /*printf("%6.3f ", cp[k][2]);*/
      k++;
    }
  }  

  getAffine(cp, 25, 2, c);
  /*printf("\n");
  for(i=-2; i<=2; i++) {
    printf("\n");
    for(j=-2; j<=2; j++) {
      tmp = c[0] + c[1]*j + c[2]*i + c[3]*j*j + c[4]*i*j + c[5]*i*i; 
    printf("%6.3f ", tmp);
    }
    }*/
  /* return: z =  c[0] + c[1]*x + c[2]*y + c[3]*x*x + c[4]*x*y + c[5]*y*y */
  /* or in:  z =  c[3]*x*x + c[5]*y*y +  c[4]*x*y + c[1]*x + c[2]*y +c[0] */

  temp = 4.0*c[3]*c[5]-c[4]*c[4];

  /* function f(x,y)=Ax^2+By^2+Cx+Dy+Exy+F has a maximum only if 4AB-E^2>0 and A<0 */
  if(temp > 0.0 && c[3] < 0.0) {
    maxc_loc[0] = (-2.0*c[5]*c[1]+c[4]*c[2])/temp; 
    maxc_loc[1] = (-2.0*c[3]*c[2]+c[4]*c[1])/temp;
  }
  else {
    maxc_loc[0] = 0.0;
    maxc_loc[1] = 0.0;
  }

  /* check sub-pixel adjustment and make sure they are not larger than half pixel */
  for(i=0; i<2; i++) {
    if(fabs(maxc_loc[i])>0.5) {
      if(maxc_loc[i]>0) maxc_loc[i] = 0.5;
      else maxc_loc[i] = -0.5;
    }
  }

  free_2dim_contig((void **) cp);

}



/*
 * do resmapling according to requirement (NN, BI or CC) 
 * ! Input
 *     px: pixel information including location and values 
 *     warp: input landsat metadata
 *     inRows: input landsat data block 
 *     current_row: current loaded row in the data block (middle row of inRows)
 * ! Output
 *     outRow: one output resampled row
 */
void resampleData(PIXEL *px, WARP_LANDSAT *warp, int current_row, int16 ***inRows, int16 *outRow)
{
  int arr_irow, iband, i, j, trow, tcol, tarr;
  int v[4][4][MAX_NBANDS];
  float D, dx, dy, weight[4][4], tref, tsum;

  /* return if it outside input image */
  if(px->urow<0 || px->urow >= warp->lnd.nrows || px->ucol<0 || px->ucol >= warp->lnd.ncols) 
    return;
  
  /* default use nearest neighbor sampling method */
  /* convert from input image space to data block space */ 
  arr_irow = px->urow - current_row + BLOCK_NROWS/2;
  /* get results from data block (memory) for fast access */ 
  if( arr_irow>=0 && arr_irow<BLOCK_NROWS && px->ucol>=0 && px->ucol<warp->lnd.ncols) {
    for(iband=0; iband<warp->nbands; iband++)
      outRow[iband] = inRows[arr_irow][px->ucol][iband];
    /*printf("irow=%d icol=%d in_urow=%d in_ucol=%d arr_irow=%d\n", i, j, px->urow, px->ucol, arr_irow);*/
  }
  
  /* if not in data block, then read input from given location (px->ucol, px->urow) */
  else if(readPixel(warp, px)==SUCCESS)
    {
      for(iband=0; iband<warp->nbands; iband++)
	outRow[iband] = px->ref[iband];
    }

  /*if(outRow[0]!=warp->lnd.fillValue) {
    readPixel(warp, px);
    printf("readPixel (%d,%d) %d %d\n", px->ucol, px->urow, outRow[0], px->ref[0]);
    }*/

  for(i=0; i<4; i++)
    for(j=0; j<4; j++) {
      for(iband=0; iband<warp->nbands; iband++)
	v[i][j][iband] = warp->lnd.fillValue;
    }

  /* process bilinear interpolation sampling method */
  if(strcmp(warp->resampleMethod, "BI")==0) {
    /*if(outRow[0]!=warp->lnd.fillValue)
      printf("NN (%d, %d) %d\n", px->ucol, px->urow, outRow[0]);*/
    /* compute nearest upleft pixel */
    px->urow = (int)((warp->lnd.uly - (px->uy + 0.5 * warp->lnd.res))/warp->lnd.res);
    px->ucol = (int)(((px->ux - 0.5 * warp->lnd.res) - warp->lnd.ulx)/warp->lnd.res);
    dx = (px->ux - 0.5 * warp->lnd.res) - warp->lnd.ulx - px->ucol * warp->lnd.res; 
    dy = warp->lnd.uly - (px->uy + 0.5 * warp->lnd.res) - px->urow * warp->lnd.res;
    D = warp->lnd.res;
    arr_irow = px->urow - current_row + BLOCK_NROWS/2;
    tcol = px->ucol;
    trow = px->urow;
    tarr = arr_irow;
    for(i=0; i<2; i++)
      for(j=0; j<2; j++) {
	arr_irow = tarr + i; px->ucol = tcol + j; px->urow = trow + i;
	if( arr_irow>=0 && arr_irow<BLOCK_NROWS && px->ucol>=0 && px->ucol<warp->lnd.ncols) {
	  for(iband=0; iband<warp->nbands; iband++)
	    v[i][j][iband] = inRows[arr_irow][px->ucol][iband];
	}
	else if(readPixel(warp, px)==SUCCESS) {
	  for(iband=0; iband<warp->nbands; iband++)
	    v[i][j][iband] = px->ref[iband];
	}

	for(iband=0; iband<warp->nbands; iband++)
	  if(v[i][j][iband] == warp->lnd.fillValue ||
	     v[i][j][iband] == warp->satValue[iband] ) return;  /* use NN value */

	/*printf("(%d,%d) %d\n", px->ucol, px->urow, v[i][j][0]);*/
      }
    /*printf("%f %f ", dx, dy);*/
    for(iband=0; iband<warp->nbands; iband++)
      outRow[iband] = (int)((v[0][0][iband]*(D-dx)*(D-dy) + v[0][1][iband]*dx*(D-dy) +
				  v[1][0][iband]*(D-dx)*dy + v[1][1][iband]*dx*dy)/(D*D) + 0.5);
    /*printf("%d %d %d %d %d\n", v[0][0][0], v[0][1][0],v[1][0][0],v[1][1][0],outRow[0]);*/
  }

  /* process cubic convolution sampling method */
  if(strcmp(warp->resampleMethod, "CC")==0) {
    /* if(outRow[0]!=warp->lnd.fillValue)
       printf("NN (%d, %d) %d\n", px->ucol, px->urow, outRow[0]);*/
    /* compute nearest upleft pixel */
    px->urow = (int)((warp->lnd.uly - (px->uy + 0.5 * warp->lnd.res))/warp->lnd.res);
    px->ucol = (int)(((px->ux - 0.5 * warp->lnd.res) - warp->lnd.ulx)/warp->lnd.res);
    arr_irow = px->urow - current_row + BLOCK_NROWS/2;
    tcol = px->ucol;
    trow = px->urow;
    tarr = arr_irow;
 
    for(i=0; i<4; i++)       
      for(j=0; j<4; j++) {
	arr_irow = tarr+i-1; px->ucol = tcol+j-1; px->urow = trow+i-1;
	if( arr_irow>=0 && arr_irow<BLOCK_NROWS && px->ucol>=0 && px->ucol<warp->lnd.ncols) {
	  for(iband=0; iband<warp->nbands; iband++)
	    v[i][j][iband] = inRows[arr_irow][px->ucol][iband];
	}
	else if(readPixel(warp, px)==SUCCESS) {
	  for(iband=0; iband<warp->nbands; iband++)
	    v[i][j][iband] = px->ref[iband];
	  /*printf("%d ", px->ref[0]);*/
	}

	for(iband=0; iband<warp->nbands; iband++)
	  if(v[i][j][iband] == warp->lnd.fillValue ||
	     v[i][j][iband] == warp->satValue[iband]) return;  /* use NN value */

	dx = (px->ux - 0.5 * warp->lnd.res) - warp->lnd.ulx - px->ucol * warp->lnd.res; 
	dy = warp->lnd.uly - (px->uy + 0.5 * warp->lnd.res) - px->urow * warp->lnd.res;
	/* approach as used in the MRT code */
	weight [i][j] = cubic(dx/warp->lnd.res)*cubic(dy/warp->lnd.res);
	/*printf("(%d %d) (%d,%d) %7.3f %7.3f %7.3f %d\n", j, i, px->ucol, px->urow, dx, dy, weight[i][j], v[i][j][0]);*/
      }
    for(iband=0; iband<warp->nbands; iband++) {
      tref = 0.0;
      tsum = 0.0;
      for(i=0; i<4; i++) 
	for(j=0; j<4; j++) { 
	  tref += v[i][j][iband] * weight[i][j];
	  tsum += weight[i][j];
	}
      /*printf("%6.3f ", tsum);*/
      outRow[iband] = (int)(tref/tsum+0.5);
    } 
    /* printf("%d %d %d %d %d\n", v[1][1][0], v[1][2][0], v[2][1][0], v[2][2][0], outRow[0]);*/
  } /* end of cubic convolution */

}




 

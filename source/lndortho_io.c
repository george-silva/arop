/**
 * ! Description
 *   read/write subroutines for Landsat orthorectification program 
 *
 * ! Credits
 *   Feng Gao (Feng.Gao@nasa.gov), Developer
 *   Jeff Masek (Jeffrey.G.Masek@nasa.gov)
 *   Robert Wolfe (Robert.E.Wolfe@nasa.gov)
 *
 * ! Revision 2.0    12/20/2007
 */

#include "lndortho.h"

/* read DEM height (one line) for a given location in output image (only depends on y or irow) */
int readHeightLine(PIXEL *px, DEM *sdem)
{
  int iline; 
  long offset; 
  double y;;
  
  /* px->irow and px->icol used center of pixel */
  y= px->uly - px->irow * px->res;
  iline = (int)((sdem->uly - y) / sdem->res);
  sdem->irow = iline; 
  if(iline>=0 && iline<sdem->nrows) {
    if(strcasecmp(sdem->fileType, "GEOTIFF")==0) {
      if (!TIFFReadScanline(sdem->fp_tiff, sdem->buf, iline, 0)) {
	fprintf(stderr, "Read line %d error\n", iline);
	return FAILURE;
      } 
    }
    else {
      offset = (long) iline * sdem->ncols * sizeof(int16);
      fseek(sdem->fp, offset, 0);
      fread((sdem->buf), sizeof(int16), sdem->ncols, sdem->fp);
    }
  }

  return SUCCESS;
}



/* read one line data from Landsat SR input file */
int readLandsatLine(WARP_LANDSAT *sr, int irow, int16 **oneRow)
{
  uint8 buf8[MAX_NCOLS];
  int16 buf[MAX_NCOLS];
  int icol, iband;
  long int offset;

  /* initialize data */  
  for(icol=0; icol<sr->lnd.ncols; icol++)   {
    for(iband=0; iband<sr->nbands; iband++)
      oneRow[icol][iband] = sr->lnd.fillValue;
  }
  
  for(iband=0; iband<sr->nbands; iband++) {

    if(strcasecmp(sr->lnd.fileType, "GEOTIFF")==0) {
      /* read data array from GeoTiff file */
      if(sr->nbyte[iband]==1) {
	if (!TIFFReadScanline(sr->fp_tiff[iband], buf8, irow, 0)) {
	  fprintf(stderr, "Read line %d error\n", irow);
	  return FAILURE;
	} 
	for(icol=0; icol<sr->lnd.ncols; icol++)
	  buf[icol] = buf8[icol];
      }
      else
	if (!TIFFReadScanline(sr->fp_tiff[iband], buf, irow, 0)) {
	  fprintf(stderr, "Read line %d error\n", irow);
	  return FAILURE;
	} 
    }
    else {
      offset = (long) (irow * sr->lnd.ncols * sr->nbyte[iband]);
      /* read data array from binary file */
      fseek(sr->fp[iband], offset, 0);
      if(sr->nbyte[iband]==1) {
	fread(buf8, sr->nbyte[iband], sr->lnd.ncols, sr->fp[iband]);
	for(icol=0; icol<sr->lnd.ncols; icol++)
	  buf[icol] = buf8[icol];
      }
      else
	fread(buf, sr->nbyte[iband], sr->lnd.ncols, sr->fp[iband]);
    }
    for(icol=0; icol<sr->lnd.ncols; icol++)
      oneRow[icol][iband] = buf[icol];
  }

  return SUCCESS;
}



/**
 * create data block in memory for fast data access in resampling process
 * the data block is in the size of BLOCK_NROWS * ncols
 * By keeping current input line as the central line of data block, data block
 * needs to be updated for each line loop
 * ! Input
 *   sr: Landsat data
 * ! Output
 *   inRows: input data block  
 *   outRow: one row for output
 */
int loadDataBlock(WARP_LANDSAT *srIn, int current_row, int16 ***inRows)
{
  int i, j, iband;
  int16 **inTemp;

  alloc_2dim_contig((void ***) (&inTemp), srIn->lnd.ncols, srIn->nbands, sizeof(int16));

  /* load first sevaral lines to memory */
  for(i=current_row-BLOCK_NROWS/2; i<=current_row+BLOCK_NROWS/2; i++) {
    if(i>=0 && i<srIn->lnd.nrows) {
      if(readLandsatLine(srIn, i, inTemp)==FAILURE) exit(1);
      for(j=0; j<srIn->lnd.ncols; j++) {  
	for(iband=0; iband<srIn->nbands; iband++)
	  inRows[i-current_row+BLOCK_NROWS/2][j][iband] = inTemp[j][iband];
      }
    }
    else {
      for(j=0; j<srIn->lnd.ncols; j++) {  
	for(iband=0; iband<srIn->nbands; iband++)
	  inRows[i-current_row+BLOCK_NROWS/2][j][iband] = srIn->lnd.fillValue;
      }
    }
  }

  free_2dim_contig((void **) inTemp);
  return SUCCESS;
}



/* read data from input Landsat image for a given pixel location */
int readPixel(WARP_LANDSAT *sr, PIXEL *px)
{
  int iband;
  long int offset;
  int16 data, buf[MAX_NCOLS];
  uint8 buf8[MAX_NCOLS];

  /* initialize data */  
  for(iband=0; iband<sr->nbands; iband++)
    px->ref[iband] = sr->lnd.fillValue;

  if(px->urow<0||px->urow>=sr->lnd.nrows||px->ucol<0||px->ucol>=sr->lnd.ncols) return FAILURE;

  for(iband=0; iband<sr->nbands; iband++) {
    if(strcasecmp(sr->lnd.fileType, "GEOTIFF")==0) {
      /* read data array from GeoTiff file */
      if(sr->nbyte[iband]==1) {
	if (!TIFFReadScanline(sr->fp_tiff[iband], buf8, px->urow, 0)) {
	  fprintf(stderr, "Read line %d error\n", px->urow);
	  return FAILURE;
	} 
	data = buf8[px->ucol];
      }
      else {
	if (!TIFFReadScanline(sr->fp_tiff[iband], buf, px->urow, 0)) {
	  fprintf(stderr, "Read line %d error\n", px->urow);
	  return FAILURE;
	} 
	data = buf[px->ucol];
      }
    }
    else {
      offset = (long) (px->urow * sr->lnd.ncols + px->ucol) * sr->nbyte[iband];
      fseek(sr->fp[iband], offset, 0);
      if(sr->nbyte[iband]==1) { 
	fread(&(buf8[px->ucol]), sr->nbyte[iband], 1, sr->fp[iband]);
	data = buf8[px->ucol];
      }
      else { 
	fread(&(buf[px->ucol]), sr->nbyte[iband], 1, sr->fp[iband]);
	data = buf[px->ucol];
      }
    }
    
    px->ref[iband] = data;
  }
  
  return SUCCESS;
}


/* write a line to output Landsat image */
int writeOneRow(OUT_LANDSAT *out, int16 **outRow)
{
  int iband, icol;
  int16 buf;

  for(iband=0; iband<out->nbands; iband++) 
    for(icol=0; icol<out->lnd.ncols; icol++) {
      buf = outRow[icol][iband];
      fwrite(&buf, out->nbyte[iband], 1, out->fp[iband]);
    }
  return SUCCESS;
}


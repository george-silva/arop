/**
 * ! Description
 *   header file for Landsat orthorectification program 
 *
 * ! Credits
 *   Feng Gao (Feng.Gao@nasa.gov), Developer
 *   Jeff Masek (Jeffrey.G.Masek@nasa.gov)
 *   Robert Wolfe (Robert.E.Wolfe@nasa.gov)
 *
 * ! Revision 2.2.6  8/22/2011
 * ! Revision 2.2.5  1/6/2010
 * ! Revision 2.2.3  2/27/2009
 * ! Revision 2.0    2/12/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>

#include "geotiffio.h"
#include "xtiffio.h"
#include "tiffio.h"

#include "cproj.h"
#include "proj.h"

#define MAX_STRLEN  1000      /* maximum allowed string length */
#define MAX_NBANDS  20        /* maximum number of bands */
#define MAX_NCOLS   100000    /* maximum number of column of imagery */
#define BLOCK_NROWS 71        /* number of lines to hold on a data block in memory */
#define MIN_COARSE_NROWS 700  /* minimum nrows for coarse image in pyramid registration */

#define MAX_NUM     1e20
#define FAILURE     -1
#define SUCCESS     0
#ifndef PI
#define PI          3.14159265
#endif
#define DEM_FILL    -32768


/* define parameters for control point searching */
typedef struct {
  int   CHIP_SIZE;                /* searching control point chip window size */
  int   CP_SEED_WIN;              /* one control point per CP_REG*CP_REG area */
  int   MAX_SIZE;                 /* maximum searching zone size */
  int   MAX_CPS;                  /* maximum control point candidates under current setting */
  int   MIN_STDEV;                /* minimum stdev for a valid control chip */
  int   MAX_SHIFT;                /* define possible maximum shift around (+-)MAX_SHIFT */
  int   MAX_NUM_HIGH_CORR;        /* acceptable maximum number of high correlation */
  int   MIN_ACCEPTABLE_NCP;       /* minimum acceptable number of control point */
  int   MAX_NUM_ITER;             /* maximum iterations allowed for re-ortho */ 
  float ACCEPTABLE_CORR;          /* acceptable control point must have corr > ACCEPTABLE_CORR */
  float MAX_AVE_ERROR;            /* acceptable maximum average error in pixel */
  int   PRELIMINARY_REGISTRATION; /* flag of preliminary registration */
  int   COARSE_SCALE;             /* reduce scale to estimate shift for preliminary registration */
  int   COARSE_CP_SEED_WIN;       /* seed window size for coarse resolution image */
  int   COARSE_MAX_SHIFT;         /* maximum shift searhing for coarse resolution image */
  float MAX_RMSE;                 /* maximum rmse allowed for control points in pixels */
  float DIFF_RMSE;                /* acceptable RMSE difference between two consecutive iteration */
} CP_PARAMETERS; 
CP_PARAMETERS CP_KEYS;

/* define semi-major axis for satellites in meters */
enum Sat_Major_Axis {Landsat1 = 7285438, 
		     Landsat2 = 7285989,
		     Landsat3 = 7285776,
		     Landsat4 = 7083465,
		     Landsat5 = 7083465,
		     Landsat7 = 7077700,
		     CBERS1   = 7148860,
		     CBERS2   = 7115798,
		     TERRA    = 7085000,
		     AWIFS    = 7195110,
		     HJ1      = 7020097,
		     MySensor = 7100000};

/* define user's choice */
enum User_Choice {REG_ONLY=0, ORTHO_ONLY=1, DO_BOTH=2, VERIFY_ONLY=3};

/* define forward and inverse projection */
enum Projection_Type {FORWARD_PROJ=0, BACKWARD_PROJ=1};
 
#define round(x) ((x)>=0? (long)((x)+0.5):(long)((x)-0.5))

typedef struct {

  char  fileName[MAX_STRLEN];  /* the band selected for control point searching */ 
  char  fileType[MAX_STRLEN];
  int   nrows;           
  int   ncols;
  float ulx;              /* x value of up-left corner */
  float uly;              /* y value of up-left corner */
  float res;              /* spatial resolution */
  uint8 fillValue;        
  int   nbyte;            /* number of byte for data */
  int   satValue;         /* DN saturation value */

  /* terrain displacement and nadir track related parameters */
  double centerLonLat[2]; /* latitude and longitude for the central pixel */
  double lpara[2];        /* saves nadir view linear equation in: irow=lpara[0]*icol+lpara[1] */
  double radius;          /* Earth radius from geocenrtic to sea level */
  double satMajorAxis;    /* satellite semi-major axis */
  double altitude;        /* Satellite altitude above sea level */
  int    porder;          /* polynomial order used for output and base precise registration */    
  float  a[6], b[6];      /* derived 1 or 2 order polynomial coefficient between output and base */

  TIFF *fp_tiff;          /* if it is a GEOTIFF file */
  FILE *fp;               /* if it is a binary file */

} LANDSAT;


typedef struct {
  int  utm_zone;          /* all GeoCover data are saved in UTM project */
  long datum;             /* default in WGS84 (=12) */
  LANDSAT lnd;
} BASE_LANDSAT;           /* define BASE Landsat structure for consistent uses */


typedef struct {
  /* define projection */
  long sys;
  int  utm_zone;
  double param[15];
  long unit;
  long datum;
} PROJECTION;

 
typedef struct {
  LANDSAT lnd;
  int  nbands;
  char fileName[MAX_NBANDS][MAX_STRLEN]; 
  int  nbyte[MAX_NBANDS];      /* number of bytes for each input */
  int  satValue[MAX_NBANDS];   /* DN saturation values for each band */ 
  char tempRotateFileName[MAX_NBANDS][MAX_STRLEN]; 
  char tempProjectFileName[MAX_NBANDS][MAX_STRLEN]; 
  char resampleMethod[MAX_STRLEN];
  /* orientation angle */
  float orientationAngle;              
  /* sensor pointing angles */
  float pointingAngle;
  /* define projection */
  PROJECTION proj;
  TIFF *fp_tiff[MAX_NBANDS];
  FILE *fp[MAX_NBANDS];
} WARP_LANDSAT;               /* structure to store WARP image info */


typedef struct {
  LANDSAT lnd;  
  int  utm_zone;
  int  nbands;
  char fileName[MAX_NBANDS][MAX_STRLEN];
  int  nbyte[MAX_NBANDS];  
  FILE *fp[MAX_NBANDS];
  int  num_redo;
  int  cps_passed;
  int  checking_passed;
} OUT_LANDSAT;           /* structure to store OUTPUT image info */



typedef struct {

  char base_fileName[MAX_STRLEN];
  char base_fileType[MAX_STRLEN];
  char base_satellite[MAX_STRLEN];   
  int  base_ncols;
  int  base_nrows;
  float base_res;
  float base_ulx;
  float base_uly;
  int  utm_zone;
  long datum;

  char warp_fileType[MAX_STRLEN];
  int  warp_ncols;
  int  warp_nrows;
  float warp_res;
  float warp_ulx;
  float warp_uly;
  float warp_ulx_degree;
  float warp_uly_degree;
  float warp_angle;
  int  warp_nbands;
  char warp_fileName[MAX_NBANDS][MAX_STRLEN];
  int  warp_nbyte[MAX_NBANDS];
  char warp_match_fileName[MAX_STRLEN];
  /* satellite type to determin major axis */
  char warp_satellite[MAX_STRLEN];   
  float warp_satellite_pangle;
  /* projection for warp image if different */
  long warp_sys;
  int  warp_utm_zone;
  double warp_param[15];
  long warp_unit;
  long warp_datum;

  char dem_fileName[MAX_STRLEN];
  /* projection for SRTM DEM data if different */
  long dem_sys;
  int  dem_utm_zone;
  double dem_param[15];
  long dem_unit;
  long dem_datum;
  
  float out_res;
  char  resampleMethod[MAX_STRLEN];
  char  extentSource[MAX_STRLEN];
  int   bkValue;
  int   satValue, satValue2;
  int   data_type;
  int   porder;
  float out_ulx, out_uly;
  float out_lrx, out_lry;

  char  out_fileName[MAX_NBANDS][MAX_STRLEN];
  char  out_match_fileName[MAX_STRLEN];
  char  cp_parameters_file[MAX_STRLEN];

  int   option;

  int   rotation_flag;            /* 1 = rotation needs; 0 = no */
  int   projection_flag;          /* 1 = reprojection needs; 0 = no */
  int   projection_dem_flag;      /* 1 = reprojection needs; 0 = no */

} IN_PARAMETERS;


typedef struct {

  double irow;           /* location in orthorectified space (output image) */
  double icol;
  double height;         /* associated height from DEM */
  double off_nadir_dis;  /* distance to nadir view (in meters) */
  double ds;             /* displacement distance caused by terrain */
  double dx;             /* shift in x direction (sample) in input image */
  double dy;             /* shift in y direction (row) in input image */
  double res;            /* pixel resolution */
  double ulx, uly;       /* Landsat upperleft coordinate in output image */ 
  int    sign;           /* left side = -1; right side = 1 */

  double ux;             /* location in un-orthorectified space */
  double uy;
  int    urow;           /* position in un-orthorectified space (input image) */
  int    ucol;
  int16  ref[MAX_NBANDS];    /* data in [iband] (from input) */

} PIXEL;

typedef struct {

  char   fileName[MAX_STRLEN];
  char   fileType[MAX_STRLEN];
  TIFF   *fp_tiff;
  FILE   *fp;
  int16  *buf;
  int    nrows;
  int    ncols;  
  double res;
  double ulx;
  double uly;
  int    irow;
  int    icol;
  PROJECTION proj;

} DEM;


/* in lndortho_compute.c */
int   extractNadirPath(LANDSAT *sr);
void  computeNadirViewfromWarp(LANDSAT *out, WARP_LANDSAT *warp);
void  computeRadius(LANDSAT *sr);
void  computeOffDis(LANDSAT *sr, PIXEL *px);
void  computeDisplacement(LANDSAT *sr, PIXEL *px);
void  getAffine(float **cp, int num_cps, int n, float a[]);
void  solveEquation(float s[][10], float ss[], int m, float a[]);
float cubic(float x);

/* in lndortho_io.c */
int readHeightLine(PIXEL *px, DEM *sdem);
int readLandsatLine(WARP_LANDSAT *sr, int irow, int16 **oneRow);
int loadDataBlock(WARP_LANDSAT *srIn, int current_row, int16 ***inRows);
int readPixel(WARP_LANDSAT *sr, PIXEL *px);
int writeOneRow(OUT_LANDSAT *sr, int16 **outRow);

/* in lndortho_gchk.c */
int  preliminaryRegistration(BASE_LANDSAT  *wbase, WARP_LANDSAT  *wwarp);
int  preciseRegistration (BASE_LANDSAT  *wbase, WARP_LANDSAT  *wwarp, DEM *sdem);
int  initialOrthoCheck(BASE_LANDSAT  *wbase, WARP_LANDSAT  *wwarp, OUT_LANDSAT *wout, DEM *sdem);
int  finalOrtho(BASE_LANDSAT  *wbase, WARP_LANDSAT  *warp, WARP_LANDSAT  *temp_warp, WARP_LANDSAT *wwarp, OUT_LANDSAT *wout, OUT_LANDSAT *out, DEM *sdem, IN_PARAMETERS *inPars);
void traceLocBackward(PIXEL *px, BASE_LANDSAT  *wbase, WARP_LANDSAT  *warp, WARP_LANDSAT *temp_warp, WARP_LANDSAT *wwarp, OUT_LANDSAT *wout, IN_PARAMETERS *inPars);
int  preciseRegistrationWithSampling (BASE_LANDSAT  *wbase, WARP_LANDSAT  *warp, WARP_LANDSAT  *temp_warp, WARP_LANDSAT  *wwarp, OUT_LANDSAT *out, OUT_LANDSAT *wout, IN_PARAMETERS *inPars);
int  matchingVerification (BASE_LANDSAT  *wbase, WARP_LANDSAT  *wwarp, OUT_LANDSAT *wout);
int  searchCPs(LANDSAT *base, LANDSAT *warp, float **match_p); 
int  computeOrthoShifts(LANDSAT *base, DEM *sdem, float **match_p, int num_total);
int  filterCandidates(LANDSAT *base, LANDSAT *warp, float **match_p, int num_total);
int  checkingCPs(LANDSAT *base, LANDSAT *warp, float **match_p, int num_total);
int  checkingTPs(LANDSAT *base, LANDSAT *out, LANDSAT *warp, float **match_p, int num_total); 
void initializePoly(LANDSAT *base, LANDSAT *out); 
void redoPreciseGeoCorrection(float **match_p, int tnum, LANDSAT *warp);
void seekSubPixelLoc(float **neighboxr_cc, float *maxc_loc);
void resampleData(PIXEL *px, WARP_LANDSAT *warp, int crow, int16 ***inRows, int16 *outRow);

/* in lndortho_util.c */
int  getInParameter(IN_PARAMETERS *inPars, int argc, char *argv[]);
int  getSRTMMetaInfo(DEM *sdem, IN_PARAMETERS *inPars);
int  getLandsatMetaInfo(IN_PARAMETERS *inPars, BASE_LANDSAT *base, WARP_LANDSAT *warp, OUT_LANDSAT *out);
void getSatMajorAxis(char *satellite, LANDSAT *lnd);
int  makeCoarse(LANDSAT *fine, LANDSAT *coarse);
int  rotateWarp(WARP_LANDSAT *warp);
int  reproject(WARP_LANDSAT *warp, BASE_LANDSAT *base);
int  reprojectDEM(DEM *sdem, BASE_LANDSAT *base);
int  convertCoor(BASE_LANDSAT *base, double *bx, double *by, PROJECTION *proj, double *wx, double *wy, int ptype);
int  convertSpace(BASE_LANDSAT *base, WARP_LANDSAT *warp, OUT_LANDSAT *out, 
		  BASE_LANDSAT *wbase, WARP_LANDSAT *wwarp, OUT_LANDSAT *wout, 
		  IN_PARAMETERS *inPars);
int  aggregateImage(LANDSAT *fine, LANDSAT *coarse);
int  writeENVIheader(WARP_LANDSAT *warp, OUT_LANDSAT *out);
void rotation(double x, double y, double *x1, double *y1, double angle);
void copyLandsat(LANDSAT *dest, LANDSAT *src);
int  createTempWarp(WARP_LANDSAT  *twarp, WARP_LANDSAT  *warp);
void updateOutCoor(LANDSAT *src, LANDSAT *dest);
void usage(char *argv[]);
void close_LANDSAT_file(LANDSAT *sr);
void clean_coarse(LANDSAT *sr);
void alloc_1dim_contig (void **, int, int);
void alloc_2dim_contig (void ***, int, int, int);
void alloc_3dim_contig (void ****, int, int, int, int);
void free_2dim_contig  (void **);
void free_3dim_contig  (void ***);

/*#define REG_DEBUG
  #define DEBUG*/
#ifdef DEBUG
#define DEBUG_icol 453
#define DEBUG_irow 1325
#endif

/**
 * ! Description
 *   georeference computation subroutinea for Landsat orthorectification program 
 *
 * ! Credits
 *   Feng Gao (Feng.Gao@nasa.gov), Developer
 *   Jeff Masek (Jeffrey.G.Masek@nasa.gov)
 *   Robert Wolfe (Robert.E.Wolfe@nasa.gov)
 *
 * ! Revision 2.2    5/20/2008
 */

#include "lndortho.h"


/** extract nadir path track from image 
    NOTE: it only works for smooth edge image
*/
int extractNadirPath(LANDSAT *sr)
{
  uint8 buf[MAX_NCOLS];
  int i, j, num, k, trend;
  int start[MAX_NCOLS][2];
  double sum_x, sum_y, sum_x2, sum_xy;
  long int offset;

  k = 0;
  for(i=0; i<sr->nrows; i++) {
    if(strcasecmp(sr->fileType, "GEOTIFF")==0) {
      /* read data array from GeoTiff file */
      if (!TIFFReadScanline(sr->fp_tiff, buf, i, 0)) {
	fprintf(stderr, "Read line %d error\n", i);
	return FAILURE;
      } 
    }
    else {
      offset = (long) (i * sr->ncols);  
      /* read data array from binary file */
      fseek(sr->fp, offset, 0);
      fread(buf, sizeof(char), sr->ncols, sr->fp);
    }
    
    for(j=0; j<sr->ncols; j++)
      if(buf[j] != sr->fillValue) {
	start[k][0] = j;
	start[k][1] = i;
	k++;
	break;
      }
  } /* end of row loop */

  /* analyze line trends using first 100 points */
  trend = 0;
  for(i=0; i<100; i++) {
    if((start[i+1][0]-start[i][0])>0) trend++;
    if((start[i+1][0]-start[i][0])<0) trend--;
  }

  /* use least square fit approah to find slope */
  num = k;
  k = 0;
  sum_x = 0.0;
  sum_y =0.0;
  sum_x2 = 0.0;
  sum_xy = 0.0;
  for(i=0; i<num-1; i++) {
    /* if direction changed */
    if((start[i+1][0]-start[i][0])*trend <= 0) continue;
    sum_x += start[i][0];
    sum_y += start[i][1];
    sum_x2 += start[i][0]*start[i][0];
    sum_xy += start[i][0]*start[i][1];
    k++;
  } 

  sr->lpara[0] = (k*sum_xy - sum_x*sum_y)/(k*sum_x2 - sum_x*sum_x);
  sr->lpara[1] = (sr->nrows/2.0) - (sr->ncols/2.0) * sr->lpara[0];

  /* return to the start of file */   
  if(strcasecmp(sr->fileType, "BINARY")==0)
    fseek(sr->fp, 0, 0);

  return SUCCESS;
}


/** 
 * use the slope of track function from warp image to determine the nadir track in out
 */ 
void computeNadirViewfromWarp(LANDSAT *out, WARP_LANDSAT *warp)
{
  float wcx, wcy, cx, cy;
  float dis, alfa, delta_x, delta_y;

  /* use slope from warp image */
  out->lpara[0] = warp->lnd.lpara[0];

  /* consider pointing angle */
  dis = warp->lnd.altitude*tan(warp->pointingAngle*PI/180.0);
  alfa = atan(-1.0*out->lpara[0]);
  delta_x = dis * sin(alfa);
  delta_y = dis * cos(alfa);
  
  wcx = warp->lnd.ncols / 2.0*warp->lnd.res - delta_x; 
  wcy = warp->lnd.nrows / 2.0*warp->lnd.res + delta_y;
  
  /* compute central point location of warp image in the output image */
  cx = ((warp->lnd.ulx + wcx) - out->ulx)/out->res;
  cy = (out->uly - (warp->lnd.uly - wcy))/out->res;

  /* then compute intercept based on the central location in output image */
  out->lpara[1] = cy - cx * out->lpara[0];
}


/**
 * Compute earth radius at any given location based on the World geodetic System 1984 (WGS84) datum
 * Note that earth radius needs local adjustment to consistent with current SRTM DEM data (also in WGS84) 
 *      SRTM_height = height_above_ellipsoid_surface - geoid_height  (above sea level)
 *      Earth_radius = radius_of_ellipsoid + geoid-height (up to sea level)
 *
 * ! Input
 *     latitude: earth geodetic latitude (normally used latitude term)
 *     longitude: earth longitude 
 * ! Output
 *     sr->radius: Earth radius in WGS84 geodetic system
 */
void computeRadius(LANDSAT *sr)
{
  /* use WGS84 as geodetic datum */
  double e_radius = 6378137.0;    /* equatorial radius  */
  double e = 0.081819190842622;   /* first eccentricity */
  /* ten by ten degrees WGS84 Geoid heights from (90, -180) to (-90, 170) */
  int geoid_height[19][36] = {   
    {13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13},
    {3,1,-2,-3,-3,-3,-1,3,1,5,9,11,19,27,31,34,33,34,33,34,28,23,17,13,9,4,4,1,-2,-2,0,2,3,2,1,1},
    {2,2,1,-1,-3,-7,-14,-24,-27,-25,-19,3,24,37,47,60,61,58,51,43,29,20,12,5,-2,-10,-14,-12,-10,-14,-12,-6,-2,3,6,4},
    {2,9,17,10,13,1,-14,-30,-39,-46,-42,-21,6,29,49,65,60,57,47,41,21,18,14,7,-3,-22,-29,-32,-32,-26,-15,-2,13,17,19,6},
    {-8,8,8,1,-11,-19,-16,-18,-22,-35,-40,-26,-12,24,45,63,62,59,47,48,42,28,12,-10,-19,-33,-43,-42,-43,-29,-2,17,23,22,6,2},
    {-12,-10,-13,-20,-31,-34,-21,-16,-26,-34,-33,-35,-26,2,33,59,52,51,52,48,35,40,33,-9,-28,-39,-48,-59,-50,-28,3,23,37,18,-1,-11},
    {-7,-5,-8,-15,-28,-40,-42,-29,-22,-26,-32,-51,-40,-17,17,31,34,44,36,28,29,17,12,-20,-15,-40,-33,-34,-34,-28,7,29,43,20,4,-6},
    {5,10,7,-7,-23,-39,-47,-34,-9,-10,-20,-45,-48,-32,-9,17,25,31,31,26,15,6,1,-29,-44,-61,-67,-59,-36,-11,21,39,49,39,22,10},
    {13,12,11,2,-11,-28,-38,-29,-10,3,1,-11,-41,-42,-16,3,17,33,22,23,2,-3,-7,-36,-59,-90,-95,-63,-24,12,53,60,58,46,36,26},
    {22,16,17,13,1,-12,-23,-20,-14,-3,14,10,-15,-27,-18,3,12,20,18,12,-13,-9,-28,-49,-62,-89,-102,-63,-9,33,58,73,74,63,50,32},
    {36,22,11,6,-1,-8,-10,-8,-11,-9,1,32,4,-18,-13,-9,4,14,12,13,-2,-14,-25,-32,-38,-60,-75,-63,-26,0,35,52,68,76,64,52},
    {51,27,10,0,-9,-11,-5,-2,-3,-1,9,35,20,-5,-6,-5,0,13,17,23,21,8,-9,-10,-11,-20,-40,-47,-45,-25,5,23,45,58,57,63},
    {46,22,5,-2,-8,-13,-10,-7,-4,1,9,32,16,4,-8,4,12,15,22,27,34,29,14,15,15,7,-9,-25,-37,-39,-23,-14,15,33,34,45},
    {21,6,1,-7,-12,-12,-12,-10,-7,-1,8,23,15,-2,-6,6,21,24,18,26,31,33,39,41,30,24,13,-2,-20,-32,-33,-27,-14,-2,5,20},
    {-15,-18,-18,-16,-17,-15,-10,-10,-8,-2,6,14,13,3,3,10,20,27,25,26,34,39,45,45,38,39,28,13,-1,-15,-22,-22,-18,-15,-14,-10},
    {-45,-43,-37,-32,-30,-26,-23,-22,-16,-10,-2,10,20,20,21,24,22,17,16,19,25,30,35,35,33,30,27,10,-2,-14,-23,-30,-33,-29,-35,-43},
    {-61,-60,-61,-55,-49,-44,-38,-31,-25,-16,-6,1,4,5,4,2,6,12,16,16,17,21,20,26,26,22,16,10,-1,-16,-29,-36,-46,-55,-54,-59},
    {-53,-54,-55,-52,-48,-42,-38,-38,-29,-26,-26,-24,-23,-21,-19,-16,-12,-8,-4,-1,1,4,4,6,5,4,2,-6,-15,-24,-33,-40,-48,-50,-53,-52},
    {-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30}
  };
  
  int step = 10;   /* geoid table step (10 degrees) */
  double latitude, longitude;
  double lrad, radius, x, y, dx, dy;
  double data[4], w[4], sum_weight, geoh;
  int i, ilat, ilon, jlat, jlon;
  
  latitude = sr->centerLonLat[1];
  longitude = sr->centerLonLat[0];

  /* compute ellipsoid radius for geodetic latitude */ 
  lrad = latitude/180.0*PI;
  x = e_radius*cos(lrad)/sqrt(1.0-e*e*sin(lrad)*sin(lrad));
  y = e_radius*(1-e*e)*sin(lrad)/sqrt(1.0-e*e*sin(lrad)*sin(lrad));
  radius = sqrt(x*x+y*y);

  /* interoplate Geoid height from ten by ten degree table */
  ilat = (int)((90-latitude)/step);
  ilon = (int)((longitude+180)/step);

  /* find lower right sample position */   
  jlat = ilat+1;
  jlon = ilon+1;
  if(jlon == 36) jlon = 0;
  

  /* use four neighbor point to interoplate geoid height */
  dy = fabs(latitude-(90-ilat*step));
  dx = fabs(longitude-(ilon*step-180));
  if(dx==0&&dy==0)
    w[0] = MAX_NUM;
  else
    /* use reversed distance to compute weights */
    w[0] = 1.0/sqrt(dx*dx+dy*dy);
  data[0] = geoid_height[ilat][ilon];
  if(step-dx==0 && dy==0)
    w[1] = MAX_NUM;
  else
    w[1] = 1.0/sqrt((step-dx)*(step-dx)+dy*dy);
  data[1] = geoid_height[ilat][jlon];
  if(step-dx==0 && step-dy==0)
    w[2] = MAX_NUM;
  else
    w[2] = 1.0/sqrt((step-dx)*(step-dx)+(step-dy)*(step-dy));
  data[2] = geoid_height[jlat][jlon];
  if(dx==0 && step-dy==0)
    w[3] = MAX_NUM;
  else
    w[3] = 1.0/sqrt(dx*dx+(step-dy)*(step-dy));
  data[3] = geoid_height[jlat][ilon];
  
  sum_weight = 0.0;
  for(i=0; i<4; i++)
    sum_weight += w[i];
  
  /* compute interoplated geoid height using reversed distance as weight */
  geoh = 0.0;
  for(i=0; i<4; i++)
    geoh += data[i]*w[i]/sum_weight;

  /*for(i=0; i<4; i++)
    printf("%f %f \n", data[i], w[i]);
    printf("%f %f %f %f\n", dx, dy, radius, geoh);*/

  /* to consistent with SRTM DEM (used WGS84), geoid height need to be added to ellipsoid radius */
  sr->radius = radius+geoh;

   /* compute altitude (assume circle since eccentricity is very small) */
  sr->altitude = sr->satMajorAxis - sr->radius;
}



/**
 * compute distance from a given pixel to nadir view line
 * ! Input:
 *     px->(irow, icol): location of a given pixel 
 *     sr->lpara: nadir view linear equation
 * ! Output
 *     px->off_nadir_dis: distance from given pixel to nadir view in meters
 *     px->sign: sign of shiftment (plus for right side to nadir and minus for left side to nadir)
 *               (nadir track start from North to South)
 */
void computeOffDis(LANDSAT *sr, PIXEL *px)
{
  double x0, y0, a, b, c;
  double dis;
  
  /* center of pixel */
  x0 = (px->icol+0.5);
  y0 = (px->irow+0.5);
  a = sr->lpara[0];
  b = -1.0;
  c = sr->lpara[1];

  /* distance from point (x0, y0) to ax+by+c=0 */
  dis = fabs(a*x0+b*y0+c)/sqrt(a*a+b*b);
  px->off_nadir_dis = dis*sr->res;

  /* decides sign of displacement (right to nadir track use -; left to nadir track use +) */
  if(px->icol * a + px->irow * b + c > 0)
    px->sign = -1;
  else
    px->sign = 1;

}


/**
 * compute displacement observed from satelllite caused by terrain effect 
 * ! Input
 *     sr->Re: earth radius (based on WGS84 geoid datum (i.e. from geocentric to sea level)
 *     px->off_nadir_dis: distance between observation location and satellite nadir location
 *     px->height: terrain height (start from sea level)
 * ! Output 
 *     px->ds: displacement in meters
 *     px->dx: shift in map x direction from orthorectified image to un-orthorectified image (use orthorectified image as base)
 *     py->dy: shift in map y direction from orthorectified image to un-orthorectified image
 */
void computeDisplacement(LANDSAT *sr, PIXEL *px)
{
  double Re, off_nadir_dis, height, Alt;
  double angle_s, angle_d, angle_dd, angle_tt, angle_zz, angle_ds;
  double LOS, AF, BF, cita;

  /* Earth radius */
  Re = sr->radius;
  /* Landsat 7 altitude (assume circle since eccentricity = 0.00117) */
  Alt = sr->altitude;
  /* distance from nadir */
  off_nadir_dis = px->off_nadir_dis;  
  /* terrain height from DEM */
  height = px->height;
  
  /* Earth centered angle */
  angle_s = off_nadir_dis / Re;
  
  /* compute line of sight */
  LOS = sqrt(Re*Re+(Re+Alt)*(Re+Alt)-2.0*Re*(Re+Alt)*cos(angle_s));
  angle_d = asin(Re*sin(angle_s)/LOS);

  /* compute temporary angles */  
  angle_tt = acos((Re+Alt)*sin(angle_d)/Re);

  /* compute displacement angle */
  AF = (Re+Alt)*cos(angle_d)-(Re+height)*sin(angle_tt);
  BF = (Re+Alt)*sin(angle_d)*height/Re;
  angle_dd = atan(BF/AF);
  angle_zz = asin((Re+Alt)*sin(angle_d+angle_dd)/Re);
  angle_ds = angle_zz - angle_s - angle_d - angle_dd;
  
  /* comput shift */
  px->ds = Re * angle_ds;
  /*printf("angle=%f ", (angle_d+angle_dd)/PI*180.0); */	     

  /* compute shift in x(icol) and y (irow) directions */
  cita = PI/2.0 - atan(-1.0*(sr->lpara[0]));
  px->dx = px->ds * cos(cita) * px->sign;
  px->dy = px->ds * sin(cita) * px->sign;

  /* convert to map space (reverse row direction to map Y-axis direction) */
  px->dy = -1.0 * px->dy;

}



/* the cubic convolution function calculate cubic convolution 
   source: ERDAS field guide (see Atkinson, IJRS, 1985) */
float cubic(float x) 
{
  float a = -0.5;
  float ax, ret;

  ax = fabs(x);
  if(ax<=1) 
    ret = (a+2)*ax*ax*ax-(a+3)*ax*ax+1;
  else if(ax<=2)
    ret = a*ax*ax*ax-5*a*ax*ax+8*a*ax-4*a;
  else
    ret = 0.0;
  return ret;
}



/**
 * create coefficient matrix and retrieve affine parameters
 * ! Input
 *     float cp[][3]  control point array 0=x, 1=y (in master image), 2= x or y (in slave image)
 *     int num_cps: number of control points
 *     int n: order of equation used 
 * ! Output
 *     float a[]: equation solution 
 */
void getAffine(float **cp, int num_cps, int n, float a[])       
{
  int i0, i1, i2, i3, j, j1, m, jf;
  float s[10][10], ss[10];

  m=(n+1)*(n+2)/2;
  for(i1=0;i1<m;i1++)
    for(i2=0;i2<m;i2++)
      s[i1][i2]=0.0;
 
  for(i0=0;i0<=n;i0++)
    for(i1=0;i1<=i0;i1++) {
      jf=i0*(i0+1)/2+i1;
      for(i2=0;i2<=n;i2++)
	for(i3=0;i3<=i2;i3++)
	  {
	    j=(i2+1)*i2/2+i3;
	    if(jf>j)
	      continue;
	    s[jf][j]=0.0;
	    for(j1=0;j1<num_cps;j1++)
	      s[jf][j]+=pow(cp[j1][0],i2-i3)*pow(cp[j1][1],i3)*pow(cp[j1][0],i0-i1)*pow(cp[j1][1],i1);
	  }
    }

  for(i1=1;i1<m;i1++)
    for(i2=0;i2<i1;i2++)
      s[i1][i2]=s[i2][i1];
  for(i1=0;i1<=n;i1++)
    for(i2=0;i2<=i1;i2++) {
      j=i1*(i1+1)/2+i2;   
      ss[j]=0.0;
      for(j1=0;j1<num_cps;j1++)
	ss[j]+=pow(cp[j1][0],i1-i2)*pow(cp[j1][1],i2)*cp[j1][2];
    }

  solveEquation(s, ss, m, a);
}

 
/**
 * find solution of a equation 
 */
void solveEquation(float s[][10], float ss[], int m, float a[])  
{
  int i1, i2, i3, ii;
  float sum, t[10][10], c;
 
  for(i1=0;i1<m;i1++) {
    for(i2=0;i2<m;i2++)
      t[i1][i2]=s[i1][i2];
    t[i1][m]=ss[i1];
  }
  
  for(i1=0;i1<m-1;i1++) {   
    c=t[i1][i1];
    ii=i1;
    for(i2=i1+1;i2<m;i2++) {
      if(t[i2][i1]>t[i1][i1])
	ii=i2; 
    }
    
    for(i3=i1;i3<=m;i3++) {
      c=t[i1][i3];
      t[i1][i3]=t[ii][i3];
      t[ii][i3]=c;
    }
    
    for(i2=i1+1;i2<m;i2++)
      for(i3=i1+1;i3<=m;i3++)
	t[i2][i3]=t[i2][i3]-t[i1][i3]*t[i2][i1]/t[i1][i1];
  }
 
  a[m-1]=t[m-1][m]/t[m-1][m-1];
  for(i1=m-2;i1>=0;i1--) {
    sum=t[i1][m];
    for(i2=m-1;i2>=i1+1;i2--)
      sum-=t[i1][i2]*a[i2];
    a[i1]=sum/t[i1][i1];
  } 
}




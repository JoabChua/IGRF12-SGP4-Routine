/*
 ============================================================================
 Name        : CProgTest.c
 Author      : Joab
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include<stdio.h>
#include<stdlib.h>
#include "math.h"
#include "string.h"
#include "sgp4ext.cpp"
#include "sgp4unit.cpp"
#include "sgp4io.cpp"
#define PI 3.141592654
#ifdef _WIN32
#include <io.h>
#endif
//
/* Private function prototypes -----------------------------------------------*/
void igrf(double magf[3], double dyear);
void getshc(double gh1[196], double gh2[196]);
double julday(int month, int day, int year);
void juliandate(int yr, int mon, int day, int hr, int min, double sec, double jdut);
void extrapsh(double dyear, double gh1[196], double gh2[196], double gha[196], double ghb[196]);
void shval3(double dlat, double dlon, int igdgc, double alt, double gha[196], double ghb[196], double ans[6], double dyear);
void dihf (double ans[6]);
void NED2Body (double yaw, double pitch, double roll, double X_NED, double Y_NED, double Z_NED);
void ECEF2Orb (double Joab_Lat, double Joab_Lon, double Joab_Angu, double Joab_X_ECEF, double Joab_Y_ECEF, double Joab_Z_ECEF);
void thetaGMST (int I, int J, int K, double UT);
void ECEF2NED (double Joab_Lat, double Joab_Lon, double Joab_X_ECEF, double Joab_Y_ECEF, double Joab_Z_ECEF);
void ECEF2ECI (double Joab_GMST, double Joab_X_ECEF, double Joab_Y_ECEF, double Joab_Z_ECEF);
void NED2ECEF (double Joab_Lat, double Joab_Lon, double ans[6]);
void ijk2lla(double r[3], double lla[3]);
void vectorprdt(double a[][3], double b[], double c[]);
void teme2ecef(double jdut[], double allr[][3]);
void ecef2llarout(double allr[][3], double lla[][3]);
void sgp4rout(double jdut[], double allr[][3], double allv[][3], int dd[], int mm[], int yy[]);
//
int main(void)
{

	  /* model and propagator codes */
	// 10 days simulation, hence the 11 rows array, 3 columns due to x, y and z axes.
	  double allr[11][3]={0}, allv[11][3]={0}, jdut[11]={0}, lla[11][3] = {0};
		int dd[11] = {0}, mm[11] = {0}, yy[11] = {0}, dyear[11]={0};
		double mfd[11][3] = {0}, magf[3] ={0};
		// sgp4 routine that runs the predefined TLE and output julian date, position and velocity vectors
		// and day, month, year for decimal year calculation needed for magnetic field model.
		sgp4rout(jdut,allr, allv,dd,mm,yy);
		// calculation for decimal year for 10 days period
		for (int c=0;c<11;c++)
		{
		dyear[c] = julday(mm[c],dd[c],yy[c]);
		}
		// transform TEME position vectors to ECEF
		teme2ecef(jdut,allr);
		// transform ECEF position vectors to LLA
		ecef2llarout(allr,lla);
	// LLA --> magnetic field data (NED) --> magnetic field data (ECEF) //
		for (int b=0;b<11;b++)
		{
		magf[0] = lla[b][0];
		magf[1] = lla[b][1];
		magf[2] = lla[b][2];
		// IGRF routine includes the transformation from NED -> ECEF frame output
		igrf(magf,dyear[b]);
		mfd[b][0] = magf[0];
		mfd[b][1] = magf[1];
		mfd[b][2] = magf[2];
		}
		// example of printing the LLA output for 10 days period.
	for (int ggg=0;ggg<11;ggg++)
		{
			for (int lll=0;lll<3;lll++)
			{
				printf("%f\t", lla[ggg][lll]);
			}
			printf("\n");
		}
	return 0;
}


/****************************************************************************/
/*                                                                          */
/*           Subroutine run all the magnetic field model routines           */
/*                                                                          */
/****************************************************************************/
void igrf(double magf[3], double dyear)
{
	double gha[196] = {0}, ghb[196] = {0}, ans[6] = {0};
	double gh1[196] = {0}, gh2[196] = {0}, dlat=magf[0],dlon=magf[1],alt=magf[2];
	int igdgc = 1;

	// load spherical harmoinc coefficients into array gh
	getshc(gh1, gh2);
	// extrapolate to get new g and h coefficients
	extrapsh(dyear, gh1, gh2, gha, ghb);
	// get x y z coordinates
	shval3(dlat, dlon, igdgc, alt, gha, ghb, ans, dyear);
	//printf("%f\t", ans[0]);printf("%f\t", ans[1]);printf("%f\t",ans[2]);printf("\n");
	NED2ECEF(dlat,dlon,ans);
	magf[0] = ans[0];magf[1] = ans[1];magf[2] = ans[2];

	// get total and horizontal intensities, declination and inclination
	dihf(ans);
	return;
}

/****************************************************************************/
/*                                                                          */
/*             Subroutine get spherical harmonic coordinates                */
/*                                                                          */
/****************************************************************************/
void getshc(double gh1[196], double gh2[196])
{
	double gg[196][2] = {
						{0,0,},
						{-29442.0,10.3,},
						{-1501.0,18.1,},
						{4797.1,-26.6,},
						{-2445.1,-8.7,},
						{3012.9,-3.3,},
						{-2845.6,-27.4,},
						{1676.7,2.1,},
						{-641.9,-14.1,},
						{1350.7,3.4,},
						{-2352.3,-5.5,},
						{-115.3,8.2,},
						{1225.6,-0.7,},
						{244.9,-0.4,},
						{582.0,-10.1,},
						{-538.4,1.8,},
						{907.6,-0.7,},
						{813.7,0.2,},
						{283.3,-1.3,},
						{120.4,-9.1,},
						{-188.7,5.3,},
						{-334.9,4.1,},
						{180.9,2.9,},
						{70.4,-4.3,},
						{-329.5,-5.2,},
						{-232.6,-0.2,},
						{360.1,0.5,},
						{47.3,0.6,},
						{192.4,-1.3,},
						{197.0,1.7,},
						{-140.9,-0.1,},
						{-119.3,-1.2,},
						{-157.5,1.4,},
						{16.0,3.4,},
						{4.1,3.9,},
						{100.2,0.0,},
						{70.0,-0.3,},
						{67.7,-0.1,},
						{-20.8,0.0,},
						{72.7,-0.7,},
						{33.2,-2.1,},
						{-129.9,2.1,},
						{58.9,-0.7,},
						{-28.9,-1.2,},
						{-66.7,0.2,},
						{13.2,0.3,},
						{7.3,0.9,},
						{-70.9,1.6,},
						{62.6,1.0,},
						{81.6,0.3,},
						{-76.1,-0.2,},
						{-54.1,0.8,},
						{-6.8,-0.5,},
						{-19.5,0.4,},
						{51.8,1.3,},
						{5.7,-0.2,},
						{15.0,0.1,},
						{24.4,-0.3,},
						{9.4,-0.6,},
						{3.4,-0.6,},
						{-2.8,-0.8,},
						{-27.4,0.1,},
						{6.8,0.2,},
						{-2.2,-0.2,},
						{24.2,0.2,},
						{8.8,0.0,},
						{10.1,-0.3,},
						{-16.9,-0.6,},
						{-18.3,0.3,},
						{-3.2,0.5,},
						{13.3,0.1,},
						{-20.6,-0.2,},
						{-14.6,0.5,},
						{13.4,0.4,},
						{16.2,-0.2,},
						{11.7,0.1,},
						{5.7,-0.3,},
						{-15.9,-0.4,},
						{-9.1,0.3,},
						{-2.0,0.3,},
						{2.1,0.0,},
						{5.4,0.0,},
						{8.8,0.0,},
						{-21.6,0.0,},
						{3.1,0.0,},
						{10.8,0.0,},
						{-3.3,0.0,},
						{11.8,0.0,},
						{0.7,0.0,},
						{-6.8,0.0,},
						{-13.3,0.0,},
						{-6.9,0.0,},
						{-0.1,0.0,},
						{7.8,0.0,},
						{8.7,0.0,},
						{1.0,0.0,},
						{-9.1,0.0,},
						{-4.0,0.0,},
						{-10.5,0.0,},
						{8.4,0.0,},
						{-1.9,0.0,},
						{-6.3,0.0,},
						{3.2,0.0,},
						{0.1,0.0,},
						{-0.4,0.0,},
						{0.5,0.0,},
						{4.6,0.0,},
						{-0.5,0.0,},
						{4.4,0.0,},
						{1.8,0.0,},
						{-7.9,0.0,},
						{-0.7,0.0,},
						{-0.6,0.0,},
						{2.1,0.0,},
						{-4.2,0.0,},
						{2.4,0.0,},
						{-2.8,0.0,},
						{-1.8,0.0,},
						{-1.2,0.0,},
						{-3.6,0.0,},
						{-8.7,0.0,},
						{3.1,0.0,},
						{-1.5,0.0,},
						{-0.1,0.0,},
						{-2.3,0.0,},
						{2.0,0.0,},
						{2.0,0.0,},
						{-0.7,0.0,},
						{-0.8,0.0,},
						{-1.1,0.0,},
						{0.6,0.0,},
						{0.8,0.0,},
						{-0.7,0.0,},
						{-0.2,0.0,},
						{0.2,0.0,},
						{-2.2,0.0,},
						{1.7,0.0,},
						{-1.4,0.0,},
						{-0.2,0.0,},
						{-2.5,0.0,},
						{0.4,0.0,},
						{-2.0,0.0,},
						{3.5,0.0,},
						{-2.4,0.0,},
						{-1.9,0.0,},
						{-0.2,0.0,},
						{-1.1,0.0,},
						{0.4,0.0,},
						{0.4,0.0,},
						{1.2,0.0,},
						{1.9,0.0,},
						{-0.8,0.0,},
						{-2.2,0.0,},
						{0.9,0.0,},
						{0.3,0.0,},
						{0.1,0.0,},
						{0.7,0.0,},
						{0.5,0.0,},
						{-0.1,0.0,},
						{-0.3,0.0,},
						{0.3,0.0,},
						{-0.4,0.0,},
						{0.2,0.0,},
						{0.2,0.0,},
						{-0.9,0.0,},
						{-0.9,0.0,},
						{-0.1,0.0,},
						{0.0,0.0,},
						{0.7,0.0,},
						{0.0,0.0,},
						{-0.9,0.0,},
						{-0.9,0.0,},
						{0.4,0.0,},
						{0.4,0.0,},
						{0.5,0.0,},
						{1.6,0.0,},
						{-0.5,0.0,},
						{-0.5,0.0,},
						{1.0,0.0,},
						{-1.2,0.0,},
						{-0.2,0.0,},
						{-0.1,0.0,},
						{0.8,0.0,},
						{0.4,0.0,},
						{-0.1,0.0,},
						{-0.1,0.0,},
						{0.3,0.0,},
						{0.4,0.0,},
						{0.1,0.0,},
						{0.5,0.0,},
						{0.5,0.0,},
						{-0.3,0.0,},
						{-0.4,0.0,},
						{-0.4,0.0,},
						{-0.3,0.0,},
						{-0.8,0.0,}};

			  int ii=0,jj=0;

			      for ( jj = 0; jj < 2; jj++)
			        {
			          for (ii = 0; ii <= 195; ii++)
			          {
			        	  if (jj==0)
			        	  {
			        	  gh1[ii] = gg[ii][jj];
			        	  }
			        	  else
			        	  {
			        		  gh2[ii] = gg[ii][jj];
			        	  }
			          }
			        }
			      return;
}

/****************************************************************************/
/*                                                                          */
/*                    Subroutine dec day of year                            */
/*                                                                          */
/****************************************************************************/

double julday(int month, int day, int year)
{

			  int days[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

			  int leap_year = (((year % 4) == 0) &&
			                   (((year % 100) != 0) || ((year % 400) == 0)));

			  double day_in_year = (days[month - 1] + day + (month > 2 ? leap_year : 0));

			  return ((double)year + (day_in_year / (365.0 + leap_year)));
}

/****************************************************************************/
/*                                                                          */
/*                         Subroutine julian day                            */
/*                                                                          */
/****************************************************************************/

void juliandate(int yr, int mon, int day, int hr, int min, double sec, double jdut)
{
	double jd, jdfrac;
	jd = 367.0 * yr - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )
			+ floor( 275 * mon / 9.0 ) + day + 1721013.5;
	//% use - 678987.0 to go to mjd directly
	jdfrac = (sec + min * 60.0 + hr *3600.0) / 86400.0;

	// check jdfrac
	if (jdfrac > 1.0)
{
    jd = jd + floor(jdfrac);
    jdfrac = jdfrac - floor(jdfrac);
}
	jdut = jd + jdfrac;
	return;
}

/****************************************************************************/
/*                                                                          */
/*              Subroutine extrapolate spherical harmonics                  */
/*                                                                          */
/****************************************************************************/

void extrapsh(double dyear, double gh1[196], double gh2[196], double gha[196], double ghb[196])
{
			double dte1 = 2015;
			int   nmax1=13;
			int   nmax2=13;
			  int   k=0;
			  int   i=0;
			  double factor=0, factor2=0;

			  factor = dyear - dte1;
			  factor2 = (dyear+1) - dte1;
			  k = nmax1 * (nmax2 + 2);
			  for ( i = 1; i <= k; i++)
			      {
			       gha[i] = gh1[i] + (factor * gh2[i]);
			      }
			  for ( i = 1; i <= k; i++)
			      {
			      ghb[i] = gh1[i] + (factor2 * gh2[i]);
			      }
			  return;
}

/****************************************************************************/
/*                                                                          */
/*                  Subroutine mag_field calculation                        */
/*                                                                          */
/****************************************************************************/

void shval3(double dlat, double dlon, int igdgc, double alt, double gha[196], double ghb[196],double ans[6], double dyear)
{
			int   nmax=13;
			// WGS84
			  double Re = 6371.2; // radius of earth
			  double a2 = 40680631.59;            /* WGS84 */
			  double b2 = 40408299.98;            /* WGS84 */
			// preallocate vectors for calculations
			  double sl[14] = {0};
			  double cl[14] = {0};
			  double p[119] = {0};
			  double q[119] = {0};
			//  double gtg = 0.9933;
			  double dtr = 0.01745329; // PI / 180 == change from deg to radian
			  double slat; // sine of latitude
			  double clat; // cosine of latitude
			  double ratio=0;
			  double aa=0, bb=0, cc=0, dd=0;
			  double sd=0.0;
			  double cd=1.0;
			  double r=0;
			  double rr=0;
			  double fm=0,fn=0;
			  int ii,j,k,l,m,n;
			  int npq=0;
			  double argument=0;
			  double power=0;

			  // initialize outputs
			  double x = 0, y = 0, z = 0, xtemp = 0, ytemp = 0, ztemp = 0;

			  // initialize counters
			  l = 1;
			  n = 0;
			  m = 1;

			// convert geodetic latitude to geocentric latitude
			//  if (igdgc == 2)
			//  {
			//	  dlat = atan2(gtg*tan(dlat*dtr),1);
			//	  dlat = dlat / dtr;
			//  }
			  r = alt;
			  argument = dlat * dtr;
			  slat = sin( argument );
			  if ((90.0 - dlat) < 0.001)
			    {
			      aa = 89.999;            /*  300 ft. from North pole  */
			    }
			  else
			    {
			      if ((90.0 + dlat) < 0.001)
			        {
			          aa = -89.999;        /*  300 ft. from South pole  */
			        }
			      else
			        {
			          aa = dlat;
			        }
			    }
			  argument = aa * dtr;
			  clat = cos( argument );
			  argument = dlon * dtr;
			  sl[1] = sin( argument );
			  cl[1] = cos( argument );
			  npq = (nmax * (nmax + 3)) / 2; // equals to 104.

//  convert geodetic to geocentric
			  if (igdgc ==1)
			  {
			  aa = a2 * clat * clat;
			  bb = b2 * slat * slat;
			  cc = aa + bb;
			  dd = sqrt( cc );
			  argument = alt * (alt + 2.0 * dd) + (a2 * aa + b2 * bb) / cc;
			  r = sqrt( argument );
			  cd = (alt + dd) / r;
			  sd = (a2 - b2) / dd * slat * clat / r;
			  argument = slat;
			  slat = slat * cd - clat * sd;
			  clat = clat * cd + argument * sd;
			  }

			  ratio = Re / r;
			  argument = sqrt( 3.0 );

			  p[1] = 2.0 * slat;
			  p[2] = 2.0 * clat;
			  p[3] = 4.5 * slat * slat - 1.5;
			  p[4] = 3.0 * argument * clat * slat;
			  q[1] = -clat;
			  q[2] = slat;
			  q[3] = -3.0 * clat * slat;
			  q[4] = argument * (slat * slat - clat * clat);

			  for ( k = 1; k <= npq; k++)
			    {
			      if (n < m)
			        {
			          m = 0;
			          n = n + 1;
			          power =  n + 2;
			          rr = pow(ratio,power);
			          fn = n;
			        }
			      fm = m;
			      if (k >= 5)
			        {
			          if (m == n)
			            {
			              aa = sqrt(1.0 - 0.5 / fm);
			              j = k - n - 1;
			              p[k] = (1.0 + 1.0/fm) * aa * clat * p[j];
			              q[k] = aa * (clat * q[j] + slat/fm * p[j]);
			              sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
			              cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
			            }
			          else
			            {
			              aa = sqrt( fn*fn - fm*fm );
			              bb = sqrt( ((fn - 1.0)*(fn-1.0)) - (fm * fm) ) / aa;
			              cc = (2.0 * fn - 1.0)/aa;
			              ii = k - n;
			              j = k - 2 * n + 1;
			              p[k] = (fn + 1.0) * (cc * slat/fn * p[ii] - bb/(fn - 1.0) * p[j]);
			              q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j];
			            }
			        }
			      aa = rr * gha[l];
			      double aatemp = rr * ghb[l];

			      if (m == 0)
			        {
			  	  	  x = x + aa * q[k];
			          z = z - aa * p[k];
			          xtemp = xtemp + aatemp * q[k];
			          ztemp = ztemp - aatemp * p[k];
			          l = l + 1;
			        }
			      else
			        {
			    	  bb = rr * gha[l+1];
			          cc = aa * cl[m] + bb * sl[m];
			          x = x + cc * q[k];
			          z = z - cc * p[k];

			      double bbtemp = rr * ghb[l+1];
			      double cctemp = aatemp * cl[m] + bbtemp * sl[m];
			      xtemp = xtemp + cctemp * q[k];
			      ztemp = ztemp - cctemp * p[k];

			      if (clat > 0)
			       {
			        y = y + (aa * sl[m] - bb * cl[m]) * fm * p[k]/((fn + 1.0) * clat);
			        ytemp = ytemp + (aatemp * sl[m] - bbtemp *cl[m]) * fm * p[k] / ((fn + 1.0) * clat);
			       }
			       else
			        {
			    	  y = y + (aa * sl[m] - bb * cl[m]) * q[k] * slat;
			          ytemp = ytemp + (aatemp * sl[m] - bbtemp * cl[m]) * q[k] * slat;
			        }
			              l = l + 2;
			        }
			      m = m + 1;
			    }

			  	  double xold = x;
			  	  x = x * cd + z * sd;
			  	  z = z * cd - xold * sd;

			  	  double xtempold = xtemp;
			  	  xtemp = xtemp * cd + ztemp * sd;
			  	  ztemp = ztemp * cd - xtempold * sd;

			  	  ans[0] = x; ans[1] = y; ans[2] = z;
			  	  ans[3] = xtemp; ans[4] = ytemp; ans[5] = ztemp;
			  	  // to push magnetic field answer back to igrf //
			  	  dlat = x; dlon = y; alt = z;
			  return;
}

/****************************************************************************/
/*                                                                          */
/*      Subroutine declination inclination horizontal total-force           */
/*                                                                          */
/****************************************************************************/

void dihf (double ans[6])
{
			  double d=0,i=0,dtemp=0,htemp=0,itemp=0,ftemp=0;
			  double h2,h,f;
			  double argument, argument2;
			  double sn = 0.0001;
			  double x = ans[0], y = ans[1], z = ans[2], xtemp = ans[3], ytemp = ans[4], ztemp = ans[5];

			  h2 = x*x + y*y;
			  h = sqrt( h2 );       /* calculate horizontal intensity */
			  f = sqrt(h2 + z*z);      /* calculate total intensity */
			  if (f < sn)
			     {
			      d = 0;        /* If d and i cannot be determined, */
			      i = 0;        /*       set equal to NaN         */
			     }
			  else
			  {
			      i = atan2(z,h);
			   if (h < sn)
			     {
			      d = 0;
			     }
			   else
			    {
			   if ((h+x) < sn)
			     {
			      d = PI;
			     }
			   else
			     {
			      d = 2.0 * atan2(y,(h+x));
			     }
			    }
			  }

			  h2 = xtemp*xtemp + ytemp*ytemp;
			  htemp = sqrt(h2);
			  ftemp = sqrt(h2 + ztemp*ztemp);
			  if (ftemp < sn)
			  {
			  dtemp = 0;    /* If d and i cannot be determined, */
			  itemp = 0;    /*       set equal to 999.0         */
			  }
			  else
			  {
			  argument = ztemp;
			  argument2 = htemp;
			  itemp = atan2(argument,argument2);
			  if (htemp < sn)
			    {
			    dtemp = 0;
			    }
			   else
			    {
			   if ((htemp + xtemp) < sn)
			    {
			    dtemp = PI;
			    }
			     else
			    {
			    argument = ytemp;
			    argument2 = htemp + xtemp;
			    dtemp = 2.0 * atan2(argument,argument2);
			    }
			    }
			  }
//			  printf("%.2f       " ,h);
//			  printf("%.2f         " ,f);
//			  printf("%.2f         " ,(i*180/PI) + itemp - itemp);
//			  printf("%.2f         \n" ,(d*180/PI) + dtemp - dtemp);
//  cout << "Horizontal Intensity 2: " << htemp << endl;
//  cout << "Total Intensity 2: " << ftemp << endl;
//  cout << "Inclination 2: " << itemp << endl;
//  cout << "declination 2: " << dtemp << endl;
			  return;
}

void ECEF2ECI (double Joab_GMST, double Joab_X_ECEF, double Joab_Y_ECEF, double Joab_Z_ECEF)
{
	int i,j,k;
	double eci[1][3]={0};

	Joab_GMST = Joab_GMST * PI / 180.0;
	double ecef[1][3] = {Joab_X_ECEF, Joab_Y_ECEF, Joab_Z_ECEF};
	double ecef2eci[3][3] = { {cos(Joab_GMST), -1*sin(Joab_GMST), 0},
							{sin(Joab_GMST), cos(Joab_GMST), 0},
							{0, 0, 1} };

// multiplying ECEF to transformation matrix
	for(i=0; i<1; ++i)
	    for(j=0; j<3; j++)
	    for(k=0; k<3; k++)
	    {
	        eci[i][j]+=ecef[i][k]*ecef2eci[k][j];
	    }

// print out the ECI coordinates
//		for(i=0;i<1;i++)
//		{
//	    for(j=0; j<3; j++)
//	    {
//	        printf("%.2f	\n", eci[i][j]);
//	    }
//		}
		return;
}

void NED2ECEF (double Joab_Lat, double Joab_Lon, double ans[6])
{
	int i,j,k;
	double ned[1][3]={0};
	double Joab_X_ECEF=0,Joab_Y_ECEF=0,Joab_Z_ECEF=0;
	Joab_X_ECEF=ans[0];Joab_Y_ECEF=ans[1];Joab_Z_ECEF=ans[2];
	Joab_Lat = Joab_Lat * PI / 180;
	Joab_Lon = Joab_Lon * PI / 180;
	double ecef[1][3] = {Joab_X_ECEF, Joab_Y_ECEF, Joab_Z_ECEF};
	double ecef2ned[3][3] = {{-1*cos (Joab_Lon)*sin(Joab_Lat), -1*sin(Joab_Lon)*sin(Joab_Lat), cos(Joab_Lat)},
							{-1*sin (Joab_Lon), cos(Joab_Lon), 0},
							{-1*cos(Joab_Lon)*cos(Joab_Lat), -1*sin(Joab_Lon)*cos(Joab_Lat), -1*sin(Joab_Lat)}};

// multiplying ECEF to transformation matrix
	for(i=0; i<1; ++i)
	{
	    for(j=0; j<3; j++)
	    for(k=0; k<3; k++)
	    {
	        ned[i][j]+=ecef[i][k]*ecef2ned[k][j];
	    }
	}
// print out the NED coordinates
//		for(i=0;i<1;i++)
//		{
//	    for(j=0; j<3; j++)
//	    {
//	    	printf("%.2f   ", ned[i][j]);
//	    }
//		}
//		printf("\n");
		ans[0] = ned[0][0];ans[1] = ned[0][1];ans[2] = ned[0][2];
		return;
}

void thetaGMST (int I, int J, int K, double UT)
{
	double JDO,JD, JD1, JD11, JD12, JD13, JD2, JD21, TT, GMST;

// Calculation for Julian Date
   	JDO = (100 * I) + J - 190002.5;
   	if (JDO < 0)
   		{
   		JDO = -1;
	   	}
	   	else JDO = 1;
	JD11 = (J+9)/12;
	JD12 = (long)JD11;
   	JD13 = ((JD12 + I) * 7) / 4;
   	JD1 = (long)JD13;
  	JD21 = (275*J) / 9;
   	JD2 = (long)JD21;
   	JD = (367 * I) - JD1 + JD2 + K + 1721013.5 + (UT/24) - (0.5* JDO) + 0.5;
   	TT = (JD - 2451545.0) / 36525.0;
   	GMST = (24110.54841 + (8640184.812866 * TT) + (0.093104 * pow(TT,2)) - (0.0000062 * pow(TT,3)) + (1.002737909350795 * 3600 * UT))/240;
//
//   	printf("%.2f	\n", JD);
//   	printf("%.2f	\n", GMST);

   	return;
}

void ECEF2Orb (double Joab_Lat, double Joab_Lon, double Joab_Angu, double Joab_X_ECEF, double Joab_Y_ECEF, double Joab_Z_ECEF)
{
	int i,j,k; // for looping use
	double orb[1][3]={0};
	double ecef[1][3] = {Joab_X_ECEF, Joab_Y_ECEF, Joab_Z_ECEF};

	// transformation matrix
	double ecef2orb[3][3] = {{(-1*cos (Joab_Lon)*sin(Joab_Lat)*cos(Joab_Angu)) - (sin(Joab_Lon)*sin(Joab_Angu)), (cos(Joab_Lon)*sin(Joab_Angu)) - (sin (Joab_Lon)*sin(Joab_Lat)*cos(Joab_Angu)), cos(Joab_Lat)*cos(Joab_Angu)},
							{(cos(Joab_Lon)*sin(Joab_Angu)*sin(Joab_Lat)) - (sin(Joab_Lon)*cos(Joab_Angu)), (cos(Joab_Lon)*cos(Joab_Angu)) + (sin(Joab_Lon)*sin(Joab_Lat)*sin(Joab_Angu)), -1*sin(Joab_Angu)*cos(Joab_Lat)},
							{-1*cos(Joab_Lat)*cos(Joab_Lon), -1*sin(Joab_Lon)*cos(Joab_Lat), -1*sin(Joab_Lat)}};

// multiplying ECEF to transformation matrix
	for(i=0; i<1; ++i)
	    for(j=0; j<3; j++)
	    for(k=0; k<3; k++)
	    {
	        orb[i][j]+=ecef[i][k]*ecef2orb[k][j];
	    }

// print out the Orbit frame coordinates
//		for(i=0;i<1;i++)
//		{
//	    for(j=0; j<3; j++)
//	    {
//	    	printf("%.2f	\n", orb[i][j]);
//	    }
//		}
		return;
}

void NED2Body (double yaw, double pitch, double roll, double X_NED, double Y_NED, double Z_NED)
{
	int i,j,k;
	double body[1][3]={0};

	double ned[1][3] = {X_NED, Y_NED, Z_NED};
	double ned2Body[3][3] = {{cos(pitch)*cos(yaw), cos(pitch)*sin(yaw), -1*sin(pitch)},
							{(sin(roll)*sin(pitch)*cos(yaw)) - (cos(roll)*sin(yaw)), (sin(roll)*sin(pitch)*sin(yaw)) + (cos(roll)*cos(roll)), sin(roll)*cos(pitch)},
							{(cos(roll)*sin(pitch)*cos(yaw)) + (sin(roll)*sin(yaw)), (cos(roll)*sin(pitch)*sin(yaw)) - (sin(roll)*cos(yaw)), cos(roll)*cos(pitch)}};

// multiplying NED to transformation matrix
	for(i=0; i<1; ++i)
	    for(j=0; j<3; j++)
	    for(k=0; k<3; k++)
	    {
	        body[i][j]+=ned[i][k]*ned2Body[k][j];
	    }

// print out the Body Frame coordinates
//		for(i=0;i<1;i++)
//		{
//	    for(j=0; j<3; j++)
//	    {
//	    	printf("%.2f	\n", body[i][j]);
//	    }
//		}
		return;
}

void ECEF2NED (double Joab_Lat, double Joab_Lon, double Joab_X_ECEF, double Joab_Y_ECEF, double Joab_Z_ECEF)
{
	int i,j,k;
	double ned[1][3]={0};

	Joab_Lat = Joab_Lat * PI / 180;
	Joab_Lon = Joab_Lon * PI / 180;
	double ecef[1][3] = {Joab_X_ECEF, Joab_Y_ECEF, Joab_Z_ECEF};
	double ecef2ned[3][3] = {{-1*cos (Joab_Lon)*sin(Joab_Lat), -1*sin (Joab_Lon), -1*cos(Joab_Lon)*cos(Joab_Lat)},
							{-1*sin(Joab_Lon)*sin(Joab_Lat), cos(Joab_Lon), -1*sin(Joab_Lon)*cos(Joab_Lat)},
							{cos(Joab_Lat), 0, -1*sin(Joab_Lat)}};

// multiplying NED to transformation matrix
	for(i=0; i<1; ++i)
	    for(j=0; j<3; j++)
	    for(k=0; k<3; k++)
	    {
	        ned[i][j]+=ecef[i][k]*ecef2ned[k][j];
	    }

// print out the ECEF coordinates
//		for(i=0;i<1;i++)
//		{
//	    for(j=0; j<3; j++)
//	    {
//	    	printf("%.2f	\t", ned[i][j]);
//	    }
//		}
//		printf("\n");
		return;
}
/****************************************************************************/
/*                                                                          */
/*          Subroutine for ECEF to LLA Operations                           */
/*                                                                          */
/****************************************************************************/
void ecef2llarout(double allr[][3], double lla[][3])
{
	int i=0,j=0;
	double r[3] = {0}, temp[3] = {0};
	for (i=0;i<12;i++)
	{
		for (j=0;j<3;j++)
		{
			r[j] = allr[i][j];
		}
		j=0;
		ijk2lla(r,temp);
		//printf("%f\t", lla[0]);printf("%f\t", lla[1]);printf("%f\t", lla[2]);printf("\n");
		for (j=0;j<3;j++)
		{
			lla[i][j] = temp[j];
		}
		j=0;
	}
	return;
}

/****************************************************************************/
/*                                                                          */
/*           Subroutine transform TEME vector to ECEF vectors               */
/*                                                                          */
/****************************************************************************/
void teme2ecef(double jdut[], double allr[][3])
{
	double gmst,st[3][3]={0},pm[3][3]={0},r[3]={0};
    //thetasa    = 7.29211514670698e-05 * (1.0  - lod/86400.0 );
    //omegaearth = thetasa;
	double cosxp,sinxp,cosyp,sinyp;
	cosxp = cos(0);
	sinxp = sin(0);
	cosyp = cos(0);
	sinyp = sin(0);
    pm[0][0] =  cosxp;
    pm[0][1] =  sinxp * sinyp;
    pm[0][2] =  sinxp * cosyp;
    pm[1][0] =  0.0;
    pm[1][1] =  cosyp;
    pm[1][2] = -sinyp;
    pm[2][0] = -sinxp;
    pm[2][1] =  cosxp * sinyp;
    pm[2][2] =  cosxp * cosyp;
    for (int j=0;j<12;j++)
    {
    for (int i=0;i<3;i++)
    {
    	r[i] = allr[j][i];
    }
    gmst= gstime(jdut[j]);
    //printf("%f\t", gmst);
    //ttt = (jdut[j] - 2451545.0) / 36525;
    st[0][0] =  cos(gmst);
    st[0][1] =  sin(gmst);
    st[0][2] =  0.0;
    st[1][0] = -sin(gmst);
    st[1][1] =  cos(gmst);
    st[1][2] =  0.0;
    st[2][0] =  0.0;
    st[2][1] =  0.0;
    st[2][2] =  1.0;
    //thetasa= 7.29211514670698e-05 * (1.0  - lod/86400.0 );
    //omegaearth[3] = {0, 0, thetasa};
    double rpef[3] = {0};
    vectorprdt(st,r,rpef);
    double rr[3] = {0};
    vectorprdt(pm,rpef,rr);
    for (int k=0;k<3;k++)
    {
    	allr[j][k] = rr[k];
    }
    }
//    rpef[][]  = st*rteme;
//    recef[][] = pm*rpef;
//    vpef[][]  = st*vteme - cross( omegaearth,rpef );
//    vecef[][] = pm*vpef;
//    cross(omegaearth,rpef,temp);
//    aecef = pm*(st*ateme - cross(omegaearth,temp) - 2.0*cross(omegaearth,vpef));
	return;
}

void vectorprdt(double a[][3], double b[], double c[])
{
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
	c[i] = c[i] + (a[i][j] * b[j]);
		}
	}
	return;
}

/****************************************************************************/
/*                                                                          */
/*             Subroutine to transform ECEF to Lat Lon Alt                  */
/*                                                                          */
/****************************************************************************/
void ijk2lla(double r[3], double lla[3])
{
	double twopi,small,re,eesqrd,temp,magr,rtasc,lon,decl,latgd,latgc;
    twopi      =     2.0*PI;
    small      =     0.00000001;         // small value for tolerances
    re         =     6378.137;
    eesqrd     =     0.006694385000;     // eccentricity of earth sqrd

    // Implementation
    magr = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    if (abs( magr ) >= 1.0e-16)
    {
    	magr= sqrt( magr );
    }
    else
    {
    	magr= 0.0;
    }

    // find longitude value
    temp = sqrt( r[0]*r[0] + r[1]*r[1] );
    if ( abs( temp ) < small )
    {
    	rtasc= signbit(r[2])*PI*0.5;
    }
    else
    {
    	rtasc= atan2( r[1], r[0] );
    }
    lon  = rtasc;
    if ( abs(lon) >= PI )
    {
    	if ( lon < 0.0  )
    	{
    		lon= twopi + lon;
    	}
        else
        {
        	lon= lon - twopi;
        }
    }
    decl = asin( r[2] / magr );
    latgd= decl;

    // iterate to find geodetic latitiude
    int i= 1;
    double olddelta = latgd + 10.0, c,sintemp,s,hellp;

    while ((abs(olddelta-latgd)>=small) && (i<10))
    {
    	olddelta= latgd;
        sintemp = sin( latgd );
        c = re  / (sqrt( 1.0 -eesqrd*sintemp*sintemp ));
        latgd= atan( (r[2]+c*eesqrd*sintemp)/temp );
        i = i + 1;
    }

    // find altitude
    if ((PI*0.5 - abs(latgd)) > PI/180.0)  // 1 deg
    {
    	hellp   = (temp/cos(latgd)) - c;
    }
    else
    {
    	s = c * (1.0 - eesqrd);
        hellp   = r[2]/sin(latgd) - s;
    }

    // convert geodetic to geocentric
    latgc = atan((1.0-eesqrd) * tan(latgd));

    // convert rad to deg for lat and lon
    latgc = latgc / PI * 180;
    latgd = latgd / PI * 180;
    lon = lon / (twopi) * 360;
    lla[0] = latgd;
    lla[1] = lon;
    lla[2] = hellp;
    return;
}

/****************************************************************************/
/*                                                                          */
/*             Subroutine get position and velocity vectors                 */
/*                                                                          */
/****************************************************************************/
void sgp4rout(double jdut[], double allr[][3], double allv[][3], int dd[], int mm[], int yy[])
{
	double ro[3], vo[3];  int jj=1;
        char typerun, opsmode;
        gravconsttype  whichconst;


// ----------------------------  locals  -------------------------------
        double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;
	double sec,  jd, tsince, stopmfe, deltamin;
        double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;
	int  year; int mon; int day; int hr; int min;
	// pre-defined TLE for the simulation work
	char longstr1[130]="1 41171U 15077E   17056.21577948  .00003135  00000-0  13988-3 0  9999";
        typedef char str3[4];
        str3 monstr[13];
	char longstr2[130]="2 41171  14.9885 264.2121 0013786 281.2478  78.6293 15.11686817 66107";
	elsetrec satrec;

// ------------------------  implementation   --------------------------
strcpy(monstr[1], "Jan");
strcpy(monstr[2], "Feb");
strcpy(monstr[3], "Mar");
strcpy(monstr[4], "Apr");
strcpy(monstr[5], "May");
strcpy(monstr[6], "Jun");
strcpy(monstr[7], "Jul");
strcpy(monstr[8], "Aug");
strcpy(monstr[9], "Sep");
strcpy(monstr[10], "Oct");
strcpy(monstr[11], "Nov");
strcpy(monstr[12], "Dec");

        opsmode = 'a';
        typerun = 'v';
        whichconst = wgs72;
        getgravconst( whichconst, &tumin, &mu, &radiusearthkm, &xke, &j2, &j3, &j4, &j3oj2 );

int rr=0;
        // ----------------- test simple propagation -------------------
        while (rr==0)//feof(infile) == 0)
          {

            if (rr == 0)
              {
                // convert the char string to sgp4 elements
                // includes initialization of sgp4
                twoline2rv( longstr1, longstr2, typerun, opsmode, whichconst,
                            &satrec );
                //fprintf(outfile, "%ld xx\n", satrec.satnum);
                // call the propagator to get the initial state vector value
                sgp4 (whichconst, &satrec,  0.0, ro,  vo);
                jd = satrec.jdsatepoch;
                invjday( jd, &year,&mon,&day,&hr,&min, &sec );
                rv2coe(ro, vo, mu, &p, &a, &ecc, &incl, &node, &argp, &nu, &m, &arglat, &truelon, &lonper );

               //feeding ro and vo into the first row of allr and allv 2D arrays
               int ii=0;
               for (ii=0;ii<3;ii++)
               {
            	   allr[0][ii] = ro[ii];
            	   allv[0][ii] = vo[ii];
               }
               jdut[0] = jd;
               dd[0] = day;mm[0] = mon;yy[0]= year;
               // set interval and period of the propagation
                tsince = 0.0;
                stopmfe = 14400.0; // 10 days
                deltamin = 1440.0; // 1 day - 24 hr - 1440 minutes

                // check so the first value isn't written twice
                if ( fabs(tsince) > 1.0e-8 )
                    tsince = tsince - deltamin;
                // ----------------- loop to perform the propagation ----------------
                while ((tsince < stopmfe) && (satrec.error == 0))
                  {
                   tsince = tsince + deltamin;

                   if(tsince > stopmfe)
                       tsince = stopmfe;

                   sgp4 (whichconst, &satrec,  tsince, ro,  vo);


                   if (satrec.error == 0) // for manual insert of timeline.
                     {
                       if ((typerun != 'v') && (typerun != 'c'))
                         {
                            jd = satrec.jdsatepoch + tsince/1440.0;
                            invjday( jd, &year,&mon,&day,&hr,&min, &sec );
                         }
                       else
                         {
                            jd = satrec.jdsatepoch + tsince/1440.0;
                            invjday( jd, &year,&mon,&day,&hr,&min, &sec );

                            for (int kk=0;kk<3;kk++)
                            {
                         	   allr[jj][kk] = ro[kk];
                         	   allv[jj][kk] = vo[kk];
                            }
                           dd[jj] = day;mm[jj] = mon;yy[jj]= year;
                           jdut[jj] = jd;
                           jj = jj + 1;
                            rv2coe(ro, vo, mu, &p, &a, &ecc, &incl, &node, &argp, &nu, &m, &arglat, &truelon, &lonper );
                         }
                     } // if satrec.error == 0
                   rr=2;
                  } // while propagating the orbit
              } // if not eof
          } // while through the input file
  return;
}  // end

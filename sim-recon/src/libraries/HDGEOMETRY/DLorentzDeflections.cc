#include "DLorentzDeflections.h"
#include <math.h>
#include <stdlib.h>

// Locate a position in array xx given x
static void locate(const double *xx,int n,double x,int *j){
  int ju,jm,jl;
  int ascnd;
  
  jl=-1;
  ju=n;
  ascnd=(xx[n-1]>=xx[0]);
  while(ju-jl>1){
    jm=(ju+jl)>>1;
    if ( (x>=xx[jm])==ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x==xx[0]) *j=0;
  else if (x==xx[n-1]) *j=n-2;
  else *j=jl; 
}

// Polynomial interpolation on a grid.
// Adapted from Numerical Recipes in C (2nd Edition), pp. 121-122.
static void polint(const double *xa, const double *ya,int n,double x, double *y,
	    double *dy){
  int i,m,ns=0;
  double den,dif,dift,ho,hp,w;

  double *c=(double *)calloc(n,sizeof(double));
  double *d=(double *)calloc(n,sizeof(double));

  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++){
    if ((dift=fabs(x-xa[i]))<dif){
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];

  for (m=1;m<n;m++){
    for (i=1;i<=n-m;i++){
      ho=xa[i-1]-x;
      hp=xa[i+m-1]-x;
      w=c[i+1-1]-d[i-1];
      if ((den=ho-hp)==0.0){
	free(c);
	free(d);
	return;
      }
      den=w/den;
      d[i-1]=hp*den;
      c[i-1]=ho*den;
      
    }
    
    *y+=(*dy=(2*ns<(n-m) ?c[ns+1]:d[ns--]));
  }
  free(c);
  free(d);
}

// Obtain slope parameters describing Lorentz deflection by interpolating 
// on the Lorentz deflection table
jerror_t DLorentzDeflections::GetLorentzCorrectionParameters(double x,
							     double y,double z,
							     double &tanz, 
							     double &tanr) const{
  int imin,ind,ind2;
  double r=sqrt(x*x+y*y);
 
   // Locate positions in x and z arrays given r and z
  locate(lorentz_x,LORENTZ_X_POINTS,r,&ind);
  locate(lorentz_z,LORENTZ_Z_POINTS,z,&ind2);
  
  // First do interpolation in z direction 
  imin=PACKAGE_Z_POINTS*(ind2/PACKAGE_Z_POINTS); // Integer division...
  double ytemp[LORENTZ_X_POINTS],ytemp2[LORENTZ_X_POINTS],dy;
  for (int j=0;j<LORENTZ_X_POINTS;j++){
    polint(&lorentz_z[imin],&lorentz_nx[j][imin],PACKAGE_Z_POINTS,z,
	   &ytemp[j],&dy);
    polint(&lorentz_z[imin],&lorentz_nz[j][imin],PACKAGE_Z_POINTS,z,
	   &ytemp2[j],&dy);
  }
  // Then do final interpolation in x direction 
  if (ind>0){
    imin=((ind+3)>LORENTZ_X_POINTS)?(LORENTZ_X_POINTS-4):(ind-1);
  }
  else imin=0;
  polint(&lorentz_x[imin],&ytemp[imin],4,r,&tanr,&dy);
  polint(&lorentz_x[imin],&ytemp2[imin],4,r,&tanz,&dy);
  
  return NOERROR;
}

// Obtain deflection of avalanche along wire due to Lorentz force
double DLorentzDeflections::GetLorentzCorrection(double x,double y,double z,
						 double alpha,double dx) const
{
  double tanz=0.,tanr=0.;
  GetLorentzCorrectionParameters(x,y,z,tanz,tanr);

  // Deflection along wire	
  double phi=atan2(y,x);
  return (-tanz*dx*sin(alpha)*cos(phi)+tanr*dx*cos(alpha));
}

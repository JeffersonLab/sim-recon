#include "DMaterialMap.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Locate a position in array xx given x
static void locate(const double *xx,int n,double x,int *j){
  int ju,jm,jl;
  int ascnd;
  
  jl=-1;
  ju=n;
  ascnd=(xx[n-1]>=xx[0]);
  while(ju-jl>1){
    jm=(ju+jl)>>1;
    if (x>=xx[jm]==ascnd)
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

// Routine to obtain the average radiation length of the material at the 
// position pos from a map of material obtained from the simulation
double DMaterialMap::GetRadLen(double x, double y, double z) const{
  double r=sqrt(x*x+y*y);
  int ind,ind2,imin;
  
  // locate positions in r and z arrays
  locate(material_z,NUM_Z_POINTS,z,&ind);
  locate(material_x,NUM_X_POINTS,r,&ind2);

  // First do the interpolation in the z direction
  double temp[NUM_X_POINTS],err;
  if (ind>0){
    imin=((ind+3)>NUM_Z_POINTS)?NUM_Z_POINTS-4:(ind-1);
  }
  else imin=0;
  for (int j=0;j<NUM_X_POINTS;j++){
    polint(&material_z[imin],&radlen[imin][j],4,z,&temp[j],&err);
  }
  // next do interpolation in r direction
  if (ind2>0){
    imin=((ind2+3)>NUM_X_POINTS)?NUM_X_POINTS-4:(ind2-1);
  }
  else imin=0;
  double fX0;
  polint(&material_x[imin],&temp[imin],4,r,&fX0,&err);

  return 1.e5/fX0;

}

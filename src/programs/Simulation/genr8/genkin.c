
 /*********************
 *                        
 *  kinematics.c            
 *                    
 **********************
 *
 * */

#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include <genkin.h>

/*
#include <itypes.h>
#include <hdr.h>
#include <param.h>
#include <disData.h>
#include <disIO.h>
#include <dataIO.h>
#include <utility.h>
#include <esr.h>
#include <omega_system.h>
#include <particleType.h>
#include <eventType.h>

*/

#define RESTFRAME -1
#define PARENTFRAME +1
 
double SQ(double x){
  double z;

  z = (x)*(x);
  return (z);
}

/**************
 * DotProduct.c
 **************
 */
double DotProduct3(const vector3_t *p1, const vector3_t *p2)
{
    return(p1->x * p2->x + p1->y * p2->y + p1->z * p2->z);
}
  
/****************
 * CrossProduct.c
 ****************
 */
vector3_t CrossProduct3(const vector3_t *p1,const vector3_t *p2)
{
  vector3_t c;

  c.x = p1->y*p2->z - p1->z*p2->y;
  c.y = -(p1->x*p2->z - p1->z*p2->x);
  c.z = p1->x*p2->y - p1->y*p2->x;

  return c;
}


/*******************
 *   get_beta()
 *
 *******************/
vector4_t get_beta(vector4_t *boost,int sign){
  
  /* find beta 4vector where beta->t = gamma
   * 
   */
  vector4_t beta;
 

  beta.space.x = sign*(boost->space.x / boost->t);
  beta.space.y = sign*(boost->space.y / boost->t);
  beta.space.z = sign*(boost->space.z / boost->t);
  /* gamma = E/m */
  
  beta.t = (boost->t) / sqrt( (boost->t) * (boost->t) - 
			  ( (boost->space.x) * (boost->space.x) +
			    (boost->space.y) * (boost->space.y) +
			   (boost->space.z) * (boost->space.z) ));
  return beta;
}

/***********
 * lorentz.c
 ***********
 */

vector4_t lorentz(const vector4_t *beta,const vector4_t *pin)
{
    vector4_t ret;
    double d,c,c2;
    d = DotProduct3(&(beta->space),&(pin->space));
    c = beta->t/(beta->t + 1.0);
    c2 = c * d + pin->t;
    ret.space.x = pin->space.x + beta->space.x * beta->t * c2;
    ret.space.y = pin->space.y + beta->space.y * beta->t * c2;
    ret.space.z = pin->space.z + beta->space.z * beta->t * c2;
    ret.t = beta->t * (pin->t + d);
    return(ret);
}
/**************
 * CMmomentum.c
 **************/

double CMmomentum(double cm_engy, double m1, double m2)
{
  double A,B,C,D,E;

  A = cm_engy*cm_engy;
  B = (m1 + m2)*(m1 + m2);
  C = (m1 - m2)*(m1-m2);
  D = (A - B) * (A - C);
  E = sqrt (D) / (2.0 * cm_engy);
  return (E);

} 
/***********
 * energy.c
 **********/

double energy(double mass, const vector3_t *p)
{
  double E,pmagsq;
  pmagsq = SQ(p->x)+SQ(p->y)+SQ(p->z);
  E = sqrt( mass*mass + pmagsq);
  return(E);
} 


/********
 * v3mag.c
 ********/

double v3mag(const vector3_t *p)
{
  double mag;
  mag = sqrt( (p->x)*(p->x) +
	     (p->y)*(p->y)+
	     (p->z)*(p->z) );
  return(mag);
} 

/*
 ***********************
 *                     *
 *  Sum4vec()   *
 *                     *
 ***********************
 */


vector4_t Sum4vec(vector4_t *vec4, int nvec4)
{
  int i;
  vector4_t temp4;
  
  temp4.space.x=0;
  temp4.space.y=0;
  temp4.space.z=0;
  temp4.t=0;

  for(i=0;i<nvec4;i++){
    temp4.space.x += (vec4 +i)->space.x;
    temp4.space.y += (vec4 +i)->space.y;
    temp4.space.z += (vec4 +i)->space.z;
    temp4.t += (vec4 +i)->t;
  }
  return temp4;
}

/************************
 *
 *  eff_mass()
 *
 ************************/
double eff_mass(vector4_t *v, int nparticles)
{
  int i;
  double mass=-1; /* set to -1 for debugging */
  vector4_t vsum;

  /*
   * initilize vsum to zero
   */

  vsum.t=0; vsum.space.x =0;vsum.space.y =0;vsum.space.z =0;

  /*
   *Sum the nparticles four vectors
   */
  for(i=0;i<nparticles; i++){
    vsum.t += v[i].t;
    vsum.space.x += v[i].space.x;
    vsum.space.y += v[i].space.y;
    vsum.space.z += v[i].space.z;
  }
  /*
   *Calculate the effective mass
   */
  mass = sqrt(vsum.t*vsum.t - (vsum.space.x*vsum.space.x +
			       vsum.space.y*vsum.space.y +
			       vsum.space.z*vsum.space.z ));
  return mass;
}




/*
 ***********************
 * lambda3pi()         *
 ***********************
 */
int lambda3pi(vector4_t *vec,int nvec, double *lambda)
{
  int i;
  vector4_t parent, beta, vecp[3];
  vector3_t analyzer;

  if(nvec !=3){
    fprintf(stderr,"Error:: lambda3pi(): you do not have 3 particles!\n");
    *lambda = -1.0;
  }
  else {
    /* define parent 
     */ 
    parent = Sum4vec(vec,nvec);
    
    /*
     *  Boost to the parent restframe
     */
    
    beta = get_beta(&parent,RESTFRAME);
    for(i=0;i<nvec;i++)
      vecp[i] = lorentz(&beta,&vec[i]);
    
    /* get the normal to the decay plane
     */
    analyzer =  CrossProduct3(&(vecp->space), &((vecp+1)->space));
    
    *lambda = DotProduct3(&analyzer,&analyzer) /
      ( (3.0/4.0)* SQ(  SQ(eff_mass(vec,nvec))/9.0  - SQ(PIMASS)));
  }  /* end of else */  
  return 1;
}
/*
 ***********************
 *                     *
 *  helicityAngles()   *
 *                     *
 ***********************
 */

int helicityAngles(vector4_t *vec,int nvec,
		   double *theta, double *phi)
{

  int i;
  vector3_t z,xhel,yhel,zhel,analyzer;
  vector4_t vecp[3],beta,parent;
  
  /* define parent 
   */
  parent = Sum4vec(vec,nvec);
 
  /* define lab frame */
//  x.x=1; x.y=0; x.z=0;
//  y.x=0; y.y=1; y.z=0;
  z.x=0; z.y=0; z.z=1;

  /* define helicity frame */
  zhel = parent.space;
  yhel = CrossProduct3(&z, &(parent.space));
  xhel =  CrossProduct3(&yhel, &zhel);/* right handed */

  /* Note that the helicity frame is invariant to the boost
   * since zhel is along and yhel is normal to the boost.
   */
  
  /*
   *  Boost to the parent restframe
   */
  
  beta = get_beta(&parent,RESTFRAME);
  for(i=0;i<nvec;i++)
    vecp[i] = lorentz(&beta,&vec[i]);
  
  /*
   * get the helicity angles of the analyzer
   *  (the particle that defines the angles)
   */
  
  switch(nvec){/* define the analyzer */
  case 2:
    analyzer = vecp->space; /* 3momentum of the 1st particle */
    break;
  case 3:
    /* use the normal to the decay plane */
    analyzer =  CrossProduct3(&(vecp->space), &((vecp+1)->space));
    break;
  default:
    fprintf(stderr,"Error(helicityAngles()):: to many decay particles\n");
    exit(-1);
  }
  


  *theta = acos(DotProduct3(&analyzer,&zhel)/
		(v3mag(&analyzer) *v3mag(&zhel)));
  *phi = atan2(DotProduct3(&analyzer,&yhel),
	      DotProduct3(&analyzer,&xhel));

  return 1;

}



/*
 ***********************
 *                     *
 *  END OF FILE        *
 *                     *
 ***********************
 */



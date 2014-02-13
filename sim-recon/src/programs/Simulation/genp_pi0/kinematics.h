
/****************************************************************/
/* kinematics.h                                                 */
/*                                                              */
/* Some basic kinematics routines and tools                     */
/*                                                              */
/* D. Lawrence                                                  */
/* 3/18/99                                                      */
/****************************************************************/


#include <math.h>


#ifndef __KINEMATICS_H__
#define __KINEMATICS_H__

/* vertex array */
typedef struct{
   double x,y,z;
}myvertex_t;
extern myvertex_t myvertex_;

/* four vector */
typedef struct{
   double E;
   double x,y,z;
}vect4;

/* four vector with particle type and position */
typedef struct{
	vect4 v;
	double x,y,z;	/* position */
	int type;     	/* GEANT MC particle type */
}mcparticle_t;


/* Routines */
#ifdef __cplusplus
extern "C" {
#endif
vect4 vect4_add(vect4 v1,vect4 v2);
vect4 vect4_sub(vect4 v1,vect4 v2);
double vect4_mul(vect4 v1,vect4 v2);
double vect4_sq(vect4 v);
vect4 vect4_boost(vect4 p,double beta);
vect4 vect3_cross(vect4 *a, vect4 *b);
vect4 vect3_cross_normalized(vect4 *a, vect4 *b);
char* part_type_str(int type);
int chargeof(int type);

double vect4_mag2(vect4 v);
double vect4_mag(vect4 v);
double vect4_theta(vect4 v);
double vect4_phi(vect4 v);
#ifdef __cplusplus
}
#endif

#endif /* __KINEMATICS_H__ */



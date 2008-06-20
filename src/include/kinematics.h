/*
*
*  kinematics.h
*
*/
/* orignally written by Dave Thompson, pointer madness and
   routine name changes for use in the CLAS by Joe Manak
   the source is currently located in c_bos_io in the files
   vector3.c, vector4.c, lorentz.c and gottjack.c
typedef struct { float32 x,y,z; }               vector3_t;
typedef struct { float32 t; vector3_t space; }  vector4_t;
*/
#include <ntypes.h>
#include <stdio.h>

/* constants */
#ifndef PI
#define PI 3.1415927
#endif

/* 3-vector routines */
float v3mag(vector3_t);
float v3magsq(vector3_t);
float v3radius(vector3_t);
float v3dot(vector3_t,vector3_t);
vector3_t v3cross(vector3_t,vector3_t);
vector3_t v3norm(vector3_t);
void v3dir(vector3_t,float *theta,float *phi);
void v3dir_deg(vector3_t,float *theta,float *phi);
float v3cos_angle(vector3_t ,vector3_t );
float v3angle(vector3_t,vector3_t);
float v3angle_deg(vector3_t, vector3_t);
vector3_t v3add(vector3_t ,vector3_t );
vector3_t v3sub(vector3_t ,vector3_t);
vector3_t v3mult(float,vector3_t);
vector3_t v3div(float,vector3_t);
vector3_t xyz2v3(void *xyz); /*convert x, y, z into a three vector*/
void v3print(FILE *stream, vector3_t vec);

/* 4-vector routines */
float v4dot(vector4_t p1,vector4_t p2);
float v4mass(vector4_t p);
float v4magsq(vector4_t p);
vector4_t v4add(vector4_t p1,vector4_t p2);
vector4_t v4sub(vector4_t p1,vector4_t p2);
vector4_t v4mult(float factor,vector4_t p);
vector4_t v4div(float divisor,vector4_t p);
vector4_t v3_to_v4(vector3_t *vec3, float mass);
vector4_t txyz2v4(void *txyz); /* convert t, x, y, z into a 4-vector*/
void v4print(FILE *stream, vector4_t vec);

/* lorentz routines */
float beta2gamma(float);
float gamma2beta(float);
float p2gamma(vector4_t);
float p2beta(vector4_t);
float beta_p_2mass(float p, float beta);
vector4_t v4boost(vector3_t beta,vector4_t vec);
vector4_t Lab2ResonanceFrame(float BeamEnergy, vector4_t ScatteredElectron, vector4_t V);

/* gottfried-jackson routines */
float costheta_gj(vector4_t decay, vector4_t beam, vector4_t product, vector4_t target, vector4_t recoil);
float phi_ty(vector4_t decay, vector4_t beam, vector4_t product, vector4_t target, vector4_t recoil);
float phi_ty_deg(vector4_t decay, vector4_t beam, vector4_t product, vector4_t target, vector4_t recoil);






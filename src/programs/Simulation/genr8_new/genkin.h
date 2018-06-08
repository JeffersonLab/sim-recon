/*
*
*  kinematics.h
*
*/

#ifndef kinematicsH
#define kinematicsH


/* type declarations */

typedef struct { float x,y,z; } vector3_t;
typedef struct { float t; vector3_t space; } vector4_t;
typedef struct { float rho,phi,z; } vector3cyl_t;

typedef struct { vector3_t e,b; } fields_t;

/* 3-vector routines */
float vec3mag(vector3_t*);
float vec3magsq(vector3_t*);
float vec3dot(vector3_t*,vector3_t*);
vector3_t *vec3cross(vector3_t*,vector3_t*);
vector3_t *vec3norm(vector3_t*);
void vec3dir(vector3_t*,float*,float*);
void vec3dir_deg(vector3_t*,float*,float*);
float vec3cos_angle(vector3_t*,vector3_t*);
float vec3angle(vector3_t*,vector3_t*);
float vec3angle_deg(vector3_t*,vector3_t*);
vector3_t *vec3sum(int,vector3_t*[]);
vector3_t *vec3add(vector3_t*,vector3_t*);
vector3_t *vec3sub(int,vector3_t*[]);
vector3_t *vec3diff(vector3_t*,vector3_t*);
vector3_t *vec3mult(float,vector3_t*);
vector3_t *vec3div(float,vector3_t*);
vector3_t *vec3make(float,float,float);
vector3_t *vec3rotate(vector3_t*,float,int);

/* 4-vector routines */
float vec4dot(vector4_t*,vector4_t*);
float vec4mag(vector4_t*);
float vec4magsq(vector4_t*);
vector4_t *vec4sum(int,vector4_t*[]);
vector4_t *vec4add(vector4_t*,vector4_t*);
vector4_t *vec4sub(int,vector4_t*[]);
vector4_t *vec4diff(vector4_t*,vector4_t*);
vector4_t *vec4mult(float,vector4_t*);
vector4_t *vec4div(float,vector4_t*);
float effmass(int,vector4_t*[]);
void pairmass(int,vector4_t*[],float[]);
float mandel_s(int,vector4_t*[]);
float mandel_t(vector4_t*,vector4_t*);
float mandel_q(vector4_t*,vector4_t*);
vector4_t *vec4make(vector3_t*,float);

/* lorentz routines */
float beta2gamma(float);
float gamma2beta(float);
float p2gamma(vector4_t*);
float p2beta(vector4_t*);
vector4_t *vec4boost(vector3_t*,vector4_t*);
fields_t *fieldboost(vector3_t*,fields_t*);

/* gottfried-jackson routines */
float costheta_gj(vector4_t* v4ptr[]);
float costheta_gj_A(vector4_t v4[]);
float phi_ty(vector4_t*[]);
float phi_ty_deg(vector4_t*[]);
float phi_ty_A(vector4_t[]);
float phi_ty_deg_A(vector4_t[]);
void gottJackGuardFlag(int flag);

/* constants */
#ifndef PI
#define PI 3.1415927
#endif

#define PIMASS 139.6

/* kinematic routines */
double SQ(double x);
double v3mag(const vector3_t *p);
double eff_mass(vector4_t *v, int nparticles);
int helicityAngles(vector4_t *vec,int nvec,
		   double *theta, double *phi);
int lambda3pi(vector4_t *vec,int nvec, double *lam);
vector4_t get_beta(vector4_t *boost,int sign);
vector4_t Sum4vec(vector4_t *vec4, int nvec4);
vector4_t lorentz(const vector4_t *beta,const vector4_t *pin);
double CMmomentum(double cm_engy, double m1, double m2);
double energy(double mass, const vector3_t *p);
#endif
/* end file */

/*
 * esr.h  --- Event Summary Records
*/


#ifndef include_esrH
#define include_esrH

typedef struct {
  vector4_t p;
  int charge;
  float xi2;
} esr_track_t;


typedef struct {
  vector3_t production;
  vector4_t beam;
  vector4_t miss;
  float mass,mmsq;
  int nprongs;
  int charge;
  esr_track_t pi[1];
}esr_nprong_t;

typedef struct {
  int groupID,index;
} particleID_t;

typedef struct {
  vector4_t p;
  int charge;
  float xi2;
  int flags;
  Particle_t particleType;
  particleID_t parent,child;
} esr_particle_t;

typedef struct {
  vector3_t production;
  vector4_t beam;
  vector4_t miss;
  float mass,mmsq;
  int nparticles;
  int ngammas;
  int npplus;
  int npminus;
  float chsq;/* full fit Xi2*/
  float cl; /* full fit confidence level */
  int charge;
  esr_particle_t p[1];
} esr_nparticle_t;

typedef struct {
  int nvertices;
  vector3_t v[1];
} esr_vertices_t; /* To store secondary vertices */

typedef struct {
  int32 nweights;
  int32 pad;   /*  This is here as a placeholder for alignment with double. */
  double weight[1];
} EventWeights_t;   /* Group of weights (for PWA) */

typedef struct {
  double real;
  double imag;
} amplitude_t;

typedef struct {
  int namps;
  amplitude_t amps[1];
} amplitudes_t;    /* Group of PWA amplitudes */


typedef struct {
  int id;
  int nwords;
  int data[1];
} generic_t;


/* Compressed ESR group */

typedef struct
{
  vector4_t  momentum;
  Particle_t type;

} esr_compressedParticle_t;

typedef struct
{
  int32     nparticles;   /* number of particles in this event  */
  int32     user1;        /* user defined */
  float32   accepted;     /* =0 if not accepted, =!0 if accepted */
  float32   weight;
  float32   user2;
  float32   user3;
  vector3_t vertex;

  esr_compressedParticle_t beam;
  esr_compressedParticle_t particles[1];

} esr_compressed_t;


/*
 * produce an ESR:
 *   take data from 'eventIn'
 *   put the esr into 'eventOut' (can be same as 'eventIn')
 *   'triggerMask' selects triggers for which esr's are
 *   produced, see tmask() for details.
 * returns: 0 - Okey, (-1) - no esr created
*/ 

int esr_processEvent(const void*eventIn,void*eventOut,int eventOutSize,int triggerMask);
int esrNparticle(const void*eventIn,void*eventOut,int eventOutSize,int triggerMask,float MinGammaEnergy);
int esrNparticleSquawMV(void *tmp,int bufsize,void *tmp2);


#endif
/* end file */

/*
 * particleType.h
*/

#ifndef particleTypeH_INCLUDED
#define particleTypeH_INCLUDED

static const char sccsid_particleTypeH[] = "@(#)particleType.h\t5.3\tCreated 8/20/97 17:13:12, \tcompiled "__DATE__;

typedef enum {

  /*
   * These constants are defined to be
   * same as GEANT. see http://wwwcn.cern.ch/asdoc/geant/H2GEANTCONS300.html
   * for more details.
  */

  Unknown        =  0,
  Gamma          =  1,
  Positron       =  2,
  Electron       =  3,
  Neutrino       =  4,
  MuonPlus       =  5,
  MuonMinus      =  6,
  Pi0            =  7,
  PiPlus         =  8,
  PiMinus        =  9,
  KLong          = 10,
  KPlus          = 11,
  KMinus         = 12,
  Neutron        = 13,
  Proton         = 14,
  AntiProton     = 15,
  KShort         = 16,
  Eta            = 17,
  Lambda         = 18,
  SigmaPlus      = 19,
  Sigma0         = 20,
  SigmaMinus     = 21,
  Xi0            = 22,
  XiMinus        = 23,
  OmegaMinus     = 24,
  AntiNeutron    = 25,
  AntiLambda     = 26,
  AntiSigmaMinus = 27,
  AntiSigma0     = 28,
  AntiSigmaPlus  = 29,
  AntiXi0        = 30,
  AntiXiPlus     = 31,
  AntiOmegaPlus  = 32,

  /* the constants defined by GEANT end here */

  /* These are E852-defined constants */

  Rho0           = 57,
  RhoPlus        = 58,
  RhoMinus       = 59,
  omega          = 60,
  EtaPrime       = 61,
  phiMeson       = 62

} Particle_t;

char *ParticleType(Particle_t);
float ParticleMass(Particle_t);
int   ParticleCharge(Particle_t);

#endif
/* end file */

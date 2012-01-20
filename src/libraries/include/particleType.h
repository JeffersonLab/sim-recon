/*
 * particleType.h
*/

#ifndef particleTypeH_INCLUDED
#define particleTypeH_INCLUDED

#include <math.h>
#include <stdio.h>
#include <string.h>

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
  Geantino       = 48,

  /* the constants defined by GEANT end here */

  /* These are E852-defined constants */

  Rho0           = 57,
  RhoPlus        = 58,
  RhoMinus       = 59,
  omega          = 60,
  EtaPrime       = 61,
  phiMeson       = 62,
  a0_980	 = 63,
  f0_980	 = 64

} Particle_t;

inline static char *ParticleType(Particle_t p)
{
  static char ret[20];
  switch (p) {
  case Unknown:
    strcpy(ret,"unknown");
    break;
  case Gamma:
    strcpy(ret,"gamma");
    break;
  case Positron:
    strcpy(ret,"positron");
    break;
  case Electron:
    strcpy(ret,"electron");
    break;
  case Neutrino:
    strcpy(ret,"neutrino");
    break;
  case MuonPlus:
    strcpy(ret,"mu+");
    break;
  case MuonMinus:
    strcpy(ret,"mu-");
    break;
  case Pi0:
    strcpy(ret,"pi0");
    break;
  case PiPlus:
    strcpy(ret,"pi+");
    break;
  case PiMinus:
    strcpy(ret,"pi-");
    break;
  case KLong:
    strcpy(ret,"kL");
    break;
  case KPlus:
    strcpy(ret,"k+");
    break;
  case KMinus:
    strcpy(ret,"k-");
    break;
  case Neutron:
    strcpy(ret,"neutron");
    break;
  case Proton:
    strcpy(ret,"proton");
    break;
  case AntiProton:
    strcpy(ret,"antiProton");
    break;
  case KShort:
    strcpy(ret,"kS");
    break;
  case Eta:
    strcpy(ret,"eta");
    break;
  case Lambda:
    strcpy(ret,"Lambda");
    break;
  case SigmaPlus:
    strcpy(ret,"Sigma+");
    break;
  case Sigma0:
    strcpy(ret,"Sigma0");
    break;
  case SigmaMinus:
    strcpy(ret,"Sigma-");
    break;
  case Xi0:
    strcpy(ret,"Xi0");
    break;
  case XiMinus:
    strcpy(ret,"Xi-");
    break;
  case OmegaMinus:
    strcpy(ret,"Omega-");
    break;
  case AntiNeutron:
    strcpy(ret,"antiNeutron");
    break;
  case AntiLambda:
    strcpy(ret,"antiLambda");
    break;
  case AntiSigmaMinus:
    strcpy(ret,"antiSigma-");
    break;
  case AntiSigma0:
    strcpy(ret,"antiSigma0");
    break;
  case AntiSigmaPlus:
    strcpy(ret,"antiSigma+");
    break;
  case AntiXi0:
    strcpy(ret,"antiXi0");
    break;
  case AntiXiPlus:
    strcpy(ret,"antiXi+");
    break;
  case AntiOmegaPlus:
    strcpy(ret,"antiOmega+");
    break;
  case Geantino:
    strcpy(ret,"geantino");
    break;
  case Rho0:
    strcpy(ret,"rho0");
    break;  
  case RhoPlus:
    strcpy(ret,"rho+");
    break;
  case RhoMinus:
    strcpy(ret,"rho-");
    break;
  case omega:
    strcpy(ret,"omega");
    break;
  case EtaPrime:
    strcpy(ret,"etaPrime");
    break;
  case phiMeson:
    strcpy(ret,"phi");
    break;
  case a0_980:
    strcpy(ret,"a0(980)");
    break;
  case f0_980:
    strcpy(ret,"f0(980)");
    break;
  default:
    sprintf(ret,"type(%d)",(int)p);
    break;
  }
  return(ret);
}

inline static double ParticleMass(Particle_t p)
{
  switch (p) {
  default:
    fprintf(stderr,"ParticleMass: Error: Unknown particle type %d,",p);
    fprintf(stderr," returning HUGE_VAL...\n");
    return HUGE_VAL;
  case Unknown:		return HUGE_VAL;
  case Gamma:		return 0;
  case Positron:	return 0.0005101;
  case Electron:	return 0.0005101;
  case Neutrino:	return 0;
  case MuonPlus:	return 0.105658;
  case MuonMinus:	return 0.105658;
  case Pi0:		return 0.13497;
  case PiPlus:		return 0.139568;
  case PiMinus:		return 0.139568;
  case KShort:		return 0.497671;
  case KLong:		return 0.497671;
  case KPlus:		return 0.49364;
  case KMinus:		return 0.49364;
  case Neutron:		return 0.93956;
  case Proton:		return 0.93827;
  case AntiProton:	return 0.93827;
  case Eta:		return 0.54745;
  case Lambda:		return 1.11568;
  case SigmaPlus:	return 1.18937;
  case Sigma0:		return 1.19264;
  case SigmaMinus:	return 1.18937;
  case Xi0:		return 1.31483;
  case XiMinus:		return 1.32131;
  case OmegaMinus:	return 1.67245;
  case AntiNeutron:	return 0.93956;
  case AntiLambda:	return 1.11568;
  case AntiSigmaMinus:	return 1.18937;
  case AntiSigma0:	return 1.19264;
  case AntiSigmaPlus:	return 1.18937;
  case AntiXi0:		return 1.31483;
  case AntiXiPlus:	return 1.32131;
  case AntiOmegaPlus:	return 1.67245;
  case Geantino:		return 0.0;
  case Rho0:		return 0.7693;
  case RhoPlus:		return 0.7693;
  case RhoMinus:	return 0.7693;
  case omega:		return 0.78257;
  case EtaPrime:	return 0.95778;
  case phiMeson:	return 1.01942;
  case a0_980:		return 0.980;
  case f0_980:		return 0.980;
  }
}

inline static int ParticleCharge(Particle_t p)
{
  switch (p) {
  default:
    fprintf(stderr,"ParticleCharge: Error: Unknown particle type %d,",p);
    fprintf(stderr," returning 0...\n");
    return 0;
  case Unknown:		return  0;
  case Gamma:		return  0;
  case Positron:	return +1;
  case Electron:	return -1;
  case Neutrino:	return  0;
  case MuonPlus:	return +1;
  case MuonMinus:	return -1;
  case Pi0:		return  0;
  case PiPlus:		return +1;
  case PiMinus:		return -1;
  case KShort:		return  0;
  case KLong:		return  0;
  case KPlus:		return +1;
  case KMinus:		return -1;
  case Neutron:		return  0;
  case Proton:		return +1;
  case AntiProton:	return -1;
  case Eta:		return  0;
  case Lambda:		return  0;
  case SigmaPlus:	return +1;
  case Sigma0:		return  0;
  case SigmaMinus:	return -1;
  case Xi0:		return  0;
  case XiMinus:		return -1;
  case OmegaMinus:	return -1;
  case AntiNeutron:	return  0;
  case AntiLambda:	return  0;
  case AntiSigmaMinus:	return -1;
  case AntiSigma0:	return  0;
  case AntiSigmaPlus:	return +1;
  case AntiXi0:		return  0;
  case AntiXiPlus:	return +1;
  case AntiOmegaPlus:	return +1;
  case Geantino:		return  0;
  case Rho0:		return  0;
  case RhoPlus:		return +1;
  case RhoMinus:	return -1;
  case omega:		return  0;
  case EtaPrime:	return  0;
  case phiMeson:	return  0;
  case a0_980:		return  0;
  case f0_980:		return  0;
  }
}

#endif

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

float ParticleMass(Particle_t);
int   ParticleCharge(Particle_t);

char *ParticleType(Particle_t p)
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
    strcpy(ret,"KL");
    break;
  case KPlus:
    strcpy(ret,"K+");
    break;
  case KMinus:
    strcpy(ret,"K-");
    break;
  case Neutron:
    strcpy(ret,"neutron");
    break;
  case Proton:
    strcpy(ret,"proton");
    break;
  case AntiProton:
    strcpy(ret,"pbar");
    break;
  case KShort:
    strcpy(ret,"Ks");
    break;
  case Eta:
    strcpy(ret,"eta");
    break;
  case Lambda:
    strcpy(ret,"lambda");
    break;
  case SigmaPlus:
    strcpy(ret,"sigma+");
    break;
  case Sigma0:
    strcpy(ret,"sigma0");
    break;
  case SigmaMinus:
    strcpy(ret,"sigma-");
    break;
  case Xi0:
    strcpy(ret,"Xi0");
    break;
  case XiMinus:
    strcpy(ret,"Xi-");
    break;
  case OmegaMinus:
    strcpy(ret,"omega-");
    break;
  case AntiNeutron:
    strcpy(ret,"nbar");
    break;

  case AntiLambda:
    strcpy(ret,"lambdabar");
    break;
  case AntiSigmaMinus:
    strcpy(ret,"sigmabar-");
    break;
  case AntiSigma0:
    strcpy(ret,"sigmabar0");
    break;
  case AntiSigmaPlus:
    strcpy(ret,"sigmabar+");
    break;
  case AntiXi0:
    strcpy(ret,"Xibar0");
    break;
  case AntiXiPlus:
    strcpy(ret,"Xibar+");
    break;
  case AntiOmegaPlus:
    strcpy(ret,"omegabar+");
    break;
  case Rho0:
    strcpy(ret,"rho0");
    break;  
  case RhoPlus:
    strcpy(ret,"rho+");
    break;
  case RhoMinus:
    strcpy(ret,"rho;");
    break;
  case omega:
    strcpy(ret,"omega");
    break;
  case EtaPrime:
    strcpy(ret,"etaprime");
    break;
  case phiMeson:
    strcpy(ret,"phi");
    break;
  default:
    sprintf(ret,"type(%d)",(int)p);
    break;
  }
  return(ret);
}




#endif
/* end file */

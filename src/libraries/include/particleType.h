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
  f0_980	 = 64,

  KStar_892_0 = 65,
  KStar_892_Plus = 66,
  KStar_892_Minus = 67,
  AntiKStar_892_0 = 68,

  K1_1400_Plus = 69,
  K1_1400_Minus = 70,

  b1_1235_Plus = 71

} Particle_t;

inline static char *ParticleType(Particle_t p)
{
  static char ret[20];
  switch (p) {
  case Unknown:
    strcpy(ret,"Unknown");
    break;
  case Gamma:
    strcpy(ret,"Gamma");
    break;
  case Positron:
    strcpy(ret,"Positron");
    break;
  case Electron:
    strcpy(ret,"Electron");
    break;
  case Neutrino:
    strcpy(ret,"Neutrino");
    break;
  case MuonPlus:
    strcpy(ret,"Muon+");
    break;
  case MuonMinus:
    strcpy(ret,"Muon-");
    break;
  case Pi0:
    strcpy(ret,"Pi0");
    break;
  case PiPlus:
    strcpy(ret,"Pi+");
    break;
  case PiMinus:
    strcpy(ret,"Pi-");
    break;
  case KLong:
    strcpy(ret,"KLong");
    break;
  case KPlus:
    strcpy(ret,"K+");
    break;
  case KMinus:
    strcpy(ret,"K-");
    break;
  case Neutron:
    strcpy(ret,"Neutron");
    break;
  case Proton:
    strcpy(ret,"Proton");
    break;
  case AntiProton:
    strcpy(ret,"AntiProton");
    break;
  case KShort:
    strcpy(ret,"KShort");
    break;
  case Eta:
    strcpy(ret,"Eta");
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
    strcpy(ret,"AntiNeutron");
    break;
  case AntiLambda:
    strcpy(ret,"AntiLambda");
    break;
  case AntiSigmaMinus:
    strcpy(ret,"AntiSigma-");
    break;
  case AntiSigma0:
    strcpy(ret,"AntiSigma0");
    break;
  case AntiSigmaPlus:
    strcpy(ret,"AntiSigma+");
    break;
  case AntiXi0:
    strcpy(ret,"AntiXi0");
    break;
  case AntiXiPlus:
    strcpy(ret,"AntiXi+");
    break;
  case AntiOmegaPlus:
    strcpy(ret,"AntiOmega+");
    break;
  case Geantino:
    strcpy(ret,"Geantino");
    break;
  case Rho0:
    strcpy(ret,"Rho0");
    break;  
  case RhoPlus:
    strcpy(ret,"Rho+");
    break;
  case RhoMinus:
    strcpy(ret,"Rho-");
    break;
  case omega:
    strcpy(ret,"omega");
    break;
  case EtaPrime:
    strcpy(ret,"EtaPrime");
    break;
  case phiMeson:
    strcpy(ret,"phiMeson");
    break;
  case a0_980:
    strcpy(ret,"a0(980)");
    break;
  case f0_980:
    strcpy(ret,"f0(980)");
    break;
  case KStar_892_0:
    strcpy(ret,"K*(892)0");
    break;
  case KStar_892_Plus:
    strcpy(ret,"K*(892)+");
    break;
  case KStar_892_Minus:
    strcpy(ret,"K*(892)-");
    break;
  case AntiKStar_892_0:
    strcpy(ret,"antiK*(892)0");
    break;
  case K1_1400_Plus:
    strcpy(ret,"K1(1400)+");
    break;
  case K1_1400_Minus:
    strcpy(ret,"K1(1400)-");
    break;
  case b1_1235_Plus:
    strcpy(ret,"b1(1235)+");
    break;
  default:
    sprintf(ret,"type(%d)",(int)p);
    break;
  }
  return(ret);
}

inline static unsigned short int IsFixedMass(Particle_t p)
{
  switch (p)
  {
  case Gamma:		return 1;
  case Positron:	return 1;
  case Electron:	return 1;
  case Neutrino:	return 1;
  case MuonPlus:	return 1;
  case MuonMinus:	return 1;
  case Pi0:      	return 1;
  case PiPlus:		return 1;
  case PiMinus:		return 1;
  case KShort:		return 1;
  case KLong:		return 1;
  case KPlus:		return 1;
  case KMinus:		return 1;
  case Neutron:		return 1;
  case Proton:		return 1;
  case AntiProton:	return 1;
  case Eta:			return 1;
  case Lambda:		return 1;
  case SigmaPlus:	return 1;
  case Sigma0:		return 1;
  case SigmaMinus:	return 1;
  case Xi0:			return 1;
  case XiMinus:		return 1;
  case OmegaMinus:	return 1;
  case AntiNeutron:	return 1;
  case AntiLambda:	return 1;
  case AntiSigmaMinus:	return 1;
  case AntiSigma0:	return 1;
  case AntiSigmaPlus:	return 1;
  case AntiXi0:		return 1;
  case AntiXiPlus:	return 1;
  case AntiOmegaPlus:	return 1;
  case Geantino:	return 1;
  case EtaPrime:	return 1;
  default: return 0;
  }
}

inline static unsigned short int IsDetachedVertex(Particle_t p)
{
  switch (p)
  {
  case MuonPlus:	return 1;
  case MuonMinus:	return 1;
  case PiPlus:		return 1;
  case PiMinus:		return 1;
  case KShort:		return 1;
  case KLong:		return 1;
  case KPlus:		return 1;
  case KMinus:		return 1;
  case Neutron:		return 1;
  case Lambda:		return 1;
  case SigmaPlus:	return 1;
  case SigmaMinus:	return 1;
  case Xi0:			return 1;
  case XiMinus:		return 1;
  case OmegaMinus:	return 1;
  case AntiNeutron:	return 1;
  case AntiLambda:	return 1;
  case AntiSigmaMinus:	return 1;
  case AntiSigmaPlus:	return 1;
  case AntiXi0:		return 1;
  case AntiXiPlus:	return 1;
  case AntiOmegaPlus:	return 1;
  default: return 0;
  }
}

inline static char* ParticleName_ROOT(Particle_t p)
{
  static char ret[40];
  switch (p) {
  case Unknown:
    strcpy(ret, "#it{X}");
    break;
  case Gamma:
    strcpy(ret, "#it{#gamma}");
    break;
  case Positron:
    strcpy(ret, "#it{e}^{+}");
    break;
  case Electron:
    strcpy(ret, "#it{e}^{-}");
    break;
  case Neutrino:
    strcpy(ret, "#it{#nu}");
    break;
  case MuonPlus:
    strcpy(ret, "#it{#mu}^{+}");
    break;
  case MuonMinus:
    strcpy(ret, "#it{#mu}^{-}");
    break;
  case Pi0:
    strcpy(ret, "#it{#pi}^{0}");
    break;
  case PiPlus:
    strcpy(ret, "#it{#pi}^{+}");
    break;
  case PiMinus:
    strcpy(ret, "#it{#pi}^{-}");
    break;
  case KLong:
    strcpy(ret, "#it{K}^{0}_{L}");
    break;
  case KPlus:
    strcpy(ret, "#it{K}^{+}");
    break;
  case KMinus:
    strcpy(ret, "#it{K}^{-}");
    break;
  case Neutron:
    strcpy(ret, "#it{n}");
    break;
  case Proton:
    strcpy(ret, "#it{p}");
    break;
  case AntiProton:
    strcpy(ret, "#it{#bar{p}}");
    break;
  case KShort:
    strcpy(ret, "#it{K}^{0}_{S}");
    break;
  case Eta:
    strcpy(ret, "#it{#eta}");
    break;
  case Lambda:
    strcpy(ret, "#it{#Lambda}");
    break;
  case SigmaPlus:
    strcpy(ret, "#it{#Sigma}^{+}");
    break;
  case Sigma0:
    strcpy(ret, "#it{#Sigma}^{0}");
    break;
  case SigmaMinus:
    strcpy(ret, "#it{#Sigma}^{-}");
    break;
  case Xi0:
    strcpy(ret, "#it{#Xi}^{0}");
    break;
  case XiMinus:
    strcpy(ret, "#it{#Xi}^{-}");
    break;
  case OmegaMinus:
    strcpy(ret, "#it{#Omega}^{-}");
    break;
  case AntiNeutron:
    strcpy(ret, "#it{#bar^{n}}");
    break;
  case AntiLambda:
    strcpy(ret, "#it{#bar^{#Lambda}}");
    break;
  case AntiSigmaMinus:
    strcpy(ret, "#it{#bar{#Sigma}}^{-}");
    break;
  case AntiSigma0:
    strcpy(ret, "#it{#bar{#Sigma}}^{0}");
    break;
  case AntiSigmaPlus:
    strcpy(ret, "#it{#bar{#Sigma}}^{+}");
    break;
  case AntiXi0:
    strcpy(ret, "#it{#bar{#Xi}}^{0}");
    break;
  case AntiXiPlus:
    strcpy(ret, "#it{#bar{#Xi}}^{+}");
    break;
  case AntiOmegaPlus:
    strcpy(ret, "#it{#bar{#Omega}}^{+}");
    break;
  case Geantino:
    strcpy(ret, "geantino");
    break;
  case Rho0:
    strcpy(ret, "#it{#rho}^{0}");
    break;
  case RhoPlus:
    strcpy(ret, "#it{#rho}^{+}");
    break;
  case RhoMinus:
    strcpy(ret, "#it{#rho}^{-}");
    break;
  case omega:
    strcpy(ret, "#it{#omega}");
    break;
  case EtaPrime:
    strcpy(ret, "#it{#eta'}");
    break;
  case phiMeson:
    strcpy(ret, "#it{#phi}");
    break;
  case a0_980:
    strcpy(ret, "#it{a}_{0}(980)");
    break;
  case f0_980:
    strcpy(ret, "#it{f}_{0}(980)");
    break;
  case KStar_892_0:
    strcpy(ret, "#it{K}*(892)^{0}");
    break;
  case KStar_892_Plus:
    strcpy(ret, "#it{K}*(892)^{+}");
    break;
  case KStar_892_Minus:
    strcpy(ret, "#it{K}*(892)^{-}");
    break;
  case AntiKStar_892_0:
    strcpy(ret, "#it{#bar{K*}}(892)^{0}");
    break;
  case K1_1400_Plus:
    strcpy(ret, "#it{K}_{1}(1400)^{+}");
    break;
  case K1_1400_Minus:
    strcpy(ret, "#it{K}_{1}(1400)^{-}");
    break;
  case b1_1235_Plus:
    strcpy(ret, "#it{b}_{1}(1235)^{+}");
    break;
  default:
    strcpy(ret, "X");
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
  case PiPlus:		return 0.139570;
  case PiMinus:		return 0.139570;
  case KShort:		return 0.497671;
  case KLong:		return 0.497671;
  case KPlus:		return 0.493677;
  case KMinus:		return 0.493677;
  case Neutron:		return 0.93956;
  case Proton:		return 0.938272;
  case AntiProton:	return 0.938272;
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
  case KStar_892_0: return 0.89594;
  case KStar_892_Plus: return 0.89166;
  case KStar_892_Minus: return 0.89166;
  case AntiKStar_892_0: return 0.89594;
  case K1_1400_Plus: return 1.403;
  case K1_1400_Minus: return 1.403;
  case b1_1235_Plus: return 1.2295;
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
  case Geantino:	return  0;
  case Rho0:		return  0;
  case RhoPlus:		return +1;
  case RhoMinus:	return -1;
  case omega:		return  0;
  case EtaPrime:	return  0;
  case phiMeson:	return  0;
  case a0_980:		return  0;
  case f0_980:		return  0;
  case KStar_892_0: return  0;
  case KStar_892_Plus: return  1;
  case KStar_892_Minus: return -1;
  case AntiKStar_892_0: return  0;
  case K1_1400_Plus: return  1;
  case K1_1400_Minus: return -1;
  case b1_1235_Plus: return 1;
  }
}

inline static int PDGtype(Particle_t p)
{
  switch (p) {
  case Unknown:		return  0;
  case Gamma:		return  22;
  case Positron:	return -11;
  case Electron:	return  11;
  case Neutrino:	return  121416;
  case MuonPlus:	return -13;
  case MuonMinus:	return  13;
  case Pi0:		return  111;
  case PiPlus:		return  211;
  case PiMinus:		return -211;
  case KShort:		return  310;
  case KLong:		return  130;
  case KPlus:		return  321;
  case KMinus:		return -321;
  case Neutron:		return  2112;
  case Proton:		return  2212;
  case AntiProton:	return -2212;
  case Eta:		return  221;
  case Lambda:		return  3122;
  case SigmaPlus:	return  3222;
  case Sigma0:		return  3212;
  case SigmaMinus:	return  3112;
  case Xi0:		return  3322;
  case XiMinus:		return  3312;
  case OmegaMinus:	return  3332;
  case AntiNeutron:	return -2112;
  case AntiLambda:	return -3122;
  case AntiSigmaMinus:	return -3112;
  case AntiSigma0:	return -3212;
  case AntiSigmaPlus:	return -3222;
  case AntiXi0:		return -3322;
  case AntiXiPlus:	return -3312;
  case AntiOmegaPlus:	return -3332;
  case Geantino:	return  0;
  case Rho0:		return  113;
  case RhoPlus:		return  213;
  case RhoMinus:	return -213;
  case omega:		return  223;
  case EtaPrime:	return  331;
  case phiMeson:	return  333;
  case a0_980:		return  9000110;
  case f0_980:		return  9010221;
  case KStar_892_0: return  313;
  case AntiKStar_892_0: return  -313;
  case KStar_892_Plus: return  323;
  case KStar_892_Minus: return -323;
  case K1_1400_Plus: return  20323;
  case K1_1400_Minus: return  -20323;
  case b1_1235_Plus: return  10213;
  default:		return  0;
  }
}

#endif

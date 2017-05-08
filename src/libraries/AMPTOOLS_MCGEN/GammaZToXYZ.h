#if !defined(GAMMAZTOXYZ)
#define GAMMAZTOXYZ

/*
 *  GammaZToXYZ.h
 *  GlueXTools
 *
 *  Modified GammaPToXYP.h, replacing proton with heavy Z target
 *  Elton 4/14/2017
 *
 */

#include "TLorentzVector.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "TH1.h"

class Kinematics;
class AmpVecs;

class GammaZToXYZ {
  
public:
  
  GammaZToXYZ( float lowMassXY, float highMassXY, float massX, float massY,
               float beamMaxE, float beamPeakE, float beamLowE, float beamHigh, ProductionMechanism::Type type );
  
  Kinematics* generate();
  
  void addResonance( float mass, float width, float bf );
  
private:

  ProductionMechanism m_prodMech;
  
  TLorentzVector m_beam;
  TLorentzVector m_target;
  
  vector< double > m_childMass;

  TH1D *cobrem_vs_E;
  TH1D *Primakoff_tdist;
  
};

#endif

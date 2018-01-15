#if !defined(BEAMPROPERTIES)
#define BEAMPROPERTIES

/*
 *  BeamProperties.cc
 *
 *  Contains histograms for beam properties to be used in event generation and fitting.  Source
 *  of beam properties is from CombremsGeneration, external ROOT file or CCDB (to be implemented). 
 *
 *  Created by Justin Stevens on 12/29/17
 */

#include <string>
#include <map>

#include "TH1.h"

class BeamProperties {
  
public:
  
  BeamProperties( TString configFile);

  inline TH1D* GetFlux() { return fluxVsEgamma; };
  inline TH1D* GetPolFrac() { return polFracVsEgamma; };
  double GetPolAngle();

private:

  void createHistograms( TString configFile );
  bool parseConfig();
  void generateCobrems();
  void fillFromROOT();
  void fillFromCCDB();

  TString mConfigFile;
  std::map<std::string,double> mBeamParametersMap;
  std::map<std::string,std::string> mBeamHistNameMap;

  bool mIsROOT;
  bool mIsCCDB;

  TH1D *fluxVsEgamma;
  TH1D *polFracVsEgamma;

};

#endif

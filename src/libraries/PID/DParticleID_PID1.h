// $Id$
//
//    File: DParticleID_PID1.h
// Created: Mon Feb 28 15:25:35 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleID_PID1_
#define _DParticleID_PID1_

#include <JANA/jerror.h>
#include "DParticleID.h"

class DParticleID_PID1:public DParticleID{
 public:
  DParticleID_PID1(JEventLoop *loop); // require JEventLoop in constructor;
  ~DParticleID_PID1();

	jerror_t GetdEdxMean_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx, Particle_t locPIDHypothesis) const;
	jerror_t GetdEdxSigma_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx, Particle_t locPIDHypothesis) const;
	jerror_t GetdEdxMean_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx, Particle_t locPIDHypothesis) const;
	jerror_t GetdEdxSigma_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx, Particle_t locPIDHypothesis) const;
	jerror_t CalcDCdEdxChiSq(DChargedTrackHypothesis *locChargedTrackHypothesis) const;
	inline double Function_dEdx(double locBetaGamma, const vector<float> &locParams) const{return locParams[0]/(locBetaGamma*locBetaGamma)+locParams[1]/locBetaGamma + locParams[2] + locParams[3]*locBetaGamma;}
	inline double Function_dEdxSigma(double locBetaGamma, const vector<float> &locParams) const{return locParams[0]/(locBetaGamma*locBetaGamma)+locParams[1]/locBetaGamma + locParams[2];}


 protected:
	vector<float> ddEdxMeanParams_FDC_Proton;
	vector<float> ddEdxMeanParams_FDC_KPlus;
	vector<float> ddEdxMeanParams_FDC_PiPlus;
	vector<float> ddEdxMeanParams_CDC_Proton;
	vector<float> ddEdxMeanParams_CDC_KPlus;
	vector<float> ddEdxMeanParams_CDC_PiPlus;

	vector<float> ddEdxSigmaParams_FDC_Proton;
	vector<float> ddEdxSigmaParams_FDC_KPlus;
	vector<float> ddEdxSigmaParams_FDC_PiPlus;
	vector<float> ddEdxSigmaParams_CDC_Proton;
	vector<float> ddEdxSigmaParams_CDC_KPlus;
	vector<float> ddEdxSigmaParams_CDC_PiPlus;


 private:
  int DEBUG_LEVEL;
  // Prohibit default constructor
  DParticleID_PID1();	

};

#endif // _DParticleID_PID1_


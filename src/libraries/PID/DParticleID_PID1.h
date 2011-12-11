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
#include "TF1.h"

class DParticleID_PID1:public DParticleID{
 public:
  DParticleID_PID1(JEventLoop *loop); // require JEventLoop in constructor;
  ~DParticleID_PID1();

	jerror_t GetdEdxMean_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx, Particle_t locPIDHypothesis);
	jerror_t GetdEdxSigma_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx, Particle_t locPIDHypothesis);
	jerror_t GetdEdxMean_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx, Particle_t locPIDHypothesis);
	jerror_t GetdEdxSigma_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx, Particle_t locPIDHypothesis);
	jerror_t CalcDCdEdxChiSq(const DChargedTrackHypothesis *locChargedTrackHypothesis, double &locChiSq, unsigned int& locNDF);

 protected:
	TF1* ddEdxMeanFunc_FDC_Proton;
	TF1* ddEdxMeanFunc_FDC_KPlus;
	TF1* ddEdxMeanFunc_CDC_Proton;
	TF1* ddEdxMeanFunc_CDC_KPlus;

	vector<TF1*> ddEdxSigmaFuncVector_CDC_Proton;
	vector<TF1*> ddEdxSigmaFuncVector_FDC_Proton;
	vector<unsigned int> ddEdxSigmaNumHitsVector_FDC_Proton;
	vector<unsigned int> ddEdxSigmaNumHitsVector_CDC_Proton;

	vector<TF1*> ddEdxSigmaFuncVector_FDC_PiPlus;
	vector<TF1*> ddEdxSigmaFuncVector_CDC_PiPlus;
	vector<unsigned int> ddEdxSigmaNumHitsVector_CDC_PiPlus;
	vector<unsigned int> ddEdxSigmaNumHitsVector_FDC_PiPlus;

	vector<TF1*> ddEdxSigmaFuncVector_FDC_KPlus;
	vector<TF1*> ddEdxSigmaFuncVector_CDC_KPlus;
	vector<unsigned int> ddEdxSigmaNumHitsVector_CDC_KPlus;
	vector<unsigned int> ddEdxSigmaNumHitsVector_FDC_KPlus;

 private:
  int DEBUG_LEVEL;
  // Prohibit default constructor
  DParticleID_PID1();	

};

#endif // _DParticleID_PID1_


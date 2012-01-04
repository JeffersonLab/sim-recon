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

	jerror_t GetdEdxMean_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx, Particle_t locPIDHypothesis);
	jerror_t GetdEdxSigma_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx, Particle_t locPIDHypothesis);
	jerror_t GetdEdxMean_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx, Particle_t locPIDHypothesis);
	jerror_t GetdEdxSigma_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx, Particle_t locPIDHypothesis);
	jerror_t CalcDCdEdxChiSq(const DChargedTrackHypothesis *locChargedTrackHypothesis, double &locChiSq, unsigned int& locNDF);
	inline double Function_dEdx(double locBetaGamma, const vector<float> &locParams){return locParams[0]*exp(locParams[1]*locBetaGamma) + locParams[2] + locParams[3]*locBetaGamma;}
	void Set_dEdxParams(vector<float>& locParamVector, float locParam1, float locParam2, float locParam3, float locParam4);


 protected:
	vector<float> ddEdxMeanParams_FDC_Proton;
	vector<float> ddEdxMeanParams_FDC_KPlus;
	vector<float> ddEdxMeanParams_CDC_Proton;
	vector<float> ddEdxMeanParams_CDC_KPlus;

	vector<vector<float> > ddEdxSigmaParamVector_CDC_Proton;
	vector<vector<float> > ddEdxSigmaParamVector_FDC_Proton;
	vector<unsigned int> ddEdxSigmaNumHitsVector_FDC_Proton;
	vector<unsigned int> ddEdxSigmaNumHitsVector_CDC_Proton;

	vector<vector<float> > ddEdxSigmaParamVector_CDC_PiPlus;
	vector<vector<float> > ddEdxSigmaParamVector_FDC_PiPlus;
	vector<unsigned int> ddEdxSigmaNumHitsVector_CDC_PiPlus;
	vector<unsigned int> ddEdxSigmaNumHitsVector_FDC_PiPlus;

	vector<vector<float> > ddEdxSigmaParamVector_CDC_KPlus;
	vector<vector<float> > ddEdxSigmaParamVector_FDC_KPlus;
	vector<unsigned int> ddEdxSigmaNumHitsVector_CDC_KPlus;
	vector<unsigned int> ddEdxSigmaNumHitsVector_FDC_KPlus;

	vector<float> dBetaGamma_PiMinus_CDC;
	vector<float> ddEdxMean_PiMinus_CDC;
	vector<float> dBetaGamma_PiMinus_FDC;
	vector<float> ddEdxMean_PiMinus_FDC;

 private:
  int DEBUG_LEVEL;
  // Prohibit default constructor
  DParticleID_PID1();	

};

#endif // _DParticleID_PID1_


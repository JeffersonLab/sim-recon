// $Id$
//
//    File: DParticleID_PID1.cc
// Created: Mon Feb 28 15:25:35 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DParticleID_PID1.h"

//---------------------------------
// DParticleID_PID1    (Constructor)
//---------------------------------
DParticleID_PID1::DParticleID_PID1(JEventLoop *loop):DParticleID(loop)
{
	TF1 *locFunc;

	//Means
   ddEdxMeanFunc_FDC_Proton = new TF1("dPID_dEdxMeanFunc_FDC_Proton", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   ddEdxMeanFunc_FDC_Proton->SetParameters(47.5481, -4.70207, 2.22134, -0.357701);
	
   ddEdxMeanFunc_CDC_Proton = new TF1("dPID_dEdxMeanFunc_CDC_Proton", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   ddEdxMeanFunc_CDC_Proton->SetParameters(40.3217, -4.33611, 2.36034, -0.381628);

   ddEdxMeanFunc_FDC_KPlus = new TF1("dPID_dEdxMeanFunc_FDC_KPlus", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   ddEdxMeanFunc_FDC_KPlus->SetParameters(11.1221, -2.65063, 1.50373, -0.0202496);
	
   ddEdxMeanFunc_CDC_KPlus = new TF1("dPID_dEdxMeanFunc_CDC_KPlus", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   ddEdxMeanFunc_CDC_KPlus->SetParameters(12.3025, -2.51467, 1.53035, -0.0120788);

	//sigmas, protons
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(14.263, -3.91119, 0.471623, -0.0867709);
	ddEdxSigmaFuncVector_CDC_Proton.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(10.5091, -4.2701, 0.638759, -0.191789);
	ddEdxSigmaFuncVector_CDC_Proton.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(5.74486, -3.46166, 0.353456, -0.0594251);
	ddEdxSigmaFuncVector_CDC_Proton.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(6.93564, -4.55521, 0.372744, -0.0809587);
	ddEdxSigmaFuncVector_CDC_Proton.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton5", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(11.3186, -6.17783, 0.269247, -0.0525885);
	ddEdxSigmaFuncVector_CDC_Proton.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton6", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(8.41072, -5.90604, 0.23806, -0.0406181);
	ddEdxSigmaFuncVector_CDC_Proton.push_back(locFunc);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_Proton.push_back(4 + 2*loc_i);

   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_Proton1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(15.0407, -4.60489, 0.501318, -0.0913598);
	ddEdxSigmaFuncVector_FDC_Proton.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_Proton2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(26.6243, -5.92033, 0.384732, -0.0820965);
	ddEdxSigmaFuncVector_FDC_Proton.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_Proton3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(10.0415, -4.956, 0.29667, -0.0551989);
	ddEdxSigmaFuncVector_FDC_Proton.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_Proton4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(14.5991, -5.9726, 0.274239, -0.0546458);
	ddEdxSigmaFuncVector_FDC_Proton.push_back(locFunc);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_Proton.push_back(3 + 3*loc_i);

	//sigmas, pi+
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.8644, -1.92051, 0.401207, -0.00510045);
	ddEdxSigmaFuncVector_CDC_PiPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(3.08596, -2.55591, 0.326127, -0.00439305);
	ddEdxSigmaFuncVector_CDC_PiPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.63183, -2.18999, 0.276187, -0.00277214);
	ddEdxSigmaFuncVector_CDC_PiPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.42353, -1.12885, 0.218266, 0.000475549);
	ddEdxSigmaFuncVector_CDC_PiPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus5", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.599551, -1.27373, 0.15483, 0.00125788);
	ddEdxSigmaFuncVector_CDC_PiPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus6", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.02319, -2.15814, 0.146108, 0.00126295);
	ddEdxSigmaFuncVector_CDC_PiPlus.push_back(locFunc);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_PiPlus.push_back(4 + 2*loc_i);

   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_PiPlus1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.972821, -1.13778, 0.28177, 0.00414764);
	ddEdxSigmaFuncVector_FDC_PiPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_PiPlus2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.72011, -1.13622, 0.205834, 0.00248226);
	ddEdxSigmaFuncVector_FDC_PiPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_PiPlus3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.783602, -1.58954, 0.175932, 0.00139313);
	ddEdxSigmaFuncVector_FDC_PiPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_PiPlus4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.666107, -1.50611, 0.144197, 0.0019479);
	ddEdxSigmaFuncVector_FDC_PiPlus.push_back(locFunc);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_PiPlus.push_back(3 + 3*loc_i);


	//sigmas, k+
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(10.3149, -3.48486, 0.590082, -0.0757512);
	ddEdxSigmaFuncVector_CDC_KPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(10.9306, -3.49101, 0.443805, -0.0527794);
	ddEdxSigmaFuncVector_CDC_KPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(7.71149, -3.19237, 0.319949, -0.0248138);
	ddEdxSigmaFuncVector_CDC_KPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(2.70214, -2.51235, 0.224439, -0.00564491);
	ddEdxSigmaFuncVector_CDC_KPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus5", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.12165, -2.16051, 0.155717, 0.000195505);
	ddEdxSigmaFuncVector_CDC_KPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus6", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.60267, -2.62908, 0.153127, -0.000775554);
	ddEdxSigmaFuncVector_CDC_KPlus.push_back(locFunc);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_KPlus.push_back(4 + 2*loc_i);

   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_KPlus1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(3.4624, -2.51387, 0.328651, -0.00906703);
	ddEdxSigmaFuncVector_FDC_KPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_KPlus2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(8.77765, -3.16076, 0.20651, 0.00175279);
	ddEdxSigmaFuncVector_FDC_KPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_KPlus3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(6.44965, -3.01797, 0.175224, -0.000185835);
	ddEdxSigmaFuncVector_FDC_KPlus.push_back(locFunc);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_KPlus4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.37552, -2.21704, 0.148445, 0.0014621);
	ddEdxSigmaFuncVector_FDC_KPlus.push_back(locFunc);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_KPlus.push_back(3 + 3*loc_i);
}

//---------------------------------
// ~DParticleID_PID1    (Destructor)
//---------------------------------
DParticleID_PID1::~DParticleID_PID1()
{

}

jerror_t DParticleID_PID1::GetdEdxMean_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx, Particle_t locPIDHypothesis){
	double locBetaGammaValue = locBeta/sqrt(1.0 - locBeta*locBeta);
	if((locPIDHypothesis == Proton) || (locPIDHypothesis == AntiProton)){
		if(locBetaGammaValue < 0.3)
			return RESOURCE_UNAVAILABLE;
		if(locBetaGammaValue > 2.15)
			return RESOURCE_UNAVAILABLE;

		locMeandEdx = (ddEdxMeanFunc_CDC_Proton->Eval(locBetaGammaValue))/1000000.0;
		return NOERROR;
	}
	if((locPIDHypothesis == KPlus) || (locPIDHypothesis == KMinus)){
		if(locBetaGammaValue < 0.65)
			return RESOURCE_UNAVAILABLE; //K+ very likely to decay, don't compute chisq!!
		if(locBetaGammaValue > 4.1)
			return RESOURCE_UNAVAILABLE;

		locMeandEdx = (ddEdxMeanFunc_CDC_KPlus->Eval(locBetaGammaValue))/1000000.0;
		return NOERROR;
	}

	if((locPIDHypothesis == PiPlus) || (locPIDHypothesis == PiMinus)){
		double locBetaGamma[44];
		double locdEdxMean[44];
		locBetaGamma[0] = 0.48;  locBetaGamma[1] = 0.8;  locBetaGamma[2] = 1.12;  locBetaGamma[3] = 1.44;  locBetaGamma[4] = 1.76;  locBetaGamma[5] = 2.08;  locBetaGamma[6] = 2.4;  locBetaGamma[7] = 2.72;  locBetaGamma[8] = 3.04;  locBetaGamma[9] = 3.36;  locBetaGamma[10] = 3.68;  locBetaGamma[11] = 4;  locBetaGamma[12] = 4.32;  locBetaGamma[13] = 4.64;  locBetaGamma[14] = 4.96;  locBetaGamma[15] = 5.28;  locBetaGamma[16] = 5.6;  locBetaGamma[17] = 5.92;  locBetaGamma[18] = 6.24;  locBetaGamma[19] = 6.56;  locBetaGamma[20] = 6.88;  locBetaGamma[21] = 7.2;  locBetaGamma[22] = 7.52;  locBetaGamma[23] = 7.84;  locBetaGamma[24] = 8.16;  locBetaGamma[25] = 8.48;  locBetaGamma[26] = 8.8;  locBetaGamma[27] = 9.12;  locBetaGamma[28] = 9.44;  locBetaGamma[29] = 9.76;  locBetaGamma[30] = 10.08;  locBetaGamma[31] = 10.4;  locBetaGamma[32] = 10.72;  locBetaGamma[33] = 11.04;  locBetaGamma[34] = 11.36;  locBetaGamma[35] = 11.68;  locBetaGamma[36] = 12;  locBetaGamma[37] = 12.32;  locBetaGamma[38] = 12.64;  locBetaGamma[39] = 12.96;  locBetaGamma[40] = 13.28;  locBetaGamma[41] = 13.6;  locBetaGamma[42] = 13.92;  locBetaGamma[43] = 14.24;
		locdEdxMean[0] = 2.40705;  locdEdxMean[1] = 2.64016;  locdEdxMean[2] = 1.97719;  locdEdxMean[3] = 1.74434;  locdEdxMean[4] = 1.60673;  locdEdxMean[5] = 1.52877;  locdEdxMean[6] = 1.47792;  locdEdxMean[7] = 1.44914;  locdEdxMean[8] = 1.44456;  locdEdxMean[9] = 1.44089;  locdEdxMean[10] = 1.44077;  locdEdxMean[11] = 1.44643;  locdEdxMean[12] = 1.46221;  locdEdxMean[13] = 1.47146;  locdEdxMean[14] = 1.48451;  locdEdxMean[15] = 1.49671;  locdEdxMean[16] = 1.50718;  locdEdxMean[17] = 1.51892;  locdEdxMean[18] = 1.53035;  locdEdxMean[19] = 1.54224;  locdEdxMean[20] = 1.5531;  locdEdxMean[21] = 1.56325;  locdEdxMean[22] = 1.57347;  locdEdxMean[23] = 1.58422;  locdEdxMean[24] = 1.59239;  locdEdxMean[25] = 1.60057;  locdEdxMean[26] = 1.61224;  locdEdxMean[27] = 1.62025;  locdEdxMean[28] = 1.6288;  locdEdxMean[29] = 1.63748;  locdEdxMean[30] = 1.64582;  locdEdxMean[31] = 1.65254;  locdEdxMean[32] = 1.65864;  locdEdxMean[33] = 1.66693;  locdEdxMean[34] = 1.67503;  locdEdxMean[35] = 1.68328;  locdEdxMean[36] = 1.68964;  locdEdxMean[37] = 1.69732;  locdEdxMean[38] = 1.70258;  locdEdxMean[39] = 1.70849;  locdEdxMean[40] = 1.71469;  locdEdxMean[41] = 1.72359;  locdEdxMean[42] = 1.73032;  locdEdxMean[43] = 1.73672;

		if(locBetaGammaValue < 0.81)
			return RESOURCE_UNAVAILABLE;
		if(locBetaGammaValue > 14.2)
			return RESOURCE_UNAVAILABLE;

		double locSlope, locIntercept;
		for(unsigned int loc_i = 0; loc_i < 43; loc_i++){
			if(locBetaGammaValue > locBetaGamma[loc_i + 1])
				continue;
			locSlope = (locdEdxMean[loc_i + 1] - locdEdxMean[loc_i])/(locBetaGamma[loc_i + 1] - locBetaGamma[loc_i]);
			locIntercept = locdEdxMean[loc_i] - locSlope*locBetaGamma[loc_i]; //y = mx + b, b = y - mx
			locMeandEdx = (locSlope*locBetaGammaValue + locIntercept)/1000000.0;
			return NOERROR;
		}
	}

	return RESOURCE_UNAVAILABLE;
}

jerror_t DParticleID_PID1::GetdEdxSigma_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx, Particle_t locPIDHypothesis){
	double locBetaGamma = locBeta/sqrt(1.0 - locBeta*locBeta);
	double locSigmadEdx_LowSide, locSigmadEdx_HighSide; //for linear interpolation/extrapolation

	vector<unsigned int> locNumHitsVector;
	vector<TF1*> locdEdxSigmaVector;
	if((locPIDHypothesis == Proton) || (locPIDHypothesis == AntiProton)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_Proton;
		locdEdxSigmaVector = ddEdxSigmaFuncVector_CDC_Proton;
	}
	if((locPIDHypothesis == KPlus) || (locPIDHypothesis == KMinus)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_KPlus;
		locdEdxSigmaVector = ddEdxSigmaFuncVector_CDC_KPlus;
	}
	if((locPIDHypothesis == PiPlus) || (locPIDHypothesis == PiMinus)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_PiPlus;
		locdEdxSigmaVector = ddEdxSigmaFuncVector_CDC_PiPlus;
	}

	double locSlope, locIntercept;

	for(unsigned int loc_i = 0; loc_i < locNumHitsVector.size(); loc_i++){
		if(locNumHitsUsedFordEdx == locNumHitsVector[loc_i]){
			locSigmadEdx = ((locdEdxSigmaVector[loc_i])->Eval(locBetaGamma))/1000000.0;
			return NOERROR;
		}
		if(loc_i == (locNumHitsVector.size() - 1)){
			locSigmadEdx_LowSide = (locdEdxSigmaVector[loc_i - 1])->Eval(locBetaGamma);
			locSigmadEdx_HighSide = (locdEdxSigmaVector[loc_i])->Eval(locBetaGamma);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i] - locNumHitsVector[loc_i - 1]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return NOERROR;
		}
		if((locNumHitsUsedFordEdx > locNumHitsVector[loc_i]) && (locNumHitsUsedFordEdx < locNumHitsVector[loc_i + 1])){
			locSigmadEdx_LowSide = (locdEdxSigmaVector[loc_i])->Eval(locBetaGamma);
			locSigmadEdx_HighSide = (locdEdxSigmaVector[loc_i + 1])->Eval(locBetaGamma);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i + 1] - locNumHitsVector[loc_i]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i + 1]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return NOERROR;
		}
	}

	return NOERROR;
}

jerror_t DParticleID_PID1::GetdEdxMean_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx, Particle_t locPIDHypothesis){
	double locBetaGammaValue = locBeta/sqrt(1.0 - locBeta*locBeta);

	if((locPIDHypothesis == Proton) || (locPIDHypothesis == AntiProton)){
		if(locBetaGammaValue < 0.3)
			return RESOURCE_UNAVAILABLE;
		if(locBetaGammaValue > 2.15)
			return RESOURCE_UNAVAILABLE;

		locMeandEdx = (ddEdxMeanFunc_FDC_Proton->Eval(locBetaGammaValue))/1000000.0;
		return NOERROR;
	}
	if((locPIDHypothesis == KPlus) || (locPIDHypothesis == KMinus)){
		if(locBetaGammaValue < 0.9)
			return RESOURCE_UNAVAILABLE; //K+ very likely to decay, don't compute chisq!!
		if(locBetaGammaValue > 4.1)
			return RESOURCE_UNAVAILABLE;

		locMeandEdx = (ddEdxMeanFunc_FDC_KPlus->Eval(locBetaGammaValue))/1000000.0;
		return NOERROR;
	}

	if((locPIDHypothesis == PiPlus) || (locPIDHypothesis == PiMinus)){
		double locBetaGamma[44];
		double locdEdxMean[44];
		locBetaGamma[0] = 0.48;  locBetaGamma[1] = 0.8;  locBetaGamma[2] = 1.12;  locBetaGamma[3] = 1.44;  locBetaGamma[4] = 1.76;  locBetaGamma[5] = 2.08;  locBetaGamma[6] = 2.4;  locBetaGamma[7] = 2.72;  locBetaGamma[8] = 3.04;  locBetaGamma[9] = 3.36;  locBetaGamma[10] = 3.68;  locBetaGamma[11] = 4;  locBetaGamma[12] = 4.32;  locBetaGamma[13] = 4.64;  locBetaGamma[14] = 4.96;  locBetaGamma[15] = 5.28;  locBetaGamma[16] = 5.6;  locBetaGamma[17] = 5.92;  locBetaGamma[18] = 6.24;  locBetaGamma[19] = 6.56;  locBetaGamma[20] = 6.88;  locBetaGamma[21] = 7.2;  locBetaGamma[22] = 7.52;  locBetaGamma[23] = 7.84;  locBetaGamma[24] = 8.16;  locBetaGamma[25] = 8.48;  locBetaGamma[26] = 8.8;  locBetaGamma[27] = 9.12;  locBetaGamma[28] = 9.44;  locBetaGamma[29] = 9.76;  locBetaGamma[30] = 10.08;  locBetaGamma[31] = 10.4;  locBetaGamma[32] = 10.72;  locBetaGamma[33] = 11.04;  locBetaGamma[34] = 11.36;  locBetaGamma[35] = 11.68;  locBetaGamma[36] = 12;  locBetaGamma[37] = 12.32;  locBetaGamma[38] = 12.64;  locBetaGamma[39] = 12.96;  locBetaGamma[40] = 13.28;  locBetaGamma[41] = 13.6;  locBetaGamma[42] = 13.92;  locBetaGamma[43] = 14.24;
		locdEdxMean[0] = 1.32105;  locdEdxMean[1] = 1.50904;  locdEdxMean[2] = 1.50873;  locdEdxMean[3] = 1.57426;  locdEdxMean[4] = 1.52213;  locdEdxMean[5] = 1.45914;  locdEdxMean[6] = 1.4271;  locdEdxMean[7] = 1.42623;  locdEdxMean[8] = 1.42268;  locdEdxMean[9] = 1.42632;  locdEdxMean[10] = 1.43233;  locdEdxMean[11] = 1.43742;  locdEdxMean[12] = 1.4467;  locdEdxMean[13] = 1.45642;  locdEdxMean[14] = 1.46337;  locdEdxMean[15] = 1.4719;  locdEdxMean[16] = 1.48252;  locdEdxMean[17] = 1.49443;  locdEdxMean[18] = 1.50234;  locdEdxMean[19] = 1.513;  locdEdxMean[20] = 1.52394;  locdEdxMean[21] = 1.5289;  locdEdxMean[22] = 1.53878;  locdEdxMean[23] = 1.54396;  locdEdxMean[24] = 1.55427;  locdEdxMean[25] = 1.56105;  locdEdxMean[26] = 1.5674;  locdEdxMean[27] = 1.57221;  locdEdxMean[28] = 1.58027;  locdEdxMean[29] = 1.5865;  locdEdxMean[30] = 1.59519;  locdEdxMean[31] = 1.60211;  locdEdxMean[32] = 1.60751;  locdEdxMean[33] = 1.61522;  locdEdxMean[34] = 1.62066;  locdEdxMean[35] = 1.6277;  locdEdxMean[36] = 1.63227;  locdEdxMean[37] = 1.63661;  locdEdxMean[38] = 1.64457;  locdEdxMean[39] = 1.65043;  locdEdxMean[40] = 1.6557;  locdEdxMean[41] = 1.66085;  locdEdxMean[42] = 1.66486;  locdEdxMean[43] = 1.66916;

		if(locBetaGammaValue < 1.77)
			return RESOURCE_UNAVAILABLE;
		if(locBetaGammaValue > 14.2)
			return RESOURCE_UNAVAILABLE;

		double locSlope, locIntercept;
		for(unsigned int loc_i = 0; loc_i < 43; loc_i++){
			if(locBetaGammaValue > locBetaGamma[loc_i + 1])
				continue;
			locSlope = (locdEdxMean[loc_i + 1] - locdEdxMean[loc_i])/(locBetaGamma[loc_i + 1] - locBetaGamma[loc_i]);
			locIntercept = locdEdxMean[loc_i] - locSlope*locBetaGamma[loc_i]; //y = mx + b, b = y - mx
			locMeandEdx = (locSlope*locBetaGammaValue + locIntercept)/1000000.0;
			return NOERROR;
		}
	}

	return RESOURCE_UNAVAILABLE;;
}

jerror_t DParticleID_PID1::GetdEdxSigma_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx, Particle_t locPIDHypothesis){
	double locBetaGamma = locBeta/sqrt(1.0 - locBeta*locBeta);
	double locSigmadEdx_LowSide, locSigmadEdx_HighSide; //for linear interpolation/extrapolation

	vector<unsigned int> locNumHitsVector;
	vector<TF1*> locdEdxSigmaVector;
	if((locPIDHypothesis == Proton) || (locPIDHypothesis == AntiProton)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_Proton;
		locdEdxSigmaVector = ddEdxSigmaFuncVector_FDC_Proton;
	}
	if((locPIDHypothesis == KPlus) || (locPIDHypothesis == KMinus)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_KPlus;
		locdEdxSigmaVector = ddEdxSigmaFuncVector_FDC_KPlus;
	}
	if((locPIDHypothesis == PiPlus) || (locPIDHypothesis == PiMinus)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_PiPlus;
		locdEdxSigmaVector = ddEdxSigmaFuncVector_FDC_PiPlus;
	}

	double locSlope, locIntercept;

	for(unsigned int loc_i = 0; loc_i < locNumHitsVector.size(); loc_i++){
		if(locNumHitsUsedFordEdx == locNumHitsVector[loc_i]){
			locSigmadEdx = ((locdEdxSigmaVector[loc_i])->Eval(locBetaGamma))/1000000.0;
			return NOERROR;
		}
		if(loc_i == (locNumHitsVector.size() - 1)){
			locSigmadEdx_LowSide = (locdEdxSigmaVector[loc_i - 1])->Eval(locBetaGamma);
			locSigmadEdx_HighSide = (locdEdxSigmaVector[loc_i])->Eval(locBetaGamma);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i] - locNumHitsVector[loc_i - 1]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return NOERROR;
		}
		if((locNumHitsUsedFordEdx > locNumHitsVector[loc_i]) && (locNumHitsUsedFordEdx < locNumHitsVector[loc_i + 1])){
			locSigmadEdx_LowSide = (locdEdxSigmaVector[loc_i])->Eval(locBetaGamma);
			locSigmadEdx_HighSide = (locdEdxSigmaVector[loc_i + 1])->Eval(locBetaGamma);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i + 1] - locNumHitsVector[loc_i]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i + 1]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return NOERROR;
		}
	}
	return NOERROR;
}

jerror_t DParticleID_PID1::CalcDCdEdxChiSq(const DChargedTrackHypothesis *locChargedTrackHypothesis, double &locChiSq, unsigned int& locNDF){
	locNDF = 0;
	locChiSq = NaN;
	unsigned int locMinimumNumberUsedHitsForConfidence = 3; //dE/dx is landau-distributed, so to approximate Gaussian must remove hits with largest dE/dx //3 means 6 or more hits originally
	Particle_t locPID = locChargedTrackHypothesis->dPID;

	const DTrackTimeBased *locTrackTimeBased = locChargedTrackHypothesis->dTrackTimeBased;
	unsigned int locNumHitsUsedFordEdx_CDC = locTrackTimeBased->dNumHitsUsedFordEdx_CDC;
	unsigned int locNumHitsUsedFordEdx_FDC = locTrackTimeBased->dNumHitsUsedFordEdx_FDC;

	bool locUseCDCHitsFlag = (locNumHitsUsedFordEdx_CDC > locMinimumNumberUsedHitsForConfidence) ? true : false;
	bool locUseFDCHitsFlag = (locNumHitsUsedFordEdx_FDC > locMinimumNumberUsedHitsForConfidence) ? true : false;

	if((locUseCDCHitsFlag == false) && (locUseFDCHitsFlag == false))
		return RESOURCE_UNAVAILABLE; //not enough hits, use other sources of information for PID

	double locDCdEdx = locTrackTimeBased->dEdx();
	double locMeandEdx_FDC, locMeandEdx_CDC, locSigmadEdx_FDC, locSigmadEdx_CDC;
	double locDeltadEdx_CDC = 0.0, locDeltadEdx_FDC = 0.0;

	double locBeta = locChargedTrackHypothesis->beta();

	if(locNumHitsUsedFordEdx_CDC >= locNumHitsUsedFordEdx_FDC){
		if(GetdEdxMean_CDC(locBeta, locNumHitsUsedFordEdx_CDC, locMeandEdx_CDC, locPID) != NOERROR)
			return RESOURCE_UNAVAILABLE;
		if(GetdEdxSigma_CDC(locBeta, locNumHitsUsedFordEdx_CDC, locSigmadEdx_CDC, locPID) != NOERROR)
			return RESOURCE_UNAVAILABLE;
		locDeltadEdx_CDC = locDCdEdx - locMeandEdx_CDC;
		locChiSq = locDeltadEdx_CDC/locSigmadEdx_CDC;
		locChiSq *= locChiSq;
		locNDF = 1;
		return NOERROR;
	}

	if(locNumHitsUsedFordEdx_FDC > locNumHitsUsedFordEdx_CDC){
		if(GetdEdxMean_FDC(locBeta, locNumHitsUsedFordEdx_FDC, locMeandEdx_FDC, locPID) != NOERROR)
			return RESOURCE_UNAVAILABLE;
		if(GetdEdxSigma_FDC(locBeta, locNumHitsUsedFordEdx_FDC, locSigmadEdx_FDC, locPID) != NOERROR)
			return RESOURCE_UNAVAILABLE;
		locDeltadEdx_FDC = locDCdEdx - locMeandEdx_FDC;
		locChiSq = locDeltadEdx_FDC/locSigmadEdx_FDC;
		locChiSq *= locChiSq;
		locNDF = 1;
		return NOERROR;
	}

	return NOERROR;
}


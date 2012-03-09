// $Id$
//
//    File: DParticleID_PID1.cc
// Created: Mon Feb 28 15:25:35 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DParticleID_PID1.h"

void DParticleID_PID1::Set_dEdxParams(vector<float>& locParamVector, float locParam1, float locParam2, float locParam3, float locParam4)
{
	locParamVector.resize(4);
	locParamVector[0] = locParam1;
	locParamVector[1] = locParam2;
	locParamVector[2] = locParam3;
	locParamVector[3] = locParam4;
}

//---------------------------------
// DParticleID_PID1    (Constructor)
//---------------------------------
DParticleID_PID1::DParticleID_PID1(JEventLoop *loop):DParticleID(loop)
{

	vector<float> locParamVector;

	//Means
	Set_dEdxParams(ddEdxMeanParams_FDC_Proton, 41.6301, -4.51759, 1.98461, -0.29211);
	Set_dEdxParams(ddEdxMeanParams_CDC_Proton, 46.9116, -4.66336, 2.32353, -0.444754);
	Set_dEdxParams(ddEdxMeanParams_FDC_KPlus, 13.4017, -2.92603, 1.41205, -0.0283139);
	Set_dEdxParams(ddEdxMeanParams_CDC_KPlus, 11.3948, -2.47605, 1.40824, -0.0309194);

	//sigmas, protons, CDC
	Set_dEdxParams(locParamVector, 9.14756, -3.25359, 0.367987, -0.0290751);
	ddEdxSigmaParamVector_CDC_Proton.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 8.39024, -3.98806, 0.608273, -0.166997);
	ddEdxSigmaParamVector_CDC_Proton.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 7.90078, -4.12579, 0.534348, -0.146632);
	ddEdxSigmaParamVector_CDC_Proton.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 9.25127, -5.00809, 0.440321, -0.114403);
	ddEdxSigmaParamVector_CDC_Proton.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 11.0353, -6.00159, 0.290354, -0.0621569);
	ddEdxSigmaParamVector_CDC_Proton.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 8.83402, -5.93664, 0.245129, -0.0439422);
	ddEdxSigmaParamVector_CDC_Proton.push_back(locParamVector);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_Proton.push_back(4 + 2*loc_i);

	//sigmas, protons, FDC
	Set_dEdxParams(locParamVector, 28.0219, -5.54908, 0.501065, -0.104623);
	ddEdxSigmaParamVector_FDC_Proton.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 8.95707, -4.22215, 0.300281, -0.0436824);
	ddEdxSigmaParamVector_FDC_Proton.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 11.9735, -4.92211, 0.274352, -0.0472131);
	ddEdxSigmaParamVector_FDC_Proton.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 14.2884, -5.91157, 0.274532, -0.0602249);
	ddEdxSigmaParamVector_FDC_Proton.push_back(locParamVector);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_Proton.push_back(3 + 3*loc_i);

	//sigmas, pi+, CDC
	Set_dEdxParams(locParamVector, 2.22871, -2.3159, 0.379665, -0.00114537);
	ddEdxSigmaParamVector_CDC_PiPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 1.05578, -2.02632, 0.287603, 0.000438264);
	ddEdxSigmaParamVector_CDC_PiPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 0.786585, -1.59372, 0.237528, 0.0013845);
	ddEdxSigmaParamVector_CDC_PiPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 0.413633, -1.23667, 0.192888, 0.00222538);
	ddEdxSigmaParamVector_CDC_PiPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 0.488919, -1.33126, 0.154777, 0.000816022);
	ddEdxSigmaParamVector_CDC_PiPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 0.827528, -2.07871, 0.142669, 0.00105942);
	ddEdxSigmaParamVector_CDC_PiPlus.push_back(locParamVector);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_PiPlus.push_back(4 + 2*loc_i);

	//sigmas, pi+, FDC
	Set_dEdxParams(locParamVector, 0.655229, -0.898369, 0.256156, 0.002824);
	ddEdxSigmaParamVector_FDC_PiPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 0.51879, -0.817423, 0.16759, 0.00399836);
	ddEdxSigmaParamVector_FDC_PiPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 0.49766, -1.16389, 0.155798, 0.0011712);
	ddEdxSigmaParamVector_FDC_PiPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 0.585499, -1.32707, 0.130717, 0.00185828);
	ddEdxSigmaParamVector_FDC_PiPlus.push_back(locParamVector);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_PiPlus.push_back(3 + 3*loc_i);


	//sigmas, k+, CDC
	Set_dEdxParams(locParamVector, 1.08389, -0.727214, 0.101523, 0.0410826);
	ddEdxSigmaParamVector_CDC_KPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 8.3022, -3.45934, 0.399704, -0.0367945);
	ddEdxSigmaParamVector_CDC_KPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 6.15751, -3.06247, 0.312122, -0.0232223);
	ddEdxSigmaParamVector_CDC_KPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 2.34946, -2.25875, 0.213915, -0.0027382);
	ddEdxSigmaParamVector_CDC_KPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 1.66756, -2.43701, 0.167655, -0.00349099);
	ddEdxSigmaParamVector_CDC_KPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 1.5591, -2.59352, 0.158265, -0.00298591);
	ddEdxSigmaParamVector_CDC_KPlus.push_back(locParamVector);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_KPlus.push_back(4 + 2*loc_i);

	//sigmas, k+, FDC
	Set_dEdxParams(locParamVector, 3.15915, -2.61119, 0.293622, -0.00531579);
	ddEdxSigmaParamVector_FDC_KPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 19.4517, -4.05641, 0.226242, -0.00897102);
	ddEdxSigmaParamVector_FDC_KPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 21.846, -4.08992, 0.188548, -0.00754031);
	ddEdxSigmaParamVector_FDC_KPlus.push_back(locParamVector);
	Set_dEdxParams(locParamVector, 2.81826, -2.90351, 0.150354, -0.00140537);
	ddEdxSigmaParamVector_FDC_KPlus.push_back(locParamVector);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_KPlus.push_back(3 + 3*loc_i);

	dBetaGamma_PiMinus_CDC.resize(43);
	ddEdxMean_PiMinus_CDC.resize(43);
	dBetaGamma_PiMinus_CDC[0] = 0.8;  dBetaGamma_PiMinus_CDC[1] = 1.12;  dBetaGamma_PiMinus_CDC[2] = 1.44;  dBetaGamma_PiMinus_CDC[3] = 1.76;  dBetaGamma_PiMinus_CDC[4] = 2.08;  dBetaGamma_PiMinus_CDC[5] = 2.4;  dBetaGamma_PiMinus_CDC[6] = 2.72;  dBetaGamma_PiMinus_CDC[7] = 3.04;  dBetaGamma_PiMinus_CDC[8] = 3.36;  dBetaGamma_PiMinus_CDC[9] = 3.68;  dBetaGamma_PiMinus_CDC[10] = 4;  dBetaGamma_PiMinus_CDC[11] = 4.32;  dBetaGamma_PiMinus_CDC[12] = 4.64;  dBetaGamma_PiMinus_CDC[13] = 4.96;  dBetaGamma_PiMinus_CDC[14] = 5.28;  dBetaGamma_PiMinus_CDC[15] = 5.6;  dBetaGamma_PiMinus_CDC[16] = 5.92;  dBetaGamma_PiMinus_CDC[17] = 6.24;  dBetaGamma_PiMinus_CDC[18] = 6.56;  dBetaGamma_PiMinus_CDC[19] = 6.88;  dBetaGamma_PiMinus_CDC[20] = 7.2;  dBetaGamma_PiMinus_CDC[21] = 7.52;  dBetaGamma_PiMinus_CDC[22] = 7.84;  dBetaGamma_PiMinus_CDC[23] = 8.16;  dBetaGamma_PiMinus_CDC[24] = 8.48;  dBetaGamma_PiMinus_CDC[25] = 8.8;  dBetaGamma_PiMinus_CDC[26] = 9.12;  dBetaGamma_PiMinus_CDC[27] = 9.44;  dBetaGamma_PiMinus_CDC[28] = 9.76;  dBetaGamma_PiMinus_CDC[29] = 10.08;  dBetaGamma_PiMinus_CDC[30] = 10.4;  dBetaGamma_PiMinus_CDC[31] = 10.72;  dBetaGamma_PiMinus_CDC[32] = 11.04;  dBetaGamma_PiMinus_CDC[33] = 11.36;  dBetaGamma_PiMinus_CDC[34] = 11.68;  dBetaGamma_PiMinus_CDC[35] = 12;  dBetaGamma_PiMinus_CDC[36] = 12.32;  dBetaGamma_PiMinus_CDC[37] = 12.64;  dBetaGamma_PiMinus_CDC[38] = 12.96;  dBetaGamma_PiMinus_CDC[39] = 13.28;  dBetaGamma_PiMinus_CDC[40] = 13.6;  dBetaGamma_PiMinus_CDC[41] = 13.92;  dBetaGamma_PiMinus_CDC[42] = 14.24;
	ddEdxMean_PiMinus_CDC[0] = 2.37646;  ddEdxMean_PiMinus_CDC[1] = 1.83021;  ddEdxMean_PiMinus_CDC[2] = 1.59472;  ddEdxMean_PiMinus_CDC[3] = 1.46196;  ddEdxMean_PiMinus_CDC[4] = 1.36924;  ddEdxMean_PiMinus_CDC[5] = 1.32853;  ddEdxMean_PiMinus_CDC[6] = 1.30056;  ddEdxMean_PiMinus_CDC[7] = 1.27378;  ddEdxMean_PiMinus_CDC[8] = 1.25816;  ddEdxMean_PiMinus_CDC[9] = 1.25344;  ddEdxMean_PiMinus_CDC[10] = 1.25623;  ddEdxMean_PiMinus_CDC[11] = 1.26464;  ddEdxMean_PiMinus_CDC[12] = 1.26824;  ddEdxMean_PiMinus_CDC[13] = 1.27456;  ddEdxMean_PiMinus_CDC[14] = 1.27943;  ddEdxMean_PiMinus_CDC[15] = 1.284;  ddEdxMean_PiMinus_CDC[16] = 1.28983;  ddEdxMean_PiMinus_CDC[17] = 1.2961;  ddEdxMean_PiMinus_CDC[18] = 1.30345;  ddEdxMean_PiMinus_CDC[19] = 1.30808;  ddEdxMean_PiMinus_CDC[20] = 1.31422;  ddEdxMean_PiMinus_CDC[21] = 1.31972;  ddEdxMean_PiMinus_CDC[22] = 1.32548;  ddEdxMean_PiMinus_CDC[23] = 1.32913;  ddEdxMean_PiMinus_CDC[24] = 1.3367;  ddEdxMean_PiMinus_CDC[25] = 1.34033;  ddEdxMean_PiMinus_CDC[26] = 1.34922;  ddEdxMean_PiMinus_CDC[27] = 1.35131;  ddEdxMean_PiMinus_CDC[28] = 1.35885;  ddEdxMean_PiMinus_CDC[29] = 1.36266;  ddEdxMean_PiMinus_CDC[30] = 1.36663;  ddEdxMean_PiMinus_CDC[31] = 1.37086;  ddEdxMean_PiMinus_CDC[32] = 1.37426;  ddEdxMean_PiMinus_CDC[33] = 1.37922;  ddEdxMean_PiMinus_CDC[34] = 1.3839;  ddEdxMean_PiMinus_CDC[35] = 1.38708;  ddEdxMean_PiMinus_CDC[36] = 1.39205;  ddEdxMean_PiMinus_CDC[37] = 1.39559;  ddEdxMean_PiMinus_CDC[38] = 1.40163;  ddEdxMean_PiMinus_CDC[39] = 1.40517;  ddEdxMean_PiMinus_CDC[40] = 1.40787;  ddEdxMean_PiMinus_CDC[41] = 1.41131;  ddEdxMean_PiMinus_CDC[42] = 1.41184;

	dBetaGamma_PiMinus_FDC.resize(44);
	ddEdxMean_PiMinus_FDC.resize(44);
	dBetaGamma_PiMinus_FDC[0] = 0.48;  dBetaGamma_PiMinus_FDC[1] = 0.8;  dBetaGamma_PiMinus_FDC[2] = 1.12;  dBetaGamma_PiMinus_FDC[3] = 1.44;  dBetaGamma_PiMinus_FDC[4] = 1.76;  dBetaGamma_PiMinus_FDC[5] = 2.08;  dBetaGamma_PiMinus_FDC[6] = 2.4;  dBetaGamma_PiMinus_FDC[7] = 2.72;  dBetaGamma_PiMinus_FDC[8] = 3.04;  dBetaGamma_PiMinus_FDC[9] = 3.36;  dBetaGamma_PiMinus_FDC[10] = 3.68;  dBetaGamma_PiMinus_FDC[11] = 4;  dBetaGamma_PiMinus_FDC[12] = 4.32;  dBetaGamma_PiMinus_FDC[13] = 4.64;  dBetaGamma_PiMinus_FDC[14] = 4.96;  dBetaGamma_PiMinus_FDC[15] = 5.28;  dBetaGamma_PiMinus_FDC[16] = 5.6;  dBetaGamma_PiMinus_FDC[17] = 5.92;  dBetaGamma_PiMinus_FDC[18] = 6.24;  dBetaGamma_PiMinus_FDC[19] = 6.56;  dBetaGamma_PiMinus_FDC[20] = 6.88;  dBetaGamma_PiMinus_FDC[21] = 7.2;  dBetaGamma_PiMinus_FDC[22] = 7.52;  dBetaGamma_PiMinus_FDC[23] = 7.84;  dBetaGamma_PiMinus_FDC[24] = 8.16;  dBetaGamma_PiMinus_FDC[25] = 8.48;  dBetaGamma_PiMinus_FDC[26] = 8.8;  dBetaGamma_PiMinus_FDC[27] = 9.12;  dBetaGamma_PiMinus_FDC[28] = 9.44;  dBetaGamma_PiMinus_FDC[29] = 9.76;  dBetaGamma_PiMinus_FDC[30] = 10.08;  dBetaGamma_PiMinus_FDC[31] = 10.4;  dBetaGamma_PiMinus_FDC[32] = 10.72;  dBetaGamma_PiMinus_FDC[33] = 11.04;  dBetaGamma_PiMinus_FDC[34] = 11.36;  dBetaGamma_PiMinus_FDC[35] = 11.68;  dBetaGamma_PiMinus_FDC[36] = 12;  dBetaGamma_PiMinus_FDC[37] = 12.32;  dBetaGamma_PiMinus_FDC[38] = 12.64;  dBetaGamma_PiMinus_FDC[39] = 12.96;  dBetaGamma_PiMinus_FDC[40] = 13.28;  dBetaGamma_PiMinus_FDC[41] = 13.6;  dBetaGamma_PiMinus_FDC[42] = 13.92;  dBetaGamma_PiMinus_FDC[43] = 14.24;
	ddEdxMean_PiMinus_FDC[0] = 1.37523;  ddEdxMean_PiMinus_FDC[1] = 1.47618;  ddEdxMean_PiMinus_FDC[2] = 1.59352;  ddEdxMean_PiMinus_FDC[3] = 1.47423;  ddEdxMean_PiMinus_FDC[4] = 1.40155;  ddEdxMean_PiMinus_FDC[5] = 1.34673;  ddEdxMean_PiMinus_FDC[6] = 1.31879;  ddEdxMean_PiMinus_FDC[7] = 1.31466;  ddEdxMean_PiMinus_FDC[8] = 1.31048;  ddEdxMean_PiMinus_FDC[9] = 1.3059;  ddEdxMean_PiMinus_FDC[10] = 1.30657;  ddEdxMean_PiMinus_FDC[11] = 1.3071;  ddEdxMean_PiMinus_FDC[12] = 1.30768;  ddEdxMean_PiMinus_FDC[13] = 1.31311;  ddEdxMean_PiMinus_FDC[14] = 1.3117;  ddEdxMean_PiMinus_FDC[15] = 1.32219;  ddEdxMean_PiMinus_FDC[16] = 1.32418;  ddEdxMean_PiMinus_FDC[17] = 1.33116;  ddEdxMean_PiMinus_FDC[18] = 1.33479;  ddEdxMean_PiMinus_FDC[19] = 1.33449;  ddEdxMean_PiMinus_FDC[20] = 1.34254;  ddEdxMean_PiMinus_FDC[21] = 1.3424;  ddEdxMean_PiMinus_FDC[22] = 1.3449;  ddEdxMean_PiMinus_FDC[23] = 1.34964;  ddEdxMean_PiMinus_FDC[24] = 1.35254;  ddEdxMean_PiMinus_FDC[25] = 1.35634;  ddEdxMean_PiMinus_FDC[26] = 1.36226;  ddEdxMean_PiMinus_FDC[27] = 1.36161;  ddEdxMean_PiMinus_FDC[28] = 1.36486;  ddEdxMean_PiMinus_FDC[29] = 1.36825;  ddEdxMean_PiMinus_FDC[30] = 1.37269;  ddEdxMean_PiMinus_FDC[31] = 1.37568;  ddEdxMean_PiMinus_FDC[32] = 1.38051;  ddEdxMean_PiMinus_FDC[33] = 1.38383;  ddEdxMean_PiMinus_FDC[34] = 1.38167;  ddEdxMean_PiMinus_FDC[35] = 1.3852;  ddEdxMean_PiMinus_FDC[36] = 1.39147;  ddEdxMean_PiMinus_FDC[37] = 1.39506;  ddEdxMean_PiMinus_FDC[38] = 1.39391;  ddEdxMean_PiMinus_FDC[39] = 1.39701;  ddEdxMean_PiMinus_FDC[40] = 1.39823;  ddEdxMean_PiMinus_FDC[41] = 1.40486;  ddEdxMean_PiMinus_FDC[42] = 1.40764;  ddEdxMean_PiMinus_FDC[43] = 1.40711;
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

		locMeandEdx = Function_dEdx(locBetaGammaValue, ddEdxMeanParams_CDC_Proton)/1000000.0;
		return NOERROR;
	}
	if((locPIDHypothesis == KPlus) || (locPIDHypothesis == KMinus)){
		if(locBetaGammaValue < 0.65)
			return RESOURCE_UNAVAILABLE; //K+ very likely to decay, don't compute chisq!!
		if(locBetaGammaValue > 4.1)
			return RESOURCE_UNAVAILABLE;

		locMeandEdx = Function_dEdx(locBetaGammaValue, ddEdxMeanParams_CDC_KPlus)/1000000.0;
		return NOERROR;
	}

	if((locPIDHypothesis == PiPlus) || (locPIDHypothesis == PiMinus)){

		if(locBetaGammaValue < 0.81)
			return RESOURCE_UNAVAILABLE;
		if(locBetaGammaValue > 14.2)
			return RESOURCE_UNAVAILABLE;

		double locSlope, locIntercept;
		for(unsigned int loc_i = 0; loc_i < 43; loc_i++){
			if(locBetaGammaValue > dBetaGamma_PiMinus_CDC[loc_i + 1])
				continue;
			locSlope = (ddEdxMean_PiMinus_CDC[loc_i + 1] - ddEdxMean_PiMinus_CDC[loc_i])/(dBetaGamma_PiMinus_CDC[loc_i + 1] - dBetaGamma_PiMinus_CDC[loc_i]);
			locIntercept = ddEdxMean_PiMinus_CDC[loc_i] - locSlope*dBetaGamma_PiMinus_CDC[loc_i]; //y = mx + b, b = y - mx
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
	vector<vector<float> > locdEdxSigmaParamVector;
	if((locPIDHypothesis == Proton) || (locPIDHypothesis == AntiProton)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_Proton;
		locdEdxSigmaParamVector = ddEdxSigmaParamVector_CDC_Proton;
	}
	if((locPIDHypothesis == KPlus) || (locPIDHypothesis == KMinus)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_KPlus;
		locdEdxSigmaParamVector = ddEdxSigmaParamVector_CDC_KPlus;
	}
	if((locPIDHypothesis == PiPlus) || (locPIDHypothesis == PiMinus)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_PiPlus;
		locdEdxSigmaParamVector = ddEdxSigmaParamVector_CDC_PiPlus;
	}

	double locSlope, locIntercept;

	for(unsigned int loc_i = 0; loc_i < locNumHitsVector.size(); loc_i++){
		if(locNumHitsUsedFordEdx == locNumHitsVector[loc_i]){
			locSigmadEdx = (Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i]))/1000000.0;
			return NOERROR;
		}
		if(loc_i == (locNumHitsVector.size() - 1)){
			locSigmadEdx_LowSide = Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i - 1]);
			locSigmadEdx_HighSide = Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i]);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i] - locNumHitsVector[loc_i - 1]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return NOERROR;
		}
		if((locNumHitsUsedFordEdx > locNumHitsVector[loc_i]) && (locNumHitsUsedFordEdx < locNumHitsVector[loc_i + 1])){
			locSigmadEdx_LowSide = Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i]);
			locSigmadEdx_HighSide = Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i + 1]);

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

		locMeandEdx = Function_dEdx(locBetaGammaValue, ddEdxMeanParams_FDC_Proton)/1000000.0;
		return NOERROR;
	}
	if((locPIDHypothesis == KPlus) || (locPIDHypothesis == KMinus)){
		if(locBetaGammaValue < 1.1)
			return RESOURCE_UNAVAILABLE; //K+ very likely to decay, don't compute chisq!!
		if(locBetaGammaValue > 4.1)
			return RESOURCE_UNAVAILABLE;

		locMeandEdx = Function_dEdx(locBetaGammaValue, ddEdxMeanParams_FDC_KPlus)/1000000.0;
		return NOERROR;
	}

	if((locPIDHypothesis == PiPlus) || (locPIDHypothesis == PiMinus)){

		if(locBetaGammaValue < 1.77)
			return RESOURCE_UNAVAILABLE;
		if(locBetaGammaValue > 14.2)
			return RESOURCE_UNAVAILABLE;

		double locSlope, locIntercept;
		for(unsigned int loc_i = 0; loc_i < 43; loc_i++){
			if(locBetaGammaValue > dBetaGamma_PiMinus_FDC[loc_i + 1])
				continue;
			locSlope = (ddEdxMean_PiMinus_FDC[loc_i + 1] - ddEdxMean_PiMinus_FDC[loc_i])/(dBetaGamma_PiMinus_FDC[loc_i + 1] - dBetaGamma_PiMinus_FDC[loc_i]);
			locIntercept = ddEdxMean_PiMinus_FDC[loc_i] - locSlope*dBetaGamma_PiMinus_FDC[loc_i]; //y = mx + b, b = y - mx
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
	vector<vector<float> > locdEdxSigmaParamVector;
	if((locPIDHypothesis == Proton) || (locPIDHypothesis == AntiProton)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_Proton;
		locdEdxSigmaParamVector = ddEdxSigmaParamVector_FDC_Proton;
	}
	if((locPIDHypothesis == KPlus) || (locPIDHypothesis == KMinus)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_KPlus;
		locdEdxSigmaParamVector = ddEdxSigmaParamVector_FDC_KPlus;
	}
	if((locPIDHypothesis == PiPlus) || (locPIDHypothesis == PiMinus)){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_PiPlus;
		locdEdxSigmaParamVector = ddEdxSigmaParamVector_FDC_PiPlus;
	}

	double locSlope, locIntercept;

	for(unsigned int loc_i = 0; loc_i < locNumHitsVector.size(); loc_i++){
		if(locNumHitsUsedFordEdx == locNumHitsVector[loc_i]){
			locSigmadEdx = (Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i]))/1000000.0;
			return NOERROR;
		}
		if(loc_i == (locNumHitsVector.size() - 1)){
			locSigmadEdx_LowSide = Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i - 1]);
			locSigmadEdx_HighSide = Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i]);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i] - locNumHitsVector[loc_i - 1]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return NOERROR;
		}
		if((locNumHitsUsedFordEdx > locNumHitsVector[loc_i]) && (locNumHitsUsedFordEdx < locNumHitsVector[loc_i + 1])){
			locSigmadEdx_LowSide = Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i]);
			locSigmadEdx_HighSide = Function_dEdx(locBetaGamma, locdEdxSigmaParamVector[loc_i + 1]);

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

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locChargedTrackHypothesis->GetT(locTrackTimeBasedVector);
	const DTrackTimeBased *locTrackTimeBased = locTrackTimeBasedVector[0];
	unsigned int locNumHitsUsedFordEdx_CDC = locTrackTimeBased->dNumHitsUsedFordEdx_CDC;
	unsigned int locNumHitsUsedFordEdx_FDC = locTrackTimeBased->dNumHitsUsedFordEdx_FDC;

	bool locUseCDCHitsFlag = (locNumHitsUsedFordEdx_CDC >= locMinimumNumberUsedHitsForConfidence) ? true : false;
	bool locUseFDCHitsFlag = (locNumHitsUsedFordEdx_FDC >= locMinimumNumberUsedHitsForConfidence) ? true : false;

	if((locUseCDCHitsFlag == false) && (locUseFDCHitsFlag == false))
		return RESOURCE_UNAVAILABLE; //not enough hits, use other sources of information for PID

	double locDCdEdx = locTrackTimeBased->dEdx();
	double locMeandEdx_FDC, locMeandEdx_CDC, locSigmadEdx_FDC, locSigmadEdx_CDC;
	double locDeltadEdx_CDC = 0.0, locDeltadEdx_FDC = 0.0;

	double locBeta = locChargedTrackHypothesis->momentum().Mag()/locChargedTrackHypothesis->energy();
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


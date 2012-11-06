// $Id$
//
//    File: DNeutralParticleHypothesis_factory.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <TMath.h>

#include "DNeutralParticleHypothesis_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DNeutralParticleHypothesis_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticleHypothesis_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	// Get Target parameters from XML
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(runnumber):NULL;
	dTargetLength = 30.0;
	double locTargetCenterZ = 65.0;
	dTargetRadius = 1.5; //FIX: grab from database!!!
	if(locGeometry)
	{
		locGeometry->GetTargetZ(locTargetCenterZ);
		locGeometry->GetTargetLength(dTargetLength);
	}
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);


  // Get the particle ID algorithms
	vector<const DParticleID *> locPIDAlgorithms;
	locEventLoop->Get(locPIDAlgorithms);
	if(locPIDAlgorithms.size() < 1){
		_DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	// Drop the const qualifier from the DParticleID pointer (I'm surely going to hell for this!)
	dPIDAlgorithm = const_cast<DParticleID*>(locPIDAlgorithms[0]);
  
	// Warn user if something happened that caused us NOT to get a dPIDAlgorithm object pointer
	if(!dPIDAlgorithm){
		_DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticleHypothesis_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	unsigned int loc_i, loc_j, loc_k, locNDF = 1;
	bool locShowerMatchFlag;
	float locMass, locMomentum, locShowerEnergy, locParticleEnergy, locPathLength, locHitTime, locFlightTime, locProjectedTime, locTimeDifference;
	float locParticleEnergyUncertainty, locShowerEnergyUncertainty, locTimeDifferenceVariance, locChiSq, locFOM;
	DVector3 locPathVector, locHitPoint;

	const DNeutralShower *locNeutralShower;
	const DChargedTrackHypothesis *locChargedTrackHypothesis;
	DNeutralParticleHypothesis *locNeutralParticleHypothesis;
	DMatrixDSym locVariances, locErrorMatrix;
	vector<const DBCALShower*> locAssociatedBCALShowers_NeutralShower;
	vector<const DFCALShower*> locAssociatedFCALShowers_NeutralShower;
	vector<const DBCALShower*> locAssociatedBCALShowers_ChargedTrack;
	vector<const DFCALShower*> locAssociatedFCALShowers_ChargedTrack;

	vector<const DChargedTrack*> locChargedTracks;
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locChargedTracks);
	locEventLoop->Get(locNeutralShowers);

	vector<Particle_t> locPIDHypotheses;
	locPIDHypotheses.push_back(Gamma);
	locPIDHypotheses.push_back(Neutron);

	double locStartTime = 0.0; //fix when rf info available!!
	DLorentzVector locSpacetimeVertex(dTargetCenter, locStartTime);
	double locVertexTimeUncertainty = 0.0;

	// Loop over DNeutralShowers
	for (loc_i = 0; loc_i < locNeutralShowers.size(); loc_i++){
		locNeutralShower = locNeutralShowers[loc_i];
		locNeutralShower->GetT(locAssociatedBCALShowers_NeutralShower);
		locNeutralShower->GetT(locAssociatedFCALShowers_NeutralShower);

		// If the DNeutralShower is matched to the DChargedTrackHypothesis with the highest FOM of ANY DChargedTrack, skip it
		locShowerMatchFlag = false;
		for (loc_j = 0; loc_j < locChargedTracks.size(); loc_j++){
			locChargedTrackHypothesis = locChargedTracks[loc_j]->dChargedTrackHypotheses[0];
			locChargedTrackHypothesis->GetT(locAssociatedBCALShowers_ChargedTrack);
			locChargedTrackHypothesis->GetT(locAssociatedFCALShowers_ChargedTrack);
			if ((locAssociatedBCALShowers_ChargedTrack.size() > 0) && (locAssociatedBCALShowers_NeutralShower.size() > 0)){
				if (locAssociatedBCALShowers_ChargedTrack[0]->id == locAssociatedBCALShowers_NeutralShower[0]->id){
					locShowerMatchFlag = true;
					break;
				}
			}
			if ((locShowerMatchFlag == false) && (locAssociatedFCALShowers_ChargedTrack.size() > 0) && (locAssociatedFCALShowers_NeutralShower.size() > 0)){
				if (locAssociatedFCALShowers_ChargedTrack[0]->id == locAssociatedFCALShowers_NeutralShower[0]->id){
					locShowerMatchFlag = true;
					break;
				}
			}
		}
		if (locShowerMatchFlag == true)
			continue; //shower matched to a DChargedTrackHypothesis with the highest FOM, not a neutral

		locHitTime = locNeutralShower->dSpacetimeVertex.T();
		locShowerEnergy = locNeutralShower->dEnergy;
		locShowerEnergyUncertainty = sqrt(locNeutralShower->dCovarianceMatrix(0, 0));
		locHitPoint = locNeutralShower->dSpacetimeVertex.Vect();

		// Loop over vertices and PID hypotheses & create DNeutralParticleHypotheses for each combination
		for (loc_k = 0; loc_k < locPIDHypotheses.size(); loc_k++)
		{
			// Calculate DNeutralParticleHypothesis Quantities (projected time at vertex for given id, etc.)
			locMass = ParticleMass(locPIDHypotheses[loc_k]);
			locParticleEnergy = locShowerEnergy; //need to correct this for neutrons!
			if (locParticleEnergy < locMass)
				continue; //not enough energy for PID hypothesis

			locParticleEnergyUncertainty = locShowerEnergyUncertainty; //need to correct this for neutrons!

			locPathVector = locHitPoint - locSpacetimeVertex.Vect();
			locPathLength = locPathVector.Mag();
			if(!(locPathLength > 0.0))
				continue; //invalid, will divide by zero when creating error matrix, so skip!
			locMomentum = sqrt(locParticleEnergy*locParticleEnergy - locMass*locMass);
			locFlightTime = locPathLength*locParticleEnergy/(locMomentum*SPEED_OF_LIGHT);
			locProjectedTime = locHitTime - locFlightTime;

			// Calculate DNeutralParticleHypothesis FOM
			locTimeDifference = locSpacetimeVertex.T() - locProjectedTime;
			locTimeDifferenceVariance = 1.0; //completely random, ok because ID disabled for neutrons anyway
			locChiSq = locTimeDifference*locTimeDifference/locTimeDifferenceVariance;
			locFOM = TMath::Prob(locChiSq, locNDF);
			if(locPIDHypotheses[loc_k] == Neutron)
				locFOM = -1.0; //disables neutron ID until the neutron energy is calculated correctly from the deposited energy in the shower

			// Build DNeutralParticleHypothesis // dEdx not set
			locNeutralParticleHypothesis = new DNeutralParticleHypothesis;
			locNeutralParticleHypothesis->AddAssociatedObject(locNeutralShower);

			locNeutralParticleHypothesis->setMass(locMass);
			locNeutralParticleHypothesis->setCharge(0.0);

			Calc_Variances(locNeutralShower, locTimeDifferenceVariance, locVariances);
			Build_ErrorMatrix(locPathVector, locParticleEnergy, locVariances, locErrorMatrix);

			locNeutralParticleHypothesis->setErrorMatrix(locErrorMatrix);
			locNeutralParticleHypothesis->clearTrackingErrorMatrix();
			locPathVector.SetMag(locMomentum);
			locNeutralParticleHypothesis->setMomentum(locPathVector);
			locNeutralParticleHypothesis->setPosition(locSpacetimeVertex.Vect());
			locNeutralParticleHypothesis->setT0(locSpacetimeVertex.T(), locVertexTimeUncertainty, SYS_NULL);
			locNeutralParticleHypothesis->setTime(locProjectedTime);
			locNeutralParticleHypothesis->setT1(locNeutralShower->dSpacetimeVertex.T(), sqrt(locNeutralShower->dCovarianceMatrix(4, 4)), locNeutralShower->dDetectorSystem);
			locNeutralParticleHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)

			locNeutralParticleHypothesis->setPID(locPIDHypotheses[loc_k]);
			locNeutralParticleHypothesis->dChiSq = locChiSq;
			locNeutralParticleHypothesis->dNDF = locNDF;
			locNeutralParticleHypothesis->dFOM = locFOM;

			_data.push_back(locNeutralParticleHypothesis);	
		} //end PID loop
	} //end DNeutralShower loop

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralParticleHypothesis_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralParticleHypothesis_factory::fini(void)
{
	return NOERROR;
}

#define DELTA(i,j) ((i==j) ? 1 : 0)
void DNeutralParticleHypothesis_factory::Calc_Variances(const DNeutralShower *locNeutralShower, double locTimeDifferenceVariance, DMatrixDSym &locVariances)
{
	// create the simplest error matrix:
	// At this point, it is assumed that error matrix of measured quantities is diagonal,
	// with elements like: sigma_Z_t = L/sqrt(12) sigma_X_t = sigma_Y_t = r0/2 
	// L=target length, r0 = target radius...
	// This means that energy-depth-polar angle relation  is neglected.
	// the order of sigmas is:  x_c, y_c, z_c, E, x_t, y_t, z_t

	locVariances.Clear();
	locVariances.ResizeTo(7, 7);

	locVariances(0, 0) = locNeutralShower->dCovarianceMatrix(1, 1);
	locVariances(1, 1) = locNeutralShower->dCovarianceMatrix(2, 2);
	locVariances(2, 2) = locNeutralShower->dCovarianceMatrix(3, 3);

	locVariances(3, 3) = pow(0.5*dTargetRadius, 2.0) ; // x_t, y_t
	locVariances(4, 4) = pow(0.5*dTargetRadius, 2.0) ; // x_t, y_t
	locVariances(5, 5) = pow(dTargetLength/sqrt(12.0), 2.0) ; // z_t
	locVariances(6, 6) = locTimeDifferenceVariance;
}

void DNeutralParticleHypothesis_factory::Build_ErrorMatrix(const DVector3 &locPathVector, double locEnergy, const DMatrixDSym& locVariances, DMatrixDSym& locErrorMatrix)
{
	unsigned int loc_i, loc_j, loc_ik, loc_jk;
   double R = locPathVector.Mag(); //path length
   double R2= locPathVector.Mag2(); //path length ^ 2
   double R3 = R*R2; //path length ^ 3
   double E_R3 = locEnergy/R3; 

	// init and fill rotation matrix, first with momentum derivatives
	DMatrix A(7, 7);
	for (loc_i = 0; loc_i < 3; loc_i++) {
		for (loc_j = 0; loc_j <3; loc_j++) {
			A[loc_i][loc_j] = E_R3*( R2*DELTA(loc_i, loc_j) - locPathVector(loc_i) * locPathVector(loc_j) );
			A[loc_i][loc_j + 4] = - A[loc_i][loc_j];
		}
	}

	// fill energy part and remember: relation between energy and photon position in calorimeter is neglected!
	A[3][3] = 1.;
	for (loc_j = 0; loc_j < 3; loc_j++)
		A[loc_j][3] = locPathVector(loc_j)/R;

	// fill spatial part where: dp_r_x/dp_x_c = - dp_r_x/dp_x_v ....
	for (loc_i = 0; loc_i < 3; loc_i++) {
		for ( loc_j = 0; loc_j < 3; loc_j++) {
			loc_ik = loc_i + 4;
			loc_jk = loc_j + 4;
			A[loc_ik][loc_j] = DELTA(loc_i, loc_j);
			A[loc_ik][loc_jk] = - A[loc_ik][loc_j];
		}
	}

	locErrorMatrix.ResizeTo(7, 7);
   locErrorMatrix = locVariances;
   locErrorMatrix = locErrorMatrix.Similarity(A);
}


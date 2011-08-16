// $Id$
//
//    File: DNeutralTrackHypothesis_factory.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DNeutralTrackHypothesis_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DNeutralTrackHypothesis_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralTrackHypothesis_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{

	// Get Target parameters from XML
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(runnumber):NULL;
	dTargetLength = 30.0;
	dTargetRadius = 1.5; //FIX: grab from database!!!
	if(locGeometry)
		locGeometry->GetTargetLength(dTargetLength);

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
jerror_t DNeutralTrackHypothesis_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	unsigned int loc_i, loc_j, loc_k, locNDF = 1;
	bool locShowerMatchFlag;
	float locMass, locMomentum, locShowerEnergy, locParticleEnergy, locPathLength, locFlightTime, locProjectedTime, locTimeDifference;
	float locParticleEnergyUncertainty, locShowerEnergyUncertainty, locTimeDifferenceVariance, locChiSq, locFOM;
	DVector3 locPathVector;

	const DNeutralShowerCandidate *locNeutralShowerCandidate;
	const DChargedTrackHypothesis *locChargedTrackHypothesis;
	const DVertex *locVertex;
	DNeutralTrackHypothesis *locNeutralTrackHypothesis;
	DKinematicData *locKinematicData;
	DMatrixDSym locVariances, locErrorMatrix;
	vector<const DBCALShower*> locAssociatedBCALShowers_NeutralShowerCandidate;
	vector<const DFCALShower*> locAssociatedFCALShowers_NeutralShowerCandidate;
	vector<const DBCALShower*> locAssociatedBCALShowers_ChargedTrack;
	vector<const DFCALShower*> locAssociatedFCALShowers_ChargedTrack;

	vector<const DChargedTrack*> locChargedTracks;
	vector<const DNeutralShowerCandidate*> locNeutralShowerCandidates;
	vector<const DVertex*> locVertices;
	locEventLoop->Get(locChargedTracks);
	locEventLoop->Get(locNeutralShowerCandidates);
	locEventLoop->Get(locVertices);

	vector<Particle_t> locPIDHypotheses;
	locPIDHypotheses.push_back(Gamma);
	locPIDHypotheses.push_back(Neutron);

	// Loop over DNeutralShowerCandidates
	for (loc_i = 0; loc_i < locNeutralShowerCandidates.size(); loc_i++){
		locNeutralShowerCandidate = locNeutralShowerCandidates[loc_i];
		locNeutralShowerCandidate->GetT(locAssociatedBCALShowers_NeutralShowerCandidate);
		locNeutralShowerCandidate->GetT(locAssociatedFCALShowers_NeutralShowerCandidate);

		// If the DNeutralShowerCandidate is matched to the DChargedTrackHypothesis with the highest FOM of ANY DChargedTrack, skip it
		locShowerMatchFlag = false;
		for (loc_j = 0; loc_j < locChargedTracks.size(); loc_j++){
			locChargedTrackHypothesis = locChargedTracks[loc_j]->dChargedTrackHypotheses[0];
			locChargedTrackHypothesis->GetT(locAssociatedBCALShowers_ChargedTrack);
			locChargedTrackHypothesis->GetT(locAssociatedFCALShowers_ChargedTrack);
			if ((locAssociatedBCALShowers_ChargedTrack.size() > 0) && (locAssociatedBCALShowers_NeutralShowerCandidate.size() > 0)){
				if (locAssociatedBCALShowers_ChargedTrack[0]->id == locAssociatedBCALShowers_NeutralShowerCandidate[0]->id){
					locShowerMatchFlag = true;
					break;
				}
			}
			if ((locShowerMatchFlag == false) && (locAssociatedFCALShowers_ChargedTrack.size() > 0) && (locAssociatedFCALShowers_NeutralShowerCandidate.size() > 0)){
				if (locAssociatedFCALShowers_ChargedTrack[0]->id == locAssociatedFCALShowers_NeutralShowerCandidate[0]->id){
					locShowerMatchFlag = true;
					break;
				}
			}
		}
		if (locShowerMatchFlag == true)
			continue; //shower matched to a DChargedTrackHypothesis with the highest FOM, not a neutral

		// Loop over vertices and PID hypotheses & create DNeutralTrackHypotheses for each combination
		for (loc_j = 0; loc_j < locVertices.size(); loc_j++){
			locVertex = locVertices[loc_j];
			for (loc_k = 0; loc_k < locPIDHypotheses.size(); loc_k++){

				// Calculate DNeutralTrackHypothesis Quantities (projected time at vertex for given id, etc.)
				locMass = ParticleMass(locPIDHypotheses[loc_k]);
				locShowerEnergy = locNeutralShowerCandidate->dEnergy;
				locParticleEnergy = locShowerEnergy; //need to correct this for neutrons!
				if (locParticleEnergy < locMass)
					continue; //not enough energy for PID hypothesis

				locShowerEnergyUncertainty = locNeutralShowerCandidate->dEnergyUncertainty;
				locParticleEnergyUncertainty = locShowerEnergyUncertainty; //need to correct this for neutrons!

				locPathVector = locNeutralShowerCandidate->dSpacetimeVertex.Vect() - locVertex->dSpacetimeVertex.Vect();
				locPathLength = locPathVector.Mag();
				if(!(locPathLength > 0.0))
					continue; //invalid, will divide by zero when creating error matrix, so skip!
				locMomentum = sqrt(locParticleEnergy*locParticleEnergy - locMass*locMass);
				locFlightTime = locPathLength*locParticleEnergy/(locMomentum*SPEED_OF_LIGHT);
				locProjectedTime = locNeutralShowerCandidate->dSpacetimeVertex.T() - locFlightTime;

				// Calculate DNeutralTrackHypothesis FOM
				locTimeDifference = locVertex->dSpacetimeVertex.T() - locProjectedTime;
				locTimeDifferenceVariance = 1.0; //completely random, ok because ID disabled for neutrons anyway
				locChiSq = locTimeDifference*locTimeDifference/locTimeDifferenceVariance;
				locFOM = TMath::Prob(locChiSq, locNDF);
				if(locPIDHypotheses[loc_k] == Neutron)
					locFOM = -1.0; //disables neutron ID until the neutron energy is calculated correctly from the deposited energy in the shower

				// Build DKinematicData //dEdx not set
				locKinematicData = new DKinematicData;
				locKinematicData->setMass(locMass);
				locKinematicData->setCharge(0.0);

				Calc_Variances(locNeutralShowerCandidate, locParticleEnergyUncertainty, locVariances);
				Build_ErrorMatrix(locPathVector, locParticleEnergy, locVariances, locErrorMatrix);

				locKinematicData->setErrorMatrix(locErrorMatrix);
				locKinematicData->clearTrackingErrorMatrix();
				locPathVector.SetMag(locMomentum);
				locKinematicData->setMomentum(locPathVector);
				locKinematicData->setPosition(locVertex->dSpacetimeVertex.Vect());
				locKinematicData->setT0(locVertex->dSpacetimeVertex.T(), locVertex->dTimeUncertainty, SYS_NULL);
				locKinematicData->setT1(locNeutralShowerCandidate->dSpacetimeVertex.T(), locNeutralShowerCandidate->dSpacetimeVertexUncertainties.T(), locNeutralShowerCandidate->dDetectorSystem);
				locKinematicData->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)

				// Build DNeutralTrackHypothesis
				locNeutralTrackHypothesis = new DNeutralTrackHypothesis;
				locNeutralTrackHypothesis->AddAssociatedObject(locVertex);
				locNeutralTrackHypothesis->AddAssociatedObject(locNeutralShowerCandidate);
				locNeutralTrackHypothesis->dKinematicData = locKinematicData;
				locNeutralTrackHypothesis->dPID = locPIDHypotheses[loc_k];
				locNeutralTrackHypothesis->dPathLength = locPathLength;
				locNeutralTrackHypothesis->dFlightTime = locFlightTime;
				locNeutralTrackHypothesis->dProjectedTime = locProjectedTime;
				locNeutralTrackHypothesis->dChiSq = locChiSq;
				locNeutralTrackHypothesis->dNDF = locNDF;
				locNeutralTrackHypothesis->dFOM = locFOM;

				_data.push_back(locNeutralTrackHypothesis);	
			} //end PID loop
		} //end DVertex loop
	} //end DNeutralShowerCandidate loop

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralTrackHypothesis_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralTrackHypothesis_factory::fini(void)
{
	return NOERROR;
}

#define DELTA(i,j) ((i==j) ? 1 : 0)
void DNeutralTrackHypothesis_factory::Calc_Variances(const DNeutralShowerCandidate *locNeutralShowerCandidate, double locParticleEnergyUncertainty, DMatrixDSym &locVariances){
	DLorentzVector locShowerPositionUncertainties = locNeutralShowerCandidate->dSpacetimeVertexUncertainties;

	// create the simplest error matrix:
	// At this point, it is assumed that error matrix of measured quantities is diagonal,
	// with elements like: sigma_Z_t = L/sqrt(12) sigma_X_t = sigma_Y_t = r0/2 
	// L=target length, r0 = target radius...
	// This means that energy-depth-polar angle relation  is neglected.
	// the order of sigmas is:  x_c, y_c, z_c, E, x_t, y_t, z_t

	locVariances.Clear();
	locVariances.ResizeTo(7, 7);

	locVariances(0,0) = pow(locShowerPositionUncertainties.X(), 2.0);
	locVariances(1,1) = pow(locShowerPositionUncertainties.Y(), 2.0);
	locVariances(2,2) = pow(locShowerPositionUncertainties.Z(), 2.0);

	locVariances[3][3] = pow(locParticleEnergyUncertainty, 2.0);

	locVariances[4][4] = pow(0.5*dTargetRadius, 2.0) ; // x_t, y_t
	locVariances[5][5] = pow(0.5*dTargetRadius, 2.0) ; // x_t, y_t
	locVariances[6][6] = pow(dTargetLength/sqrt(12.0), 2.0) ; // z_t
}

void DNeutralTrackHypothesis_factory::Build_ErrorMatrix(const DVector3 &locPathVector, double locEnergy, const DMatrixDSym& locVariances, DMatrixDSym& locErrorMatrix)
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


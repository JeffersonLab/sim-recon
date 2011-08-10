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
	float locTimeDifferenceVariance, locChiSq, locFOM;
	DVector3 locPathVector;

	const DNeutralShowerCandidate *locNeutralShowerCandidate;
	const DChargedTrackHypothesis *locChargedTrackHypothesis;
	const DVertex *locVertex;
	DNeutralTrackHypothesis *locNeutralTrackHypothesis;
	DKinematicData *locKinematicData;
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
				locPathVector = locNeutralShowerCandidate->dSpacetimeVertex.Vect() - locVertex->dSpacetimeVertex.Vect();
				locPathLength = locPathVector.Mag();
				locMomentum = sqrt(locParticleEnergy*locParticleEnergy - locMass*locMass);
				locFlightTime = locPathLength*locParticleEnergy/(locMomentum*SPEED_OF_LIGHT);
				locProjectedTime = locNeutralShowerCandidate->dSpacetimeVertex.T() - locFlightTime;

				// Calculate DNeutralTrackHypothesis FOM
				locTimeDifference = locVertex->dSpacetimeVertex.T() - locProjectedTime;
				locTimeDifferenceVariance = 1.0; //completely random, ok because ID disabled for neutrons anyway
				locChiSq = locTimeDifference*locTimeDifference/locTimeDifferenceVariance;
				locFOM = TMath::Prob(locChiSq, locNDF);
				if(locPIDHypotheses[loc_k] == Neutron)
					locFOM = 0.0; //disables neutron ID until the neutron energy is calculated correctly from the deposited energy in the shower

				// Build DKinematicData
				locKinematicData = new DKinematicData;
				locKinematicData->setMass(locMass);
				locKinematicData->setCharge(0.0);
				locKinematicData->clearErrorMatrix(); // FIXME!!!
				locKinematicData->clearTrackingErrorMatrix();
				locPathVector.SetMag(locMomentum);
				locKinematicData->setMomentum(locPathVector);
				locKinematicData->setPosition(locVertex->dSpacetimeVertex.Vect());
				locKinematicData->setT1(locNeutralShowerCandidate->dSpacetimeVertex.T(), locNeutralShowerCandidate->dSpacetimeVertexUncertainties.T(), locNeutralShowerCandidate->dDetectorSystem);
				locKinematicData->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
				//error matrix & dEdx not set

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



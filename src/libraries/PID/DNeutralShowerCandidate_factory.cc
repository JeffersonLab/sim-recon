// $Id$
//
//    File: DNeutralShowerCandidate_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DNeutralShowerCandidate_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DNeutralShowerCandidate_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralShowerCandidate_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralShowerCandidate_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	unsigned int loc_i, loc_j, loc_k, loc_l;
	const DChargedTrackHypothesis* locChargedTrackHypothesis;
	const DChargedTrack *locChargedTrack;
	const DBCALShower *locBCALShower;
	const DFCALShower *locFCALShower;
	DNeutralShowerCandidate *locNeutralShowerCandidate;
	vector<const DBCALShower*> locAssociatedBCALShowers;
	vector<const DFCALShower*> locAssociatedFCALShowers;
	vector<const DChargedTrackHypothesis*> locAssociatedChargedTrackHypotheses;
	JObject::oid_t locObjectID;
	bool locChargedTrackHypothesisMatchFlag, locChargedTrackMatchFlag, locAnyChargedTrackMatchFlag;

	vector<const DChargedTrack*> locChargedTracks;
	vector<const DBCALShower*> locBCALShowers;
	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locChargedTracks);
	locEventLoop->Get(locBCALShowers);
	locEventLoop->Get(locFCALShowers);

	// Loop over all DBCALShowers, see if they match to DChargedTrackHypotheses.  
	// Create DNeutralShowerCandidate for each shower that is not matched to each DChargedTrackHypothesis of any DChargedTrack
	// Each DChargedTrackHypothesis that a shower is matched to is added as an associated object of the DNeutralShowerCandidate
	// If one of these DChargedTrackHypothesis objects is chosen for the actual track, then the DNeutralShowerCandidate should not be considered a neutral track
		// This choice is outside of the scope of this factory (in order to return a PID-independent list of showers), but is implemented during creation of the DNeutralTrackHypothesis objects
	for(loc_i = 0; loc_i < locBCALShowers.size(); loc_i++){
		locBCALShower = locBCALShowers[loc_i];
		locObjectID = locBCALShower->id;
		locAssociatedChargedTrackHypotheses.resize(0);
		locAnyChargedTrackMatchFlag = false;
		for(loc_j = 0; loc_j < locChargedTracks.size(); loc_j++){
			locChargedTrack = locChargedTracks[loc_j];
			locChargedTrackMatchFlag = false; //false in case no dChargedTrackHypotheses (not really possible...) 
			for(loc_k = 0; loc_k < locChargedTrack->dChargedTrackHypotheses.size(); loc_k++){
				locChargedTrackMatchFlag = true; //will set to false if no match for any of the DChargedTrackHypotheses
				locChargedTrackHypothesis = locChargedTrack->dChargedTrackHypotheses[loc_k];
				// Get the list of BCALShowers associated with this locChargedTrackHypothesis
				locChargedTrackHypothesis->GetT(locAssociatedBCALShowers);
				locChargedTrackHypothesisMatchFlag = false;
				for(loc_l = 0; loc_l < locAssociatedBCALShowers.size(); loc_l++){
					if(locAssociatedBCALShowers[loc_l]->id == locObjectID){ //check to see if any of the associated showers is the same as the one from the main loop
						locChargedTrackHypothesisMatchFlag = true;
						locAssociatedChargedTrackHypotheses.push_back(locChargedTrackHypothesis); //this charged track associated with this shower
						break;
					}
				}
				if (locChargedTrackHypothesisMatchFlag == false){
					locChargedTrackMatchFlag = false; //not matched to each DChargedTrackHypothesis for this DChargedTrack (ok for DNeutralShowerCandidate so far)
					break;
				}
			} //end of DChargedTrackHypothesis loop
			if (locChargedTrackMatchFlag == true){
				locAnyChargedTrackMatchFlag = true;
				break; // matched to every DChargedTrackHypothesis for this DChargedTrack: not a neutral, stop search
			}
		} //end of DChargedTrack loop
		if (locAnyChargedTrackMatchFlag == true)
			continue; // matched to every DChargedTrackHypothesis for a DChargedTrack: not a DNeutralShowerCandidate

		// create DNeutralShowerCandidate
		locNeutralShowerCandidate = new DNeutralShowerCandidate(locBCALShower);
		locNeutralShowerCandidate->AddAssociatedObject(locBCALShower);
		for (loc_j = 0; loc_j < locAssociatedChargedTrackHypotheses.size(); loc_j++) //if you later choose one of these DChargedTrackHypotheses for the DChargedTrack, this DNeutralShowerCandidate won't be a neutral
			locNeutralShowerCandidate->AddAssociatedObject(locAssociatedChargedTrackHypotheses[loc_j]);
		_data.push_back(locNeutralShowerCandidate);
	}

	// Loop over all DFCALShowers, see if they match to DChargedTrackHypotheses.  
	// Create DNeutralShowerCandidate for each shower that is not matched to each DChargedTrackHypothesis of any DChargedTrack
	// Each DChargedTrackHypothesis that a shower is matched to is added as an associated object of the DNeutralShowerCandidate
	// If one of these DChargedTrackHypothesis objects is chosen for the actual track, then the DNeutralShowerCandidate should not be considered a neutral track
		// This choice is outside of the scope of this factory (in order to return a PID-independent list of showers), but is implemented during creation of the DNeutralTrackHypothesis objects
	for(loc_i = 0; loc_i < locFCALShowers.size(); loc_i++){
		locFCALShower = locFCALShowers[loc_i];
		locObjectID = locFCALShower->id;
		locAssociatedChargedTrackHypotheses.resize(0);
		locAnyChargedTrackMatchFlag = false;
		for(loc_j = 0; loc_j < locChargedTracks.size(); loc_j++){
			locChargedTrack = locChargedTracks[loc_j];
			locChargedTrackMatchFlag = false; //false in case no dChargedTrackHypotheses (not really possible...) 
			for(loc_k = 0; loc_k < locChargedTrack->dChargedTrackHypotheses.size(); loc_k++){
				locChargedTrackMatchFlag = true; //will set to false if no match for any of the DChargedTrackHypotheses
				locChargedTrackHypothesis = locChargedTrack->dChargedTrackHypotheses[loc_k];
				// Get the list of FCALShowers associated with this locChargedTrackHypothesis
				locChargedTrackHypothesis->GetT(locAssociatedFCALShowers);
				locChargedTrackHypothesisMatchFlag = false;
				for(loc_l = 0; loc_l < locAssociatedFCALShowers.size(); loc_l++){
					if(locAssociatedFCALShowers[loc_l]->id == locObjectID){ //check to see if any of the associated showers is the same as the one from the main loop
						locChargedTrackHypothesisMatchFlag = true;
						locAssociatedChargedTrackHypotheses.push_back(locChargedTrackHypothesis); //this charged track associated with this shower
						break;
					}
				}
				if (locChargedTrackHypothesisMatchFlag == false){
					locChargedTrackMatchFlag = false; //not matched to each DChargedTrackHypothesis for this DChargedTrack (ok for DNeutralShowerCandidate so far)
					break;
				}
			} //end of DChargedTrackHypothesis loop
			if (locChargedTrackMatchFlag == true){
				locAnyChargedTrackMatchFlag = true;
				break; // matched to every DChargedTrackHypothesis for this DChargedTrack: not a neutral, stop search
			}
		} //end of DChargedTrack loop
		if (locAnyChargedTrackMatchFlag == true)
			continue; // matched to every DChargedTrackHypothesis for a DChargedTrack: not a DNeutralShowerCandidate

		// create DNeutralShowerCandidate
		locNeutralShowerCandidate = new DNeutralShowerCandidate(locFCALShower);
		locNeutralShowerCandidate->AddAssociatedObject(locFCALShower);
		for (loc_j = 0; loc_j < locAssociatedChargedTrackHypotheses.size(); loc_j++) //if you later choose one of these DChargedTrackHypotheses for the DChargedTrack, this DNeutralShowerCandidate won't be a neutral
			locNeutralShowerCandidate->AddAssociatedObject(locAssociatedChargedTrackHypotheses[loc_j]);
		_data.push_back(locNeutralShowerCandidate);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralShowerCandidate_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralShowerCandidate_factory::fini(void)
{
	return NOERROR;
}



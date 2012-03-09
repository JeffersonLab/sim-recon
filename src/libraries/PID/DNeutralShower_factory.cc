// $Id$
//
//    File: DNeutralShower_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DNeutralShower_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DNeutralShower_factory::init(void)
{
	//this parameter controls what BCAL reconstruction algorithm to use
 	USE_KLOE = 1;
	gPARMS->SetDefaultParameter("BCALRECON:USE_KLOE", USE_KLOE);


	if(USE_KLOE){    
		cout << "Using KLOE BCAL clustering." << endl;
	}
	else{   
		cout << "Using alternative (experimental) BCAL clustering." << endl;
	}

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralShower_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralShower_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	const DChargedTrackHypothesis* locChargedTrackHypothesis;
	const DBCALShower *locBCALShower;
	const DFCALShower *locFCALShower;
	DNeutralShower *locNeutralShower;
	vector<const DBCALShower*> locAssociatedBCALShowers;
	vector<const DFCALShower*> locAssociatedFCALShowers;
	vector<const DChargedTrackHypothesis*> locAssociatedChargedTrackHypotheses;
	bool locShowerMatchedToTrackFlag;

	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	vector<const DBCALShower*> locBCALShowers;
	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locChargedTrackHypotheses);
	if (USE_KLOE) {
		locEventLoop->Get(locBCALShowers, "KLOE");
	} else {	  
		locEventLoop->Get(locBCALShowers);
	}
	locEventLoop->Get(locFCALShowers);

	// Loop over all DBCALShowers, see if they match to DChargedTrackHypotheses.  
	// Create DNeutralShower for each shower that is not matched to any DChargedTrackHypothesis
	// The chance of an actual neutral shower matching to a bogus DChargedTrackHypothesis is very small
	for(unsigned int loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i){
		locBCALShower = locBCALShowers[loc_i];
		locShowerMatchedToTrackFlag = false;
		for(unsigned int loc_j = 0; loc_j < locChargedTrackHypotheses.size(); ++loc_j){
			locChargedTrackHypothesis = locChargedTrackHypotheses[loc_j];
			locChargedTrackHypothesis->GetT(locAssociatedBCALShowers);
			for(unsigned int loc_k = 0; loc_k < locAssociatedBCALShowers.size(); ++loc_k){
				if(locBCALShower == locAssociatedBCALShowers[loc_k]){
					locShowerMatchedToTrackFlag = true;
					break;
				}
			}
			if(locShowerMatchedToTrackFlag)
				break;
		}
		if(locShowerMatchedToTrackFlag)
			continue;

		// create DNeutralShower
		locNeutralShower = new DNeutralShower(locBCALShower);
		locNeutralShower->AddAssociatedObject(locBCALShower);
		_data.push_back(locNeutralShower);
	}

	// Loop over all DFCALShowers, see if they match to DChargedTrackHypotheses.  
	// Create DNeutralShower for each shower that is not matched to any DChargedTrackHypothesis
	// The chance of an actual neutral shower matching to a bogus DChargedTrackHypothesis is very small
	for(unsigned int loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i){
		locFCALShower = locFCALShowers[loc_i];
		locShowerMatchedToTrackFlag = false;
		for(unsigned int loc_j = 0; loc_j < locChargedTrackHypotheses.size(); ++loc_j){
			locChargedTrackHypothesis = locChargedTrackHypotheses[loc_j];
			locChargedTrackHypothesis->GetT(locAssociatedFCALShowers);
			for(unsigned int loc_k = 0; loc_k < locAssociatedFCALShowers.size(); ++loc_k){
				if(locFCALShower == locAssociatedFCALShowers[loc_k]){
					locShowerMatchedToTrackFlag = true;
					break;
				}
			}
			if(locShowerMatchedToTrackFlag)
				break;
		}
		if(locShowerMatchedToTrackFlag)
			continue;

		// create DNeutralShower
		locNeutralShower = new DNeutralShower(locFCALShower);
		locNeutralShower->AddAssociatedObject(locFCALShower);
		_data.push_back(locNeutralShower);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralShower_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralShower_factory::fini(void)
{
	return NOERROR;
}



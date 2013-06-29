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

inline bool DNeutralShower_SortByEnergy(const DNeutralShower* locNeutralShower1, const DNeutralShower* locNeutralShower2)
{
	// sort by increasing energy in the 1's and 0.1s digits (MeV): pseudo-random
	return int(locNeutralShower1->dEnergy*10000.0)%100 < int(locNeutralShower2->dEnergy*10000.0)%100;
}

//------------------
// init
//------------------
jerror_t DNeutralShower_factory::init(void)
{
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
	locEventLoop->Get(locBCALShowers);
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
		locNeutralShower = new DNeutralShower();
		locNeutralShower->dDetectorSystem = SYS_BCAL;
		locNeutralShower->dEnergy = locBCALShower->E;
		double locEnergyUncertainty = (locBCALShower->E >= 0.0) ? locBCALShower->E*sqrt( 0.0598*0.0598/locBCALShower->E + 0.0094*0.0094 ) : 1e-3; //last updated at svn revision 9242
		locNeutralShower->dSpacetimeVertex.SetXYZT(locBCALShower->x, locBCALShower->y, locBCALShower->z, locBCALShower->t);
		locNeutralShower->dCovarianceMatrix.ResizeTo(5, 5);
		locNeutralShower->dCovarianceMatrix(0, 0) = locEnergyUncertainty*locEnergyUncertainty;
		locNeutralShower->dCovarianceMatrix(1, 1) = locBCALShower->xErr*locBCALShower->xErr;
		locNeutralShower->dCovarianceMatrix(2, 2) = locBCALShower->yErr*locBCALShower->yErr;
		locNeutralShower->dCovarianceMatrix(3, 3) = locBCALShower->zErr*locBCALShower->zErr;
		locNeutralShower->dCovarianceMatrix(4, 4) = locBCALShower->tErr*locBCALShower->tErr;
		//NEED CORRELATIONS!
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
		locNeutralShower = new DNeutralShower();
		locNeutralShower->dDetectorSystem = SYS_FCAL;
		locNeutralShower->dEnergy = locFCALShower->getEnergy();
		double locEnergyUncertainty = (locFCALShower->getEnergy() >= 0.0) ? 0.042*sqrt(locFCALShower->getEnergy()) + 0.0001 : 1e-3; //from old DPhoton_factory::makeFCalPhoton() function
		locNeutralShower->dSpacetimeVertex.SetVect(locFCALShower->getPosition());
		locNeutralShower->dSpacetimeVertex.SetT(locFCALShower->getTime());
		locNeutralShower->dCovarianceMatrix.ResizeTo(5, 5);
		locNeutralShower->dCovarianceMatrix(0, 0) = locEnergyUncertainty*locEnergyUncertainty;
		locNeutralShower->dCovarianceMatrix(1, 1) = locFCALShower->getPositionError().X()*locFCALShower->getPositionError().X();
		locNeutralShower->dCovarianceMatrix(2, 2) = locFCALShower->getPositionError().Y()*locFCALShower->getPositionError().Y();
		locNeutralShower->dCovarianceMatrix(3, 3) = locFCALShower->getPositionError().Z()*locFCALShower->getPositionError().Z();
		locNeutralShower->dCovarianceMatrix(4, 4) = 0.0; //not stored in DFCALShower
		//NEED CORRELATIONS!
		locNeutralShower->AddAssociatedObject(locFCALShower);

		_data.push_back(locNeutralShower);
	}

	sort(_data.begin(), _data.end(), DNeutralShower_SortByEnergy);

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



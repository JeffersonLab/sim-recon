// $Id$
//
//    File: DNeutralShower_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DNeutralShower_factory.h"

inline bool DNeutralShower_SortByEnergy(const DNeutralShower* locNeutralShower1, const DNeutralShower* locNeutralShower2)
{
	// truncate the shower energies: in units of MeV, ignore all digits that are 10s-place and above
	// then sort by increasing energy: pseudo-random

	//guard against NaN: necessary since casting to int
	bool locFirstIsNaN = (!(locNeutralShower1->dEnergy > -1.0) && !(locNeutralShower1->dEnergy < 1.0));
	bool locSecondIsNaN = (!(locNeutralShower2->dEnergy > -1.0) && !(locNeutralShower2->dEnergy < 1.0));
	if(locFirstIsNaN)
		return false;
	if(locSecondIsNaN)
		return true;
	double locE1 = locNeutralShower1->dEnergy - double(int(locNeutralShower1->dEnergy*100.0))/100.0;
	double locE2 = locNeutralShower2->dEnergy - double(int(locNeutralShower2->dEnergy*100.0))/100.0;
	return (locE1 < locE2);
}

//------------------
// init
//------------------
jerror_t DNeutralShower_factory::init(void)
{
	dResourcePool_TMatrixFSym->Set_ControlParams(20, 20, 20);
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralShower_factory::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralShower_factory::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	// Loop over all DBCALShowers, create DNeutralShower if didn't match to any tracks
		// The chance of an actual neutral shower matching to a bogus track is very small
	JObject::oid_t locShowerID = 0;
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
	{
		if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[loc_i]))
			continue;

		// create DNeutralShower
		DNeutralShower* locNeutralShower = new DNeutralShower();
		locNeutralShower->dBCALFCALShower = static_cast<const JObject*>(locBCALShowers[loc_i]);
		locNeutralShower->dDetectorSystem = SYS_BCAL;
		locNeutralShower->dShowerID = locShowerID;
		++locShowerID;

		locNeutralShower->dEnergy = locBCALShowers[loc_i]->E;
		locNeutralShower->dSpacetimeVertex.SetXYZT(locBCALShowers[loc_i]->x, locBCALShowers[loc_i]->y, locBCALShowers[loc_i]->z, locBCALShowers[loc_i]->t);
		auto locCovMatrix = dResourcePool_TMatrixFSym->Get_SharedResource();
		locCovMatrix->ResizeTo(5, 5);
		*locCovMatrix = locBCALShowers[loc_i]->ExyztCovariance;
		locNeutralShower->dCovarianceMatrix = locCovMatrix;

		locNeutralShower->AddAssociatedObject(locBCALShowers[loc_i]);

		_data.push_back(locNeutralShower);
	}

	// Loop over all DFCALShowers, create DNeutralShower if didn't match to any tracks
		// The chance of an actual neutral shower matching to a bogus track is very small
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
	{
		if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[loc_i]))
			continue;

		// create DNeutralShower
		DNeutralShower* locNeutralShower = new DNeutralShower();
		locNeutralShower->dBCALFCALShower = static_cast<const JObject*>(locFCALShowers[loc_i]);
		locNeutralShower->dDetectorSystem = SYS_FCAL;
		locNeutralShower->dShowerID = locShowerID;
		++locShowerID;

		locNeutralShower->dEnergy = locFCALShowers[loc_i]->getEnergy();
		locNeutralShower->dSpacetimeVertex.SetVect(locFCALShowers[loc_i]->getPosition());
		locNeutralShower->dSpacetimeVertex.SetT(locFCALShowers[loc_i]->getTime());

		auto locCovMatrix = dResourcePool_TMatrixFSym->Get_SharedResource();
		locCovMatrix->ResizeTo(5, 5);
		*locCovMatrix = locFCALShowers[loc_i]->ExyztCovariance;
		locNeutralShower->dCovarianceMatrix = locCovMatrix;

		locNeutralShower->AddAssociatedObject(locFCALShowers[loc_i]);

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


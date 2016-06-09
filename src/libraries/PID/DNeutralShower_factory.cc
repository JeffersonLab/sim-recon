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
		locNeutralShower->dDetectorSystem = SYS_BCAL;
		locNeutralShower->dShowerID = locShowerID;
		++locShowerID;

		locNeutralShower->dEnergy = locBCALShowers[loc_i]->E;
		double locEnergyUncertainty = (locBCALShowers[loc_i]->E >= 0.0) ? locBCALShowers[loc_i]->E*sqrt( 0.0598*0.0598/locBCALShowers[loc_i]->E + 0.0094*0.0094 ) : 1e-3; //last updated at svn revision 9242
		locNeutralShower->dSpacetimeVertex.SetXYZT(locBCALShowers[loc_i]->x, locBCALShowers[loc_i]->y, locBCALShowers[loc_i]->z, locBCALShowers[loc_i]->t);

		locNeutralShower->dCovarianceMatrix.ResizeTo(5, 5);
		locNeutralShower->dCovarianceMatrix(0, 0) = locEnergyUncertainty*locEnergyUncertainty;
		locNeutralShower->dCovarianceMatrix(1, 1) = locBCALShowers[loc_i]->xErr*locBCALShowers[loc_i]->xErr;
		locNeutralShower->dCovarianceMatrix(2, 2) = locBCALShowers[loc_i]->yErr*locBCALShowers[loc_i]->yErr;
		locNeutralShower->dCovarianceMatrix(3, 3) = locBCALShowers[loc_i]->zErr*locBCALShowers[loc_i]->zErr;
		locNeutralShower->dCovarianceMatrix(4, 4) = locBCALShowers[loc_i]->tErr*locBCALShowers[loc_i]->tErr;
		//NEED CORRELATIONS!
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
		locNeutralShower->dDetectorSystem = SYS_FCAL;
		locNeutralShower->dShowerID = locShowerID;
		++locShowerID;

		locNeutralShower->dEnergy = locFCALShowers[loc_i]->getEnergy();
		double locEnergyUncertainty = (locFCALShowers[loc_i]->getEnergy() >= 0.0) ? 0.042*sqrt(locFCALShowers[loc_i]->getEnergy()) + 0.0001 : 1e-3; //from old DPhoton_factory::makeFCalPhoton() function
		locNeutralShower->dSpacetimeVertex.SetVect(locFCALShowers[loc_i]->getPosition());
		locNeutralShower->dSpacetimeVertex.SetT(locFCALShowers[loc_i]->getTime());

		locNeutralShower->dCovarianceMatrix.ResizeTo(5, 5);
		locNeutralShower->dCovarianceMatrix(0, 0) = locEnergyUncertainty*locEnergyUncertainty;
		locNeutralShower->dCovarianceMatrix(1, 1) = locFCALShowers[loc_i]->getPositionError().X()*locFCALShowers[loc_i]->getPositionError().X();
		locNeutralShower->dCovarianceMatrix(2, 2) = locFCALShowers[loc_i]->getPositionError().Y()*locFCALShowers[loc_i]->getPositionError().Y();
		locNeutralShower->dCovarianceMatrix(3, 3) = locFCALShowers[loc_i]->getPositionError().Z()*locFCALShowers[loc_i]->getPositionError().Z();
		locNeutralShower->dCovarianceMatrix(4, 4) = 0.0; //not stored in DFCALShower
		//NEED CORRELATIONS!
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


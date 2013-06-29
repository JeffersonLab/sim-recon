// $Id$
//
//    File: DChargedTrackHypothesis_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DChargedTrackHypothesis_factory.h"
using namespace jana;

inline bool DChargedTrackHypothesis_SortByEnergy(const DChargedTrackHypothesis* locChargedTrackHypothesis1, const DChargedTrackHypothesis* locChargedTrackHypothesis2)
{
	// sort by increasing energy in the 1's and 0.1s digits (MeV): pseudo-random
	return int(locChargedTrackHypothesis1->energy()*10000.0)%100 < int(locChargedTrackHypothesis2->energy()*10000.0)%100;
}

//------------------
// init
//------------------
jerror_t DChargedTrackHypothesis_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrackHypothesis_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
  // Get the particle ID algorithms
	vector<const DParticleID *> locPIDAlgorithms;
	locEventLoop->Get(locPIDAlgorithms);
	if(locPIDAlgorithms.size() < 1){
		_DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	// Drop the const qualifier from the DParticleID pointer (I'm surely going to hell for this!)
	dPIDAlgorithm = locPIDAlgorithms[0];
  
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
jerror_t DChargedTrackHypothesis_factory::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);

	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); loc_i++)
	{
		DChargedTrackHypothesis* locChargedTrackHypothesis = Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBasedVector[loc_i], locEventRFBunches[0], false);
		if(locChargedTrackHypothesis != NULL)
			_data.push_back(locChargedTrackHypothesis);
	}
	sort(_data.begin(), _data.end(), DChargedTrackHypothesis_SortByEnergy);

	return NOERROR;
}

DChargedTrackHypothesis* DChargedTrackHypothesis_factory::Create_ChargedTrackHypothesis(JEventLoop* locEventLoop, const DTrackTimeBased* locTrackTimeBased, const DEventRFBunch* locEventRFBunch, bool locRFTimeFixedFlag)
{
	DMatrixDSym locCovarianceMatrix(7);

	unsigned int locTOFIndex, locSCIndex;
	double locTempProjectedTime = 0.0, locPathLength = 0.0, locFlightTime = 0.0;

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	DChargedTrackHypothesis* locChargedTrackHypothesis = new DChargedTrackHypothesis();
	locChargedTrackHypothesis->AddAssociatedObject(locTrackTimeBased);
	locChargedTrackHypothesis->candidateid = locTrackTimeBased->candidateid;
	locChargedTrackHypothesis->dRT = locTrackTimeBased->rt;

	// Chi square and degree-of-freedom data from the track fit
	locChargedTrackHypothesis->dChiSq_Track=locTrackTimeBased->chisq;
	locChargedTrackHypothesis->dNDF_Track=locTrackTimeBased->Ndof;
	
	//Set DKinematicData Members
	DKinematicData *locKinematicData = locChargedTrackHypothesis;
	*locKinematicData = *(static_cast<const DKinematicData*>(locTrackTimeBased));
	locCovarianceMatrix = locTrackTimeBased->errorMatrix();

	//useful kinematic variable
	double betagamma=locTrackTimeBased->momentum().Mag()/locTrackTimeBased->mass();

	// Use time-based tracking time as initial guess
	double locInitialStartTime = locChargedTrackHypothesis->t0(); // used to reject hits that are not in time with the track
	double locInitialStartTimeUncertainty = locChargedTrackHypothesis->t0_err();
	locChargedTrackHypothesis->setTime(locChargedTrackHypothesis->t0());
	locCovarianceMatrix(6, 6) = locInitialStartTimeUncertainty*locInitialStartTimeUncertainty;
	locChargedTrackHypothesis->setT0(locInitialStartTime, locInitialStartTimeUncertainty, SYS_CDC);
	locChargedTrackHypothesis->setT1(locInitialStartTime, locInitialStartTimeUncertainty, SYS_CDC);
	locChargedTrackHypothesis->setPathLength(NaN, 0.0); //zero uncertainty (for now)

	// Match to the start counter using the result of the time-based fit
	locTempProjectedTime = locInitialStartTime; // to reject hits that are not in time with the track
	pair<double,double>locTempSCdEdx;
	locChargedTrackHypothesis->dStartCounterdEdx=0.; 
	if (dPIDAlgorithm->MatchToSC(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locSCHits, locTempProjectedTime, locSCIndex, locPathLength, locFlightTime,&locTempSCdEdx) == NOERROR)
	{
		locChargedTrackHypothesis->setT0(locTempProjectedTime, 0.3, SYS_START); //uncertainty guess for now
		locCovarianceMatrix(6, 6) = 0.3*0.3; // guess for now //will be overriden by other detector systems if hit match
		locChargedTrackHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now) //will be overriden by other detector systems if hit match
		locChargedTrackHypothesis->dStartCounterdEdx=locTempSCdEdx.first;
		locChargedTrackHypothesis->dStartCounterdEdx_norm_residual=(locTempSCdEdx.second-2.204+1.893/betagamma)/1.858;

		locChargedTrackHypothesis->AddAssociatedObject(locSCHits[locSCIndex]);
	}

	// Try matching the track with hits in the outer detectors
	locTempProjectedTime = locInitialStartTime; // to reject hits that are not in time with the track
	deque<const DBCALShower*> locMatchedBCALShowers;
	if (dPIDAlgorithm->MatchToBCAL(locTrackTimeBased->rt, locBCALShowers, locMatchedBCALShowers, locTempProjectedTime, locPathLength, locFlightTime) == NOERROR)
	{
		for(unsigned int loc_j = 0; loc_j < locMatchedBCALShowers.size(); ++loc_j)
			locChargedTrackHypothesis->AddAssociatedObject(locMatchedBCALShowers[loc_j]);
		locChargedTrackHypothesis->setTime(locTempProjectedTime);
		locCovarianceMatrix(6, 6) = 0.00255*pow(locChargedTrackHypothesis->momentum().Mag(), -2.52) + 0.220;
		locCovarianceMatrix(6, 6) = locCovarianceMatrix(6, 6)*locCovarianceMatrix(6, 6);
		locChargedTrackHypothesis->setT1(locMatchedBCALShowers[0]->t, locMatchedBCALShowers[0]->tErr, SYS_BCAL);
		locChargedTrackHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
	}
	locTempProjectedTime = locInitialStartTime; // to reject hits that are not in time with the track
	pair<double,double>locTOFdEdx;
	if (dPIDAlgorithm->MatchToTOF(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locTOFPoints, locTempProjectedTime, locTOFIndex, locPathLength, locFlightTime,&locTOFdEdx) == NOERROR)
	{
	  locChargedTrackHypothesis->AddAssociatedObject(locTOFPoints[locTOFIndex]);
	  locChargedTrackHypothesis->dTOFdEdx=locTOFdEdx.first;
	  locChargedTrackHypothesis->dTOFdEdx_norm_residual
	    =(locTOFdEdx.second-2.462+1.488/betagamma)/1.273;
		if(locChargedTrackHypothesis->t1_detector() == SYS_CDC)
		{
			locChargedTrackHypothesis->setTime(locTempProjectedTime);
			locCovarianceMatrix(6, 6) = 0.08*0.08;
			locChargedTrackHypothesis->setT1(locTOFPoints[locTOFIndex]->t, 0.08, SYS_TOF);
			locChargedTrackHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
		}
	}
	locTempProjectedTime = locInitialStartTime; // to reject hits that are not in time with the track
	double locFCALdEdx=0.;
	deque<const DFCALShower*> locMatchedFCALShowers;
	if (dPIDAlgorithm->MatchToFCAL(locTrackTimeBased->rt, locFCALShowers, locMatchedFCALShowers, locTempProjectedTime, locPathLength, locFlightTime,&locFCALdEdx) == NOERROR)
	{
		for(unsigned int loc_j = 0; loc_j < locMatchedFCALShowers.size(); ++loc_j)
			locChargedTrackHypothesis->AddAssociatedObject(locMatchedFCALShowers[loc_j]);
		locChargedTrackHypothesis->dFCALdEdx=locFCALdEdx;
		if(locChargedTrackHypothesis->t1_detector() == SYS_CDC)
		{
			locChargedTrackHypothesis->setTime(locTempProjectedTime);
			locCovarianceMatrix(6, 6) = 0.6*0.6; // straight-line fit to high momentum data
			locChargedTrackHypothesis->setT1(locMatchedFCALShowers[0]->getTime(), 0.6, SYS_FCAL);
			locChargedTrackHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
		}
	}
	locChargedTrackHypothesis->setErrorMatrix(locCovarianceMatrix);

	//Calculate PID ChiSq, NDF, FOM
	dPIDAlgorithm->Calc_ChargedPIDFOM(locChargedTrackHypothesis, locEventRFBunch, locRFTimeFixedFlag);

	return locChargedTrackHypothesis;
}

//------------------
// erun
//------------------
jerror_t DChargedTrackHypothesis_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DChargedTrackHypothesis_factory::fini(void)
{
	return NOERROR;
}



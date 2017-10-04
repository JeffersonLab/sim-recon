// $Id$
//
//    File: DChargedTrackHypothesis_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DChargedTrackHypothesis_factory.h"

inline bool DChargedTrackHypothesis_SortByEnergy(const DChargedTrackHypothesis* locChargedTrackHypothesis1, const DChargedTrackHypothesis* locChargedTrackHypothesis2)
{
	// truncate the track energies: in units of MeV, ignore all digits that are 10s-place and above
	// then sort by increasing energy: pseudo-random

	//guard against NaN: necessary since casting to int
	bool locFirstIsNaN = (!(locChargedTrackHypothesis1->energy() > -1.0) && !(locChargedTrackHypothesis1->energy() < 1.0));
	bool locSecondIsNaN = (!(locChargedTrackHypothesis2->energy() > -1.0) && !(locChargedTrackHypothesis2->energy() < 1.0));
	if(locFirstIsNaN)
		return false;
	if(locSecondIsNaN)
		return true;
	double locE1 = locChargedTrackHypothesis1->energy() - double(int(locChargedTrackHypothesis1->energy()*100.0))/100.0;
	double locE2 = locChargedTrackHypothesis2->energy() - double(int(locChargedTrackHypothesis2->energy()*100.0))/100.0;
	return (locE1 < locE2);
}

//------------------
// init
//------------------
jerror_t DChargedTrackHypothesis_factory::init(void)
{
	//Setting this flag makes it so that JANA does not delete the objects in _data.  This factory will manage this memory. 
	SetFactoryFlag(NOT_OBJECT_OWNER);
	dResourcePool_ChargedTrackHypothesis = new DResourcePool<DChargedTrackHypothesis>();
	dResourcePool_ChargedTrackHypothesis->Set_ControlParams(30, 20, 200, 2000, 0);
	dResourcePool_TMatrixFSym = std::make_shared<DResourcePool<TMatrixFSym>>();
	dResourcePool_TMatrixFSym->Set_ControlParams(20, 20, 50);
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrackHypothesis_factory::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	locEventLoop->GetSingle(dPIDAlgorithm);
 	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrackHypothesis_factory::evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber)
{
	//Recycle
	dResourcePool_ChargedTrackHypothesis->Recycle(dCreated);
	dCreated.clear();
	_data.clear();

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);
	if (locEventRFBunches.size() == 0)
	   return NOERROR;

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	map<JObject::oid_t, vector<DChargedTrackHypothesis*> > locChargedTrackHypotheses;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); loc_i++)
	{
		DChargedTrackHypothesis* locChargedTrackHypothesis = Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBasedVector[loc_i], locDetectorMatches, locEventRFBunches[0]);
		locChargedTrackHypotheses[locChargedTrackHypothesis->Get_TrackTimeBased()->candidateid].push_back(locChargedTrackHypothesis);
	}

	//choose the first hypothesis from each track, and sort by increasing energy in the 1's and 0.1s digits (MeV): pseudo-random
	vector<DChargedTrackHypothesis*> locTracksToSort;
	map<JObject::oid_t, vector<DChargedTrackHypothesis*> >::iterator locIterator = locChargedTrackHypotheses.begin();
	for(; locIterator != locChargedTrackHypotheses.end(); ++locIterator)
		locTracksToSort.push_back(locIterator->second[0]);
	sort(locTracksToSort.begin(), locTracksToSort.end(), DChargedTrackHypothesis_SortByEnergy);

	//now loop through the sorted vector, grab all of the hypotheses for each of those tracks from the map, and save them all in _data
	for(size_t loc_i = 0; loc_i < locTracksToSort.size(); loc_i++)
	{
		JObject::oid_t locTrackID = locTracksToSort[loc_i]->Get_TrackTimeBased()->candidateid;
		_data.insert(_data.end(), locChargedTrackHypotheses[locTrackID].begin(), locChargedTrackHypotheses[locTrackID].end());
	}

	dCreated = _data;
	return NOERROR;
}

DChargedTrackHypothesis* DChargedTrackHypothesis_factory::Create_ChargedTrackHypothesis(JEventLoop* locEventLoop, const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, const DEventRFBunch* locEventRFBunch)
{
	DChargedTrackHypothesis* locChargedTrackHypothesis = Get_Resource();
	locChargedTrackHypothesis->Share_FromInput_Kinematics(static_cast<const DKinematicData*>(locTrackTimeBased));
	locChargedTrackHypothesis->Set_TrackTimeBased(locTrackTimeBased);

	auto locCovarianceMatrix = dResourcePool_TMatrixFSym->Get_SharedResource();
	locCovarianceMatrix->ResizeTo(7, 7);
	if(locChargedTrackHypothesis->errorMatrix() != nullptr)
		*locCovarianceMatrix = *(locChargedTrackHypothesis->errorMatrix());

	// RF Time
	if(locEventRFBunch->dTimeSource != SYS_NULL)
	{
		double locPropagatedRFTime = dPIDAlgorithm->Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch);
		locChargedTrackHypothesis->Set_T0(locPropagatedRFTime, locEventRFBunch->dTimeVariance, locEventRFBunch->dTimeSource);
	}

	// Start Counter
	shared_ptr<const DSCHitMatchParams> locSCHitMatchParams;
	if(dPIDAlgorithm->Get_BestSCMatchParams(locTrackTimeBased, locDetectorMatches, locSCHitMatchParams))
	{
		locChargedTrackHypothesis->Set_SCHitMatchParams(locSCHitMatchParams);

		double locPropagatedTime = locSCHitMatchParams->dHitTime - locSCHitMatchParams->dFlightTime;
		locChargedTrackHypothesis->setTime(locPropagatedTime);

//		double locPropagatedTimeUncertainty = sqrt(locSCHitMatchParams->dHitTimeVariance + locSCHitMatchParams->dFlightTimeVariance);
		double locPropagatedTimeUncertainty = 0.3;
//		double locFlightTimePCorrelation = locDetectorMatches->Get_FlightTimePCorrelation(locTrackTimeBased, SYS_START); //uncomment when ready!!
//		Add_TimeToTrackingMatrix(locChargedTrackHypothesis, locCovarianceMatrix, locSCHitMatchParams->dFlightTimeVariance, locSCHitMatchParams->dHitTimeVariance, locFlightTimePCorrelation); //uncomment when ready!!
		(*locCovarianceMatrix)(6, 6) = 0.3*0.3+locSCHitMatchParams->dFlightTimeVariance;

		if(locEventRFBunch->dTimeSource == SYS_NULL)
			locChargedTrackHypothesis->Set_T0(locPropagatedTime, locPropagatedTimeUncertainty, SYS_START); //update when ready
	}

	// MATCHES
	shared_ptr<const DBCALShowerMatchParams> locBCALShowerMatchParams;
	shared_ptr<const DTOFHitMatchParams> locTOFHitMatchParams;
	shared_ptr<const DFCALShowerMatchParams> locFCALShowerMatchParams;
	if(dPIDAlgorithm->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
		locChargedTrackHypothesis->Set_BCALShowerMatchParams(locBCALShowerMatchParams);
	if(dPIDAlgorithm->Get_BestTOFMatchParams(locTrackTimeBased, locDetectorMatches, locTOFHitMatchParams))
		locChargedTrackHypothesis->Set_TOFHitMatchParams(locTOFHitMatchParams);
	if(dPIDAlgorithm->Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams))
		locChargedTrackHypothesis->Set_FCALShowerMatchParams(locFCALShowerMatchParams);

	//PID
	if(locChargedTrackHypothesis->t1_detector() == SYS_BCAL)
	{
		auto locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
		const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
		locChargedTrackHypothesis->setTime(locBCALShower->t - locBCALShowerMatchParams->dFlightTime);
//		double locFlightTimePCorrelation = locDetectorMatches->Get_FlightTimePCorrelation(locTrackTimeBased, SYS_BCAL); //uncomment when ready!!
//		Add_TimeToTrackingMatrix(locChargedTrackHypothesis, locCovarianceMatrix.get(), locBCALShowerMatchParams->dFlightTimeVariance, locBCALShower->tErr()*locBCALShower->tErr(), locFlightTimePCorrelation); //uncomment when ready!!
		(*locCovarianceMatrix)(6 , 6) = 0.3*0.3+locBCALShowerMatchParams->dFlightTimeVariance;
	}

	// TOF
	if(locChargedTrackHypothesis->t1_detector() == SYS_TOF)
	{
		auto locTOFHitMatchParams = locChargedTrackHypothesis->Get_TOFHitMatchParams();
		locChargedTrackHypothesis->setTime(locTOFHitMatchParams->dHitTime - locTOFHitMatchParams->dFlightTime);
//		double locFlightTimePCorrelation = locDetectorMatches->Get_FlightTimePCorrelation(locTrackTimeBased, SYS_TOF); //uncomment when ready!!
//		Add_TimeToTrackingMatrix(locChargedTrackHypothesis, locCovarianceMatrix.get(), locTOFHitMatchParams->dFlightTimeVariance, locTOFHitMatchParams->dHitTimeVariance, locFlightTimePCorrelation); //uncomment when ready!!
		(*locCovarianceMatrix)(6, 6) = 0.1*0.1+locTOFHitMatchParams->dFlightTimeVariance;
	}

	//FCAL
	if(locChargedTrackHypothesis->t1_detector() == SYS_FCAL)
	{
		auto locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
		const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
		locChargedTrackHypothesis->setTime(locFCALShower->getTime() - locFCALShowerMatchParams->dFlightTime);
//		double locFlightTimePCorrelation = locDetectorMatches->Get_FlightTimePCorrelation(locTrackTimeBased, SYS_FCAL); //uncomment when ready!!
//		Add_TimeToTrackingMatrix(locChargedTrackHypothesis, locCovarianceMatrix.get(), locFCALShowerMatchParams->dFlightTimeVariance, locFCALShower->dCovarianceMatrix(4, 4), locFlightTimePCorrelation); //uncomment when ready!!
		(*locCovarianceMatrix)(6, 6) = 0.5*0.5+locFCALShowerMatchParams->dFlightTimeVariance;
	}

	locChargedTrackHypothesis->setErrorMatrix(locCovarianceMatrix);

	//Calculate PID ChiSq, NDF, FOM
	locChargedTrackHypothesis->Set_TimeAtPOCAToVertex(locChargedTrackHypothesis->time());
	dPIDAlgorithm->Calc_ChargedPIDFOM(locChargedTrackHypothesis);

	return locChargedTrackHypothesis;
}

void DChargedTrackHypothesis_factory::Add_TimeToTrackingMatrix(DChargedTrackHypothesis* locChargedTrackHypothesis, TMatrixFSym* locCovarianceMatrix, double locFlightTimeVariance, double locHitTimeVariance, double locFlightTimePCorrelation) const
{
	DVector3 locMomentum = locChargedTrackHypothesis->momentum();

	//extract momentum 3x3 portion of covariance matrix
	TMatrixFSym locCov_PxPyPz(3);
	for(size_t loc_i = 0; loc_i < 3; ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < 3; ++loc_j)
			locCov_PxPyPz(loc_i, loc_j) = (*locCovarianceMatrix)(loc_i, loc_j);
	}

	//convert px/py/pz to p/theta/phi
	double locPMag = locMomentum.Mag();
	double locPPerpSq = locMomentum.Px()*locMomentum.Px() + locMomentum.Py()*locMomentum.Py();
	TMatrix locJ_CartSpher(3, 3); //jacobian
	locJ_CartSpher(0, 0) = locMomentum.Px()/locPMag;
	locJ_CartSpher(0, 1) = locMomentum.Py()/locPMag;
	locJ_CartSpher(0, 2) = locMomentum.Pz()/locPMag;
	locJ_CartSpher(1, 0) = locMomentum.Px()*locMomentum.Pz()/(locPMag*locPMag*sqrt(locPPerpSq));
	locJ_CartSpher(1, 1) = locMomentum.Py()*locMomentum.Pz()/(locPMag*locPMag*sqrt(locPPerpSq));
	locJ_CartSpher(1, 2) = sqrt(locPPerpSq)/(locPMag*locPMag);
	locJ_CartSpher(2, 0) = -1.0*locMomentum.Py()/locPPerpSq;
	locJ_CartSpher(2, 1) = locMomentum.Px()/locPPerpSq;
	locJ_CartSpher(2, 2) = 0.0;
	TMatrixFSym locCov_PThetaPhi = locCov_PxPyPz.Similarity(locJ_CartSpher);

	//add flight time and hit time to p/theta/phi covariance matrix
	TMatrixFSym locCov_PThetaPhiFTimeHTime(5);
	double locSigmaPMag = locCov_PThetaPhi(0, 0);
	for(size_t loc_i = 0; loc_i < 3; ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < 3; ++loc_j)
			locCov_PThetaPhiFTimeHTime(loc_i, loc_j) = locCov_PThetaPhi(loc_i, loc_j);
	}
	locCov_PThetaPhiFTimeHTime(0, 3) = locFlightTimePCorrelation*locSigmaPMag*sqrt(locFlightTimeVariance);
	locCov_PThetaPhiFTimeHTime(0, 4) = locFlightTimePCorrelation*locSigmaPMag*sqrt(locHitTimeVariance);
	locCov_PThetaPhiFTimeHTime(1, 3) = 0.0;
	locCov_PThetaPhiFTimeHTime(1, 4) = 0.0;
	locCov_PThetaPhiFTimeHTime(2, 3) = 0.0;
	locCov_PThetaPhiFTimeHTime(2, 4) = 0.0;

	locCov_PThetaPhiFTimeHTime(3, 0) = locFlightTimePCorrelation*locSigmaPMag*sqrt(locFlightTimeVariance);
	locCov_PThetaPhiFTimeHTime(3, 1) = 0.0;
	locCov_PThetaPhiFTimeHTime(3, 2) = 0.0;
	locCov_PThetaPhiFTimeHTime(3, 3) = locFlightTimeVariance;
	locCov_PThetaPhiFTimeHTime(3, 4) = sqrt(locFlightTimeVariance*locHitTimeVariance); //correlation = 1

	locCov_PThetaPhiFTimeHTime(4, 0) = locFlightTimePCorrelation*locSigmaPMag*sqrt(locHitTimeVariance);
	locCov_PThetaPhiFTimeHTime(4, 1) = 0.0;
	locCov_PThetaPhiFTimeHTime(4, 2) = 0.0;
	locCov_PThetaPhiFTimeHTime(4, 3) = sqrt(locFlightTimeVariance*locHitTimeVariance); //correlation = 1
	locCov_PThetaPhiFTimeHTime(4, 4) = locHitTimeVariance;

	//convert p/theta/phi/flight-time/hit-time to px/py/pz/t
	TMatrix locJ_SpherCart(4, 5); //jacobian
	double locTheta = locMomentum.Theta();
	locJ_SpherCart(0, 0) = locMomentum.Px()/locPMag;
	locJ_SpherCart(0, 1) = locMomentum.Px()/tan(locTheta);
	locJ_SpherCart(0, 2) = -1.0*locMomentum.Py();
	locJ_SpherCart(0, 3) = 0.0;
	locJ_SpherCart(0, 4) = 0.0;
	locJ_SpherCart(1, 0) = locMomentum.Py()/locPMag;
	locJ_SpherCart(1, 1) = locMomentum.Py()/tan(locTheta);
	locJ_SpherCart(1, 2) = locMomentum.Px();
	locJ_SpherCart(1, 3) = 0.0;
	locJ_SpherCart(1, 4) = 0.0;
	locJ_SpherCart(2, 0) = locMomentum.Pz()/locPMag;
	locJ_SpherCart(2, 1) = -1.0*locPMag*sin(locTheta);
	locJ_SpherCart(2, 2) = 0.0;
	locJ_SpherCart(2, 3) = 0.0;
	locJ_SpherCart(2, 4) = 0.0;
	locJ_SpherCart(3, 0) = 0.0;
	locJ_SpherCart(3, 1) = 0.0;
	locJ_SpherCart(3, 2) = 0.0;
	locJ_SpherCart(3, 3) = -1.0;
	locJ_SpherCart(3, 4) = 1.0;
	TMatrixFSym locCov_PxPyPzT = locCov_PThetaPhiFTimeHTime.Similarity(locJ_SpherCart);

	//add time terms to the orignal covariance matrix
	(*locCovarianceMatrix)(0, 6) = locCov_PxPyPzT(0, 3);
	(*locCovarianceMatrix)(1, 6) = locCov_PxPyPzT(1, 3);
	(*locCovarianceMatrix)(2, 6) = locCov_PxPyPzT(2, 3);
	(*locCovarianceMatrix)(3, 6) = 0.0;
	(*locCovarianceMatrix)(4, 6) = 0.0;
	(*locCovarianceMatrix)(5, 6) = 0.0;
	(*locCovarianceMatrix)(6, 0) = locCov_PxPyPzT(3, 0);
	(*locCovarianceMatrix)(6, 1) = locCov_PxPyPzT(3, 1);
	(*locCovarianceMatrix)(6, 2) = locCov_PxPyPzT(3, 2);
	(*locCovarianceMatrix)(6, 3) = 0.0;
	(*locCovarianceMatrix)(6, 4) = 0.0;
	(*locCovarianceMatrix)(6, 5) = 0.0;
	(*locCovarianceMatrix)(6, 6) = locCov_PxPyPzT(3, 3);
}

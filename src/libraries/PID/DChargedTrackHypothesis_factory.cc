// $Id$
//
//    File: DChargedTrackHypothesis_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DChargedTrackHypothesis_factory.h"

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
	locEventLoop->GetSingle(dPIDAlgorithm);
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
	if (locEventRFBunches.size() == 0)
	   return NOERROR;

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); loc_i++)
	{
		DChargedTrackHypothesis* locChargedTrackHypothesis = Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBasedVector[loc_i], locDetectorMatches, locEventRFBunches[0]);
		if(locChargedTrackHypothesis != NULL)
			_data.push_back(locChargedTrackHypothesis);
	}
	sort(_data.begin(), _data.end(), DChargedTrackHypothesis_SortByEnergy);

	return NOERROR;
}

DChargedTrackHypothesis* DChargedTrackHypothesis_factory::Create_ChargedTrackHypothesis(JEventLoop* locEventLoop, const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, const DEventRFBunch* locEventRFBunch) const
{
	DChargedTrackHypothesis* locChargedTrackHypothesis = new DChargedTrackHypothesis();
	locChargedTrackHypothesis->AddAssociatedObject(locTrackTimeBased);
	locChargedTrackHypothesis->candidateid = locTrackTimeBased->candidateid;

	// Chi square and degree-of-freedom data from the track fit
	locChargedTrackHypothesis->dChiSq_Track = locTrackTimeBased->chisq;
	locChargedTrackHypothesis->dNDF_Track = locTrackTimeBased->Ndof;

	//Set DKinematicData Members
	DKinematicData *locKinematicData = locChargedTrackHypothesis;
	*locKinematicData = *(static_cast<const DKinematicData*>(locTrackTimeBased));

	DMatrixDSym locCovarianceMatrix = locChargedTrackHypothesis->errorMatrix();

	// CDC/FDC
	locChargedTrackHypothesis->setTime(locTrackTimeBased->t0());
	locChargedTrackHypothesis->setT0(locTrackTimeBased->t0(), locTrackTimeBased->t0_err(), SYS_CDC);
	locChargedTrackHypothesis->setT1(locTrackTimeBased->t0(), locTrackTimeBased->t0_err(), SYS_CDC);
	locChargedTrackHypothesis->setPathLength(numeric_limits<double>::quiet_NaN(), 0.0);
	locCovarianceMatrix(6, 6) = locTrackTimeBased->t0_err()*locTrackTimeBased->t0_err();

	// Start Counter
	DSCHitMatchParams locSCHitMatchParams;
	if(dPIDAlgorithm->Get_BestSCMatchParams(locTrackTimeBased, locDetectorMatches, locSCHitMatchParams))
	{
		locChargedTrackHypothesis->dSCHitMatchParams = locSCHitMatchParams;
		double locPropagatedTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime;
//		double locPropagatedTimeUncertainty = sqrt(locSCHitMatchParams.dHitTimeVariance + locSCHitMatchParams.dFlightTimeVariance);
		double locPropagatedTimeUncertainty = 0.3;
		locChargedTrackHypothesis->setT0(locPropagatedTime, locPropagatedTimeUncertainty, SYS_START);

//		double locFlightTimePCorrelation = locDetectorMatches->Get_FlightTimePCorrelation(locTrackTimeBased, SYS_START); //uncomment when ready!!
//		Add_TimeToTrackingMatrix(locChargedTrackHypothesis, locSCHitMatchParams.dFlightTimeVariance, locSCHitMatchParams.dHitTimeVariance, locFlightTimePCorrelation); //uncomment when ready!!
		locCovarianceMatrix(6, 6) = locPropagatedTimeUncertainty*locPropagatedTimeUncertainty; //delete when ready!!

		//add associated objects
		vector<DSCHitMatchParams> locSCHitMatchParams;
		locDetectorMatches->Get_SCMatchParams(locTrackTimeBased, locSCHitMatchParams);
		for(size_t loc_i = 0; loc_i < locSCHitMatchParams.size(); ++loc_i)
			locChargedTrackHypothesis->AddAssociatedObject(locSCHitMatchParams[loc_i].dSCHit);
	}

	// BCAL
	DShowerMatchParams locBCALShowerMatchParams;
	if(dPIDAlgorithm->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
	{
		locChargedTrackHypothesis->dBCALShowerMatchParams = locBCALShowerMatchParams;
		const DBCALShower* locBCALShower = dynamic_cast<const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
		locChargedTrackHypothesis->setT1(locBCALShower->t, locBCALShower->tErr, SYS_BCAL);
		locChargedTrackHypothesis->setTime(locBCALShower->t - locBCALShowerMatchParams.dFlightTime);
		locChargedTrackHypothesis->setPathLength(locBCALShowerMatchParams.dPathLength, 0.0);
//		double locFlightTimePCorrelation = locDetectorMatches->Get_FlightTimePCorrelation(locTrackTimeBased, SYS_BCAL); //uncomment when ready!!
//		Add_TimeToTrackingMatrix(locChargedTrackHypothesis, locBCALShowerMatchParams.dFlightTimeVariance, locBCALShower->tErr*locBCALShower->tErr, locFlightTimePCorrelation); //uncomment when ready!!
		locCovarianceMatrix(6, 6) = 0.00255*pow(locChargedTrackHypothesis->momentum().Mag(), -2.52) + 0.220; //delete when ready!!
		locCovarianceMatrix(6, 6) *= locCovarianceMatrix(6, 6); //delete when ready!!

		//add associated objects
		vector<DShowerMatchParams> locShowerMatchParams;
		locDetectorMatches->Get_BCALMatchParams(locTrackTimeBased, locShowerMatchParams);
		for(size_t loc_i = 0; loc_i < locShowerMatchParams.size(); ++loc_i)
			locChargedTrackHypothesis->AddAssociatedObject(locShowerMatchParams[loc_i].dShowerObject);
	}

	// TOF
	DTOFHitMatchParams locTOFHitMatchParams;
	if(dPIDAlgorithm->Get_BestTOFMatchParams(locTrackTimeBased, locDetectorMatches, locTOFHitMatchParams))
	{
		locChargedTrackHypothesis->dTOFHitMatchParams = locTOFHitMatchParams;
		if(locChargedTrackHypothesis->t1_detector() == SYS_CDC)
		{
			const DTOFPoint* locTOFPoint = locTOFHitMatchParams.dTOFPoint;
			locChargedTrackHypothesis->setT1(locTOFPoint->t, locTOFPoint->tErr, SYS_TOF);
			locChargedTrackHypothesis->setTime(locTOFPoint->t - locTOFHitMatchParams.dFlightTime);
			locChargedTrackHypothesis->setPathLength(locTOFHitMatchParams.dPathLength, 0.0);
//			double locFlightTimePCorrelation = locDetectorMatches->Get_FlightTimePCorrelation(locTrackTimeBased, SYS_TOF); //uncomment when ready!!
//			Add_TimeToTrackingMatrix(locChargedTrackHypothesis, locTOFHitMatchParams.dFlightTimeVariance, locTOFPoint->tErr*locTOFPoint->tErr, locFlightTimePCorrelation); //uncomment when ready!!
			locCovarianceMatrix(6, 6) = 0.08*0.08; //delete when ready!!
		}

		//add associated objects
		vector<DTOFHitMatchParams> locTOFHitMatchParams;
		locDetectorMatches->Get_TOFMatchParams(locTrackTimeBased, locTOFHitMatchParams);
		for(size_t loc_i = 0; loc_i < locTOFHitMatchParams.size(); ++loc_i)
			locChargedTrackHypothesis->AddAssociatedObject(locTOFHitMatchParams[loc_i].dTOFPoint);
	}

	//FCAL
	DShowerMatchParams locFCALShowerMatchParams;
	if(dPIDAlgorithm->Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams))
	{
		locChargedTrackHypothesis->dFCALShowerMatchParams = locFCALShowerMatchParams;
		if(locChargedTrackHypothesis->t1_detector() == SYS_CDC)
		{
			const DFCALShower* locFCALShower = dynamic_cast<const DFCALShower*>(locFCALShowerMatchParams.dShowerObject);
//			locChargedTrackHypothesis->setT1(locFCALShower->getTime(), sqrt(locFCALShower->dCovarianceMatrix(4, 4)), SYS_FCAL); //uncomment when ready!!
			locChargedTrackHypothesis->setT1(locFCALShower->getTime(), 0.5, SYS_FCAL);
			locChargedTrackHypothesis->setTime(locFCALShower->getTime() - locFCALShowerMatchParams.dFlightTime);
			locChargedTrackHypothesis->setPathLength(locFCALShowerMatchParams.dPathLength, 0.0);
//			double locFlightTimePCorrelation = locDetectorMatches->Get_FlightTimePCorrelation(locTrackTimeBased, SYS_FCAL); //uncomment when ready!!
//			Add_TimeToTrackingMatrix(locChargedTrackHypothesis, locFCALShowerMatchParams.dFlightTimeVariance, locFCALShower->dCovarianceMatrix(4, 4), locFlightTimePCorrelation); //uncomment when ready!!
			locCovarianceMatrix(6, 6) = 0.6*0.6; // straight-line fit to high momentum data //delete when ready!!
		}

		//add associated objects
		vector<DShowerMatchParams> locShowerMatchParams;
		locDetectorMatches->Get_FCALMatchParams(locTrackTimeBased, locShowerMatchParams);
		for(size_t loc_i = 0; loc_i < locShowerMatchParams.size(); ++loc_i)
			locChargedTrackHypothesis->AddAssociatedObject(locShowerMatchParams[loc_i].dShowerObject);
	}

	locChargedTrackHypothesis->setErrorMatrix(locCovarianceMatrix);

	//Calculate PID ChiSq, NDF, FOM
	dPIDAlgorithm->Calc_ChargedPIDFOM(locChargedTrackHypothesis, locEventRFBunch, true);

	return locChargedTrackHypothesis;
}

void DChargedTrackHypothesis_factory::Add_TimeToTrackingMatrix(DChargedTrackHypothesis* locChargedTrackHypothesis, double locFlightTimeVariance, double locHitTimeVariance, double locFlightTimePCorrelation) const
{
	DMatrixDSym locCovarianceMatrix = locChargedTrackHypothesis->errorMatrix();
	DVector3 locMomentum = locChargedTrackHypothesis->momentum();

	//extract momentum 3x3 portion of covariance matrix
	DMatrixDSym locCov_PxPyPz(3);
	for(size_t loc_i = 0; loc_i < 3; ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < 3; ++loc_j)
			locCov_PxPyPz(loc_i, loc_j) = locCovarianceMatrix(loc_i, loc_j);
	}

	//convert px/py/pz to p/theta/phi
	double locPMag = locMomentum.Mag();
	double locPPerpSq = locMomentum.Px()*locMomentum.Px() + locMomentum.Py()*locMomentum.Py();
	DMatrix locJ_CartSpher(3, 3); //jacobian
	locJ_CartSpher(0, 0) = locMomentum.Px()/locPMag;
	locJ_CartSpher(0, 1) = locMomentum.Py()/locPMag;
	locJ_CartSpher(0, 2) = locMomentum.Pz()/locPMag;
	locJ_CartSpher(1, 0) = locMomentum.Px()*locMomentum.Pz()/(locPMag*locPMag*sqrt(locPPerpSq));
	locJ_CartSpher(1, 1) = locMomentum.Py()*locMomentum.Pz()/(locPMag*locPMag*sqrt(locPPerpSq));
	locJ_CartSpher(1, 2) = sqrt(locPPerpSq)/(locPMag*locPMag);
	locJ_CartSpher(2, 0) = -1.0*locMomentum.Py()/locPPerpSq;
	locJ_CartSpher(2, 1) = locMomentum.Px()/locPPerpSq;
	locJ_CartSpher(2, 2) = 0.0;
	DMatrixDSym locCov_PThetaPhi = locCov_PxPyPz.Similarity(locJ_CartSpher);

	//add flight time and hit time to p/theta/phi covariance matrix
	DMatrixDSym locCov_PThetaPhiFTimeHTime(5);
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
	DMatrix locJ_SpherCart(4, 5); //jacobian
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
	DMatrixDSym locCov_PxPyPzT = locCov_PThetaPhiFTimeHTime.Similarity(locJ_SpherCart);

	//add time terms to the orignal covariance matrix
	locCovarianceMatrix(0, 6) = locCov_PxPyPzT(0, 3);
	locCovarianceMatrix(1, 6) = locCov_PxPyPzT(1, 3);
	locCovarianceMatrix(2, 6) = locCov_PxPyPzT(2, 3);
	locCovarianceMatrix(3, 6) = 0.0;
	locCovarianceMatrix(4, 6) = 0.0;
	locCovarianceMatrix(5, 6) = 0.0;
	locCovarianceMatrix(6, 0) = locCov_PxPyPzT(3, 0);
	locCovarianceMatrix(6, 1) = locCov_PxPyPzT(3, 1);
	locCovarianceMatrix(6, 2) = locCov_PxPyPzT(3, 2);
	locCovarianceMatrix(6, 3) = 0.0;
	locCovarianceMatrix(6, 4) = 0.0;
	locCovarianceMatrix(6, 5) = 0.0;
	locCovarianceMatrix(6, 6) = locCov_PxPyPzT(3, 3);

	locChargedTrackHypothesis->setErrorMatrix(locCovarianceMatrix);
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



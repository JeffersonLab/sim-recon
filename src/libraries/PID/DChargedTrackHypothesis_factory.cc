// $Id$
//
//    File: DChargedTrackHypothesis_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <TMath.h>

#include "DChargedTrackHypothesis_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DChargedTrackHypothesis_factory::init(void)
{
	//this parameter controls what BCAL reconstruction algorithm to use
	//the same parameter is used in DNeutralShowerCandidate_factory
	USE_KLOE = 1;
	gPARMS->SetDefaultParameter("BCALRECON:USE_KLOE", USE_KLOE);

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
	dPIDAlgorithm = const_cast<DParticleID*>(locPIDAlgorithms[0]);
  
	// Warn user if something happened that caused us NOT to get a dPIDAlgorithm object pointer
	if(!dPIDAlgorithm){
		_DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	dTargetZCenter = 0.0;
	// Get Target parameters from XML
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(runnumber):NULL;
	if(locGeometry)
		locGeometry->GetTargetZ(dTargetZCenter);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrackHypothesis_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	unsigned int loc_i;

	const DTrackTimeBased *locTrackTimeBased;
	DChargedTrackHypothesis *locChargedTrackHypothesis;

	unsigned int locTOFIndex, locSCIndex;
	double locTempProjectedTime = 0.0, locProjectedTime = 0.0, locPathLength = 0.0, locProjectedTimeUncertainty = 0.0, locFlightTime = 0.0;
	double locPropagatedRFTime, locSTHitTime, locSTRFTimeDifference, locTimeDifference, locTimingChiSq;

	bool locMatchedOuterDetectorFlag;
	double locRFTime = 0.0;
	double locRFBunchFrequency = 2.004;
	double locChiSq_Total, locChiSq_DCdEdx;
	unsigned int locNDF_Total, locNDF_DCdEdx;
	bool locUseDCdEdxForPIDFlag;

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	vector<const DTOFPoint*> locTOFPoints;
	vector<const DBCALShower*> locBCALShowers;
	vector<const DFCALShower*> locFCALShowers;
	vector<const DBCALShower*> locMatchedBCALShowers;
	vector<const DFCALShower*> locMatchedFCALShowers;
	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locTrackTimeBasedVector);
	locEventLoop->Get(locTOFPoints);
	if (USE_KLOE) {
		locEventLoop->Get(locBCALShowers, "KLOE");
	} else { 
		locEventLoop->Get(locBCALShowers);
	}
	locEventLoop->Get(locFCALShowers);

  for (loc_i = 0; loc_i < locTrackTimeBasedVector.size(); loc_i++){
		locChargedTrackHypothesis = new DChargedTrackHypothesis();
		locTrackTimeBased = locTrackTimeBasedVector[loc_i];
		locChargedTrackHypothesis->AddAssociatedObject(locTrackTimeBased);
		locChargedTrackHypothesis->candidateid = locTrackTimeBased->candidateid;
		
		// Chi square and degree-of-freedom data from the track fit
		locChargedTrackHypothesis->dChiSq_Track=locTrackTimeBased->chisq;
		locChargedTrackHypothesis->dNDF_Track=locTrackTimeBased->Ndof;	
		
		//Set DKinematicData Members
		DKinematicData *locKinematicData = locChargedTrackHypothesis;
		*locKinematicData = *locTrackTimeBased;
		locChargedTrackHypothesis->dPID = dPIDAlgorithm->IDTrack(locChargedTrackHypothesis->charge(), locChargedTrackHypothesis->mass()); //mass used in track fit to create DTrackWireBased

		// Calculate DC dE/dx ChiSq
		// Compute the dEdx for the hits on the track
		locUseDCdEdxForPIDFlag = false; //true when enabled
		if(dPIDAlgorithm->CalcDCdEdxChiSq(locChargedTrackHypothesis, locChiSq_DCdEdx, locNDF_DCdEdx) == NOERROR)
			locUseDCdEdxForPIDFlag = true;
		locChargedTrackHypothesis->dChiSq_DCdEdx = locChiSq_DCdEdx;
		locChargedTrackHypothesis->dNDF_DCdEdx = locNDF_DCdEdx;

		// Initialize projected time to estimate from track
		locProjectedTime = locChargedTrackHypothesis->t0();

		locMatchedOuterDetectorFlag = false;
		// Try matching the track with hits in the outer detectors
		locMatchedBCALShowers.resize(0);
		locMatchedFCALShowers.resize(0);
		locTempProjectedTime = locChargedTrackHypothesis->t0(); // to reject hits that are not in time with the track
		if (dPIDAlgorithm->MatchToBCAL(locTrackTimeBased->rt, locBCALShowers, locMatchedBCALShowers, locTempProjectedTime, locPathLength, locFlightTime) == NOERROR){
			for(unsigned int loc_j = 0; loc_j < locMatchedBCALShowers.size(); ++loc_j)
				locChargedTrackHypothesis->AddAssociatedObject(locMatchedBCALShowers[loc_j]);
			locProjectedTime = locTempProjectedTime;
			locProjectedTimeUncertainty = 0.00255*pow(locChargedTrackHypothesis->momentum().Mag(), -2.52) + 0.220;
			locChargedTrackHypothesis->setT1(locMatchedBCALShowers[0]->t, locMatchedBCALShowers[0]->tErr, SYS_BCAL);
			locChargedTrackHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
			locMatchedOuterDetectorFlag = true;
		}
		locTempProjectedTime = locChargedTrackHypothesis->t0(); // to reject hits that are not in time with the track
		if (dPIDAlgorithm->MatchToTOF(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locTOFPoints, locTempProjectedTime, locTOFIndex, locPathLength, locFlightTime) == NOERROR){
			locChargedTrackHypothesis->AddAssociatedObject(locTOFPoints[locTOFIndex]);
			if (locMatchedOuterDetectorFlag == false){
				locProjectedTime = locTempProjectedTime;
				locProjectedTimeUncertainty = 0.08;
				locChargedTrackHypothesis->setT1(locTOFPoints[locTOFIndex]->t, NaN, SYS_TOF);
				locChargedTrackHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
				locMatchedOuterDetectorFlag = true;
			}
		}
		locTempProjectedTime = locChargedTrackHypothesis->t0(); // to reject hits that are not in time with the track
		if (dPIDAlgorithm->MatchToFCAL(locTrackTimeBased->rt, locFCALShowers, locMatchedFCALShowers, locTempProjectedTime, locPathLength, locFlightTime) == NOERROR){
			for(unsigned int loc_j = 0; loc_j < locMatchedFCALShowers.size(); ++loc_j)
				locChargedTrackHypothesis->AddAssociatedObject(locMatchedFCALShowers[loc_j]);
			if (locMatchedOuterDetectorFlag == false){
				locMatchedOuterDetectorFlag = true;
				locProjectedTime = locTempProjectedTime;
				locProjectedTimeUncertainty = 0.6; // straight-line fit to high momentum data
				locChargedTrackHypothesis->setT1(locMatchedFCALShowers[0]->getTime(), NaN, SYS_FCAL);
				locChargedTrackHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
			}
		}
		if (locMatchedOuterDetectorFlag == false){
			locTempProjectedTime = locChargedTrackHypothesis->t0(); // to reject hits that are not in time with the track
			if (dPIDAlgorithm->MatchToSC(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locSCHits, locTempProjectedTime, locSCIndex, locPathLength, locFlightTime) == NOERROR){
				locProjectedTime = locTempProjectedTime;
				locProjectedTimeUncertainty = 1.0; // uncertainty unknown
				locChargedTrackHypothesis->setT1(locSCHits[locSCIndex]->t, NaN, SYS_START);
				locChargedTrackHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
				locChargedTrackHypothesis->AddAssociatedObject(locSCHits[locSCIndex]);
			}else{
				locProjectedTime = locChargedTrackHypothesis->t0();
				locChargedTrackHypothesis->setT1(NaN, 0.0, SYS_NULL); //initialize
				locChargedTrackHypothesis->setPathLength(NaN, 0.0); //zero uncertainty (for now)
			}
		}

		//Calculate PID ChiSq, NDF, FOM
		locNDF_Total = 0;
		locChiSq_Total = 0.0;
		//locNDF_Total=1;
		//locChiSq_Total=	locChargedTrackHypothesis->dChiSq_Track/float(locChargedTrackHypothesis->dNDF_Track);

		if((locChargedTrackHypothesis->t0_detector() == SYS_START) && (locMatchedOuterDetectorFlag == true)){ //use timing info to determine particle ID
			// Use ST hit to select RF beam bucket
			locPropagatedRFTime = locRFTime + (locChargedTrackHypothesis->z() - dTargetZCenter)/SPEED_OF_LIGHT;
			locSTHitTime = locChargedTrackHypothesis->t0();
			locSTRFTimeDifference = locSTHitTime - locPropagatedRFTime; 
			while(fabs(locSTRFTimeDifference) > locRFBunchFrequency/2.0){
				locPropagatedRFTime += (locSTRFTimeDifference > 0.0) ? locRFBunchFrequency : -1.0*locRFBunchFrequency;
				locSTRFTimeDifference = locSTHitTime - locPropagatedRFTime;
			}
			// Compare time difference between RF & TOF/BCAL/FCAL times at the vertex
			locTimeDifference = locPropagatedRFTime - locProjectedTime;
			// Calculate ChiSq, FOM
			locTimingChiSq = locTimeDifference*locTimeDifference/(locProjectedTimeUncertainty*locProjectedTimeUncertainty);
			locChargedTrackHypothesis->dChiSq_Timing = locTimingChiSq;
			locChargedTrackHypothesis->dNDF_Timing = 1;

			locChiSq_Total += locTimingChiSq;
			locNDF_Total += 1;
		}else{ //not enough timing information
			locChargedTrackHypothesis->dChiSq_Timing = 0.0;
			locChargedTrackHypothesis->dNDF_Timing = 0;
		}
		if(locUseDCdEdxForPIDFlag == true){
			locChiSq_Total += locChargedTrackHypothesis->dChiSq_DCdEdx;
			locNDF_Total += locChargedTrackHypothesis->dNDF_DCdEdx;
		}
		locChargedTrackHypothesis->dChiSq = locChiSq_Total;
		locChargedTrackHypothesis->dNDF = locNDF_Total;
		locChargedTrackHypothesis->dFOM = (locNDF_Total > 0) ? TMath::Prob(locChiSq_Total, locNDF_Total) : NaN;

		_data.push_back(locChargedTrackHypothesis);
	}

	return NOERROR;
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



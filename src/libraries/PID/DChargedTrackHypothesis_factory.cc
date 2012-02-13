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

	unsigned int locBCALIndex, locFCALIndex, locTOFIndex, locSCIndex;
	double locProjectedTime = 0.0, locPathLength = 0.0, locProjectedTimeUncertainty = 0.0, locFlightTime = 0.0;
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
		locChargedTrackHypothesis->dTrackTimeBased = locTrackTimeBased;

		locChargedTrackHypothesis->dPID = dPIDAlgorithm->IDTrack(locTrackTimeBased->charge(), locTrackTimeBased->mass()); //mass used in track fit to create DTrackWireBased

			// Calculate DC dE/dx ChiSq
			// Compute the dEdx for the hits on the track
			locUseDCdEdxForPIDFlag = false; //true when enabled
			if(dPIDAlgorithm->CalcDCdEdxChiSq(locChargedTrackHypothesis, locChiSq_DCdEdx, locNDF_DCdEdx) == NOERROR)
				locUseDCdEdxForPIDFlag = true;
			locChargedTrackHypothesis->dChiSq_DCdEdx = locChiSq_DCdEdx;
			locChargedTrackHypothesis->dNDF_DCdEdx = locNDF_DCdEdx;

		// Initialize projected time to estimate from track
		locProjectedTime = locTrackTimeBased->t0();

		locMatchedOuterDetectorFlag = false;
		// Try matching the track with hits in the outer detectors
		if (dPIDAlgorithm->MatchToBCAL(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locBCALShowers, locProjectedTime, locBCALIndex, locPathLength, locFlightTime) == NOERROR){
			locChargedTrackHypothesis->dProjectedTime = locProjectedTime;
			locProjectedTimeUncertainty = 0.00255*pow(locTrackTimeBased->momentum().Mag(), -2.52) + 0.220;
			locChargedTrackHypothesis->dFlightTime = locFlightTime;
			locChargedTrackHypothesis->dPathLength = locPathLength;
			locChargedTrackHypothesis->AddAssociatedObject(locBCALShowers[locBCALIndex]);
			locChargedTrackHypothesis->dMatchedTimeDetector = SYS_BCAL;
			locMatchedOuterDetectorFlag = true;
		}
		if (dPIDAlgorithm->MatchToTOF(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locTOFPoints, locProjectedTime, locTOFIndex, locPathLength, locFlightTime) == NOERROR){
			locChargedTrackHypothesis->AddAssociatedObject(locTOFPoints[locTOFIndex]);
			if (locMatchedOuterDetectorFlag == false){
				locChargedTrackHypothesis->dProjectedTime = locProjectedTime;
				locProjectedTimeUncertainty = 0.08;
				locChargedTrackHypothesis->dFlightTime = locFlightTime;
				locChargedTrackHypothesis->dPathLength = locPathLength;
				locChargedTrackHypothesis->dMatchedTimeDetector = SYS_TOF;
				locMatchedOuterDetectorFlag = true;
			}
		}
		if (dPIDAlgorithm->MatchToFCAL(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locFCALShowers, locProjectedTime, locFCALIndex, locPathLength, locFlightTime) == NOERROR){
			locChargedTrackHypothesis->AddAssociatedObject(locFCALShowers[locFCALIndex]);
			if (locMatchedOuterDetectorFlag == false){
				locMatchedOuterDetectorFlag = true;
				locChargedTrackHypothesis->dProjectedTime = locProjectedTime;
				locProjectedTimeUncertainty = 0.6; // straight-line fit to high momentum data
				locChargedTrackHypothesis->dFlightTime = locFlightTime;
				locChargedTrackHypothesis->dPathLength = locPathLength;
				locChargedTrackHypothesis->dMatchedTimeDetector = SYS_FCAL;
			}
		}
		if (locMatchedOuterDetectorFlag == false){
			if (dPIDAlgorithm->MatchToSC(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locSCHits, locProjectedTime, locSCIndex, locPathLength, locFlightTime) == NOERROR){
				locChargedTrackHypothesis->dProjectedTime = locProjectedTime;
				locProjectedTimeUncertainty = 1.0; // uncertainty unknown
				locChargedTrackHypothesis->dFlightTime = locFlightTime;
				locChargedTrackHypothesis->dPathLength = locPathLength;
				locChargedTrackHypothesis->dMatchedTimeDetector = SYS_START;
				locChargedTrackHypothesis->AddAssociatedObject(locSCHits[locSCIndex]);
			}else{
				locChargedTrackHypothesis->dProjectedTime = NaN; //initialize
				locChargedTrackHypothesis->dPathLength = NaN; //initialize
				locChargedTrackHypothesis->dFlightTime = NaN; //initialize
				locChargedTrackHypothesis->dMatchedTimeDetector = SYS_NULL;
			}
		}

		//Calculate PID ChiSq, NDF, FOM
		locNDF_Total = 0;
		locChiSq_Total = 0.0;
		if((locTrackTimeBased->t0_detector() == SYS_START) && (locMatchedOuterDetectorFlag == true)){ //use timing info to determine particle ID
			// Use ST hit to select RF beam bucket
			locPropagatedRFTime = locRFTime + (locTrackTimeBased->z() - dTargetZCenter)/SPEED_OF_LIGHT;
			locSTHitTime = locTrackTimeBased->t0();
			locSTRFTimeDifference = locSTHitTime - locPropagatedRFTime; 
			while(fabs(locSTRFTimeDifference) > locRFBunchFrequency/2.0){
				locPropagatedRFTime += (locSTRFTimeDifference > 0.0) ? locRFBunchFrequency : -1.0*locRFBunchFrequency;
				locSTRFTimeDifference = locSTHitTime - locPropagatedRFTime;
			}
			// Compare time difference between RF & TOF/BCAL/FCAL times at the vertex
			locTimeDifference = locPropagatedRFTime - locChargedTrackHypothesis->dProjectedTime;
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



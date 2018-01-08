// $Id$
//
//    File: DEventProcessor_coherent_peak_skim.cc
// Created: Tue Mar 28 10:57:49 EDT 2017
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-696.el6.x86_64 x86_64)
//

#include "DEventProcessor_coherent_peak_skim.h"

// Routine used to create our DEventProcessor

extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new DEventProcessor_coherent_peak_skim()); //register this plugin
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_coherent_peak_skim::init(void)
{
	// This is called once at program startup.

	dTimingCutMap[Gamma][SYS_BCAL] = 3.0;
	dTimingCutMap[Gamma][SYS_FCAL] = 5.0;
	dTimingCutMap[Proton][SYS_NULL] = -1.0;
	dTimingCutMap[Proton][SYS_TOF] = 2.5;
	dTimingCutMap[Proton][SYS_BCAL] = 2.5;
	dTimingCutMap[Proton][SYS_FCAL] = 3.0;
	dTimingCutMap[PiPlus][SYS_NULL] = -1.0;
	dTimingCutMap[PiPlus][SYS_TOF] = 2.0;
	dTimingCutMap[PiPlus][SYS_BCAL] = 2.5;
	dTimingCutMap[PiPlus][SYS_FCAL] = 3.0;
	dTimingCutMap[PiMinus][SYS_NULL] = -1.0;
	dTimingCutMap[PiMinus][SYS_TOF] = 2.0;
	dTimingCutMap[PiMinus][SYS_BCAL] = 2.5;
	dTimingCutMap[PiMinus][SYS_FCAL] = 3.0;
	dTimingCutMap[Electron][SYS_NULL] = -1.0;
	dTimingCutMap[Electron][SYS_TOF] = 2.0;
	dTimingCutMap[Electron][SYS_BCAL] = 2.5;
	dTimingCutMap[Electron][SYS_FCAL] = 3.0;
	dTimingCutMap[Positron][SYS_NULL] = -1.0;
	dTimingCutMap[Positron][SYS_TOF] = 2.0;
	dTimingCutMap[Positron][SYS_BCAL] = 2.5;
	dTimingCutMap[Positron][SYS_FCAL] = 3.0;

	dShowerEOverPCut = 0.75;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_coherent_peak_skim::brun(jana::JEventLoop* locEventLoop, int32_t locRunNumber)
{
	// This is called whenever the run number changes
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	dCoherentPeakRange = pair<double, double>(8.4, 9.0);
	map<string, double> photon_beam_param;
	if(locEventLoop->GetCalib("/ANALYSIS/beam_asymmetry/coherent_energy", photon_beam_param) == false)
		dCoherentPeakRange = pair<double, double>(photon_beam_param["cohmin_energy"], photon_beam_param["cohedge_energy"]);
	dCoherentPeakRange.first -= 0.2; //same as below
	dCoherentPeakRange.second += 0.2; //in case the high-side edge is fluctuating a lot

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locRunNumber);
	locGeometry->GetTargetZ(dTargetCenterZ);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_coherent_peak_skim::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// locEventLoop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	//
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// locEventLoop->Get(mydataclasses);
	//
	// japp->RootFillLock(this);
	//  ... fill historgrams or trees ...
	// japp->RootFillUnLock(this);

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software

	//THE BELOW ALGORITHM IS BAD.  WHY?
	//Because junk or accidental tracks could have voted on which is the correct RF bunch

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "PreSelect");

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);
/*
	bool locIsHadronicEventFlag = false;
	for(auto locTrack : locChargedTracks)
	{
		//make sure at least one track isn't a lepton!! (read: e+/-. assume all muons come from pion decays)
		auto locChargedHypo = locTrack->Get_Hypothesis(PiPlus);
		if(locChargedHypo == nullptr)
			locChargedHypo = locTrack->Get_Hypothesis(PiMinus);
		if(locChargedHypo == nullptr)
			locChargedHypo = locTrack->Get_BestFOM();
		if(locChargedHypo == nullptr)
			continue; //should be impossible!!

		//timing cut: is it consistent with an e+/-??
		auto locDetector = locChargedHypo->t1_detector();
		double locDeltaT = locChargedHypo->time() - locChargedHypo->t0();
		if(fabs(locDeltaT) > dTimingCutMap[Electron][locDetector])
		{
			locIsHadronicEventFlag = true; //not an electron!!
			break;
		}

		//compute shower-E/p, cut
		if(Cut_ShowerEOverP(locChargedHypo))
		{
			locIsHadronicEventFlag = true; //not an electron!!
			break;
		}
	}
*/
	bool locIsTrackEventFlag = !locChargedTracks.empty();

	//see if is in coherent peak
	bool locCoherentPeakFlag = false;
	if(std::isnan(locEventRFBunch->dTime))
		locCoherentPeakFlag = true;
	for(auto& locBeamPhoton : locBeamPhotons)
	{
		if((locBeamPhoton->energy() < dCoherentPeakRange.first) || (locBeamPhoton->energy() > dCoherentPeakRange.second))
			continue; //not in coherent peak

		double locBeamRFDeltaT = locBeamPhoton->time() - locEventRFBunch->dTime;
		if(fabs(locBeamRFDeltaT) > 0.5*dBeamBunchPeriod)
			continue;
		locCoherentPeakFlag = true;
		break;
	}

	/******************************************************** OPTIONAL: SKIMS *******************************************************/

	//Optional: Save event to output REST file. Use this to create physical skims.
	const DEventWriterREST* locEventWriterREST = NULL;
	locEventLoop->GetSingle(locEventWriterREST);
//	if(locIsHadronicEventFlag && locCoherentPeakFlag)
//		locEventWriterREST->Write_RESTEvent(locEventLoop, "hadronic_coherent_peak"); //string is part of output file name
	if(locCoherentPeakFlag && locIsTrackEventFlag)
		locEventWriterREST->Write_RESTEvent(locEventLoop, "coherent_peak"); //string is part of output file name

	return NOERROR;
}

bool DEventProcessor_coherent_peak_skim::Cut_ShowerEOverP(const DChargedTrackHypothesis* locChargedHypo) const
{
	//compute shower-E/p, cut
	double locP = locChargedHypo->momentum().Mag();
	double locShowerEOverP = -1.0;
	auto locFCALShowerMatchParams = locChargedHypo->Get_FCALShowerMatchParams();
	auto locBCALShowerMatchParams = locChargedHypo->Get_BCALShowerMatchParams();
	if(locFCALShowerMatchParams != NULL)
	{
		const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
		locShowerEOverP = locFCALShower->getEnergy()/locP;
	}
	else if(locBCALShowerMatchParams != NULL)
	{
		const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
		locShowerEOverP = locBCALShower->E/locP;
	}

	return (locShowerEOverP <= dShowerEOverPCut);
}

double DEventProcessor_coherent_peak_skim::Step_TimeToNearInputTime(double locTimeToStep, double locTimeToStepTo) const
{
	double locDeltaT = locTimeToStepTo - locTimeToStep;
	int locNumBucketsToShift = (locDeltaT > 0.0) ? int(locDeltaT/dBeamBunchPeriod + 0.5) : int(locDeltaT/dBeamBunchPeriod - 0.5);
	return (locTimeToStep + dBeamBunchPeriod*double(locNumBucketsToShift));
}


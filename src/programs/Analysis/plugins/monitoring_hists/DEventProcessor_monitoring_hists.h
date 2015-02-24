// $Id$
//
//    File: DEventProcessor_monitoring_hists.h
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_monitoring_hists_
#define _DEventProcessor_monitoring_hists_

#include "TFile.h"
#include "TROOT.h"

#include "JANA/JEventProcessor.h"

#include "DANA/DApplication.h"
#include "ANALYSIS/DHistogramActions.h"

using namespace jana;

class DEventProcessor_monitoring_hists : public JEventProcessor
{
	public:
		DEventProcessor_monitoring_hists(){};
		~DEventProcessor_monitoring_hists(){};
		const char* className(void){return "DEventProcessor_monitoring_hists";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DHistogramAction_TrackMultiplicity dHistogramAction_TrackMultiplicity;
		DHistogramAction_ThrownParticleKinematics dHistogramAction_ThrownParticleKinematics;
		DHistogramAction_DetectedParticleKinematics dHistogramAction_DetectedParticleKinematics;
		DHistogramAction_GenReconTrackComparison dHistogramAction_GenReconTrackComparison;
		DHistogramAction_NumReconstructedObjects dHistogramAction_NumReconstructedObjects;
		DHistogramAction_ReconnedThrownKinematics dHistogramAction_ReconnedThrownKinematics;
		DHistogramAction_EventVertex dHistogramAction_EventVertex;

		DHistogramAction_DetectorMatchParams dHistogramAction_DetectorMatchParams;
		DHistogramAction_Neutrals dHistogramAction_Neutrals;
		DHistogramAction_DetectorPID dHistogramAction_DetectorPID;
		DHistogramAction_DetectorMatching dHistogramAction_DetectorMatching;
		DHistogramAction_Reconstruction dHistogramAction_Reconstruction;
//		DHistogramAction_ObjectMemory dHistogramAction_ObjectMemory;

		TH1D* dHist_IsEvent; //counts events
};

#endif // _DEventProcessor_monitoring_hists_


// $Id$
//
//    File: JEventProcessor_dirc_truth_hists.h
// Created: Thu Mar 19 10:25:58 EDT 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_dirc_truth_hists_
#define _JEventProcessor_dirc_truth_hists_

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TVector3.h"

#include <JANA/JEventProcessor.h>

#include <TRACKING/DTrackTimeBased.h>
#include <DIRC/DDIRCTruthPoint.h>
#include <DIRC/DDIRCTruthHit.h>
#include <DIRC/DDIRCHit.h>

#include <TRACKING/DMCThrown.h>

class JEventProcessor_dirc_truth_hists:public jana::JEventProcessor{
	public:
		JEventProcessor_dirc_truth_hists();
		~JEventProcessor_dirc_truth_hists();
		const char* className(void){return "JEventProcessor_dirc_truth_hists";}
		enum { kTruthHitMax = 5000 };
		enum { kGenTrkMax = 20 };

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		TTree* dircTree;
		Int_t event;

		Int_t nGenTrk;
		Float_t m_genTrk_px[kGenTrkMax];
		Float_t m_genTrk_py[kGenTrkMax];
		Float_t m_genTrk_pz[kGenTrkMax];
		Float_t m_genTrk_E[kGenTrkMax];
		Float_t m_genTrk_x[kGenTrkMax];
		Float_t m_genTrk_y[kGenTrkMax];
		Float_t m_genTrk_z[kGenTrkMax];
		Float_t m_genTrk_t[kGenTrkMax];

		Int_t nTruthHit;
		Float_t m_hit_x[kTruthHitMax];
		Float_t m_hit_y[kTruthHitMax];
		Float_t m_hit_z[kTruthHitMax];
		Float_t m_hit_t[kTruthHitMax];
		Float_t m_hitpixel_w[kTruthHitMax];
		Float_t m_hitpixel_y[kTruthHitMax];

		TH1F *hTruthHitT, *hTruthHitDeltaT, *hTruthHitE, *hTruthHitLambda, *hTruthPointM;
		TH2F *hTruthPoint, *hTruthHitXY, *hTruthHitZX, *hTruthHitYLocW, *hTruthHitMissingYLocW;
		TH2F *hTruthHitIncidentAngleLocW, *hTruthHitIncidentAngleY;
		TH3F *hTruthHitYLocWT;
		TH2F *hNTruthHit;

};

#endif // _JEventProcessor_dirc_truth_hists_


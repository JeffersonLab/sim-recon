// $Id$
//
//    File: DEventProcessor_track_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_track_hists_
#define _DEventProcessor_track_hists_

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include "JANA/JFactory.h"
#include "JANA/JEventProcessor.h"
#include "JANA/JEventLoop.h"

class DEventProcessor_track_hists:public JEventProcessor{

	public:
		DEventProcessor_track_hists();
		~DEventProcessor_track_hists();
		void FillAll(float what, int nhits, float theta, float phi, float p, float weight=1.0);
		void EffVsX(TH1F *out, TH2F* in, int numerator=NMATCHED);

		enum{
			NHITS_THROWN =1,
			NHITS_FOUND,
			NHITS_THROWN_AND_FOUND,
			NHITS_FOUND_DIFFERENT,
			NHITS_THROWN_UNUSED,
			NTHROWN,
			NFOUND,
			NFITTABLE,
			NMATCHED,
			NHITMATCHED,

			R_FOUND_TO_FITTABLE,
			R_THROWN_AND_FOUND_TO_THROWN,
			R_THROWN_AND_FOUND_TO_FOUND,
			R_MATCHED_TO_FITTABLE,
			R_HITMATCHED_TO_FITTABLE,
			
			NBINS
		};
		
		TH1F *stats;
		TH1F *frac_from_thrown;
		TH2F *stats_vs_theta, *stats_vs_phi, *stats_vs_p, *stats_vs_nhits;
		TH2F *dp_over_p_vs_p, *dp_over_p_vs_theta, *dpcandidate_over_p_vs_theta;
		TH2F *pthrown_over_pfound_vs_p, *sinthrown_over_sinfound_vs_sin;
		TH2F *pcandidatethrown_over_pfound_vs_p;
		TH2F *phithrown_over_phifound_vs_phi;
		TH1F *eff_vs_theta, *eff_vs_phi, *eff_vs_p, *eff_vs_nhits;
		TH1F *eff_vs_theta_hm, *eff_vs_phi_hm, *eff_vs_p_hm, *eff_vs_nhits_hm;
		TH1F *dist_same, *dist_diff;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		unsigned int Nevents;
		unsigned int Ncdchits;
		unsigned int Nfdchits;
};

#endif // _DEventProcessor_track_hists_


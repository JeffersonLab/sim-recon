// $Id$
//
//    File: track2.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _track2_
#define _track2_

#include <TObject.h>
#include <TVector3.h>


class track2:public TObject{

	public:

		unsigned long event;

		float q_thrown;
		int PID_thrown;
		TVector3 pthrown;
		bool isreconstructable;

		int Nstereo;
		int Ncdc;
		int Nfdc;
		int NLR_bad_stereo;
		int NLR_bad;

		int PID_hypothesized;
		float q_hypothesized;
		float FOM_hypothesized;

		float q_timebased;
		int PID_timebased;
		bool dTrackReconstructedFlag_TimeBased;
		TVector3 pfit;
		double trk_chisq;
		int trk_Ndof;
		double delta_pt_over_pt;
		double delta_theta;	// mrad
		double delta_phi;		// mrad
		int num_timebased;

		float q_wirebased;
		int PID_wirebased;
		bool dTrackReconstructedFlag_WireBased;
		TVector3 pfit_wire;
		double trk_chisq_wb;
		int trk_Ndof_wb;
		double delta_pt_over_pt_wire;
		double delta_theta_wire;	// mrad
		double delta_phi_wire;		// mrad
		int num_wirebased;

		float q_candidate;
		int PID_candidate;
		bool dTrackReconstructedFlag_Candidate;
		TVector3 pcan;
		double delta_pt_over_pt_can;
		double delta_theta_can;	// mrad
		double delta_phi_can;		// mrad
		int num_candidates;

	private:
		ClassDef(track2,1);

};

#endif // _track2_


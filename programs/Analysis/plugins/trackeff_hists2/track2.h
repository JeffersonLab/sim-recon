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

		TVector3 pthrown;
		TVector3 pfit;
		TVector3 pfit_wire;
		TVector3 pcan;
		double trk_chisq;
		int trk_Ndof;
		double trk_chisq_wb;
		int trk_Ndof_wb;
		double delta_pt_over_pt;
		double delta_theta;	// mrad
		double delta_phi;		// mrad
		bool isreconstructable;
		int Nstereo;
		int Ncdc;
		int Nfdc;
		int NLR_bad_stereo;
		int NLR_bad;
		unsigned long event;

	private:
		ClassDef(track2,1);

};

#endif // _track2_


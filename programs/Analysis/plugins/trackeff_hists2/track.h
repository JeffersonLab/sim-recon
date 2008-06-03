// $Id$
//
//    File: DEventProcessor_trackeff_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _track_
#define _track_

#include <TObject.h>
#include <TVector3.h>


class track:public TObject{

	public:

		TVector3 pthrown;
		TVector3 pfit;
		double likelihood;
		double chisq;
		double pt_pull;
		double theta_pull;
		double phi_pull;
		bool isreconstructable;

	private:
		ClassDef(track,1);

};

#endif // _track_


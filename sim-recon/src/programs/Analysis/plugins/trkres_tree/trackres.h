// $Id$
//
//    File: DEventProcessor_trackeff_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _trackres_
#define _trackres_

#include <TObject.h>
#include <TVector3.h>


class trackres:public TObject{

	public:

		int event;
		TVector3 recon;
		TVector3 thrown;
		double deltak;
		double pt_res; // relative
		double p_res;  // relative
		double theta_res; // absolute
		double phi_res;  // absolute

	private:
		ClassDef(trackres,1);

};

#endif // _trackres_


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

#include "track_info.h"

class track:public TObject{

	public:

		TVector3 pthrown;
		int event;
		int track;
		int Ncdc;
		int Nfdc;
		int mech; // largest value of mech from DMCTrajectoryPoint with R<60.0
		double dtheta_mech;	// angle change due to mech (radians)
		double dp_mech;		// momentum change due to mech (GeV/c)
		
		track_info can;	// Track candidate
		track_info trkwb;	// wire-based fit
		track_info trktb;	// time-based fit

	private:
		ClassDef(track,1);

};

#endif // _track_


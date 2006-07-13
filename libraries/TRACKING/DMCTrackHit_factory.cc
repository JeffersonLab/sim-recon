// $Id$
//
//    File: DMCTrackHit_factory.cc
// Created: Mon Apr  4 08:18:07 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DMCTrackHit_factory.h"
#include "GlueX.h"


//------------------
// evnt
//------------------
jerror_t DMCTrackHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	/// The objects of this class are all produced by the event source.
	/// (See DEventSourceHDDM for example). The only reason we even need
	/// this class here is for the toString() method used by hd_dump.

	return NOERROR;
}


//------------------
// toString
//------------------
const string DMCTrackHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("     id:   r(cm): phi(rad):  z(cm): track: primary:    system:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DMCTrackHit *mctrackhit = _data[i];

		printnewrow();
		
		printcol("%x", mctrackhit->id);
		printcol("%3.1f", mctrackhit->r);
		printcol("%1.3f", mctrackhit->phi);
		printcol("%3.1f", mctrackhit->z);
		printcol("%d", mctrackhit->track);
		printcol(mctrackhit->primary ? "Y":"N");
		printcol(SystemName(mctrackhit->system));
		printrow();
	}

	return _table;
}

// $Id$
//
//    File: DFactory_DTrackHit_MC.cc
// Created: Tue Aug 23 05:29:23 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DFactory_DTrackHit_MC.h"
#include "DMCTrackHit.h"
#include "Dtrkhit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DTrackHit_MC::evnt(DEventLoop *loop, int eventnumber)
{
	/// Here we just copy the data from the DMCTrackHit factory.
	/// Note that we create objects of type Dtrkhit which is derived
	/// from the DTrackHit class. We store them and the factory
	/// presents them as DTrackHit objects. The reason for this is
	/// that the attributes added by Dtrkhit are used internally
	/// by the DTrackCandidate factory and so are normally "hidden".
	/// This saves creating and deleting a whole other set of objects
	/// in the DTrackCandidate factory which hold the same
	/// information as is already in DTrackHit.

	vector<const DMCTrackHit*> dmctrackhits;
	loop->Get(dmctrackhits);
	for(unsigned int i=0;i<dmctrackhits.size(); i++){
		_data.push_back(new Dtrkhit(dmctrackhits[i]));
	}

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DTrackHit_MC::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!
	
	printheader("id:    x:     y:     z:    r:    phi:   system:");

	for(unsigned int i=0; i<_data.size(); i++){
		DTrackHit *trackhit = _data[i];

		printnewrow();
		printcol("%d",	trackhit->id);
		printcol("%1.3f",	trackhit->x);
		printcol("%1.3f",	trackhit->y);
		printcol("%1.3f",	trackhit->z);
		printcol("%1.3f",	trackhit->r);
		printcol("%1.3f",	trackhit->phi);
		printcol("%s",	SystemName(trackhit->system));
		printrow();
	}

	return _table;

}

// $Id$
//
//    File: DTrackHit_factory_MC.cc
// Created: Tue Aug 23 05:29:23 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DTrackHit_factory_MC.h"
#include "DMCTrackHit.h"
#include "Dmctrk_hit.h"

//------------------
// evnt
//------------------
DTrackHit_factory_MC::DTrackHit_factory_MC()
{
	// Set defaults
	EXCLUDE_SECONDARIES = false;
	MAX_HIT_R_FDC = 56.3; // cm

	gPARMS->SetDefaultParameter("TRKMC:EXCLUDE_SECONDARIES",		EXCLUDE_SECONDARIES);
	gPARMS->SetDefaultParameter("TRKMC:MAX_HIT_R_FDC",				MAX_HIT_R_FDC);

}

//------------------
// evnt
//------------------
jerror_t DTrackHit_factory_MC::evnt(JEventLoop *loop, int eventnumber)
{
	/// Here we just copy the data from the DMCTrackHit factory.
	/// Note that we create objects of type Dmctrk_hit which is derived
	/// from the DTrackHit class. We store them and the factory
	/// presents them as DTrackHit objects. The reason for this is
	/// that the attributes added by Dmctrk_hit are used internally
	/// by the DTrackCandidate factory and so are normally "hidden".
	/// This saves creating and deleting a whole other set of objects
	/// in the DTrackCandidate factory which hold the same
	/// information as is already in DTrackHit.
	///
	/// Some filters are also applied here to limit what data is
	/// used by the tracking code. If the TRK:EXCLUDE_SECONDARIES
	/// configuration parameter is set, then hits that do not have
	/// the primary flag set are filtered out. A filter is also
	/// made of FDC hits with an R value greater than
	/// TRK:MAX_HIT_R_FDC.
	///
	/// This factory also records the track number and primary
	/// flag for each of the hits that do pass the filters.
	/// These are kept in private members which are
	/// vectors that can be accessed via the GetTrackNumbers()
	/// and GetPrimaryFlags() routines.

	vector<const DMCTrackHit*> dmctrackhits;
	loop->Get(dmctrackhits);
	
	tracknumber.clear();
	primaryflag.clear();
	for(unsigned int i=0;i<dmctrackhits.size(); i++){
		const DMCTrackHit *dmctrackhit = dmctrackhits[i];
		if(EXCLUDE_SECONDARIES)
			if(!dmctrackhit->primary)continue;
		if(dmctrackhit->system == SYS_FDC)
			if(dmctrackhit->r>MAX_HIT_R_FDC)continue;

		Dmctrk_hit *t = new Dmctrk_hit(dmctrackhit);
		//t->InitCovarianceMatrix();
		_data.push_back(t);
		tracknumber.push_back(dmctrackhit->track);
		primaryflag.push_back(dmctrackhit->primary);
	}

	return NOERROR;
}


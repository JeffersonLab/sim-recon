// $Id$
//
//    File: DTrackEfficiency_factory.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <vector>
#include <map>
using namespace std;

#include "DTrackEfficiency_factory.h"
#include "DTrackCandidate.h"
#include "DMCTrackHit.h"
#include "DTrackHit.h"
#include "DMCThrown.h"
#include "DTrack.h"

class Dthrown_found{
	public:
		oid_t id; 			///< thrown id
		oid_t trackid;	///< found id
		int Nhits_thrown_and_found;
};


//------------------
// evnt
//------------------
jerror_t DTrackEfficiency_factory::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DTrackCandidate*> trackcandidates;
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DTrackHit*> trackhits;
	vector<const DMCThrown*> mcthrowns;
	vector<const DTrack*> tracks;

	string TRACKHIT_SOURCE;
	jparms.GetParameter("TRK:TRACKHIT_SOURCE", TRACKHIT_SOURCE);
	
	JFactory<DTrackCandidate> *factory_tc = loop->Get(trackcandidates);
	JFactory<DMCTrackHit> *factory_mcth = loop->Get(mctrackhits);
	loop->Get(trackhits, TRACKHIT_SOURCE.c_str());
	loop->Get(mcthrowns);
	JFactory<DTrack> *factory_trk = loop->Get(tracks);
	
	// First, we need to loop over all found tracks and figure out
	// which thrown track (if any) they correspond to. We record the
	// thrown track and the number of hits from the thrown track
	// for each so that it can be used in the loop over thrown tracks
	// below.
	//
	// To make it easier below, this is stored in such a way that
	// it can be indexed via thrown id. This is done by
	// using the special Dthrown_found class that includes an
	// oid_t id member so we can use GetByID. (Actually,
	// since it's possible to have more than one found track
	// correponding to the same thrown track, we use GetVectorByID
	// and pick the one with the most hits)
	vector<const Dthrown_found*> thrown_founds;
	for(unsigned int i=0;i<tracks.size();i++){
		const DTrack *track = tracks[i];
		const DTrackCandidate *trackcandidate = factory_tc->GetByIDT(track->candidateid);
		if(!trackcandidate)continue;

		// Loop over hits on this track. The hit ids for DMCTrackHit
		// and DTrackHit:MC are one-to-one so just use the hit ids
		// to reference the DMCTrackHit
		map<oid_t,int> hist;
		//int hist[10] = {0,0,0,0,0,0,0,0,0,0};
		vector<oid_t> hitid = trackcandidate->hitid;
		for(unsigned int k=0; k<hitid.size(); k++){
			const DMCTrackHit *mctrackhit = factory_mcth->GetByIDT(hitid[k]);
			if(!mctrackhit)continue;
			if(!mctrackhit->primary)continue;
			if((mctrackhit->system&&(SYS_CDC|SYS_FDC)) == 0)continue;

			if(hist.find(mctrackhit->track)==hist.end())
				hist[mctrackhit->track]=1;
			else
				hist[mctrackhit->track]++;
		}
		
		// At this point the "hist" map is filled but is indexed by the GEANT
		// track number. We need to find the DANA id of the DMCThrown object
		// corresponding to the track number in hist with the most hits.
		// First, we get the GEANT track number, then we have to assume that
		// the order of the DMCThrown objects is by track number since they
		// do not have the track info in them anywhere.
		int trackno = 0;
		map<oid_t,int>::const_iterator iter;
		for(iter=hist.begin(); iter!=hist.end(); iter++){
			if(trackno){
				if(iter->second > hist[trackno])trackno = iter->first;
			}else{
				trackno = iter->first;
			}
		}
		oid_t thrownid = 0;
		if(mcthrowns.size()>=(unsigned int)trackno)
			if(trackno>0 && trackno<=(int)mcthrowns.size())thrownid = mcthrowns[trackno-1]->id;
		
		// Remember this result
		Dthrown_found *tf = new Dthrown_found;
		thrown_founds.push_back(tf);
		tf->id = thrownid; // thrown track's DANA id
		tf->Nhits_thrown_and_found = hist[trackno];
		tf->trackid = track->id;
	}
	
	// Loop over thrown tracks. 
	for(unsigned int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];
	
		DTrackEfficiency *trkeff = new DTrackEfficiency();
		_data.push_back(trkeff);
		
		// One-to-one correspondance with DMCThrown
		trkeff->id = mcthrown->id;

		// Find total hits from primary thrown track
		trkeff->Nhits_thrown = 0;
		for(unsigned int j=0;j<mctrackhits.size();j++){
			const DMCTrackHit* mctrackhit = mctrackhits[j];
			if(!mctrackhit->primary)continue;
			if((mctrackhit->system&&(SYS_CDC|SYS_FDC)) == 0)continue;

			// The following is tricky. It assumes the track number from
			// HDDM truth tags (which is what DMCTrackHit::track is
			// filled with) corresponds to the position of the thrown
			// track in HDDM. This seems to be true, though I'm not
			// sure if it's guaranteed.
			if(mctrackhit->track == (int)i+1)trkeff->Nhits_thrown++;
		}

		// Find reconstructed track (if any) corresponding to this one
		// There may be several reconstructed tracks corresponding to
		// this. Pick the one with the most hits from this track.
		// (Should we use the largest fraction of hits?)
		vector<const Dthrown_found*> my_thrown_founds;
		for(uint j=0; j<thrown_founds.size(); j++){
			if(thrown_founds[j]->id == mcthrown->id)my_thrown_founds.push_back(thrown_founds[j]);
		}		
		
		const Dthrown_found *my_thrown_found = NULL;
		for(unsigned int j=0; j<my_thrown_founds.size(); j++){
			if(my_thrown_found){
				if(my_thrown_found->Nhits_thrown_and_found < my_thrown_founds[j]->Nhits_thrown_and_found){
					my_thrown_found = my_thrown_founds[j];
				}
			}else{
				my_thrown_found = my_thrown_founds[j];
			}
		}
		
		// Fill in some values. Use defaults if no track found.
		trkeff->Nhits_thrown_and_found = 0;
		trkeff->Nhits_found = 0;
		trkeff->trackid = 0;
		if(my_thrown_found){
			const DTrack *track = factory_trk->GetByIDT(my_thrown_found->trackid); // this HAS to succeed
			const DTrackCandidate *trackcandidate = factory_tc->GetByIDT(track->candidateid);
			trkeff->Nhits_thrown_and_found = my_thrown_found->Nhits_thrown_and_found;
			trkeff->Nhits_found = trackcandidate->hitid.size();
			trkeff->trackid = track->id;
		}
		
		trkeff->Nhits_found_different = trkeff->Nhits_found - trkeff->Nhits_thrown_and_found;
		trkeff->Nhits_thrown_unused = trkeff->Nhits_thrown - trkeff->Nhits_thrown_and_found;
		trkeff->fittable = trkeff->Nhits_thrown >=10;
	}
	
	// Delete the thrown_found objects
	for(unsigned int i=0; i<thrown_founds.size(); i++){
		// have to typecast this back from being const
		delete (Dthrown_found*)thrown_founds[i];
	}

	return NOERROR;
}


//------------------
// toString
//------------------
const string DTrackEfficiency_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("id: Nthrown: Nfound: Nthrown_and_found: frac_found: frac_thrown: fittable: found: trackid:");
	
	for(unsigned int i=0; i<_data.size(); i++){

		DTrackEfficiency *trkeff = _data[i];
		float fraction_from_thrown = (float)trkeff->Nhits_thrown_and_found/(float)trkeff->Nhits_found;
		
		printnewrow();
		
		printcol("%x", trkeff->id);
		printcol("%d", trkeff->Nhits_thrown);
		printcol("%d", trkeff->Nhits_found);
		printcol("%d", trkeff->Nhits_thrown_and_found);
		printcol("%3.0f%%", 100.0*(float)trkeff->Nhits_thrown_and_found/(float)trkeff->Nhits_found);
		printcol("%3.0f%%", 100.0*(float)trkeff->Nhits_thrown_and_found/(float)trkeff->Nhits_thrown);
		printcol("%s", trkeff->fittable ? "Y":"N");
		printcol("%s", fraction_from_thrown>=0.70 ? "Y":"N");
		printcol("%x", trkeff->trackid);

		printrow();
	}

	return _table;
}

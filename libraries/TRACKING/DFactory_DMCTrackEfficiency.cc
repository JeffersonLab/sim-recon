// $Id$
//
//    File: DFactory_DMCTrackEfficiency.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <vector>
using namespace std;

#include "DFactory_DMCTrackEfficiency.h"
#include "DFactory_DMCTrackCandidate.h"
#include "DFactory_DMCCheatHit.h"
#include "DFactory_DMCThrown.h"
#include "DFactory_DMCReconstructed.h"


//------------------
// evnt
//------------------
derror_t DFactory_DMCTrackEfficiency::evnt(DEventLoop *loop, int eventnumber)
{
	vector<const DMCTrackCandidate*> mctrackcandidates;
	vector<const DMCCheatHit*> mccheathits;
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCReconstructed*> mcreconstructeds;
	
	loop->Get(mctrackcandidates, "B");
	loop->Get(mccheathits);
	loop->Get(mcthrowns);
	loop->Get(mcreconstructeds);
	
	// The vector of DMCReconstructed objects should correspond
	// exactly to the vector of DMCTrackCandidate objects.
	if(mctrackcandidates.size() != mcreconstructeds.size()){
		cerr<<__FILE__<<":"<<__LINE__<<" DMCTrackCandidate and DMCReconstructed";
		cerr<<" have different sizes ("<<mctrackcandidates.size()<<" and "<<mcreconstructeds.size()<<")!!!"<<endl;
		return VALUE_OUT_OF_RANGE;
	}

	// Loop over thrown tracks. 
	for(unsigned int i=0;i<mcthrowns.size();i++){
		DMCTrackEfficiency *trkeff = new DMCTrackEfficiency();
		_data.push_back(trkeff);
		
		trkeff->track = i+1; // I'm not sure if this guaranteed to be right

		// Find total hits from primary thrown track
		trkeff->Nhits_thrown = 0;
		for(unsigned int j=0;j<mccheathits.size();j++){
			const DMCCheatHit* mccheathit = mccheathits[j];
			if(!mccheathit->primary)continue;
			if(mccheathit->system>2)continue;
			if(mccheathit->track == trkeff->track)trkeff->Nhits_thrown++;
		}

		// Find reconstructed track (if any) corresponding to this one
		// There may be several reconstructed tracks corresponding to
		// this. Pick the one with the most hits from this track.
		// (Should we use the largest fraction of hits?)
		trkeff->Nhits_thrown_and_found = 0;
		trkeff->Nhits_found = 0;
		trkeff->index_DMCReconstructed = -1;
		int max_Ntrackhits = 0;
		for(unsigned int j=0; j<mctrackcandidates.size(); j++){
			const DMCTrackCandidate *mctc = mctrackcandidates[j];
			if(mctc->track == (int)i+1){
				if(mctc->Ntrackhits > max_Ntrackhits){
					max_Ntrackhits = mctc->Ntrackhits;
					trkeff->Nhits_thrown_and_found = max_Ntrackhits;
					trkeff->Nhits_found = mctc->Nhits;
					trkeff->index_DMCReconstructed = j;
				}
			}
		}

		trkeff->Nhits_found_different = trkeff->Nhits_found - trkeff->Nhits_thrown_and_found;
		trkeff->Nhits_thrown_unused = trkeff->Nhits_thrown - trkeff->Nhits_thrown_and_found;
		trkeff->fittable = trkeff->Nhits_thrown >=10;
	}
	
	return NOERROR;
}


//------------------
// toString
//------------------
const string DFactory_DMCTrackEfficiency::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: Nthrown: Nfound: Nthrown_and_found: fraction_found: fittable:");
	
	for(unsigned int i=0; i<_data.size(); i++){

		DMCTrackEfficiency *trkeff = _data[i];
		
		printnewrow();
		
		printcol("%d", i);
		printcol("%d", trkeff->Nhits_thrown);
		printcol("%d", trkeff->Nhits_found);
		printcol("%d", trkeff->Nhits_thrown_and_found);
		printcol("%3.0f%%", 100.0*(float)trkeff->Nhits_thrown_and_found/(float)trkeff->Nhits_thrown);
		printcol("%s", trkeff->fittable ? "Y":"N");

		printrow();
	}

	return _table;
}

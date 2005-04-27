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
derror_t DFactory_DMCTrackEfficiency::evnt(int eventnumber)
{
	vector<const DMCTrackCandidate*> mctrackcandidates;
	vector<const DMCCheatHit*> mccheathits;
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCReconstructed*> mcreconstructeds;
	
	event->Get(mctrackcandidates);
	event->Get(mccheathits);
	event->Get(mcthrowns);
	event->Get(mcreconstructeds);
	
	// The vector of DMCReconstructed objects should correspond
	// exactly to the vector of DMCTrackCandidate objects.
	if(mctrackcandidates.size() != mcreconstructeds.size()){
		cerr<<__FILE__<<":"<<__LINE__<<" DMCTrackCandidate and DMCReconstructed";
		cerr<<" have different sizes!!!"<<endl;
		return VALUE_OUT_OF_RANGE;
	}

	// Loop over thrown tracks. 
	for(int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		DMCTrackEfficiency *trkeff = new DMCTrackEfficiency();
		_data.push_back(trkeff);
		
		trkeff->track = i+1; // I'm not sure if this guaranteed to be right

		// Find total hits from primary thrown track
		trkeff->Nhits_thrown = 0;
		for(int j=0;j<mccheathits.size();j++){
			const DMCCheatHit* mccheathit = mccheathits[j];
			if(!mccheathit->primary)continue;
			if(mccheathit->track == trkeff->track)trkeff->Nhits_thrown++;
		}

		// Find reconstructed track (if any) corresponding to this one
		trkeff->index_DMCReconstructed = -1;
		for(int j=0; j<mcreconstructeds.size(); j++){
			if(mcreconstructeds[j]->thrownid == i)trkeff->index_DMCReconstructed = j;
		}
		
		
		// Loop over reconstructed tracks hits (if one was found)
		trkeff->Nhits_thrown_and_found = 0;
		trkeff->Nhits_found = 0;
		if(trkeff->index_DMCReconstructed >= 0){
			const DMCTrackCandidate *mctc = mctrackcandidates[trkeff->index_DMCReconstructed];

			trkeff->Nhits_found = mctc->Nhits;
			for(int j=0;j<mctc->Nhits;j++){
				const DMCCheatHit *mccheathit = mccheathits[mctc->ihit[j]];
				
				// Only consider primary track hits for now
				if(!mccheathit->primary)continue;
				if(mccheathit->track == trkeff->track)trkeff->Nhits_thrown_and_found++;
			}
		}

		trkeff->Nhits_found_different = trkeff->Nhits_found - trkeff->Nhits_thrown_and_found;
		trkeff->Nhits_thrown_unused = trkeff->Nhits_thrown - trkeff->Nhits_thrown_and_found;
		
	}
	
	return NOERROR;
}


//------------------
// toString
//------------------
const string DFactory_DMCTrackEfficiency::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: Nthrown: Nfound: Nthrown_and_found: fraction_found: fittable:");
	
	for(int i=0; i<_data.size(); i++){

		DMCTrackEfficiency *trkeff = _data[i];
		
		printnewrow();
		
		printcol("%d", i);
		printcol("%d", trkeff->Nhits_thrown);
		printcol("%d", trkeff->Nhits_found);
		printcol("%d", trkeff->Nhits_thrown_and_found);
		printcol("%3.0f%%", 100.0*(float)trkeff->Nhits_thrown_and_found/(float)trkeff->Nhits_thrown);
		printcol("%s", trkeff->Nhits_thrown>3 ? "Y":"N");

		printrow();
	}

	return _table;
}

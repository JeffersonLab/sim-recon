// $Id$
//
//    File: DChargedTruthMatch.h
// Created: Sun Jan 31 08:45:38 EST 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DChargedTruthMatch_
#define _DChargedTruthMatch_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include <PID/DKinematicData.h>

class DMCThrown;
class DChargedTrack;

class DChargedTruthMatch:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DChargedTruthMatch);

		int Nhits_thrown;					// Number of hits on truth track
		int Nhits_thrown_selector;		// Number of hits selector found on truth track
		int Nhits_recon;					// Number of hits selector found on reconstructed track
		int Nhits_both;					// Number of hits in both the Nhits_thrown_selector and Nhits_recon sets
		const DMCThrown* thrown;		// Pointer to thrown track
		const DChargedTrack* recon;	// Pointer to reconstructed track
		double fom;							// Figure of merit for match

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "Nhits_thrown", "%d", Nhits_thrown);
			AddString(items, "Nhits_thrown_selector", "%d", Nhits_thrown_selector);
			AddString(items, "Nhits_recon", "%d", Nhits_recon);
			AddString(items, "Nhits_both", "%d", Nhits_both);
			AddString(items, "fom", "%f", fom);
		}
};

#endif // _DChargedTruthMatch_


// $Id$
//
//    File: DMCFitStats.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCFitStats_
#define _DMCFitStats_

#include <TH1.h>
#include <TH2.h>

#include "DFactory.h"

class DEvent;

class DMCFitStats{
	public:
		HDCLASSDEF(DMCFitStats);
		
		DMCFitStats();
		~DMCFitStats();
		void AddEvent(DEvent *event);
		void DMCFitStats::Finalize(void);
		
		
		TH1F *stats, *frac, *h4_dist, *h4_dist_primary;
		TH2F *delta_p;
		TH1F *delta_p_over_p;
		TH1F *hits_per_thrown_track; 

		/// The "stats" histogram has a special meaning for each bin:
		///
		///  1 - Number of hits in found track that were actually from same track
		///  2 - Number of hits in found track that were actually from different track
		///  3 - Number of hits in actual track which corresponds to this found one
		///  4 - Same as above, but which also were included in the found track
		///  5 - Number of hits not included in any found track
		///  6 - Total number of cheat hits
		///  7 - Total number of tracks thrown (FDC+CDC had more than 3 hits)
		///  8 - Total number of tracks found
		///
		///
		/// The "frac" histogram also has special meanings for each bin:
		///
		///  1 - Ratio of found to "thrown" tracks (thrown means >3 hits in CDC+FDC)
		///  2 - Fraction of cheat hits used in at least one found track
		///  3 - Fraction of hits assigned to track they actually came from
		///  4 - Fraction of hits in thrown track assigned to found track
		
		ClassDef(DMCFitStats,1)
};

#endif // _DMCFitStats_


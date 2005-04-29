// $Id$
//
//    File: DMCTrackEfficiency.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCTrackEfficiency_
#define _DMCTrackEfficiency_

#include <Rtypes.h>

#include "DFactory.h"

class DMCTrackEfficiency{
	public:
		HDCLASSDEF(DMCTrackEfficiency);

		/// Objects in this class correspond 1 to 1 with those in DMCThrown
		int Nhits_thrown;					///< Num hits in thrown (primary) track
		int Nhits_found;					///< Num hits in reconstructed (primary) track
		int Nhits_thrown_and_found;	///< Num hits in both thrown and found (primary) tracks
		int Nhits_found_different;		///< Nfound - Nhits_thrown_and_found
		int Nhits_thrown_unused;		///< Nthrown - Nhits_thrown_and_found
		int fittable;						///< non-zero if more than 3 pimary hits

		int track;							///< track number from GEANT (same as DMCCheatHits)
		int index_DMCReconstructed;	///< index to corresponding DMCReconstructed
		
		ClassDef(DMCTrackEfficiency,1)
};

#endif // _DMCTrackEfficiency_


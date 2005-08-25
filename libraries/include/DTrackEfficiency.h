// $Id$
//
//    File: DTrackEfficiency.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackEfficiency_
#define _DTrackEfficiency_

#include "DObject.h"
#include "DFactory.h"

class DTrackEfficiency:public DObject{
	public:
		HDCLASSDEF(DTrackEfficiency);

		/// Objects in this class correspond 1 to 1 with those in DThrown
		int Nhits_thrown;					///< Num hits in thrown (primary) track
		int Nhits_found;					///< Num hits in reconstructed (primary) track
		int Nhits_thrown_and_found;	///< Num hits in both thrown and found (primary) tracks
		int Nhits_found_different;		///< Nfound - Nhits_thrown_and_found
		int Nhits_thrown_unused;		///< Nthrown - Nhits_thrown_and_found
		int fittable;						///< non-zero if more than 3 pimary hits

		identifier_t trackid;			///< id of DTrack (if any)
};

#endif // _DTrackEfficiency_


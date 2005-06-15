// $Id$
//
//    File: DMCTrackCandidate.h
// Created: Sun Apr  3 12:38:16 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCTrackCandidate_
#define _DMCTrackCandidate_

#include "DFactory.h"

#define MAX_IHITS 128

class DMCTrackCandidate{
	public:
		HDCLASSDEF(DMCTrackCandidate);
		
		int Nhits;
		int ihit[MAX_IHITS];	///< index of hits in MCCheatHits factory
		float x0,y0;			///< center of circle
		float z_vertex;		///< z coordinate of vertex
		float dphidz;			///< dphi/dz in radians per cm
		float q;					///< electric charge 
		float p, p_trans;		///< total and transverse momenta in GeV/c
		float phi, theta;		///< theta and phi in radians
};

#endif // _DMCTrackCandidate_


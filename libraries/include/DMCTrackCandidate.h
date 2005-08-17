// $Id$
//
//    File: DMCTrackCandidate.h
// Created: Sun Apr  3 12:38:16 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCTrackCandidate_
#define _DMCTrackCandidate_

#include "DObject.h"
#include "DFactory.h"

#define MAX_IHITS 256

class DMCTrackCandidate:public DObject{
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
		int track;				///< track from which most hits came from
		int Ntrackhits;		///< number of hits in this candidate that came from "track"
};

#endif // _DMCTrackCandidate_


// $Id$
//
//    File: DTrackCandidate.h
// Created: Sun Apr  3 12:38:16 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_
#define _DTrackCandidate_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

#define MAX_IHITS 256

class DTrackCandidate:public JObject{
	public:
		HDCLASSDEF(DTrackCandidate);
		
		vector<oid_t> hitid;	///< ids of hits in DTrackHit factory
		float x0,y0;			///< center of circle
		float z_vertex;		///< z coordinate of vertex
		float dzdphi;			///< dz/dphi in cm per radian
		float q;					///< electric charge 
		float p, p_trans;		///< total and transverse momenta in GeV/c
		float phi, theta;		///< theta and phi in radians
};

#endif // _DTrackCandidate_


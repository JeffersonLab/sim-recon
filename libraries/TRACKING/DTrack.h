// $Id$
//
//    File: DTrack.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrack_
#define _DTrack_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DReferenceTrajectory;

class DTrack:public JObject{
	public:
		HDCLASSDEF(DTrack);
		
		float q;				///< electric charge
		float p;				///< Total momentum in GeV/c
		float theta,phi;	///< Inital theta and phi angles in radians
		float x,y,z;		///< Vertex position in cm
		oid_t candidateid;	///< id of DTrackCandidate this came from
		float chisq;			///< reduced Chi-squared for the track
		const DReferenceTrajectory *rt; ///< pointer to reference trjectory representing this track
};

#endif // _DTrack_


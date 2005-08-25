// $Id$
//
//    File: DTrack.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrack_
#define _DTrack_

#include "DObject.h"
#include "DFactory.h"

class DMCThrown;

class DTrack:public DObject{
	public:
		HDCLASSDEF(DTrack);
		
		float q;				///< electric charge
		float p;				///< Total momentum in GeV/c
		float theta,phi;	///< Inital theta and phi angles in radians
		float x,y,z;		///< Vertex position in cm
		identifier_t candidateid;	///< id of DTrackCandidate this came from
};

#endif // _DTrack_


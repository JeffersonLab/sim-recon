// $Id$
//
//    File: DMCThrown.h
// Created: Sun Apr  3 12:22:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCThrown_
#define _DMCThrown_

#include "DFactory.h"

class DMCThrown{
	public:
		HDCLASSDEF(DMCThrown);
		
		int type;			///< GEANT particle ID
		float q;				///< electric charge
		float p;				///< Total momentum in GeV/c
		float E;				///< Total energy in GeV
		float theta,phi;	///< Inital theta and phi angles in radians
		float x,y,z;		///< Vertex position in cm
		float mass;			///< Mass in GeV/c^2
};

#endif // _DMCThrown_


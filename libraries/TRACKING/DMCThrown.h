// $Id$
//
//    File: DMCThrown.h
// Created: Sun Apr  3 12:22:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCThrown_
#define _DMCThrown_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

#include "PID/DKinematicData.h"

class DMCThrown:public DKinematicData{
	public:
		HDCLASSDEF(DMCThrown);
		
		int type;			///< GEANT particle ID
		int pdgtype;		///< PDG particle type (not used by GEANT)
		int myid;			///< id of this particle from original generator
		int parentid;		///< id of parent of this particle from original generator
		int mech;			///< production mechanism of this partcle (generator specific)
		float q;				///< electric charge
		float p;				///< Total momentum in GeV/c
		float E;				///< Total energy in GeV
		float theta,phi;	///< Inital theta and phi angles in radians
		float x,y,z;		///< Vertex position in cm
		float mass;			///< Mass in GeV/c^2

};

#endif // _DMCThrown_


// $Id$
//
//    File: DMCReconstructed.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCReconstructed_
#define _DMCReconstructed_

#include "DFactory.h"

class DMCThrown;

class DMCReconstructed{
	public:
		HDCLASSDEF(DMCReconstructed);
		
		int type;			///< GEANT particle ID
		float q;				///< electric charge
		float p;				///< Total momentum in GeV/c
		float E;				///< Total energy in GeV
		float theta,phi;	///< Inital theta and phi angles in radians
		float x,y,z;		///< Vertex position in cm
		float mass;			///< Mass in GeV/c^2
		int thrownid;		///< index to closest match in DMCThrown
		float thrown_delta_p;	///< Magnitude of momentum diff. with thrownid
		
		void FindClosestThrown(vector<const DMCThrown*> &mcthrowns);
};

#endif // _DMCReconstructed_


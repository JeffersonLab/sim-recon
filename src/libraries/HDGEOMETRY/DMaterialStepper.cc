// $Id$
//
//    File: DMaterialStepper.cc
// Created: Tue May  6 16:28:38 EDT 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.11.1 i386)
//

#include "DGeometry.h"
#include "DMaterialStepper.h"

//---------------------------------
// DMaterialStepper    (Constructor)
//---------------------------------
DMaterialStepper::DMaterialStepper(const DGeometry *dgeom)
{
	this->dgeom = dgeom;
	this->bfield = dgeom->GetBfield();

	stepper = NULL;
}

//---------------------------------
// DMaterialStepper    (Constructor)
//---------------------------------
DMaterialStepper::DMaterialStepper(const DGeometry *dgeom, double q, const DVector3 &pos, const DVector3 &mom)
{
	this->dgeom = dgeom;
	this->q = q;
	this->pos = pos;
	this->mom = mom;

	this->bfield = dgeom->GetBfield();
	stepper = NULL;
}

//---------------------------------
// ~DMaterialStepper    (Destructor)
//---------------------------------
DMaterialStepper::~DMaterialStepper()
{
	if(stepper)delete stepper;
}

//---------------------------------
// GetTraversedMaterialsZ
//---------------------------------
void DMaterialStepper::GetTraversedMaterialsZ(double z_end, vector<DMaterialStep> &materialsteps)
{
	if(!stepper)
		stepper = new DMagneticFieldStepper(bfield, q, &pos, &mom);
	else
		stepper->SetStartingParams(q, &pos, &mom);


	// The method here is to figure out which detector system we are in
	// and then call that system's SwimToEdge method to find the point
	// where the track would exit that system. We also try swimming up to
	// the given z-plane to see if that is closer. If the edge of the
	// detector system is closer, we add a DMaterialStep and repeat the
	// whole process. If the z-plane is closer, then we add the last
	// DMaterialStep.
	
	for(int i=0; i<100; i++){ // limit the number of material steps

		// Find the detector system we're in and the point at 
		// which this track would exit it.
		DVector3 pos_detector = pos;
		DVector3 mom_detector = mom;
		//double s_detector;
		//const DMaterial *mat=NULL;
		switch(WhereAmI()){
			case SYS_CDC:
				//SwimToEdgeCDC(pos_detector, mom_detector, &s_detector);
				//mat = dgeom->GetDMaterial("CDchamberGas");
				break;
			default:
				_DBG__;
		}
	}

}

//---------------------------------
// WhereAmI
//---------------------------------
DetectorSystem_t WhereAmI(void)
{
	/// Determine which detector system we are currently in according to
	/// our position vector.
	
	return SYS_NULL;
}




// $Id$
//
//    File: DMaterialStepper.h
// Created: Tue May  6 16:28:38 EDT 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.11.1 i386)
//

#ifndef _DMaterialStepper_
#define _DMaterialStepper_

#include <vector>
using std::vector;

#include <GlueX.h>
#include <HDGEOMETRY/DMaterialStep.h>

#include <TRACKING/DMagneticFieldStepper.h>
#include <JANA/jerror.h>

class DGeometry;

class DMaterialStepper{
	public:
		DMaterialStepper(const DGeometry *dgeom);
		DMaterialStepper(const DGeometry *dgeom, double q, const DVector3 &pos, const DVector3 &mom);
		virtual ~DMaterialStepper();
		
		const DVector3& GetPosition(void){return pos;}
		const DVector3& GetMomentum(void){return mom;}
		double GetCharge(void){return q;}
		const DGeometry* GetDGeometry(void){return dgeom;}
		
		void SetPosition(DVector3 pos){this->pos = pos;}
		void SetMomentum(DVector3 mom){this->mom = mom;}
		void SetCharge(double q){this->q = q;}

		void GetTraversedMaterialsZ(double z_end, vector<DMaterialStep> &materialsteps);
		DetectorSystem_t WhereAmI(void);
		
	protected:
		double q;
		DVector3 pos;
		DVector3 mom;
		const DGeometry *dgeom;
		const DMagneticFieldMap *bfield;
		DMagneticFieldStepper *stepper;
	
	private:
		DMaterialStepper(){} // Force instantiater to provide us with a DGeometry pointer

};

#endif // _DMaterialStepper_


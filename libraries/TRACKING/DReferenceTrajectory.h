// $Id$
//
//    File: DReferenceTrajectory.h
// Created: Wed Jul 19 13:42:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#ifndef _DReferenceTrajectory_
#define _DReferenceTrajectory_

#include <vector>
using std::vector;

#include <DVector3.h>

#include <JANA/jerror.h>

#include <DCoordinateSystem.h>

class DMagneticFieldMap;
class DTrackCandidate;

class DReferenceTrajectory{
	
	/// This class is a utility class used by the TRACKING package. It
	/// is used to swim a particle through the (inhomogeneous) magnetic
	/// field, accounting for energy loss in the geometry, and recording
	/// each step so that derivatives and distances can be calculated.
	/// Because this uses the coordinates defined in the reference 
	/// trajectory method (B. Mecking Nucl. Instr. and Methods. 
	/// 203 (1982) 299-305), the angles needed to rotate into the
	/// lab frame are saved as well. This used the DMagneticFieldStepper
	/// class for swimming through the field.

	public:


		class swim_step_t:public DCoordinateSystem{
			public:
				DVector3 mom;
				double Ro;
				double s; // distance along RT
		};

		DReferenceTrajectory(const DMagneticFieldMap *
									, double q=1.0
									, swim_step_t *swim_steps=NULL
									, int max_swim_steps=0
									, double step_size=-1.0);
		DReferenceTrajectory(const DReferenceTrajectory& rt);
		DReferenceTrajectory& operator=(const DReferenceTrajectory& rt);
		virtual ~DReferenceTrajectory();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DReferenceTrajectory";}
		
		double DistToRT(double x, double y, double z){return DistToRT(DVector3(x,y,z));}
		double DistToRT(DVector3 hit);
		double DistToRT(const DCoordinateSystem *wire, double *s=NULL);
		double DistToRTBruteForce(const DCoordinateSystem *wire, double *s=NULL);
		double DistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double *s=NULL);
		double DistToRTBruteForce(const DCoordinateSystem *wire, const swim_step_t *step, double *s=NULL);
		double Straw_dx(const DCoordinateSystem *wire, double radius);
		swim_step_t* FindClosestSwimStep(const DCoordinateSystem *wire);
		void Swim(const DVector3 &pos, const DVector3 &mom, double q=-1000.0);
		DVector3 GetLastDOCAPoint(void);
		void GetLastDOCAPoint(DVector3 &pos, DVector3 &mom);
		double GetLastDistAlongWire(void){return last_dist_along_wire;}
		void SetStepSize(double step_size){this->step_size=step_size;}

		const swim_step_t *GetLastSwimStep(void){return last_swim_step;}
		swim_step_t *swim_steps;
		int Nswim_steps;
		float q;

	protected:
	
		int max_swim_steps;
		bool own_swim_steps;
		double step_size;
		const DMagneticFieldMap *bfield;
		
		double last_phi;							///< last phi found in DistToRT
		const swim_step_t* last_swim_step;	///< last swim step used in DistToRT
		double last_dist_along_wire;
		double last_dz_dphi;
	
	private:
		DReferenceTrajectory(){} // force use of constructor with arguments.

};

#endif // _DReferenceTrajectory_


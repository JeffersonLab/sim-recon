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

#include <TVector3.h>

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
				TVector3 mom;
				double Ro;
				double s; // distance along RT
		};

		DReferenceTrajectory(const DMagneticFieldMap *
									, double q, const TVector3 &pos, const TVector3 &mom
									, swim_step_t *swim_steps=NULL
									, int max_swim_steps=0
									, double step_size=-1.0);

		virtual ~DReferenceTrajectory();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DReferenceTrajectory";}
		
		double DistToRT(double x, double y, double z){return DistToRT(TVector3(x,y,z));}
		double DistToRT(TVector3 hit);
		double DistToRT(const DCoordinateSystem *wire, double L, double *s=NULL);
		double DistToRTBruteForce(const DCoordinateSystem *wire, double L, double *s=NULL);
		double DistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double *s=NULL);
		double DistToRTBruteForce(const DCoordinateSystem *wire, const swim_step_t *step, double *s=NULL);
		swim_step_t* FindClosestSwimStep(const DCoordinateSystem *wire, double L);
		void Reswim(const TVector3 &pos, const TVector3 &mom);

		swim_step_t *swim_steps;
		int Nswim_steps;
		float q;

	protected:
	
		int max_swim_steps;
		bool own_swim_steps;
		double step_size;
		const DMagneticFieldMap *bfield;
	
	private:
		DReferenceTrajectory(){} // force use of constructor with arguments.

};

#endif // _DReferenceTrajectory_


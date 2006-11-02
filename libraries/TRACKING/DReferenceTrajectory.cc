// $Id$
//
//    File: DReferenceTrajectory.cc
// Created: Wed Jul 19 13:42:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#include <TVector3.h>

#include "DReferenceTrajectory.h"
#include "DTrackCandidate.h"
#include "DMagneticFieldStepper.h"

//---------------------------------
// DReferenceTrajectory    (Constructor)
//---------------------------------
DReferenceTrajectory::DReferenceTrajectory(const DMagneticFieldMap *bfield
														, const DTrackCandidate *tc
														, swim_step_t *swim_steps
														, int max_swim_steps)
{
	// Assume we start on the beamline
	TVector3 p;
	p.SetMagThetaPhi(tc->p, tc->theta, tc->phi);
	TVector3 pos(0.0, 0.0, tc->z_vertex);

	// Real work is done in SwimRT
	SwimRT(bfield, tc->q, pos, p, swim_steps, max_swim_steps);
}

//---------------------------------
// DReferenceTrajectory    (Constructor)
//---------------------------------
DReferenceTrajectory::DReferenceTrajectory(const DMagneticFieldMap *bfield
														, double q, const TVector3 &pos, const TVector3 &p
														, swim_step_t *swim_steps
														, int max_swim_steps)
{
	SwimRT(bfield, q, pos, p, swim_steps, max_swim_steps);
}

//---------------------------------
// SwimRT
//---------------------------------
void DReferenceTrajectory::SwimRT(const DMagneticFieldMap *bfield
														, double q, const TVector3 &pos, const TVector3 &p
														, swim_step_t *swim_steps
														, int max_swim_steps)
{

	// It turns out that the greatest bottleneck in speed here comes from
	// allocating/deallocating the large block of memory required to hold
	// all of the trajectory info. The preferred way of calling this is 
	// with a pointer allocated once at program startup. This code block
	// though allows it to be allocated here if necessary.
	if(!swim_steps){
		own_swim_steps = true;
		this->max_swim_steps = 50000;
		this->swim_steps = new swim_step_t[this->max_swim_steps];
	}else{
		own_swim_steps = false;
		this->max_swim_steps = max_swim_steps;
		this->swim_steps = swim_steps;
	}

	// Initialize stepper.
	this->q = q;
	DMagneticFieldStepper stepper(bfield, q, &pos, &p);
		
	// Step until we hit a boundary (don't track more than 20 meters)
	swim_step_t *swim_step = this->swim_steps;
	Nswim_steps = 0;
	for(double s=0; fabs(s)<2000.0; Nswim_steps++, swim_step++){

		if(Nswim_steps>=this->max_swim_steps){
			cerr<<__FILE__<<":"<<__LINE__<<" Too many steps in trajectory. Truncating..."<<endl;
			break;
		}

		stepper.GetDirs(swim_step->xdir, swim_step->ydir, swim_step->zdir);
		stepper.GetPosition(swim_step->pos);
		stepper.GetMomentum(swim_step->mom);
		swim_step->Ro = stepper.GetRo();
		swim_step->s = s;
//cout<<__FILE__<<":"<<__LINE__<<" s="<<s<<" ";swim_step->pos.Print();

		// Exit loop if we leave the tracking volume
		if(swim_step->pos.Perp()>65.0){break;} // ran into BCAL
		if(swim_step->pos.Z()>650.0){break;} // ran into FCAL
		if(swim_step->pos.Z()<-50.0){break;} // ran into UPV

		// Swim to next
		s += stepper.Step(NULL);
	}

	// OK. At this point the positions of the trajectory in the lab
	// frame have been recorded along with the momentum of the
	// particle and the directions of reference trajectory
	// coordinate system at each point.
}

//---------------------------------
// ~DReferenceTrajectory    (Destructor)
//---------------------------------
DReferenceTrajectory::~DReferenceTrajectory()
{
	if(own_swim_steps){
		delete[] swim_steps;
	}
}

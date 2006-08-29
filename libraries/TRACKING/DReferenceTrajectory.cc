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

	// It turns out that the greatest bottleneck in speed here comes from
	// allocating/deallocating the large block of memory required to hold
	// all of the trajectory info. The preferred way of calling this is 
	// with a pointer allocated once at program startup. This code block
	// though allows it to be allocated here if necessary.
	if(!swim_steps){
		own_swim_steps = true;
		this->max_swim_steps = 2000;
		this->swim_steps = new swim_step_t[this->max_swim_steps];
	}else{
		own_swim_steps = false;
		this->max_swim_steps = max_swim_steps;
		this->swim_steps = swim_steps;
	}

	// Initialize stepper with track candidate values. This always
	// assumes we start on the beamline.
	q= tc->q;
	TVector3 p;
	p.SetMagThetaPhi(tc->p, tc->theta, tc->phi);
	TVector3 pos(0.0, 0.0, tc->z_vertex);
	DMagneticFieldStepper stepper(bfield, q, &pos, &p);
	
	// Set step size to 2cm. The stepper currently uses a constant step size.
	double step_size = 2.0;
	stepper.SetStepSize(step_size);
		
	// First step is starting position
	swim_step_t *swim_step = &swim_steps[Nswim_steps++];
	stepper.GetDirs(swim_step->xdir, swim_step->ydir, swim_step->zdir);
	stepper.GetMomentum(swim_step->mom);

	// Step until we hit a boundary
	int max_steps = (int)(2000.0/step_size); // don't track more than 20 meters along path
	Nswim_steps = 0;
	for(int i=0; i<max_steps; i++){
		//swim_step_t swim_step;
		if(Nswim_steps>max_swim_steps){
			cerr<<__FILE__<<":"<<__LINE__<<" Too many steps in trajectory. Truncating..."<<endl;
			break;
		}
		swim_step_t *swim_step = &swim_steps[Nswim_steps++];
		stepper.Step(&swim_step->pos);
		stepper.GetDirs(swim_step->xdir, swim_step->ydir, swim_step->zdir);
		stepper.GetMomentum(swim_step->mom);
		swim_step->Ro = stepper.GetRo();
		
		// Exit loop if we leave the tracking volume
		if(swim_step->pos.Perp()>65.0){break;} // ran into BCAL
		if(swim_step->pos.Z()>650.0){break;} // ran into FCAL
		if(swim_step->pos.Z()<-50.0){break;} // ran into UPV
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

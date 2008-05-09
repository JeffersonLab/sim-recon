#include <math.h>

#include <DVector3.h>
#include <DMatrix.h>

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DANA/DApplication.h"
#include "DMagneticFieldStepper.h"
#include "DTrackCandidate.h"
#include "DTrack_factory_ALT3.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include "DKalmanFilter.h"
#include "DReferenceTrajectory.h"

//------------------
// DTrack_factory_ALT3   (Constructor)
//------------------
DTrack_factory_ALT3::DTrack_factory_ALT3(){

  TOF_MASS = 0.13957018;
  MIN_FDC_HIT_PROB=0.1;
  MIN_HITS=1;
}

//------------------
// DTrack_factory_ALT3   (Destructor)
//------------------
DTrack_factory_ALT3::~DTrack_factory_ALT3()
{

  for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
}



//------------------
jerror_t DTrack_factory_ALT3::init(void){return NOERROR;}
jerror_t DTrack_factory_ALT3::erun(void){return NOERROR;}
jerror_t DTrack_factory_ALT3::fini(void){return NOERROR;}

//------------------
// brun
//------------------
jerror_t DTrack_factory_ALT3::brun(JEventLoop *loop, int runnumber)
{
  // Get pointer to DMagneticFieldMap field object
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();
  dgeom  = dapp->GetDGeometry(runnumber);

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory_ALT3::evnt(JEventLoop *loop, int eventnumber)
{
  // Store current event number
  this->eventnumber = eventnumber;
  
  // Get input data
  trackcandidates.clear();
  cdctrackhits.clear();
  fdctrackhits.clear();
  loop->Get(trackcandidates);
  loop->Get(cdctrackhits);
  loop->Get(fdctrackhits,"CORRECTED");	

  // Allocate more DReferenceTrajectory objects if needed.
  // These each have a large enough memory footprint that
  // it causes noticable performance problems if we allocated
  // and deallocated them every event. Therefore, we allocate
  // when needed, but recycle them on the next event.
  // They are deleted in the factory deconstructor.
  while(rtv.size()<trackcandidates.size())rtv.push_back(new DReferenceTrajectory(bfield));
  
  // Loop over track candidates
  for(unsigned int i=0; i<trackcandidates.size(); i++){ 
    const DTrackCandidate *tc = trackcandidates[i];
    DReferenceTrajectory *rt = rtv[i];
    double chisq;

     // Find reference trajectory by swimming through the field
    DVector3 pos = tc->position();
    DVector3 mom = tc->momentum(); 	
    rt->Swim(pos, mom, tc->charge());

    // Initialize Kalman filter with B-field
    DKalmanFilter fit(bfield,dgeom);

    // Add hits to be put through Kalman engine
    unsigned int num_matched_hits=0;
    DVector3 last_pos,last_mom;  
    for(unsigned int j=0; j<fdctrackhits.size(); j++){
      const DFDCPseudo *hit = fdctrackhits[j];
      double variance=1.0; //guess for now

      // Find residual 
      pos(0)=hit->x;
      pos(1)=hit->y;
      pos(2)=hit->wire->origin(2);
      double resi=rt->DistToRT(pos);
    
      // Use an un-normalized gaussian so that for a residual
      // of zero, we get a probability of 1.0.
      double p = finite(resi) ? exp(-resi*resi/2./variance):0.0;
      if(p>=MIN_FDC_HIT_PROB){
	fit.AddHit(hit->x,hit->y,hit->wire->origin(2),hit->covxx,hit->covxy,
		   hit->covyy);
	const swim_step_t *last_step=rt->GetLastSwimStep();
     
	if (last_step!=NULL){
	  last_pos=last_step->origin;
	  last_mom=last_step->mom;
	}
	num_matched_hits++;
      }
    }
   
    if (num_matched_hits>=MIN_HITS){
    // Set the initial parameters from the track candidate
      fit.SetSeed(tc->charge(),last_pos,last_mom);

      // Kalman filter 
      jerror_t error=fit.KalmanLoop(TOF_MASS,chisq);

      if (error==NOERROR){
	// Create a new track object
	DTrack *track = new DTrack;
	track->q=tc->charge();
	
	DVector3 mom,pos;
	fit.GetMomentum(mom);
	fit.GetPosition(pos);
	
	track->x=pos(0);
	track->y=pos(1);
	track->z=pos(2);
	track->p=mom.Mag();
	track->phi=mom.Phi();
	if(track->phi<0.0)track->phi+=2.0*M_PI;
	track->theta=mom.Theta();
	track->chisq=fit.GetChiSq();
	track->candidateid=tc->id;
	track->rt=rt;
	
	// Fill in DKinematicData part
	track->setMass(0.0);
	track->setMomentum(mom);
	track->setPosition(pos);
	track->setCharge(track->q);

	_data.push_back(track);
      }
    }
  }
  

  return NOERROR;
}


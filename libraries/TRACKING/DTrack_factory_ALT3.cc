#include <math.h>

#include <DVector3.h>
#include <DMatrix.h>

#include <JANA/JEventLoop.h>

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
  MIN_HITS=3;
}

//------------------
// DTrack_factory_ALT3   (Destructor)
//------------------
DTrack_factory_ALT3::~DTrack_factory_ALT3(){};

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
  
  // Loop over track candidates
  for(unsigned int i=0; i<trackcandidates.size(); i++){ 
    const DTrackCandidate *tc = trackcandidates[i];
    double chisq;

    // Initialize Kalman filter with B-field
    DKalmanFilter fit(bfield);

    // Initialize the stepper 
    DMagneticFieldStepper stepper(bfield, tc->charge()); 

    // Add hits to be put through Kalman engine
    unsigned int num_matched_hits=0;
    DVector3 last_pos,last_mom;
    for(unsigned int j=0; j<fdctrackhits.size(); j++){
      const DFDCPseudo *hit = fdctrackhits[j];
      double variance=5.0; //guess for now
      DVector3 norm(0,0,1), origin(0,0,hit->wire->origin(2));
 
      // Swim from target position to measurement plane and compute residual
      DVector3 pos = tc->position();
      DVector3 mom = tc->momentum(); 
      stepper.SwimToPlane(pos,mom,origin,norm,NULL); 
   
      // compute residual
      double resi=sqrt((pos(0)-hit->x)*(pos(0)-hit->x)+(pos(1)-hit->y)*(pos(1)-hit->y));

      // Use an un-normalized gaussian so that for a residual
      // of zero, we get a probability of 1.0.
      double p = finite(resi) ? exp(-resi*resi/2./variance):0.0;
      if(p>=MIN_FDC_HIT_PROB){
	fit.AddHit(hit->x,hit->y,hit->wire->origin(2),hit->covxx,hit->covxy,
		   hit->covyy);
	last_pos=pos;
	last_mom=mom;
	num_matched_hits++;
      }
    }
   
    if (num_matched_hits>MIN_HITS){
    // Set the initial parameters from the track candidate
      fit.SetSeed(tc->charge(),last_pos,last_mom);

      // Kalman filter 
      fit.KalmanLoop(TOF_MASS,chisq);

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

      // Fill in DKinematicData part
      track->setMass(0.0);
      track->setMomentum(mom);
      track->setPosition(pos);
      track->setCharge(track->q);

      _data.push_back(track);	
    }
  }
  

  return NOERROR;
}


//------------------
const string DTrack_factory_ALT3::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: q:       p:   theta:   phi:    x:    y:    z:");
	
	for(unsigned int i=0; i<_data.size(); i++){

		DTrack *track = _data[i];

		printnewrow();
		
		printcol("%x", i);
		printcol("%+d", (int)track->q);
		printcol("%3.3f", track->p);
		printcol("%1.3f", track->theta);
		printcol("%1.3f", track->phi);
		printcol("%2.2f", track->x);
		printcol("%2.2f", track->y);
		printcol("%2.2f", track->z);

		printrow();
	}
	
	return _table;
}

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
#include "FDC/DFDCSegment.h"
#include "DKalmanFilter.h"
#include "DReferenceTrajectory.h"

//------------------
// DTrack_factory_ALT3   (Constructor)
//------------------
DTrack_factory_ALT3::DTrack_factory_ALT3(){

  TOF_MASS = 0.13957018;
  MIN_FDC_HIT_PROB=0.1;
  MIN_HITS=6;
}

//------------------
// DTrack_factory_ALT3   (Destructor)
//------------------
DTrack_factory_ALT3::~DTrack_factory_ALT3()
{
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

  // Loop over track candidates
  for(unsigned int i=0; i<trackcandidates.size(); i++){ 
    const DTrackCandidate *tc = trackcandidates[i];    
    vector<const DFDCSegment *>segments;
    vector<const DCDCTrackHit *>cdchits;
    tc->GetT(segments);
    tc->GetT(cdchits);

     // Initialize energy loss sum
    double dEsum=0.;
    unsigned int num_matched_hits=0;
    DVector3 last_pos,last_mom; 
    double R=0;

    // Initialize Kalman filter with B-field
    DKalmanFilter fit(bfield,dgeom);

    for (unsigned int k=0;k<cdchits.size();k++){
      double r=cdchits[k]->wire->origin.Perp();
      if (r>R) R=r;

      fit.AddCDCHit(cdchits[k]);
    }
    num_matched_hits+=cdchits.size();

    if (segments.size()==0){    
      // Initialize the stepper 
      DMagneticFieldStepper stepper(bfield,tc->charge());
      
      // Find reference trajectory by swimming through the field
      DVector3 pos = tc->position();
      DVector3 mom = tc->momentum(); 	
      DVector3 norm(0,0,1);

      if (cdchits.size()>0){
	stepper.SwimToRadius(pos,mom,R+0.8,NULL);
      }
      last_pos=pos;
      last_mom=mom;
      
      // Add hits to be put through Kalman engine
      for(unsigned int j=0; j<fdctrackhits.size(); j++){
	const DFDCPseudo *hit = fdctrackhits[j];
	double variance=1.0; //guess for now
	
	// Swim to FDC hit plane
	stepper.SwimToPlane(pos,mom,hit->wire->origin,norm,NULL);

	// Find residual 
	double resi=sqrt((hit->x-pos.x())*(hit->x-pos.x())
			 +(hit->y-pos.y())*(hit->y-pos.y()));
	
	// Use an un-normalized gaussian so that for a residual
	// of zero, we get a probability of 1.0.
	double p = finite(resi) ? exp(-resi*resi/2./variance):0.0;
	if(p>=MIN_FDC_HIT_PROB){
	  fit.AddHit(hit->x,hit->y,hit->wire->origin(2),hit->covxx,
		     hit->covyy,hit->covxy,hit->dE);
	  last_pos=pos;
	  last_mom=mom;
	  num_matched_hits++;
	  dEsum+=hit->dE;
	}
      }      
    }
    else{	
      const DFDCSegment *segment=NULL;
      for (unsigned m=0;m<segments.size();m++){
	segment=segments[m];
	for (unsigned n=0;n<segment->hits.size();n++){
	  const DFDCPseudo *hit=segment->hits[n];
	  fit.AddHit(hit->x,hit->y,hit->wire->origin(2),hit->covxx,
		     hit->covyy,hit->covxy,hit->dE);
	  num_matched_hits++;
	  dEsum+=hit->dE;
	}

      }
      GetPositionAndMomentum(segment,last_pos,last_mom);
     
    }

    if (num_matched_hits>=MIN_HITS){
    // Set the initial parameters from the track candidate
      fit.SetSeed(tc->charge(),last_pos,last_mom);

      // Kalman filter 
      jerror_t error=fit.KalmanLoop(TOF_MASS);

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
	//track->rt=rt;

	// Fill in DKinematicData part
	track->setMass(0.0);
	track->setMomentum(mom);
	track->setPosition(pos);
	track->setCharge(track->q);

	// dEdx
	track->dE=dEsum;
	track->ds=fit.GetActivePathLength();

	_data.push_back(track);
      }
    }
  }
  

  return NOERROR;
}


// Obtain position and momentum at the exit of a given package using the 
// helical track model.
//
jerror_t DTrack_factory_ALT3::GetPositionAndMomentum(const DFDCSegment *segment,
					      DVector3 &pos, DVector3 &mom){
  // Position of track segment at last hit plane of package
  double x=segment->xc+segment->rc*cos(segment->Phi1);
  double y=segment->yc+segment->rc*sin(segment->Phi1);
  double z=segment->hits[0]->wire->origin(2);
 
  // Track parameters
  double kappa=segment->S(0,0);
  double phi0=segment->S(1,0);
  double tanl=segment->S(3,0);
  double z0=segment->S(4,0);

  // Useful intermediate variables
  double cosp=cos(phi0);
  double sinp=sin(phi0);
  double sperp=(z-z0)/tanl;
  double sin2ks=sin(2.*kappa*sperp);
  double cos2ks=cos(2.*kappa*sperp); 
  kappa=fabs(kappa);  // magnitude of curvature

  // Get Bfield
  double Bx,By,Bz,B;
  bfield->GetField(x,y,z,Bx,By,Bz);
  B=sqrt(Bx*Bx+By*By+Bz*Bz);
  
  // Momentum
  double px=(cosp*cos2ks-sinp*sin2ks)*0.003*B/2./kappa;
  double py=(sinp*cos2ks+cosp*sin2ks)*0.003*B/2./kappa;
  double pz=0.003*B*tanl/2./kappa;

  //if (sqrt(px*px+py*py)>PT_MAX) return VALUE_OUT_OF_RANGE;

  pos.SetXYZ(x,y,z);
  mom.SetXYZ(px,py,pz);

  return NOERROR;
}

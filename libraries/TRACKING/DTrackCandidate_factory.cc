// $Id$
//
//    File: DTrackCandidate_factory.cc
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//


#include "DTrackCandidate_factory.h"
#include "DANA/DApplication.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCSegment.h"
#include "DMagneticFieldStepper.h"
#include "DRiemannFit.h"
#include <TROOT.h>
#include <TH2F.h>

#define CUT 10.
#define RADIUS_CUT 28.
#define BEAM_VAR 0.01 // cm^2
#define Z_VERTEX 65.0

bool cdc_fdc_match(double p, double dist){
  if (p<=0.25) return true;
  if (dist < 1.+7.6/(p-0.25)) return true;
  return false;
}



jerror_t DTrackCandidate_factory::brun(JEventLoop* eventLoop,int runnumber){
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  bfield = dapp->GetBfield();

  // Get the position of the exit of the CDC endplate from DGeometry
  double endplate_z,endplate_dz,endplate_rmin;
  dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  cdc_endplate(2)=endplate_z;

  if(DEBUG_HISTS){
    dapp->Lock();
    match_dist=(TH2F*)gROOT->FindObject("match_dist");
    if (!match_dist) match_dist=new TH2F("match_dist","Matching distance",
					 60,0.,60.,100,0,25.);
    match_dist_vs_p=(TH2F*)gROOT->FindObject("match_dist_vs_p");
    if (!match_dist_vs_p) 
      match_dist_vs_p=new TH2F("match_dist_vs_p","Matching distance vs p",
			       50,0.,7.,100,0,25.);
    dapp->Unlock();
  }

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory::evnt(JEventLoop *loop, int eventnumber)
{
  vector<const DTrackCandidate*> cdctrackcandidates;
  vector<unsigned int> cdc_forward_ids;
  vector<DVector3> cdc_endplate_projections;
  vector<const DTrackCandidate*> fdctrackcandidates;
  
  loop->Get(cdctrackcandidates, "CDC");
  loop->Get(fdctrackcandidates, "FDCCathodes");
	
  // Normal vector for CDC endplate
  DVector3 norm(0,0,1);

  // Loop through the list of CDC candidates, flagging those that point toward 
  // the FDC.  The others are immediately put in the final candidate list.
  for(unsigned int i=0; i<cdctrackcandidates.size(); i++){	
    const DTrackCandidate *srccan = cdctrackcandidates[i];
    DVector3 mom=srccan->momentum();
    DVector3 pos=srccan->position();
    double theta=mom.Theta();
    
    // Propagate track to CDC endplate
    if (theta<atan(endplate_rmax/(cdc_endplate(2)-pos.z()))){
      DMagneticFieldStepper stepper(bfield,srccan->charge()); 
      if (stepper.SwimToPlane(pos,mom,cdc_endplate,norm,NULL)==false){
	cdc_endplate_projections.push_back(pos);
	cdc_forward_ids.push_back(i);
      }
    }
    else{
      DTrackCandidate *can = new DTrackCandidate;
      
      can->setMass(srccan->mass());
      can->setMomentum(srccan->momentum());
      can->setPosition(srccan->position());
      can->setCharge(srccan->charge());
      
      _data.push_back(can);
    }
  }

  // Loop through the list of FDC candidates looking for matches between the
  // CDC and the FDC in the transition region.
  for(unsigned int i=0; i<fdctrackcandidates.size(); i++){
    const DTrackCandidate *srccan = fdctrackcandidates[i];
    // Get the segment data
    vector<const DFDCSegment *>segments;
    srccan->GetT(segments);

    bool got_match=false;

    if (cdc_forward_ids.size()>0){
      double diff_min=1000.; // candidate matching difference

      // Propagate FDC track candidate back to CDC endplate 
      DVector3 mom=-1.*srccan->momentum();
      DVector3 pos;
      pos(0)=segments[0]->hits[0]->x;
      pos(1)=segments[0]->hits[0]->y;
      pos(2)=segments[0]->hits[0]->wire->origin(2);

      DMagneticFieldStepper stepper(bfield,srccan->charge()); 
      if (stepper.SwimToPlane(pos,mom,cdc_endplate,norm,NULL)==false){
	unsigned int jmin=0;
	double radius=0.;
	for (unsigned int j=0;j<cdc_forward_ids.size();j++){
	  unsigned int cdc_index=cdc_forward_ids[j];
	  double diff=(cdc_endplate_projections[cdc_index]-pos).Mag();
		 
	  if (diff<diff_min){
	    diff_min=diff;
	    jmin=j;
	    radius=pos.Perp();
	  }
	}
	if (DEBUG_HISTS){
	  match_dist->Fill(radius,diff_min);
	  match_dist_vs_p->Fill(srccan->momentum().Mag(),diff_min);
	}

	// Initialize the momentum vector to the FDC candidate
	DVector3 mom=srccan->momentum();
	
	if (cdc_fdc_match(mom.Mag(),diff_min)){
	  unsigned int cdc_index=cdc_forward_ids[jmin];

	  // Check that the charges match
	  if (cdctrackcandidates[cdc_index]->charge()==srccan->charge()){
	    // Remove the CDC candidate from the list and indicate that we 
	    // found a match
	    cdc_forward_ids.erase(cdc_forward_ids.begin()+jmin);	  
	    got_match=true;

	    if (radius<RADIUS_CUT){
	      vector<const DCDCTrackHit *>cdchits;
	      cdctrackcandidates[cdc_index]->GetT(cdchits);
	      
	      DRiemannFit fit;
	      for (unsigned int k=0;k<cdchits.size();k++){
		fit.AddHit(cdchits[k]->wire->origin.x(),
			   cdchits[k]->wire->origin.y(),
			   cdchits[k]->wire->origin.z());
	      }
	      for (unsigned int k=0;k<segments.size();k++){
		for (unsigned int n=0;n<segments[k]->hits.size();n++){
		  double covxx=segments[k]->hits[n]->covxx;
		  double covyy=segments[k]->hits[n]->covyy;
		  double covxy=segments[k]->hits[n]->covxy;
		  double x=segments[k]->hits[n]->x;
		  double y=segments[k]->hits[n]->y;
		  double z=segments[k]->hits[n]->wire->origin(2);
		  fit.AddHit(x,y,z,covxx,covyy,covxy);
		}
	      }	
	      // Fake point at origin
	      fit.AddHit(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	      // Fit the points to a circle
	      fit.FitCircle();
	      
	      // Compute new transverse momentum
	      double Bx,By,Bz;
	      bfield->GetField(segments[0]->hits[0]->x,segments[0]->hits[0]->y,
			       segments[0]->hits[0]->wire->origin(2),Bx,By,Bz);
	      double B=sqrt(Bx*Bx+By*By+Bz*Bz);
	      double pt=0.003*B*fit.rc;
	      double theta=srccan->momentum().Theta(); 
	      mom.SetMagThetaPhi(pt/sin(theta),theta,srccan->momentum().Phi());
	    }

	    // Put the candidate in the combined list
	    DTrackCandidate *can = new DTrackCandidate;
	    
	    can->setMass(srccan->mass());
	    can->setMomentum(mom);
	    can->setPosition(srccan->position());
	    can->setCharge(srccan->charge());
	    
	    _data.push_back(can);	    
	  }
	}
      }
    }
 
    if (!got_match){
      // Put the fdc candidate in the combined list
      DTrackCandidate *can = new DTrackCandidate;
      
      can->setMass(srccan->mass());
      can->setMomentum(srccan->momentum());
      can->setPosition(srccan->position());
      can->setCharge(srccan->charge());
      
      _data.push_back(can);
    }
  }
  
  // Unmatched CDC track candidates
  for (unsigned int j=0;j<cdc_forward_ids.size();j++){	  
    DTrackCandidate *can = new DTrackCandidate;
    const DTrackCandidate *cdccan = cdctrackcandidates[cdc_forward_ids[j]];
    
    can->setMass(cdccan->mass());
    can->setMomentum(cdccan->momentum());
    can->setPosition(cdccan->position());
    can->setCharge(cdccan->charge());
    
    _data.push_back(can);
  }	  
  
  return NOERROR;
}


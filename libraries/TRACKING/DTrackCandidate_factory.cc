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
#include "DHelicalFit.h"
#include <TROOT.h>
#include <TH2F.h>

#define CUT 10.
#define RADIUS_CUT 50.0
#define BEAM_VAR 0.01 // cm^2
#define Z_VERTEX 65.0
#define Z_MIN 45.
#define Z_MAX 85.
#define EPS 0.001

bool cdc_fdc_match(double p, double dist){
  if (p<=0.25) return true;
  if (dist < 1.+7.6/(p-0.25)) return true;
  return false;
}

bool cdchit_cmp(const DCDCTrackHit *a, const DCDCTrackHit *b){
  double ra=a->wire->origin.Perp();
  double rb=b->wire->origin.Perp();
  return (rb>ra);
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
    if (!match_dist){
      match_dist=new TH2F("match_dist","Matching distance",
			  60,0.,60.,100,0,25.);
      match_dist->SetXTitle("r (cm)");
      match_dist->SetYTitle("#Deltar (cm)");
    }
    match_dist_vs_p=(TH2F*)gROOT->FindObject("match_dist_vs_p");
    if (!match_dist_vs_p) {
      match_dist_vs_p=new TH2F("match_dist_vs_p","Matching distance vs p",
			       50,0.,7.,100,0,25.);
      match_dist_vs_p->SetYTitle("#Deltar (cm)");
      match_dist_vs_p->SetXTitle("p (GeV/c)");
    }
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
    if (theta<M_PI/6. && fdctrackcandidates.size()>0){      
      DMagneticFieldStepper stepper(bfield,srccan->charge()); 
      if (stepper.SwimToPlane(pos,mom,cdc_endplate,norm,NULL)==false){
	cdc_endplate_projections.push_back(pos);
	cdc_forward_ids.push_back(i);
      }
    }
    else{
      vector<const DCDCTrackHit *>cdchits;
      srccan->GetT(cdchits);

      DTrackCandidate *can = new DTrackCandidate;
      
      can->setMass(srccan->mass());
      can->setMomentum(srccan->momentum());
      can->setPosition(srccan->position());
      can->setCharge(srccan->charge());

      for (unsigned int n=0;n<cdchits.size();n++)
	can->AddAssociatedObject(cdchits[n]);    

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
    // Initialize the stepper 
    DMagneticFieldStepper stepper(bfield,srccan->charge());
 
    if (cdc_forward_ids.size()>0){
      double diff_min=1000.; // candidate matching difference

      // Propagate FDC track candidate back to CDC endplate 
      DVector3 mom_back,pos;
      GetPositionAndMomentum(segments[0],pos,mom_back);
      mom_back=-1.0*mom_back;
      if (stepper.SwimToPlane(pos,mom_back,cdc_endplate,norm,NULL)==false){
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
	    vector<const DCDCTrackHit *>cdchits;
	    cdctrackcandidates[cdc_index]->GetT(cdchits);	  
	    got_match=true;

	    if (radius<RADIUS_CUT){
	      DHelicalFit fit;
	      for (unsigned int k=0;k<cdchits.size();k++){	
		if (fabs(cdchits[k]->wire->stereo)<EPS){
		  double cov=0.8*0.8/12.;  //guess
		  fit.AddHitXYZ(cdchits[k]->wire->origin.x(),
				cdchits[k]->wire->origin.y(),
				cdchits[k]->wire->origin.z(),cov,cov,0.);
		}
	      }
	      for (unsigned int k=0;k<segments.size();k++){
		for (unsigned int n=0;n<segments[k]->hits.size();n++){
		  double covxx=segments[k]->hits[n]->covxx;
		  double covyy=segments[k]->hits[n]->covyy;
		  double covxy=segments[k]->hits[n]->covxy;
		  double x=segments[k]->hits[n]->x;
		  double y=segments[k]->hits[n]->y;
		  double z=segments[k]->hits[n]->wire->origin(2);
		
		  if (segments[k]->hits[n]->status<20) printf("uncorrected?\n");
		  fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
		}
	      }	
	      // Fake point at origin
	      fit.AddHitXYZ(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	      // Fit the points to a circle
	      fit.FitCircleRiemannCorrected(segments[0]->rc);
	      
	      // Compute new transverse momentum
	      double Bz_avg=0.;
	      for (unsigned int m=0;m<segments[0]->hits.size();m++){
		double Bx,By,Bz;
		DFDCPseudo *hit=segments[0]->hits[m];
		bfield->GetField(hit->x,hit->y,hit->wire->origin.z(),Bx,By,Bz);
		Bz_avg-=Bz;
	      }

	      Bz_avg/=float(segments[0]->hits.size());
	      double pt=0.003*Bz_avg*fit.r0;
	      double theta=srccan->momentum().Theta();
	      double phi=atan2(-fit.x0,fit.y0);
	      if (srccan->charge()<0) phi+=M_PI;
	      mom.SetMagThetaPhi(pt/sin(theta),theta,phi);
	    }

	    // Put the candidate in the combined list
	    DTrackCandidate *can = new DTrackCandidate;
	      
	    can->setMass(srccan->mass());
	    can->setMomentum(mom);
	    if (srccan->position().Z()>Z_MIN && srccan->position().Z()<Z_MAX)
	      can->setPosition(srccan->position());
	    else if (cdctrackcandidates[cdc_index]->position().Z()>Z_MIN 
		     && cdctrackcandidates[cdc_index]->position().Z()<Z_MAX)
	      can->setPosition(cdctrackcandidates[cdc_index]->position());
	    else{
	      DVector3 center(0,0,Z_VERTEX);
	      can->setPosition(center);
	    }
	    can->setCharge(srccan->charge());
	    
	    for (unsigned int m=0;m<segments.size();m++)
	      for (unsigned int n=0;n<segments[m]->hits.size();n++)
		can->AddAssociatedObject(segments[m]->hits[n]);
	    for (unsigned int n=0;n<cdchits.size();n++)
	      can->AddAssociatedObject(cdchits[n]); 
	    
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
      
      for (unsigned int m=0;m<segments.size();m++)
	for (unsigned int n=0;n<segments[m]->hits.size();n++)
	  can->AddAssociatedObject(segments[m]->hits[n]);
      
      // Try to gather up stray CDC hits from candidates that were not 
      // matched with the previous algorithm.  Use axial wires for now.
      DVector3 pos=srccan->position();
      DVector3 mom=srccan->momentum();
      for (unsigned int j=0;j<cdc_forward_ids.size();j++){
	vector<const DCDCTrackHit *>cdchits;
	cdctrackcandidates[cdc_forward_ids[j]]->GetT(cdchits);	 
	sort(cdchits.begin(),cdchits.end(),cdchit_cmp);
	
	unsigned int num_match=0;
	unsigned int num_axial=0;
	for (unsigned int m=0;m<cdchits.size();m++){
	  if (fabs(cdchits[m]->wire->stereo)<EPS){
	    double R=cdchits[m]->wire->origin.Perp();
	    
	    // Swim out from the "vertex"
	    if (stepper.SwimToRadius(pos,mom,R,NULL)==false){
	      double dx=cdchits[m]->wire->origin.x()-pos.x();
	      double dy=cdchits[m]->wire->origin.y()-pos.y();
	      double dr=sqrt(dx*dx+dy*dy);
	      
	      // Use an un-normalized gaussian so that for a residual
	      // of zero, we get a probability of 1.0.
	      double variance=1.0;
	      double prob = finite(dr) ? exp(-dr*dr/2./variance):0.0;
	      if (prob>0.1) num_match++;
	    }
	    num_axial++;
	  }
	}
	if (num_match==num_axial && num_match>0){
	  // Remove the CDC candidate from the list
	  cdc_forward_ids.erase(cdc_forward_ids.begin()+j); 

	  // Add the CDC hits to the track candidate
	  for (unsigned int m=0;m<cdchits.size();m++)
		can->AddAssociatedObject(cdchits[m]);
	  break;
	}
      }

      _data.push_back(can);
    }
  }
  
  // Unmatched CDC track candidates
  for (unsigned int j=0;j<cdc_forward_ids.size();j++){	  
    DTrackCandidate *can = new DTrackCandidate;
    const DTrackCandidate *cdccan = cdctrackcandidates[cdc_forward_ids[j]]; 
    vector<const DCDCTrackHit *>cdchits;
    cdccan->GetT(cdchits);

    can->setMass(cdccan->mass());
    can->setMomentum(cdccan->momentum());
    can->setPosition(cdccan->position());
    can->setCharge(cdccan->charge());
    
    for (unsigned int n=0;n<cdchits.size();n++)
      can->AddAssociatedObject(cdchits[n]);
   
    _data.push_back(can);
  }	
  
  return NOERROR;
}

// Obtain position and momentum at the exit of a given package using the 
// helical track model.
jerror_t DTrackCandidate_factory::GetPositionAndMomentum(
					      const DFDCSegment *segment,
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
  //B=sqrt(Bx*Bx+By*By+Bz*Bz);
  B=fabs(Bz);

  // Momentum
  double px=(cosp*cos2ks-sinp*sin2ks)*0.003*B/2./kappa;
  double py=(sinp*cos2ks+cosp*sin2ks)*0.003*B/2./kappa;
  double pz=0.003*B*tanl/2./kappa;

  pos.SetXYZ(x,y,z);
  mom.SetXYZ(px,py,pz);

  return NOERROR;
}



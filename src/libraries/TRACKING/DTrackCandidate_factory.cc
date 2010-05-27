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
  if (dist<10. && dist < 1.5+1.5/p) return true;
  return false;
}

bool cdchit_cmp(const DCDCTrackHit *a, const DCDCTrackHit *b){
  if (a->wire==NULL || b->wire==NULL){
    cout << "Null pointer in CDC hit list??" << endl;
    return false;
  }
  return (b->wire->ring>a->wire->ring);
}



jerror_t DTrackCandidate_factory::brun(JEventLoop* eventLoop,int runnumber){
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  bfield = dapp->GetBfield();

  // Get the position of the exit of the CDC endplate from DGeometry
  double endplate_z,endplate_dz,endplate_rmin;
  dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  cdc_endplate.SetZ(endplate_z);

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
  vector<unsigned int> cdc_backward_ids;
  vector<int> forward_matches;
  vector<DVector3> cdc_endplate_projections;
  vector<const DTrackCandidate*> fdctrackcandidates;
  
  loop->Get(cdctrackcandidates, "CDC");
  loop->Get(fdctrackcandidates, "FDCCathodes");

  // Fill list of ids for FDC candidates
  for(unsigned int i=0; i<fdctrackcandidates.size(); i++)
    forward_matches.push_back(0);
	
  // Normal vector for CDC endplate
  DVector3 norm(0,0,1);

  //Loop through the list of CDC candidates, flagging those that point toward 
  // the FDC.  The others are put in the final candidate list.
  for(unsigned int i=0; i<cdctrackcandidates.size(); i++){	
    const DTrackCandidate *srccan = cdctrackcandidates[i];
    DVector3 mom=srccan->momentum();
    DVector3 pos=srccan->position();
    double theta=mom.Theta();
    
    // Propagate track to CDC endplate
    if (theta<M_PI/4. && fdctrackcandidates.size()>0){
      DMagneticFieldStepper stepper(bfield,srccan->charge()); 
      if (stepper.SwimToPlane(pos,mom,cdc_endplate,norm,NULL)==false){
	cdc_endplate_projections.push_back(pos);
	cdc_forward_ids.push_back(i);
      }
    }
    else{
      cdc_backward_ids.push_back(i);
    }
  }

  // Variables for candidate number accounting
  unsigned int num_forward_cdc_cands_remaining=cdc_forward_ids.size();
  unsigned int num_fdc_cands_remaining=fdctrackcandidates.size();
 
  // Loop through the list of FDC candidates looking for matches between the
  // CDC and the FDC in the transition region.
  for(unsigned int i=0; i<fdctrackcandidates.size(); i++){
    const DTrackCandidate *srccan = fdctrackcandidates[i];

    // Momentum and position vectors for the FDC candidate
    DVector3 mom=srccan->momentum();
    DVector3 pos=srccan->position();

    // Get the segment data
    vector<const DFDCSegment *>segments;
    srccan->GetT(segments);

    bool got_match=false;
    // Initialize the stepper 
    DMagneticFieldStepper stepper(bfield,srccan->charge());
 
    if (cdc_forward_ids.size()>0){
      double diff_min=1000.; // candidate matching difference

      // Propagate FDC track candidate back to CDC endplate 
      DVector3 mom_back,pos_back;
      GetPositionAndMomentum(segments[0],pos_back,mom_back);
      mom_back=-1.0*mom_back;
      if (stepper.SwimToPlane(pos_back,mom_back,cdc_endplate,norm,NULL)
	  ==false){
	unsigned int jmin=0;
	double radius=0.;

	for (unsigned int j=0;j<cdc_forward_ids.size();j++){
	  unsigned int cdc_index=cdc_forward_ids[j];
	  
	  // Check that the charges match
	  if (cdctrackcandidates[cdc_index]->charge()==srccan->charge()){
	    double diff=(cdc_endplate_projections[j]-pos_back).Mag();
 
	    if (diff<diff_min){
	      diff_min=diff;
	      jmin=j;
	      radius=pos_back.Perp();
	    }
	  }
	}
	if (DEBUG_HISTS){
	  match_dist->Fill(radius,diff_min);
	  match_dist_vs_p->Fill(srccan->momentum().Mag(),diff_min);
	}

	got_match=cdc_fdc_match(mom.Mag(),diff_min);
		
	if (got_match){
	  unsigned int cdc_index=cdc_forward_ids[jmin];

	  // Remove the CDC candidate from the id list because we 
	  // found a match
	  cdc_forward_ids.erase(cdc_forward_ids.begin()+jmin);
	  num_forward_cdc_cands_remaining--;

	  // Mark the FDC candidate as matched
	  forward_matches[i]=1;
	  num_fdc_cands_remaining--;

	  // Get the associated cdc hits
	  vector<const DCDCTrackHit *>cdchits;
	  cdctrackcandidates[cdc_index]->GetT(cdchits);

	  // Create new track candidate object 
	  DTrackCandidate *can = new DTrackCandidate;

	  // Add cdc and fdc hits to the track as associated objects
	  unsigned int num_fdc_hits=0;
	  for (unsigned int m=0;m<segments.size();m++)
	    for (unsigned int n=0;n<segments[m]->hits.size();n++){
	      can->AddAssociatedObject(segments[m]->hits[n]);
	      num_fdc_hits++;
	    }
	  for (unsigned int n=0;n<cdchits.size();n++)
	    can->AddAssociatedObject(cdchits[n]); 
	  
	  // variables for calculating average Bz
	  double Bz_avg=0.,Bx,By,Bz;
	  unsigned int num_hits=0;
	  
	  // Redo circle fit with additional hits
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
	      double z=segments[k]->hits[n]->wire->origin.z();
	      
	      fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	      bfield->GetField(x,y,z,Bx,By,Bz);
	      Bz_avg-=Bz;
	    }
	    num_hits+=segments[k]->hits.size();
	  }	
	  // Fake point at origin
	  fit.AddHitXYZ(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);

	  // if the cdc and fdc candidates do not agree as to the particle's
	  // charge, set the charge according to which detector has more hits
	  if (srccan->charge()!=cdctrackcandidates[cdc_index]->charge() 
	      && cdchits.size()>num_fdc_hits){
	    can->setCharge(cdctrackcandidates[cdc_index]->charge());
	  }
	  else
	    can->setCharge(srccan->charge());

	  // Fit the points to a circle
	  if (fit.FitCircleRiemannCorrected(segments[0]->rc)==NOERROR){	  
	    // Compute new transverse momentum
	    Bz_avg/=double(num_hits);
	    double pt=0.003*Bz_avg*fit.r0;
	    
	    // Determine the polar angle
	    double theta=srccan->momentum().Theta();
	    if (num_fdc_hits<cdchits.size())
	      theta=cdctrackcandidates[cdc_index]->momentum().Theta();

	    // Determine the azimuthal angle
	    double phi=atan2(-fit.x0,fit.y0);
	    if (can->charge()<0) phi+=M_PI;

	    // Set the momentum
	    mom.SetMagThetaPhi(pt/sin(theta),theta,phi);
	  }

	  // Set the mass and momentum
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
	  	 	  
	  _data.push_back(can);	    
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

      unsigned int num_fdc_hits=0;
      for (unsigned int m=0;m<segments.size();m++)
	for (unsigned int n=0;n<segments[m]->hits.size();n++){
	  can->AddAssociatedObject(segments[m]->hits[n]);
	  num_fdc_hits++;
	}
      
      // Try to gather up stray CDC hits from candidates that were not 
      // matched with the previous algorithm.  Use axial wires for now.
      for (unsigned int j=0;j<cdc_forward_ids.size();j++){
	vector<const DCDCTrackHit *>cdchits;
	cdctrackcandidates[cdc_forward_ids[j]]->GetT(cdchits);	 
	sort(cdchits.begin(),cdchits.end(),cdchit_cmp);
	
	unsigned int num_match=0;
	unsigned int num_axial=0;

	DVector3 pos=srccan->position();
	DVector3 mom=srccan->momentum();
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
	  num_forward_cdc_cands_remaining--;
	 
	  // Mark the FDC candidate as matched
	  forward_matches[i]=1;
	  num_fdc_cands_remaining--;

	  // variables for calculating average Bz
	  double Bz_avg=0.,Bx,By,Bz;
	  unsigned int num_hits=0;
	  
	  // Add the CDC hits to the track candidate and redo helical fit
	  DHelicalFit fit;
	  for (unsigned int m=0;m<cdchits.size();m++){
	    if (fabs(cdchits[m]->wire->stereo)<EPS){
	      double cov=0.8*0.8/12.;  //guess
	      fit.AddHitXYZ(cdchits[m]->wire->origin.x(),
			    cdchits[m]->wire->origin.y(),
			    cdchits[m]->wire->origin.z(),cov,cov,0.);
	    }
	    can->AddAssociatedObject(cdchits[m]);
	  }
	  for (unsigned int k=0;k<segments.size();k++){
	    for (unsigned int n=0;n<segments[k]->hits.size();n++){
	      double covxx=segments[k]->hits[n]->covxx;
	      double covyy=segments[k]->hits[n]->covyy;
	      double covxy=segments[k]->hits[n]->covxy;
	      double x=segments[k]->hits[n]->x;
	      double y=segments[k]->hits[n]->y;
	      double z=segments[k]->hits[n]->wire->origin.z();
	      
	      fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	      bfield->GetField(x,y,z,Bx,By,Bz);
	      Bz_avg-=Bz;
	    }
	    num_hits+=segments[k]->hits.size();
	  }	
	  // Fake point at origin
	  fit.AddHitXYZ(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);

	  // if the cdc and fdc candidates do not agree as to the particle's
	  // charge, set the charge according to which detector has more hits
	  unsigned int cdc_index=cdc_forward_ids[j];
	  if (srccan->charge()!=cdctrackcandidates[cdc_index]->charge() 
	      && cdchits.size()>num_fdc_hits){
	    can->setCharge(cdctrackcandidates[cdc_index]->charge());
	  }
	  else
	    can->setCharge(srccan->charge());


	  // Fit the points to a circle
	  if (fit.FitCircleRiemannCorrected(segments[0]->rc)==NOERROR){
	    // Compute new transverse momentum
	    Bz_avg/=double(num_hits);
	    double pt=0.003*Bz_avg*fit.r0;
	    double theta=can->momentum().Theta();
	    double phi=atan2(-fit.x0,fit.y0);
	    if (can->charge()<0) phi+=M_PI;
	    mom.SetMagThetaPhi(pt/sin(theta),theta,phi);
	    can->setMomentum(mom);
	  }
	    
	  _data.push_back(can);
	  
	  break;
	}
      } // loop over forward cdc candidates
    }
  }

  // If starting with the fdc candidates did not lead to a complete set of
  // CDC-FDC matches, try looping over the remaing CDC candidates that point
  // toward the FDC.
  vector<int>cdc_forward_matches(cdc_forward_ids.size());
  for (unsigned int i=0;i<cdc_forward_ids.size();i++){ 
    const DTrackCandidate *srccan=cdctrackcandidates[cdc_forward_ids[i]];
    cdc_forward_matches[i]=0;
    
    // Get hits already linked to this candidate from associated objects
    vector<const DCDCTrackHit *>cdchits;
    srccan->GetT(cdchits);
   
    // loop over fdc candidates
    for (unsigned int k=0;k<fdctrackcandidates.size();k++){
      if (forward_matches[k]==0){
	const DTrackCandidate *fdccan = fdctrackcandidates[k];

	// Get the segment data
	vector<const DFDCSegment *>segments;
	fdccan->GetT(segments);

	// Initialize the stepper 
	DMagneticFieldStepper stepper(bfield,srccan->charge());

	// Momentum and position vectors for the CDC candidate
	DVector3 mom=srccan->momentum();
	DVector3 pos=srccan->position();

	// Try to match unmatched fdc candidates
	int num_match=0;
	int num_hits=0;
	for (unsigned int m=0;m<segments.size();m++){
	  for (unsigned int n=0;n<segments[m]->hits.size();n++){
	    unsigned int ind=segments[m]->hits.size()-1-n;
	    if (stepper.SwimToPlane(pos,mom,
				    segments[m]->hits[ind]->wire->origin,
				    norm,NULL)
		==false){
	      double dx=segments[m]->hits[ind]->x-pos.x();
	      double dy=segments[m]->hits[ind]->y-pos.y();
	      double dr=sqrt(dx*dx+dy*dy);
	      
	      // Use an un-normalized gaussian so that for a residual
	      // of zero, we get a probability of 1.0.
	      double variance=1.0;
	      double prob = finite(dr) ? exp(-dr*dr/2./variance):0.0;
	      if (prob>0.1) num_match++;
	      num_hits++;
	    }
	  }
	}
	if (double(num_match)/double(num_hits)>0.5){
	  DTrackCandidate *can = new DTrackCandidate;

	  can->setMass(srccan->mass());
	  can->setPosition(srccan->position());
	  can->setCharge(srccan->charge());
	  can->setMomentum(srccan->momentum());
	  
	  // mark the fdc track candidate as matched
	  forward_matches[k]=1; 
	  num_fdc_cands_remaining--;

	  //Mark the cdc track candidate as matched
	  cdc_forward_matches[i]=1;
	  
	  // variables for calculating average Bz
	  double Bz_avg=0.,Bx,By,Bz;
	  unsigned int num_hits=0;

	  // Redo helical fit with all available hits
	  DHelicalFit fit; 
	  for (unsigned int m=0;m<cdchits.size();m++){
	    if (fabs(cdchits[m]->wire->stereo)<EPS){
	      double cov=0.8*0.8/12.;  //guess
	      fit.AddHitXYZ(cdchits[m]->wire->origin.x(),
			    cdchits[m]->wire->origin.y(),
			    cdchits[m]->wire->origin.z(),cov,cov,0.);
	    }
	    can->AddAssociatedObject(cdchits[m]);
	  } 	
	  for (unsigned int m=0;m<segments.size();m++){
	    for (unsigned int n=0;n<segments[m]->hits.size();n++){
	      double covxx=segments[m]->hits[n]->covxx;
	      double covyy=segments[m]->hits[n]->covyy;
	      double covxy=segments[m]->hits[n]->covxy;
	      double x=segments[m]->hits[n]->x;
	      double y=segments[m]->hits[n]->y;
	      double z=segments[m]->hits[n]->wire->origin.z();
	      
	      fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	      bfield->GetField(x,y,z,Bx,By,Bz);
	      Bz_avg-=Bz;

	      can->AddAssociatedObject(segments[m]->hits[n]);
	    }
	    num_hits+=segments[m]->hits.size();
	  }	
	  // Fake point at origin
	  fit.AddHitXYZ(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	  // Fit the points to a circle
	  if (fit.FitCircleRiemannCorrected(segments[0]->rc)==NOERROR){
	    // Compute new transverse momentum
	    Bz_avg/=double(num_hits);
	    double pt=0.003*Bz_avg*fit.r0;
	    double theta=srccan->momentum().Theta();
	    double phi=atan2(-fit.x0,fit.y0);
	    if (can->charge()<0) phi+=M_PI;
	    mom.SetMagThetaPhi(pt/sin(theta),theta,phi);
	    can->setMomentum(mom);
	  }

	  _data.push_back(can);
	  break;
	}
      }
    } // loop over fdc candidates
  }

  // If there are FDC candidates remaining, use the best track candidate 
  // knowledge we have so far to do one final check for matches
  if (num_fdc_cands_remaining>0){
    for (unsigned int j=0;j<_data.size();j++){
      const DTrackCandidate *srccan= _data[j];
      
      // Get hits already linked to this candidate from associated objects
      vector<const DCDCTrackHit *>cdchits;
      _data[j]->GetT(cdchits);
      vector<const DFDCPseudo *>fdchits;
      _data[j]->GetT(fdchits);	 	 

      for (unsigned int k=0;k<fdctrackcandidates.size();k++){
	if (forward_matches[k]==0){
	  const DTrackCandidate *fdccan = fdctrackcandidates[k];

	  // Get the segment data
	  vector<const DFDCSegment *>segments;
	  fdccan->GetT(segments);
	  
	  // Initialize the stepper 
	  DMagneticFieldStepper stepper(bfield,srccan->charge());

	  // Momentum and position vectors for the CDC candidate
	  DVector3 mom=srccan->momentum();
	  DVector3 pos=srccan->position();
		  
	  // Try to match unmatched fdc candidates
	  int num_match=0;
	  int num_hits=0;
	  for (unsigned int m=0;m<segments.size();m++){
	    for (unsigned int i=0;i<segments[m]->hits.size();i++){
	      unsigned int ind=segments[m]->hits.size()-1-i;
	      if (stepper.SwimToPlane(pos,mom,
				      segments[m]->hits[ind]->wire->origin,
				      norm,NULL)
		==false){
		double dx=segments[m]->hits[ind]->x-pos.x();
		double dy=segments[m]->hits[ind]->y-pos.y();
		double dr=sqrt(dx*dx+dy*dy);
		
		// Use an un-normalized gaussian so that for a residual
	      // of zero, we get a probability of 1.0.
		double variance=1.0;
		double prob = finite(dr) ? exp(-dr*dr/2./variance):0.0;
	      if (prob>0.1) num_match++;
	      num_hits++;
	      }
	    }
	  }   
	  if (double(num_match)/double(num_hits)>0.5){
	    forward_matches[k]=1; 
	    
	    // variables for calculating average Bz
	    double Bz_avg=0.,Bx,By,Bz;
	    unsigned int num_hits=0;
	    
	    // Redo helical fit with all available hits
	    DHelicalFit fit; 
	    for (unsigned int m=0;m<cdchits.size();m++){
	      if (fabs(cdchits[m]->wire->stereo)<EPS){
		double cov=0.8*0.8/12.;  //guess
		fit.AddHitXYZ(cdchits[m]->wire->origin.x(),
			      cdchits[m]->wire->origin.y(),
			      cdchits[m]->wire->origin.z(),cov,cov,0.);
	      }
	    }
	    for (unsigned int m=0;m<fdchits.size();m++){
	      double x=fdchits[m]->x;
	      double y=fdchits[m]->y;
	      double z=fdchits[m]->wire->origin.z();
	      double covxx=fdchits[m]->covxx;
	      double covyy=fdchits[m]->covyy;
	      double covxy=fdchits[m]->covxy;
	      fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	      
	      bfield->GetField(x,y,z,Bx,By,Bz);
	      Bz_avg-=Bz;
	      num_hits++;
	    }
	    for (unsigned int m=0;m<segments.size();m++){
	      for (unsigned int n=0;n<segments[m]->hits.size();n++){
		double covxx=segments[m]->hits[n]->covxx;
		double covyy=segments[m]->hits[n]->covyy;
		double covxy=segments[m]->hits[n]->covxy;
		double x=segments[m]->hits[n]->x;
		double y=segments[m]->hits[n]->y;
		double z=segments[m]->hits[n]->wire->origin.z();
		
		fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
		bfield->GetField(x,y,z,Bx,By,Bz);
		Bz_avg-=Bz;
		
		_data[j]->AddAssociatedObject(segments[m]->hits[n]);
	      }
	      num_hits+=segments[m]->hits.size();
	    }	

	    // Fake point at origin
	    fit.AddHitXYZ(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	    // Fit the points to a circle
	    if (fit.FitCircleRiemannCorrected(segments[0]->rc)==NOERROR){
	      // Compute new transverse momentum
	      Bz_avg/=double(num_hits);
	      double pt=0.003*Bz_avg*fit.r0;
	      double theta=_data[j]->momentum().Theta();
	      double phi=atan2(-fit.x0,fit.y0);
	      if (_data[j]->charge()<0) phi+=M_PI;
	      mom.SetMagThetaPhi(pt/sin(theta),theta,phi);
	      _data[j]->setMomentum(mom);
	    }
	  }
	}
      } // loop over fdc track candidates
    }
  }

  // Copy remaining FDC candidates to final candidate list
  for (unsigned int i=0;i<forward_matches.size();i++){
    if (forward_matches[i]==0){
      DTrackCandidate *can = new DTrackCandidate;
      const DTrackCandidate *fdccan = fdctrackcandidates[i]; 
      
      can->setMass(fdccan->mass());
      can->setMomentum(fdccan->momentum());
      can->setPosition(fdccan->position());
      can->setCharge(fdccan->charge());
      
      // Get the segment data
      vector<const DFDCSegment *>segments;
      fdccan->GetT(segments);
      for (unsigned int m=0;m<segments.size();m++){
	for (unsigned int n=0;n<segments[m]->hits.size();n++){
	  can->AddAssociatedObject(segments[m]->hits[n]);
	}
      }
      
      _data.push_back(can);
    }
  }
  
  // Unmatched CDC track candidates
  for (unsigned int j=0;j<cdc_forward_ids.size();j++){
    if (cdc_forward_matches[j]==0){
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
  }
  for (unsigned int j=0;j<cdc_backward_ids.size();j++){	  
    DTrackCandidate *can = new DTrackCandidate;
    const DTrackCandidate *cdccan = cdctrackcandidates[cdc_backward_ids[j]]; 
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
  double z=segment->hits[0]->wire->origin.z();

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



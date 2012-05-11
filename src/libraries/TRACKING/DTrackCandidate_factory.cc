// $Id$
//
//    File: DTrackCandidate_factory.cc
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//


#include "DTrackCandidate_factory.h"
#include "DANA/DApplication.h"
#include "DMagneticFieldStepper.h"

#include <TROOT.h>
#include <TH2F.h>

#define CUT 10.
#define RADIUS_CUT 50.0
#define BEAM_VAR 0.01 // cm^2
#define EPS 0.001

//------------------
// cdc_fdc_match
//------------------
inline bool cdc_fdc_match(double p_fdc,double p_cdc,double dist){
  //double frac=fabs(1.-p_cdc/p_fdc);  
  if (dist<10. && dist < 1.5+1.5/p_fdc
      // && frac<0.5
      ) return true;
  return false;
}

//------------------
// SegmentSortByLayerincreasing
//------------------
inline bool SegmentSortByLayerincreasing(const DFDCSegment* const &segment1, const DFDCSegment* const &segment2) {
	// Compare DFDCSegment->DFDCPseudo[0]->DFDCWire->layer
	int layer1 = 100; // defaults just in case there is a segment with no hits
	int layer2 = 100;
	
	if(segment1->hits.size()>0)layer1=segment1->hits[0]->wire->layer;
	if(segment2->hits.size()>0)layer2=segment2->hits[0]->wire->layer;

	return layer1 < layer2;
}

//------------------
// CDCHitSortByLayerincreasing
//------------------
inline bool CDCHitSortByLayerincreasing(const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2) {
	// Used to sort CDC hits by layer (ring) with innermost layer hits first

	// if same ring, sort by wire number
	if(hit1->wire->ring == hit2->wire->ring){
		return hit1->wire->straw < hit2->wire->straw;
	}

	return hit1->wire->ring < hit2->wire->ring;
}

//------------------
// FDCHitSortByLayerincreasing
//------------------
inline bool FDCHitSortByLayerincreasing(const DFDCPseudo* const &hit1, const DFDCPseudo* const &hit2) {
	// Used to sort CDC hits by layer (ring) with innermost layer hits first

	// if same ring, sort by wire number
	if(hit1->wire->layer == hit2->wire->layer){
		return hit1->wire->wire < hit2->wire->wire;
	}

	return hit1->wire->layer < hit2->wire->layer;
}


//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory::brun(JEventLoop* eventLoop,int runnumber){
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  bfield = dapp->GetBfield();

  // Get the position of the exit of the CDC endplate from DGeometry
  double endplate_z,endplate_dz,endplate_rmin;
  dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  cdc_endplate.SetZ(endplate_z+endplate_dz);

  dgeom->GetTargetZ(TARGET_Z);
 
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

  // List of cdc hits
  vector<const DCDCTrackHit*>mycdchits;
  loop->Get(mycdchits);

  // Vector to keep track of cdc hits used in candidates
  vector<unsigned int>used_cdc_hits(mycdchits.size());

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
    if (theta<M_PI_4 && fdctrackcandidates.size()>0){
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
	
	// JANA does not maintain the order that the segments were added
	// as associated objects. Therefore, we need to reorder them here
	// so segment[0] is the most upstream segment.
	sort(segments.begin(), segments.end(), SegmentSortByLayerincreasing);

	bool got_match=false;
    if (cdc_forward_ids.size()>0){
      double diff_min=1000.; // candidate matching difference
      unsigned int jmin=0;
      double radius=0.;
      // Magnitude of the momentum of the potential cdc match
      double p_cdc=0.;
      
      for (unsigned int j=0;j<cdc_forward_ids.size();j++){
	//unsigned int cdc_index=cdc_forward_ids[j];
	//const DTrackCandidate *cdccan=cdctrackcandidates[cdc_index];
	
	// Check that the charges match
	//if (cdccan->charge()==srccan->charge())
	  {
	  double diff=(cdc_endplate_projections[j]-pos).Mag();
	  
	  if (diff<diff_min){
	    diff_min=diff;
	    jmin=j;
	    radius=pos.Perp();
	    p_cdc=cdctrackcandidates[cdc_forward_ids[jmin]]->momentum().Mag();
	  }
	}
      }
      if (DEBUG_HISTS){
	match_dist->Fill(radius,diff_min);
	match_dist_vs_p->Fill(srccan->momentum().Mag(),diff_min);
      }

      double frac=fabs(1.-p_cdc/mom.Mag());  
      if (frac<0.5){      
	got_match=cdc_fdc_match(mom.Mag(),p_cdc,diff_min);
		
	if (got_match){
	  unsigned int cdc_index=cdc_forward_ids[jmin];
	  const DTrackCandidate *cdccan=cdctrackcandidates[cdc_index];
	  
	  // Mark the FDC candidate as matched
	  forward_matches[i]=1;
	  num_fdc_cands_remaining--;
	  
	  // Get the associated cdc hits
	  vector<const DCDCTrackHit *>cdchits;
	  cdccan->GetT(cdchits);
	  
	  // Sort CDC hits by layer
	  sort(cdchits.begin(), cdchits.end(), CDCHitSortByLayerincreasing);
	  
	  // Create new track candidate object 
	  DTrackCandidate *can = new DTrackCandidate;
	  can->setCharge(srccan->charge());
	  can->used_cdc_indexes=cdccan->used_cdc_indexes;
	  
	  // Add cdc and fdc hits to the track as associated objects
	  unsigned int num_fdc_hits=0;
	  for (unsigned int m=0;m<segments.size();m++)
	    for (unsigned int n=0;n<segments[m]->hits.size();n++){
	      can->AddAssociatedObject(segments[m]->hits[n]);
	      num_fdc_hits++;
	    }
	  for (unsigned int n=0;n<cdchits.size();n++){
	    used_cdc_hits[cdccan->used_cdc_indexes[n]]=1;
	    can->AddAssociatedObject(cdchits[n]); 
	  }
	  
	  // variables for calculating average Bz
	  double Bz_avg=0.;
	  unsigned int num_hits=0;
	  
	  // Redo circle fit with additional hits
	  DHelicalFit fit;
	  for (unsigned int k=0;k<cdchits.size();k++){	
	    if (cdchits[k]->is_stereo==false){
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
	      double x=segments[k]->hits[n]->xy.X();
	      double y=segments[k]->hits[n]->xy.Y();
	      double z=segments[k]->hits[n]->wire->origin.z();
	      
	      fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	      //bfield->GetField(x,y,z,Bx,By,Bz);
	      //Bz_avg-=Bz;
	      Bz_avg-=bfield->GetBz(x,y,z);
	    }
	    num_hits+=segments[k]->hits.size();
	  }	
	  // Fake point at origin
	  fit.AddHitXYZ(0.,0.,TARGET_Z,BEAM_VAR,BEAM_VAR,0.);
	  
	  // Fit the points to a circle
	  if (fit.FitCircleRiemannCorrected(segments[0]->rc)==NOERROR){	  
	    // Compute new transverse momentum
	    Bz_avg/=double(num_hits);
	    
	    // Determine the polar angle
	    /*
	      unsigned int num_cdc_hits=cdchits.size();
	      double theta
	      =(srccan->momentum().Theta()*num_fdc_hits
	      +cdccan->momentum().Theta()*num_cdc_hits)/
	      (num_fdc_hits+num_cdc_hits);
	    */
	    double theta=srccan->momentum().Theta();
	    
	    //if (num_fdc_hits<cdchits.size())
	    //  theta=cdccan->momentum().Theta();
	    fit.tanl=tan(M_PI_2-theta);
	    fit.z_vertex=cdccan->position().Z();
	    fit.q=can->charge();
	    
	    GetPositionAndMomentum(fit,Bz_avg,cdchits[0]->wire->origin,pos,mom);
	    can->setMomentum(mom);
	    can->setPosition(pos);
	  }
	  else{
	    can->setMomentum(srccan->momentum());
	    can->setPosition(srccan->position());
	    
	  }
	
	  // Set the mass and momentum
	  can->setMass(srccan->mass());
	  
	  _data.push_back(can);	    
	  
	  // Remove the CDC candidate from the id list because we 
	  // found a match
	  cdc_forward_ids.erase(cdc_forward_ids.begin()+jmin);
	  num_forward_cdc_cands_remaining--;
	}
	else{  
	  // Try to gather up stray CDC hits from candidates that were not 
	  // matched with the previous algorithm. 
	  for (unsigned int j=0;j<cdc_forward_ids.size();j++){
	    const DTrackCandidate *cdccan=cdctrackcandidates[cdc_forward_ids[j]];
	    
	    vector<const DCDCTrackHit *>cdchits;
	    cdccan->GetT(cdchits);	 
	    sort(cdchits.begin(),cdchits.end(),CDCHitSortByLayerincreasing);
	    
	    unsigned int num_match=0;
	    unsigned int num_cdc=0;
	    
	    DVector3 pos=srccan->position();
	    DVector3 mom=srccan->momentum();
	    double q=srccan->charge();
	    for (unsigned int m=0;m<cdchits.size();m++){
	      // Use an un-normalized gaussian so that for a residual
	      // of zero, we get a probability of 1.0.
	      double variance=1.0;
	      // Use a helical approximation to the track to match both axial and
	      // stereo wires
	      double dr2=DocaToHelix(cdchits[m],q,pos,mom);
	      double prob=finite(dr2) ? exp(-dr2/(2.*variance)):0.0;

	      if (prob>0.1) num_match++;
	      
	      num_cdc++;
	    }
	    if (num_match>0 && float(num_match)/float(num_cdc)>0.5){
	      // Mark this fdc candidate as having a match to a cdc candidate
	      forward_matches[i]=1;
	      
	      // Put the fdc candidate in the combined list
	      DTrackCandidate *can = new DTrackCandidate;
	      
	      can->setMass(srccan->mass());
	      //can->setMomentum(srccan->momentum());
	      //can->setPosition(srccan->position());
	      can->setCharge(srccan->charge());
	      
	      can->used_cdc_indexes=cdccan->used_cdc_indexes;
	      
	      unsigned int num_fdc_hits=0;
	      for (unsigned int m=0;m<segments.size();m++)
		for (unsigned int n=0;n<segments[m]->hits.size();n++){
		  can->AddAssociatedObject(segments[m]->hits[n]);
		  num_fdc_hits++;
		}
	      
	      // Mark the FDC candidate as matched
	      forward_matches[i]=1;
	      num_fdc_cands_remaining--;
	      
	      // variables for calculating average Bz
	      double Bz_avg=0.;
	      unsigned int num_hits=0;
	      
	      // Add the CDC hits to the track candidate and redo helical fit
	      DHelicalFit fit;
	      for (unsigned int m=0;m<cdchits.size();m++){
		if (cdchits[m]->is_stereo==false){
		  double cov=0.8*0.8/12.;  //guess
		  fit.AddHitXYZ(cdchits[m]->wire->origin.x(),
				cdchits[m]->wire->origin.y(),
				cdchits[m]->wire->origin.z(),cov,cov,0.);
		}
		used_cdc_hits[cdccan->used_cdc_indexes[m]]=1;
		
		can->AddAssociatedObject(cdchits[m]);
	      }
	      for (unsigned int k=0;k<segments.size();k++){
		for (unsigned int n=0;n<segments[k]->hits.size();n++){
		  double covxx=segments[k]->hits[n]->covxx;
		  double covyy=segments[k]->hits[n]->covyy;
		  double covxy=segments[k]->hits[n]->covxy;
		  double x=segments[k]->hits[n]->xy.X();
		  double y=segments[k]->hits[n]->xy.Y();
		  double z=segments[k]->hits[n]->wire->origin.z();
		
		  fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
		  //bfield->GetField(x,y,z,Bx,By,Bz);
		  //Bz_avg-=Bz;
		  Bz_avg-=bfield->GetBz(x,y,z);
		}
		num_hits+=segments[k]->hits.size();
	      }	
	      // Fake point at origin
	      fit.AddHitXYZ(0.,0.,TARGET_Z,BEAM_VAR,BEAM_VAR,0.);
	      
	      // Fit the points to a circle
	      if (fit.FitCircleRiemannCorrected(segments[0]->rc)==NOERROR){
		// Compute new transverse momentum
		Bz_avg/=double(num_hits);
		
		// Determine the polar angle
		/*
		  unsigned int num_cdc_hits=cdchits.size();
		  double theta
		  =(srccan->momentum().Theta()*num_fdc_hits
		  +cdccan->momentum().Theta()*num_cdc_hits)/
		  (num_fdc_hits+num_cdc_hits);
		*/
		double theta=srccan->momentum().Theta();
		
		fit.tanl=tan(M_PI_2-theta);
		fit.z_vertex=cdccan->position().Z();
		fit.q=GetCharge(fit,segments[0]->hits[0],cdccan->position());
		can->setCharge(fit.q);
		GetPositionAndMomentum(fit,Bz_avg,cdchits[0]->wire->origin,pos,mom);
		can->setMomentum(mom);
		can->setPosition(pos);
	      }
	      else{
		can->setMomentum(srccan->momentum());
		can->setPosition(srccan->position());	    
	      }
	      
	    
	      _data.push_back(can);
	      
	      // Remove the CDC candidate from the list
	      cdc_forward_ids.erase(cdc_forward_ids.begin()+j); 
	      num_forward_cdc_cands_remaining--;
	      
	      break;
	    }

	  } // loop over forward cdc candidates
	}
      } // check for momentum ratio
    }
  } // loop over fdc candidates

  // If starting with the fdc candidates did not lead to a complete set of
  // CDC-FDC matches, try looping over the remaining CDC candidates that point
  // toward the FDC.
  vector<int>cdc_forward_matches(cdc_forward_ids.size());
  for (unsigned int i=0;i<cdc_forward_ids.size();i++){ 
    const DTrackCandidate *srccan=cdctrackcandidates[cdc_forward_ids[i]];
    cdc_forward_matches[i]=0;
 
    // Get hits already linked to this candidate from associated objects
    vector<const DCDCTrackHit *>cdchits;
    srccan->GetT(cdchits);
    sort(cdchits.begin(), cdchits.end(), CDCHitSortByLayerincreasing);
 
    // Store the magnitude of the momentum
    double p_cdc=srccan->momentum().Mag();
   
    // loop over fdc candidates
    for (unsigned int k=0;k<fdctrackcandidates.size();k++){
      if (forward_matches[k]==0){
	const DTrackCandidate *fdccan = fdctrackcandidates[k];
	double p_fdc=fdccan->momentum().Mag();
	if (fabs(1-p_fdc/p_cdc)<0.5){

	  // Get the segment data
	  vector<const DFDCSegment *>segments;
	  fdccan->GetT(segments);
	  sort(segments.begin(), segments.end(), SegmentSortByLayerincreasing);
	  
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
		double dx=segments[m]->hits[ind]->xy.X()-pos.x();
		double dy=segments[m]->hits[ind]->xy.Y()-pos.y();
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
	  //if (double(num_match)/double(num_hits)>0.5)
	  if (num_match>3) {
	    DTrackCandidate *can = new DTrackCandidate;
	    
	    can->setMass(srccan->mass());
	    can->setCharge(srccan->charge());
	    can->used_cdc_indexes=srccan->used_cdc_indexes;
	    
	    // mark the fdc track candidate as matched
	    forward_matches[k]=1; 
	    num_fdc_cands_remaining--;
	    
	    //Mark the cdc track candidate as matched
	    cdc_forward_matches[i]=1;
	    
	    // variables for calculating average Bz
	    double Bz_avg=0.;
	    unsigned int num_hits=0;
	    
	    // Redo helical fit with all available hits
	    DHelicalFit fit; 
	    for (unsigned int m=0;m<cdchits.size();m++){
	      if (cdchits[m]->is_stereo==false){
		double cov=0.8*0.8/12.;  //guess
		fit.AddHitXYZ(cdchits[m]->wire->origin.x(),
			      cdchits[m]->wire->origin.y(),
			      cdchits[m]->wire->origin.z(),cov,cov,0.);
	      }
	      used_cdc_hits[srccan->used_cdc_indexes[m]]=1;
	      can->AddAssociatedObject(cdchits[m]);
	    } 	
	    for (unsigned int m=0;m<segments.size();m++){
	      for (unsigned int n=0;n<segments[m]->hits.size();n++){
		double covxx=segments[m]->hits[n]->covxx;
		double covyy=segments[m]->hits[n]->covyy;
		double covxy=segments[m]->hits[n]->covxy;
		double x=segments[m]->hits[n]->xy.X();
		double y=segments[m]->hits[n]->xy.Y();
		double z=segments[m]->hits[n]->wire->origin.z();
		
		fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
		//bfield->GetField(x,y,z,Bx,By,Bz);
		//Bz_avg-=Bz;
		Bz_avg-=bfield->GetBz(x,y,z);
		
		can->AddAssociatedObject(segments[m]->hits[n]);
	      }
	      num_hits+=segments[m]->hits.size();
	    }	
	    // Fake point at origin
	    fit.AddHitXYZ(0.,0.,TARGET_Z,BEAM_VAR,BEAM_VAR,0.);
	    // Fit the points to a circle
	    if (fit.FitCircleRiemannCorrected(segments[0]->rc)==NOERROR){
	      // Compute new transverse momentum
	      Bz_avg/=double(num_hits);
	      
	      // Determine the polar angle
	      /*
		unsigned int num_cdc_hits=cdchits.size();
		double theta
		=(fdccan->momentum().Theta()*num_match
		+srccan->momentum().Theta()*num_cdc_hits)/
		(num_match+num_cdc_hits);
	      */
	      double theta=fdccan->momentum().Theta();
	      
	      fit.tanl=tan(M_PI_2-theta);
	      fit.z_vertex=srccan->position().Z();
	      fit.q=GetCharge(fit,segments[0]->hits[0],srccan->position());
	      can->setCharge(fit.q);
	      GetPositionAndMomentum(fit,Bz_avg,cdchits[0]->wire->origin,pos,mom);
	      can->setMomentum(mom);
	      can->setPosition(pos);
	    }
	    else{
	      can->setMomentum(srccan->momentum());
	      can->setPosition(srccan->position());
	    }
	    
	    _data.push_back(can);
	    break;
	  }
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
      srccan->GetT(cdchits);
      vector<const DFDCPseudo *>fdchits;
      srccan->GetT(fdchits);
	  
      sort(cdchits.begin(), cdchits.end(), CDCHitSortByLayerincreasing);
	  sort(fdchits.begin(), fdchits.end(), FDCHitSortByLayerincreasing);

      for (unsigned int k=0;k<fdctrackcandidates.size();k++){
	if (forward_matches[k]==0){
	  const DTrackCandidate *fdccan = fdctrackcandidates[k];

	  // Get the segment data
	  vector<const DFDCSegment *>segments;
	  fdccan->GetT(segments);
      sort(segments.begin(), segments.end(), SegmentSortByLayerincreasing);

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
		double dx=segments[m]->hits[ind]->xy.X()-pos.x();
		double dy=segments[m]->hits[ind]->xy.Y()-pos.y();
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
	    double Bz_avg=0.;
	    unsigned int num_hits=0;

	    _data[j]->used_cdc_indexes=srccan->used_cdc_indexes;

	    // Redo helical fit with all available hits
	    DHelicalFit fit; 
	    for (unsigned int m=0;m<cdchits.size();m++){
	      used_cdc_hits[srccan->used_cdc_indexes[m]]=1;

	      if (cdchits[m]->is_stereo==false){
		double cov=0.8*0.8/12.;  //guess
		fit.AddHitXYZ(cdchits[m]->wire->origin.x(),
			      cdchits[m]->wire->origin.y(),
			      cdchits[m]->wire->origin.z(),cov,cov,0.);
	      }
	    }
	    for (unsigned int m=0;m<fdchits.size();m++){
	      double x=fdchits[m]->xy.X();
	      double y=fdchits[m]->xy.Y();
	      double z=fdchits[m]->wire->origin.z();
	      double covxx=fdchits[m]->covxx;
	      double covyy=fdchits[m]->covyy;
	      double covxy=fdchits[m]->covxy;
	      fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	      
	      //bfield->GetField(x,y,z,Bx,By,Bz);
	      //Bz_avg-=Bz;
	      Bz_avg-=bfield->GetBz(x,y,z);

	      num_hits++;
	    }
	    for (unsigned int m=0;m<segments.size();m++){
	      for (unsigned int n=0;n<segments[m]->hits.size();n++){
		double covxx=segments[m]->hits[n]->covxx;
		double covyy=segments[m]->hits[n]->covyy;
		double covxy=segments[m]->hits[n]->covxy;
		double x=segments[m]->hits[n]->xy.X();
		double y=segments[m]->hits[n]->xy.Y();
		double z=segments[m]->hits[n]->wire->origin.z();
		
		fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
		//bfield->GetField(x,y,z,Bx,By,Bz);
		//Bz_avg-=Bz;
		Bz_avg-=bfield->GetBz(x,y,z);
		
		_data[j]->AddAssociatedObject(segments[m]->hits[n]);
	      }
	      num_hits+=segments[m]->hits.size();
	    }	

	    // Fake point at origin
	    fit.AddHitXYZ(0.,0.,TARGET_Z,BEAM_VAR,BEAM_VAR,0.);
	    // Fit the points to a circle
	    if (fit.FitCircleRiemannCorrected(segments[0]->rc)==NOERROR){
	      // Compute new transverse momentum
	      Bz_avg/=double(num_hits);
	      
	      /*
	      unsigned int num_cdc_hits=cdchits.size();
	      double theta
		=(fdccan->momentum().Theta()*num_match
		  +srccan->momentum().Theta()*num_cdc_hits)/
		(num_match+num_cdc_hits);
	      */
	      double theta=fdccan->momentum().Theta();
	      
	      fit.tanl=tan(M_PI_2-theta);
	      fit.z_vertex=srccan->position().Z();
	      fit.q=GetCharge(fit,segments[0]->hits[0],_data[j]->position());
	      _data[j]->setCharge(fit.q);
	 
	      GetPositionAndMomentum(fit,Bz_avg,cdchits[0]->wire->origin,pos,mom);
	      _data[j]->setMomentum(mom);
	      _data[j]->setPosition(pos);
	    }
	    else{
	      _data[j]->setMomentum(srccan->momentum());
	      _data[j]->setPosition(srccan->position());
	    }
	  

	  }
	}
      } // loop over fdc track candidates
    }
  }

  // Unmatched CDC track candidates
  for (unsigned int j=0;j<cdc_forward_ids.size();j++){
    if (cdc_forward_matches[j]==0){
      DTrackCandidate *can = new DTrackCandidate;
      const DTrackCandidate *cdccan = cdctrackcandidates[cdc_forward_ids[j]];
      vector<const DCDCTrackHit *>cdchits;
      cdccan->GetT(cdchits);
	  sort(cdchits.begin(), cdchits.end(), CDCHitSortByLayerincreasing);
      
      can->setMass(cdccan->mass());
      can->setMomentum(cdccan->momentum());
      can->setPosition(cdccan->position());
      can->setCharge(cdccan->charge());

      for (unsigned int n=0;n<cdchits.size();n++){
	used_cdc_hits[cdccan->used_cdc_indexes[n]]=1;
	can->AddAssociatedObject(cdchits[n]);
      }
      _data.push_back(can);
    }	
  }
  for (unsigned int j=0;j<cdc_backward_ids.size();j++){	  
    DTrackCandidate *can = new DTrackCandidate;
    const DTrackCandidate *cdccan = cdctrackcandidates[cdc_backward_ids[j]]; 
    vector<const DCDCTrackHit *>cdchits;
    cdccan->GetT(cdchits);
	sort(cdchits.begin(), cdchits.end(), CDCHitSortByLayerincreasing);
    
    can->setMass(cdccan->mass());
    can->setMomentum(cdccan->momentum());
    can->setPosition(cdccan->position());
    can->setCharge(cdccan->charge());
    
    for (unsigned int n=0;n<cdchits.size();n++){
      used_cdc_hits[cdccan->used_cdc_indexes[n]]=1;
      can->AddAssociatedObject(cdchits[n]);
    }
   
    _data.push_back(can);
  }	

  unsigned int num_unmatched_cdcs=0;
  for (unsigned int i=0;i<used_cdc_hits.size();i++){
    if (!used_cdc_hits[i]) num_unmatched_cdcs++;
  }

  // Copy remaining FDC candidates to final candidate list
  for (unsigned int i=0;i<forward_matches.size();i++){
    if (forward_matches[i]==0){
      // Magnetic field
      double Bz_avg=0.;	  
      double denom=0.;

      DTrackCandidate *can = new DTrackCandidate;
      const DTrackCandidate *fdccan = fdctrackcandidates[i]; 
 // Get the segment data
      vector<const DFDCSegment *>segments;
      fdccan->GetT(segments);
      sort(segments.begin(), segments.end(), SegmentSortByLayerincreasing);
      for (unsigned int m=0;m<segments.size();m++){
	for (unsigned int n=0;n<segments[m]->hits.size();n++){
	  can->AddAssociatedObject(segments[m]->hits[n]);
	  if (num_unmatched_cdcs>0){
	     double x=segments[m]->hits[n]->xy.X();
	     double y=segments[m]->hits[n]->xy.Y();
	     double z=segments[m]->hits[n]->wire->origin.z();
	     Bz_avg-=bfield->GetBz(x,y,z);
	     denom+=1.;
	  }
	}
      }
      if (num_unmatched_cdcs>0) Bz_avg/=denom;

      
      DVector3 mom=fdccan->momentum();
      DVector3 pos=fdccan->position();
      double q=fdccan->charge();

      // Try to match cdc hits that were not linked into track candidates
      // with the fdc candidate
      if (num_unmatched_cdcs>0){
	unsigned int num_match_cdc=0;
	unsigned int num_axial=0;
	unsigned int id_for_smallest_r=0;
	double rmin=1000.;
	vector<int>matched_ids;

	for (unsigned int k=0;k<used_cdc_hits.size();k++){
	  if (!used_cdc_hits[k]){
	    // Use an un-normalized gaussian so that for a residual
	    // of zero, we get a probability of 1.0.
	    double variance=1.0;
	    // Use a helical approximation to the track to match both axial and
	    // stereo wires
	    double dr2=DocaToHelix(mycdchits[k],q,pos,mom);
	    double prob=finite(dr2) ? exp(-dr2/(2.*variance)):0.0;  

	    if (prob>0.5){
	      double r=mycdchits[k]->wire->origin.Perp();
	      if (r<rmin){
		rmin=r;
		id_for_smallest_r=k;
	      }

	      num_match_cdc++;
	      num_unmatched_cdcs--;
	      used_cdc_hits[k]=1;
	      
	      can->used_cdc_indexes.push_back(k);
	      can->AddAssociatedObject(mycdchits[k]);
	
	      if (mycdchits[k]->is_stereo==false){
		num_axial++;
		matched_ids.push_back(k);
	      }
	    }
	  }	 
	} 
	// If we matched to more than three cdc hits, try to improve the candidate
	// parameters
	if (num_match_cdc>3){

	  DHelicalFit fit;
	  fit.r0=segments[0]->rc;
	  fit.x0=segments[0]->xc;
	  fit.y0=segments[0]->yc;
	  fit.phi=segments[0]->phi0;
	
	  // If we have more than one axial hit, redo the circle fit with the 
	  // additional data
	  if (num_axial>1){	  
	    for (unsigned int m=0;m<matched_ids.size();m++){
	      double cov=0.8*0.8/12.;  //guess
	      unsigned int ind=matched_ids[m];
	      fit.AddHitXYZ(mycdchits[ind]->wire->origin.x(),
			    mycdchits[ind]->wire->origin.y(),
			    mycdchits[ind]->wire->origin.z(),cov,cov,0.);
	    }
	    for (unsigned int k=0;k<segments.size();k++){
	      for (unsigned int n=0;n<segments[k]->hits.size();n++){
		double covxx=segments[k]->hits[n]->covxx;
		double covyy=segments[k]->hits[n]->covyy;
		double covxy=segments[k]->hits[n]->covxy;
		double x=segments[k]->hits[n]->xy.X();
		double y=segments[k]->hits[n]->xy.Y();
		double z=segments[k]->hits[n]->wire->origin.z();
		
		fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);	
	      }
	    }
	    // Fake point at origin
	    fit.AddHitXYZ(0.,0.,TARGET_Z,BEAM_VAR,BEAM_VAR,0.);
	    
	    // Redo the fit
	    fit.FitCircleRiemannCorrected(fit.r0);
	  }
	  fit.tanl=tan(M_PI_2-mom.Theta());
	  fit.z_vertex=0; // this will be changed later
	  fit.q=q;
	  GetPositionAndMomentum(fit,Bz_avg,
				 mycdchits[id_for_smallest_r]->wire->origin,
				 pos,mom);
	  // Find the z-position at the new position in x and y
	  DVector2 xy0(pos.X(),pos.Y());
	  double tworc=2.*fit.r0;
	  double ratio=(segments[0]->hits[0]->xy-xy0).Mod()/tworc;
	  double sperp=(ratio<1.)?tworc*asin(ratio):tworc*M_PI_2;

	  pos.SetZ(segments[0]->hits[0]->wire->origin.z()-sperp*fit.tanl);
    	} // if we succeeded in matching at least one cdc
      } // if we have cdc hits not matched to candidates yet...
      
     
      can->setMass(fdccan->mass());
      can->setMomentum(mom);
      can->setPosition(pos);
      can->setCharge(q);
      
      _data.push_back(can);
    }
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
  //double kappa=segment->q/(2.*segment->rc);
  double phi0=segment->phi0;
  double tanl=segment->tanl;
  double z0=segment->z_vertex;

  // Useful intermediate variables
  double cosp=cos(phi0);
  double sinp=sin(phi0);
  double sperp=(z-z0)/tanl;
  //double twoks=2.*kappa*sperp;
  double twoks=segment->q*sperp/segment->rc;
  double sin2ks=sin(twoks);
  double cos2ks=cos(twoks); 

  // Get Bfield
  double B=fabs(bfield->GetBz(x,y,z));

  // Momentum
  double pt=0.003*B*segment->rc;
  double px=pt*(cosp*cos2ks-sinp*sin2ks);
  double py=pt*(sinp*cos2ks+cosp*sin2ks);
  double pz=pt*tanl;

  pos.SetXYZ(x,y,z);
  mom.SetXYZ(px,py,pz);

  return NOERROR;
}


// Get the position and momentum at a fixed radius from the beam line
jerror_t DTrackCandidate_factory::GetPositionAndMomentum(DHelicalFit &fit,
							 double Bz,
							 const DVector3 &origin,
							 DVector3 &pos,
							 DVector3 &mom){
  double r2=90.0;
  double xc=fit.x0;
  double yc=fit.y0;
  double rc=fit.r0;
  double rc2=rc*rc;
  double xc2=xc*xc;
  double yc2=yc*yc;
  double xc2_plus_yc2=xc2+yc2;
  double a=(r2-xc2_plus_yc2-rc2)/(2.*rc);
  double temp1=yc*sqrt(xc2_plus_yc2-a*a);
  double temp2=xc*a;
  double cosphi_plus=(temp2+temp1)/xc2_plus_yc2;
  double cosphi_minus=(temp2-temp1)/xc2_plus_yc2;

  // Direction tangent and transverse momentum
  double tanl=fit.tanl;
  double pt=0.003*Bz*rc;

  if(!isfinite(temp1) || !isfinite(temp2)){
    // We did not find an intersection between the two circles, so return 
    // sensible defaults for pos and mom.

    double my_phi=fit.phi;
    double sinphi=sin(my_phi);
    double cosphi=cos(my_phi);
    double D=-fit.q*rc-xc/sinphi;
    pos.SetXYZ(-D*sinphi,D*cosphi,fit.z_vertex);
    mom.SetXYZ(pt*cosphi,pt*sinphi,pt*tanl);

    return NOERROR;

  }

  double phi_plus=acos(cosphi_plus);
  double phi_minus=acos(cosphi_minus);
  double x_plus=xc+rc*cosphi_plus;
  double x_minus=xc+rc*cosphi_minus;
  double y_plus=yc+rc*sin(phi_plus);
  double y_minus=yc+rc*sin(phi_minus);

  // if the resulting radial position on the circle from the fit does not agree
  // with the radius to which we are matching, we have the wrong sign for phi+ 
  // or phi-
  double r2_plus=x_plus*x_plus+y_plus*y_plus;
  double r2_minus=x_minus*x_minus+y_minus*y_minus;  
  if (fabs(r2-r2_plus)>EPS){
    phi_plus*=-1.;
    y_plus=yc+rc*sin(phi_plus);
  }
  if (fabs(r2-r2_minus)>EPS){
    phi_minus*=-1.;
    y_minus=yc+rc*sin(phi_minus);
  }

  // Choose phi- or phi+ depending on proximity to one of the cdc hits
  double xwire=origin.x();
  double ywire=origin.y();
  double dx=x_minus-xwire;
  double dy=y_minus-ywire;
  double d2_minus=dx*dx+dy*dy;
  dx=x_plus-xwire;
  dy=y_plus-ywire;
  double d2_plus=dx*dx+dy*dy;
  if (d2_plus>d2_minus){
    phi_minus*=-1.;
    if (fit.q<0) phi_minus+=M_PI;  
    //double dphi=M_PI_2-phi_minus-fit.phi;
    //while (dphi>2.*M_PI) dphi-=2*M_PI;
    //while (dphi<-2.*M_PI) dphi+=2*M_PI;   
    //if (dphi<-M_PI) dphi+=2*M_PI;
    //if (dphi>M_PI) dphi-=2*M_PI;
    //dphi=0.;
    //pos.SetXYZ(x_minus,y_minus,fit.z_vertex+fit.q*rc*dphi*tanl);
    pos.SetXYZ(x_minus,y_minus,fit.z_vertex);
    mom.SetXYZ(pt*sin(phi_minus),pt*cos(phi_minus),pt*tanl);

  }
  else{
    phi_plus*=-1.;   
    if (fit.q<0) phi_plus+=M_PI;
    //double dphi=M_PI_2-phi_plus-fit.phi;
    //while (dphi>2.*M_PI) dphi-=2*M_PI;
    //while (dphi<-2.*M_PI) dphi+=2*M_PI;
    //if (dphi<-M_PI) dphi+=2*M_PI;
    //if (dphi>M_PI) dphi-=2*M_PI; 
    //dphi=0.;
    //pos.SetXYZ(x_plus,y_plus,fit.z_vertex+fit.q*rc*dphi*tanl); 
    pos.SetXYZ(x_plus,y_plus,fit.z_vertex); 
    mom.SetXYZ(pt*sin(phi_plus),pt*cos(phi_plus),pt*tanl);

  }

  return NOERROR;
}

// Routine to do a crude match between cdc wires and a helical approximation to
// the trajectory
double DTrackCandidate_factory::DocaToHelix(const DCDCTrackHit *hit,
					    double q,
					    const DVector3 &pos,
					    const DVector3 &mom){
  double pt=mom.Perp();
  double phi=mom.Phi();
  double sinphi=sin(phi);
  double cosphi=cos(phi);
  double tanl=tan(M_PI_2-mom.Theta());
  double x0=pos.X();
  double y0=pos.Y();
  double z0=pos.Z();
  double B=bfield->GetBz(x0,y0,z0);
  double twokappa=0.003*B*q/pt;
  
  double sperp=0;
  DVector3 origin=hit->wire->origin;
  double z0w=origin.z();
  DVector3 dir=hit->wire->udir;
  double uz=dir.z();
  
  DVector3 wirepos;
  double old_doca2=1e8;
  double doca2=old_doca2;
  double z=z0;
  double r2=x0*x0+y0*y0;
  while (z>50. && z<170. && r2<3600.){
    old_doca2=doca2;
    sperp-=1.;
    double twoks=twokappa*sperp;
    double sin2ks=sin(twoks);
    double cos2ks=cos(twoks);
    double x=x0+(cosphi*sin2ks-sinphi*(1.-cos2ks))/twokappa;
    double y=y0+(sinphi*sin2ks+cosphi*(1.-cos2ks))/twokappa;
    r2=x*x+y*y;
    z=z0+sperp*tanl;
    
    wirepos=origin+((z-z0w)/uz)*dir;
    double xw=wirepos.x();
    double yw=wirepos.y();
    
    doca2=(x-xw)*(x-xw)+(y-yw)*(y-yw);
    if (doca2>old_doca2){
      break;
    } 
  }

  return old_doca2;
}


// Find the particle charge given the helical fit result and an fdc hit 
// on the track
double DTrackCandidate_factory::GetCharge(DHelicalFit &fit,
					  const DFDCPseudo *fdchit,
					  const DVector3 &pos
					  ){
  // Get circle parameters
  double rc=fit.r0;
  double xc=fit.x0;
  double yc=fit.y0;
  // Compute phi rotation from "vertex" to fdc hit
  double dphi=(fdchit->wire->origin.Z()-pos.z())/(rc*fit.tanl);
  double Phi1=atan2(pos.Y()-yc,pos.X()-xc);
  
  // Positive and negative changes in phi
  double phiplus=Phi1+dphi;
  double phiminus=Phi1-dphi;
  DVector2 plus(xc+rc*cos(phiplus),yc+rc*sin(phiplus));
  DVector2 minus(xc+rc*cos(phiminus),yc+rc*sin(phiminus));

  // Compute differences 
  double d2plus=(plus-fdchit->xy).Mod2();
  double d2minus=(minus-fdchit->xy).Mod2();
    
  if (d2minus<d2plus)
    return (-1.0);
  else
    return (+1.0);
}

// $Id$
//
//    File: DTrackCandidate_factory_FDCCathodes.cc
// Created: Tue Nov  6 13:37:08 EST 2007
// Creator: staylor (on Linux ifarml1.jlab.org 2.4.21-47.0.1.ELsmp i686)
//
/// This factory links segments in the FDC packages into track candidates 
/// by swimming through the field from one package to the next.

#include "DTrackCandidate_factory_FDCCathodes.h"
#include "DANA/DApplication.h"
#include "DMagneticFieldStepper.h"
#include "FDC/DFDCPseudo_factory.h"
#include "FDC/DFDCSegment_factory.h"
#include "DRiemannFit.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>

#define MATCH_RADIUS(p) (2.79+2.88/(p)/(p))
#define MAX_SEGMENTS 20
#define HALF_PACKAGE 6.0
#define FDC_OUTER_RADIUS 50.0 
#define BEAM_VAR 0.01 // cm^2
#define HIT_CHI2_CUT 10.0
#define Z_VERTEX 65.0
#define PT_MAX 9.0

///
/// DTrackCandidate_factory_FDCCathodes::brun():
///
jerror_t DTrackCandidate_factory_FDCCathodes::brun(JEventLoop* eventLoop, int eventNo) {
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();
  
  if(DEBUG_HISTS){
    dapp->Lock();
    match_dist=(TH2F*)gROOT->FindObject("match_dist");
    if (!match_dist) match_dist=new TH2F("match_dist","Matching distance",
					 50,0.,7,100,0,25.);
    dapp->Unlock();
  }
    
  return NOERROR;
}

// Local routine for sorting segments by charge and curvature
bool DTrackCandidate_segment_cmp(const DFDCSegment *a, const DFDCSegment *b){
  double k1=a->S(0,0),k2=b->S(0,0);
  double q1=k1/fabs(k1),q2=k2/fabs(k2);
  if (q1!=q2) return q1<q2;
  return fabs(k1)<fabs(k2);  
}

//------------------
// evnt:  main segment linking routine
//------------------
jerror_t DTrackCandidate_factory_FDCCathodes::evnt(JEventLoop *loop, int eventnumber)
{
  vector<const DFDCSegment*>segments;
  eventLoop->Get(segments);
  std::sort(segments.begin(), segments.end(), DTrackCandidate_segment_cmp);

  // Group segments by package
  vector<DFDCSegment*>package[4]; 
  for (unsigned int i=0;i<segments.size();i++){
     const DFDCSegment *segment=segments[i];
 package[(segment->hits[0]->wire->layer-1)/6].push_back((DFDCSegment*)segment);
  }
      
  double zpackage[4];  // z-positions of entrances to FDC packages 
  zpackage[0]=DFDCGeometry::GetZpackage(1);
  zpackage[1]=DFDCGeometry::GetZpackage(7);
  zpackage[2]=DFDCGeometry::GetZpackage(13);
  zpackage[3]=DFDCGeometry::GetZpackage(19);
  DFDCSegment *match2=NULL;
  DFDCSegment *match3=NULL;
  DFDCSegment *match4=NULL;
  unsigned int match_id=0;

  // Bail if there are too many segments
  if (package[0].size()+package[1].size()+package[2].size()
      +package[3].size()>MAX_SEGMENTS)
    return UNRECOVERABLE_ERROR; 

  // First deal with tracks with segments in the first package
  if (package[0].size()>0){
    // Loop over segments in the first package, matching them to segments in 
    // the second, third, and fourth (most downstream) packages.
    for (unsigned int i=0;i<package[0].size();i++){
      DFDCSegment *segment=package[0][i];
      
      // Tracking parameters from first segment
      double tanl=segment->S(3,0);
      double phi0=segment->S(1,0);
      double kappa=segment->S(0,0);
      double zvertex=segment->S(4,0);
      // Sign of the charge
      double q=kappa/fabs(kappa);
      
      double qsum=q;
      // Start filling vector of segments belonging to current track    
      vector<DFDCSegment*>segments; 
      segments.push_back(segment);

      // Check that the tangent of the dip angle makes sense for FDC hits
      if (segment->S(3,0)<=0.0){
	continue;
      }
 
      // Try matching to package 2
      if (package[1].size()>0 && 
	  (match2=GetTrackMatch(q,zpackage[1],segment,package[1],match_id))
	  !=NULL){
	// Insert the segment from package 2 into the track 
	segments.push_back(match2);
	
	// remove the segment from the list 
	package[1].erase(package[1].begin()+match_id);

	qsum+=match2->S(0,0)/fabs(match2->S(0,0));

	// Try matching to package 3
	if (package[2].size()>0 && 
	    (match3=GetTrackMatch(q,zpackage[2],match2,package[2],match_id))
	    !=NULL){
	  // Insert the segment from package 3 into the track 
	  segments.push_back(match3);

	  // remove the segment from the list 
	  package[2].erase(package[2].begin()+match_id);

	  qsum+=match3->S(0,0)/fabs(match3->S(0,0));

	  // Try matching to package 4
	  if (package[3].size()>0 && 
	      (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	       match_id))!=NULL){
	    // Insert the segment from package 4 into the track 
	    //segments.push_back(match4);

	    // remove the segment from the list 
	    package[3].erase(package[3].begin()+match_id);

	    qsum+=match4->S(0,0)/fabs(match4->S(0,0));
	  }
	}
	// No match in package 3, try for 4
	else if(package[3].size()>0 && 
		(match4=GetTrackMatch(q,zpackage[3],match2,package[3],
		 match_id))!=NULL){
	  // Insert the segment from package 4 into the track 
	  // segments.push_back(match4);

	  // remove the segment from the list 
	  package[3].erase(package[3].begin()+match_id);

	  qsum+=match4->S(0,0)/fabs(match4->S(0,0));
	}
      }
      // No match in package 2, try for 3
      else if (package[2].size()>0 && 
	       (match3=GetTrackMatch(q,zpackage[2],segment,package[2],
		match_id))!=NULL){
	// Insert the segment from package 3 into the track
	segments.push_back(match3);

	// remove the segment from the list 
	package[2].erase(package[2].begin()+match_id);

	qsum+=match3->S(0,0)/fabs(match3->S(0,0));
	
	// Try matching to package 4
	if (package[3].size()>0 && 
	    (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	     match_id))!=NULL){
	  // Insert the segment from package 4 into the track 
	  //segments.push_back(match4);

	  // remove the segment from the list 
	  package[3].erase(package[3].begin()+match_id);

	  qsum+=match4->S(0,0)/fabs(match4->S(0,0));
	}
      }    
      // No match to package 2 or 3, try 4
      else if (package[3].size()>0 && 
	       (match4=GetTrackMatch(q,zpackage[3],segment,package[3],
		match_id))!=NULL){
	// Insert the segment from package 4 into the track 
	//segments.push_back(match4);
	
	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);

	qsum+=match4->S(0,0)/fabs(match4->S(0,0));
      }

      if (qsum>0) q=1.;
      else if (qsum<0) q=-1.;
      else q=0.;

      if (segments.size()>1){
	DRiemannFit fit;
	for (unsigned int m=0;m<segments.size();m++){
	  for (unsigned int n=0;n<segments[m]->hits.size();n++){
	    DFDCPseudo *hit=segments[m]->hits[n];
	    double covxx=hit->covxx;
	    double covyy=hit->covyy;
	    double covxy=hit->covxy;
	    double hit_chi2=segments[segments.size()-1]->track[n].chi2;
	    /*if (hit_chi2>1.){
	      covxx*=hit_chi2;      
	      covyy*=hit_chi2;
	      covxy*=hit_chi2;
	      }*/
	    fit.AddHit(hit->x,hit->y,hit->wire->origin(2),covxx,covyy,covxy);
	  }
	}
	// Fake point at origin
	fit.AddHit(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	if (fit.DoFit(segments[0]->rc)==NOERROR){      
	  // Charge
	  if (q==0) q=fit.q;
	  // Curvature
	  segments[1]->S(0,0)=kappa=q/2./fit.rc;
	  // Estimate for azimuthal angle
	  phi0=atan2(-fit.xc,fit.yc); 
	  if (q<0) phi0+=M_PI;
	  segments[1]->S(1,0)=phi0;
	  // remaining tracking parameters
	  tanl=fit.tanl;
	  zvertex=fit.zvertex;
	 
	  segments[1]->S(3,0)=tanl;
	  segments[1]->S(4,0)=zvertex;
	  segments[1]->xc=fit.xc;
	  segments[1]->yc=fit.yc;
	  segments[1]->rc=fit.rc;
	  
	  // Try to match to package 2 again.
	  if (match2==NULL && package[1].size()>0 &&
	      (match2=GetTrackMatch(q,zpackage[1],segments[1],package[1],
				    match_id))!=NULL){ 
	    // remove the segment from the list 
	    package[1].erase(package[1].begin()+match_id);

	    qsum+=match2->S(0,0)/fabs(match2->S(0,0));
	  }
	  
	  // Try to match to package 3 again.
	  if (match3==NULL && package[2].size()>0 &&
	      (match3=GetTrackMatch(q,zpackage[2],segments[1],package[2],
				    match_id))!=NULL){
	    // Insert the segment from package 3 into the track
	    segments.push_back(match3);
	    
	    // remove the segment from the list 
	    package[2].erase(package[2].begin()+match_id);
	    
	    qsum+=match3->S(0,0)/fabs(match3->S(0,0));

	    // Redo the fit with the additional hits from package 3
	    DRiemannFit fit;
	    for (unsigned int m=0;m<segments.size();m++){
	      for(unsigned int n=0;n<segments[m]->hits.size();n++){
		DFDCPseudo *hit=segments[m]->hits[n];
		double covxx=hit->covxx;
		double covyy=hit->covyy;
		double covxy=hit->covxy;
		double hit_chi2=segments[segments.size()-1]->track[n].chi2;
		/*if (hit_chi2>1.){
		  covxx*=hit_chi2;      
		  covyy*=hit_chi2;
		covxy*=hit_chi2;
		}*/
		fit.AddHit(hit->x,hit->y,hit->wire->origin(2),covxx,covyy,covxy);
	      }
	    }
	    // Fake point at origin
	    fit.AddHit(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);

	    if (qsum>0) q=1.;
	    else if (qsum<0) q=-1.;
	    else q=0.;

	    if (fit.DoFit(segments[1]->rc)==NOERROR){     
	      // Charge
	      if (q==0) q=fit.q;
	      // Curvature
	      segments[2]->S(0,0)=kappa=q/2./fit.rc;
	      // Estimate for azimuthal angle
	      phi0=atan2(-fit.xc,fit.yc); 
	      if (q<0) phi0+=M_PI;
	      segments[2]->S(1,0)=phi0;
	      // remaining tracking parameters
	      tanl=fit.tanl;
	      zvertex=fit.zvertex;
	
	      segments[2]->S(3,0)=tanl;
	      segments[2]->S(4,0)=zvertex;
	      segments[2]->xc=fit.xc;
	      segments[2]->yc=fit.yc;
	      segments[2]->rc=fit.rc;
	      
	      // If we failed to match to package 4, try again.
	      if (match4==NULL && package[3].size()>0 && 
		  (match4=GetTrackMatch(q,zpackage[3],segments[2],package[3],
					match_id))!=NULL){
		// remove the segment from the list 
		package[3].erase(package[3].begin()+match_id); 
	      }
	    }
	    // Try to match to package 4 again.
	    if (match4==NULL && package[3].size()>0 &&
		(match4=GetTrackMatch(q,zpackage[3],segments[1],package[3],
				      match_id))!=NULL){ 
	      // remove the segment from the list
	      package[3].erase(package[3].begin()+match_id);
	    }
	  }
	}
      }

      DVector3 mom,mom2,pos;
      double Bx,By,Bz;
      int middle=segment->hits.size()/2;
      
      bfield->GetField(segment->hits[middle]->x,segment->hits[middle]->y,
		       segment->hits[middle]->wire->origin(2),Bx,By,Bz);
      double B=sqrt(Bx*Bx+By*By+Bz*Bz);
      double pt=0.003*B/2./fabs(kappa);
      double theta=M_PI_2-atan(tanl);

      mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);
      pos.SetXYZ(0,0,zvertex);

      // Create new track, starting with the current segment
      DTrackCandidate *track = new DTrackCandidate;
      track->setPosition(pos);
      track->setMomentum(mom);
      track->setCharge(q);
          
      _data.push_back(track); 
    }
  }

  // Next try to link segments starting at package 2
  if (package[1].size()>0 ){
    // Loop over segments in the 2nd package, matching them to segments in 
    // the third and fourth (most downstream) packages.
    for (unsigned int i=0;i<package[1].size();i++){
      DFDCSegment *segment=package[1][i];
      
      // tracking parameters
      double tanl=segment->S(3,0);
      double phi0=segment->S(1,0);
      double kappa=segment->S(0,0); 
      double zvertex=segment->S(4,0);
      // Sign of the charge
      double q=kappa/fabs(kappa);
      
      double qsum=q;
      
      // Start filling vector of segments belonging to current track    
      vector<DFDCSegment*>segments; 
      segments.push_back(segment);

      // Check that the tangent of the dip angle makes sense for FDC hits
      if (segment->S(3,0)<=0.0) continue;

      // Try matching to package 3
      if (package[2].size()>0 && 
	  (match3=GetTrackMatch(q,zpackage[2],segment,package[2],match_id))
	  !=NULL){
	// Insert the segment from package 3 into the track
	segments.push_back(match3);

	// remove the segment from the list 
	package[2].erase(package[2].begin()+match_id);

	qsum+=match3->S(0,0)/fabs(match3->S(0,0));

	// Try matching to package 4
	if (package[3].size()>0 && 
	    (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	     match_id))!=NULL){
	  // Insert the segment from package 4 into the track 
	  //segments.push_back(match4);

	  // remove the segment from the list 
	  package[3].erase(package[3].begin()+match_id);

	  qsum+=match4->S(0,0)/fabs(match4->S(0,0));	  
	}	
      }
      // No match in 3, try for 4
      else if (package[3].size()>0 && 
	       (match4=GetTrackMatch(q,zpackage[3],segment,package[3],
		match_id))!=NULL){
	// Insert the points in the segment from package 4 into the track 
	//segments.push_back(match4);

	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);

	qsum+=match4->S(0,0)/fabs(match4->S(0,0));
      }
      
      if (qsum>0) q=1.;
      else if (qsum<0) q=-1.;
      else q=0.;

      if (segments.size()>1){
	DRiemannFit fit;
	for (unsigned int m=0;m<segments.size();m++){
	  for (unsigned int n=0;n<segments[m]->hits.size();n++){
	    DFDCPseudo *hit=segments[m]->hits[n];
	    double covxx=hit->covxx;
	    double covyy=hit->covyy;
	    double covxy=hit->covxy;
	    double hit_chi2=segments[m]->track[n].chi2;
	    /*if (hit_chi2>1.){
	      covxx*=hit_chi2;      
	      covyy*=hit_chi2;
	      covxy*=hit_chi2;
	      }*/
	    fit.AddHit(hit->x,hit->y,hit->wire->origin(2),covxx,covyy,covxy);
	  }
	}
	// Fake point at origin
	fit.AddHit(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	if (fit.DoFit(segments[0]->rc)==NOERROR){
	  // Charge
	  if (q==0) q=fit.q;
	  // Curvature
	  segments[1]->S(0,0)=kappa=q/2./fit.rc;
	  // Estimate for azimuthal angle
	  phi0=atan2(-fit.xc,fit.yc); 
	  if (q<0) phi0+=M_PI;      
	  // remaining tracking parameters
	  tanl=fit.tanl;
	  zvertex=fit.zvertex;

	  segments[1]->S(3,0)=tanl;
	  segments[1]->S(4,0)=zvertex;
	  segments[1]->xc=fit.xc;
	  segments[1]->yc=fit.yc;
	  segments[1]->rc=fit.rc;
	  
	  // Try to match to package 4 again.
	  if (match4==NULL && package[3].size()>0 &&
	      (match4=GetTrackMatch(q,zpackage[3],segments[1],package[3],match_id))
	      !=NULL){
	    // remove the segment from the list 
	    package[3].erase(package[3].begin()+match_id);
	  }
	}
      }

      DVector3 mom,pos;      
      double Bx,By,Bz;
      int middle=segment->hits.size()/2;
      
      bfield->GetField(segment->hits[middle]->x,segment->hits[middle]->y,
		       segment->hits[middle]->wire->origin(2),Bx,By,Bz);
      double B=sqrt(Bx*Bx+By*By+Bz*Bz);
      double pt=0.003*B/2./fabs(kappa);
      double theta=M_PI_2-atan(tanl);

      mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);
      pos.SetXYZ(0,0,zvertex);

      // Create new track, starting with the current segment
      DTrackCandidate *track = new DTrackCandidate;
      track->setPosition(pos);
      track->setMomentum(mom);
      track->setCharge(q);
      
      _data.push_back(track); 
    }
  }
  
  // Next try to link segments starting at package 3
  if(package[2].size()>0){
    // Loop over segments in the 3rd package, matching them to segments in 
    // the fourth (most downstream) packages.
    for (unsigned int i=0;i<package[2].size();i++){
      DFDCSegment *segment=package[2][i];
 
      // tracking parameters
      double tanl=segment->S(3,0);
      double phi0=segment->S(1,0);
      double kappa=segment->S(0,0);
      double zvertex=segment->S(4,0);
      // Sign of the charge
      double q=kappa/fabs(kappa);
      
      // Start filling vector of segments belonging to current track    
      vector<DFDCSegment*>segments; 
      segments.push_back(segment);
      double qsum=q;
	
      // Check that the tangent of the dip angle makes sense for FDC hits
      if (segment->S(3,0)<=0.0) continue;
      
      // Try matching to package 4
      if (package[3].size()>0 && 
	  (match4=GetTrackMatch(q,zpackage[3],segment,package[3],match_id))
	  !=NULL){
	// Insert the segment from package 4 into the track 
	segments.push_back(match4);
	
	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);
	
	qsum+=match4->S(0,0)/fabs(match4->S(0,0));
      }	
      
      if (qsum>0) q=1.;
      else if (qsum<0) q=-1.;
      else q=0.;   

      if (segments.size()>1){
	DRiemannFit fit;
	for (unsigned int m=0;m<segments.size();m++){
	  for (unsigned int n=0;n<segments[m]->hits.size();n++){
	    DFDCPseudo *hit=segments[m]->hits[n];
	    double covxx=hit->covxx;
	    double covyy=hit->covyy;
	    double covxy=hit->covxy;
	    double hit_chi2=segments[m]->track[n].chi2;
	    /*if (hit_chi2>1.){
	      covxx*=hit_chi2;      
	      covyy*=hit_chi2;
	      covxy*=hit_chi2;
	      }*/
	    fit.AddHit(hit->x,hit->y,hit->wire->origin(2),covxx,covyy,covxy);
	  }
	}
	// Fake point at origin
	fit.AddHit(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	if (fit.DoFit(segments[0]->rc)==NOERROR){     	
	  // Charge
	  if (q==0) q=fit.q;
	  // Curvature
	  kappa=q/2./fit.rc;
	  // Estimate for azimuthal angle
	  phi0=atan2(-fit.xc,fit.yc); 
	  if (q<0) phi0+=M_PI; 
	  // remaining tracking parameters
	  tanl=fit.tanl;
	  zvertex=fit.zvertex;
	}
      }
      
    
      DVector3 mom,pos;
      double Bx,By,Bz;
      int middle=segment->hits.size()/2;
   
      bfield->GetField(segment->hits[middle]->x,segment->hits[middle]->y,
		       segment->hits[middle]->wire->origin(2),Bx,By,Bz);
      double B=sqrt(Bx*Bx+By*By+Bz*Bz);
      double pt=0.003*B/2./fabs(kappa);
      double theta=M_PI_2-atan(tanl);

      mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);
      pos.SetXYZ(0,0,zvertex);
      
      // Create new track, starting with the current segment
      DTrackCandidate *track = new DTrackCandidate;
      track->setPosition(pos);
      track->setMomentum(mom);
      track->setCharge(q);
     
      _data.push_back(track); 
    }
  }

  // Now collect stray segments in package 4
  for (unsigned int k=0;k<package[3].size();k++){
    DFDCSegment *segment=package[3][k];
    
    DVector3 pos,mom;
    double kappa=segment->S(0,0);
    double tanl=segment->S(3,0);
    double phi0=segment->S(1,0);
    double Bx,By,Bz;
    int middle=segment->hits.size()/2;
    
    bfield->GetField(segment->hits[middle]->x,segment->hits[middle]->y,
		     segment->hits[middle]->wire->origin(2),Bx,By,Bz);
    double B=sqrt(Bx*Bx+By*By+Bz*Bz);
    double pt=0.003*B/2./fabs(kappa);
    double theta=M_PI_2-atan(tanl);
    
    mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);
    pos.SetXYZ(0,0,segment->S(4,0));
     
    // Create new track, starting with the current segment
    DTrackCandidate *track = new DTrackCandidate;
    track->setPosition(pos);
    track->setMomentum(mom);    
    track->setCharge(segment->S(0,0)/fabs(segment->S(0,0)));

    _data.push_back(track); 
  }
  

  return NOERROR;
}

// Swim track from one package to the next and look for a match to a segment
// in the new package
DFDCSegment *DTrackCandidate_factory_FDCCathodes::GetTrackMatch(double q,double z, 
						 DFDCSegment *segment,
						 vector<DFDCSegment*>package,
						 unsigned int &match_id){
  DFDCSegment *match=NULL;
  DVector3 norm;  // normal to FDC planes
  norm.SetXYZ(0.,0.,1.);

  // Initialize the stepper 
  DMagneticFieldStepper stepper(bfield, q); 

  // Get the position and momentum at the exit of the package for the 
  // current segment
  DVector3 pos,mom,origin(0.,0.,z);
  if (GetPositionAndMomentum(segment,pos,mom)!=NOERROR) return NULL;

  // Match to the next package by swimming the track through the field
  double diff_min=1000.,diff;
  double s;
  if (stepper.SwimToPlane(pos,mom,origin,norm,&s)==false){
    for (unsigned int j=0;j<package.size();j++){
      DFDCSegment *segment2=package[j];
	  
      // Skip to next segment if the sign of the charge is wrong
      //double kappa=segment2->S(0,0);
      //double q2=kappa/fabs(kappa);

      // if (q2!=q) continue;
     
      // Check that the tangent of the dip angle makes sense for FDC hits
      if (segment2->S(3,0)<=0.0) continue;
 
      double x2=segment2->hits[segment2->hits.size()-1]->x;
      double y2=segment2->hits[segment2->hits.size()-1]->y;
      diff=sqrt((pos(0)-x2)*(pos(0)-x2)+(pos(1)-y2)*(pos(1)-y2));

      if (diff<diff_min&&diff<MATCH_RADIUS(mom.Mag())){
	diff_min=diff;
	match=segment2;
	match_id=j;
      }
    }
  }
  
  // If matching in the forward direction did not work, try swimming and
 // matching backwards...
  if (match==NULL){
    diff_min=1000.;
    for (unsigned int i=0;i<package.size();i++){
      DFDCSegment *segment2=package[i];
      if (GetPositionAndMomentum(segment2,pos,mom)==NOERROR){
        mom=-1.0*mom;
        origin(2)=segment->hits[0]->wire->origin(2);
        if (stepper.SwimToPlane(pos,mom,origin,norm,&s)==false){
          double x2=segment->hits[0]->x;
          double y2=segment2->hits[0]->y;
          diff=sqrt((pos(0)-x2)*(pos(0)-x2)+(pos(1)-y2)*(pos(1)-y2));
          if (diff<diff_min&&diff<MATCH_RADIUS(mom.Mag())){
	    diff_min=diff;
	    match=segment2;
            match_id=i;
          }
        }	
      }       
    }
  }
  if(DEBUG_HISTS){
    match_dist->Fill(mom.Mag(),diff_min);
  }
  return match;
}


// Obtain position and momentum at the exit of a given package using the 
// helical track model.
//
jerror_t DTrackCandidate_factory_FDCCathodes::GetPositionAndMomentum(DFDCSegment *segment,
					      DVector3 &pos, DVector3 &mom){
  // Position of track segment at last hit plane of package
  double x=segment->xc+segment->rc*cos(segment->Phi1);
  double y=segment->yc+segment->rc*sin(segment->Phi1);
  double z=segment->hits[0]->wire->origin(2);
 
  // Make sure that the position makes sense!
  if (sqrt(x*x+y*y)>FDC_OUTER_RADIUS) return VALUE_OUT_OF_RANGE;

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


//------------------
// toString
//------------------
const string DTrackCandidate_factory_FDCCathodes::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	printheader("row:    x:     y:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		//DTrackCandidate *myDTrackCandidate = _data[i];
	
		printnewrow();
		printcol("%d",	i);
//		printcol("%1.3f",	myDTrackCandidate->x);
//		printcol("%3.2f",	myDTrackCandidate->y);
		printrow();
	}

	return _table;
}

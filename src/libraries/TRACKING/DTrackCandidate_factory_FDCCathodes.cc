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
#include "DHelicalFit.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>

//#define MATCH_RADIUS(p) (2.79+2.88/(p)/(p))
#define MATCH_RADIUS(p) (2.75+2./(p+0.25)/(p+0.25))
#define MAX_SEGMENTS 20
#define HALF_PACKAGE 6.0
#define FDC_OUTER_RADIUS 50.0 
#define BEAM_VAR 1.0 // cm^2
#define HIT_CHI2_CUT 10.0
#define Z_VERTEX 65.0
#define Z_MAX 85.0
#define Z_MIN 45.0

///
/// DTrackCandidate_factory_FDCCathodes::brun():
///
jerror_t DTrackCandidate_factory_FDCCathodes::brun(JEventLoop* eventLoop, 
						   int runnumber) {
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  dgeom->GetFDCZ(z_wires);

  if(DEBUG_HISTS){
    dapp->Lock();
    match_dist_fdc=(TH2F*)gROOT->FindObject("match_dist_fdc");
    if (!match_dist_fdc) 
      match_dist_fdc=new TH2F("match_dist_fdc",
		  "Matching distance for connecting FDC segments",
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
  // abort if there are no segments
  if (segments.size()==0.) return NOERROR;

  std::sort(segments.begin(), segments.end(), DTrackCandidate_segment_cmp);

  // Group segments by package
  vector<DFDCSegment*>package[4]; 
  for (unsigned int i=0;i<segments.size();i++){
     const DFDCSegment *segment=segments[i];
 package[(segment->hits[0]->wire->layer-1)/6].push_back((DFDCSegment*)segment);
  }
      
  double zpackage[4];  // z-positions of entrances to FDC packages 
  zpackage[0]=z_wires[0];
  zpackage[1]=z_wires[6];
  zpackage[2]=z_wires[12];
  zpackage[3]=z_wires[18];
  DFDCSegment *match2=NULL;
  DFDCSegment *match3=NULL;
  DFDCSegment *match4=NULL;
  unsigned int match_id=0;
  double xc,yc,rc;

  // Bail if there are too many segments
  if (package[0].size()+package[1].size()+package[2].size()
      +package[3].size()>MAX_SEGMENTS)
    return UNRECOVERABLE_ERROR; 

  // First deal with tracks with segments in the first package
  vector<int>pack1_matched(package[0].size());
  if (package[0].size()>0){
    // Loop over segments in the first package, matching them to segments in 
    // the second, third, and fourth (most downstream) packages.
    for (unsigned int i=0;i<package[0].size();i++){
      DFDCSegment *segment=package[0][i];
      match2=match3=match4=NULL;
      
      // Tracking parameters from first segment
      double tanl=segment->S(3,0);
      double phi0=segment->S(1,0);
      double z_vertex=segment->S(4,0);
      xc=segment->xc;
      yc=segment->yc;
      rc=segment->rc;
      // Sign of the charge
      double q=segment->S(0,0)/fabs(segment->S(0,0));
      
      double qsum=q;
      // Start filling vector of segments belonging to current track    
      vector<DFDCSegment*>mysegments; 
      mysegments.push_back(segment);
 
      // Try matching to package 2
      if (package[1].size()>0 && 
	  (match2=GetTrackMatch(q,zpackage[1],segment,package[1],match_id))
	  !=NULL){
	// Insert the segment from package 2 into the track 
	mysegments.push_back(match2);
	
	// remove the segment from the list 
	package[1].erase(package[1].begin()+match_id);

	qsum+=match2->S(0,0)/fabs(match2->S(0,0));

	// Try matching to package 3
	if (package[2].size()>0 && 
	    (match3=GetTrackMatch(q,zpackage[2],match2,package[2],match_id))
	    !=NULL){
	  // Insert the segment from package 3 into the track 
	  mysegments.push_back(match3);

	  // remove the segment from the list 
	  package[2].erase(package[2].begin()+match_id);

	  qsum+=match3->S(0,0)/fabs(match3->S(0,0));

	  // Try matching to package 4
	  if (package[3].size()>0 && 
	      (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	       match_id))!=NULL){
	    // Insert the segment from package 4 into the track 
	    mysegments.push_back(match4);

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
	  mysegments.push_back(match4);

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
	mysegments.push_back(match3);

	// remove the segment from the list 
	package[2].erase(package[2].begin()+match_id);

	qsum+=match3->S(0,0)/fabs(match3->S(0,0));
	
	// Try matching to package 4
	if (package[3].size()>0 && 
	    (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	     match_id))!=NULL){
	  // Insert the segment from package 4 into the track 
	  mysegments.push_back(match4);

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
	mysegments.push_back(match4);
	
	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);

	qsum+=match4->S(0,0)/fabs(match4->S(0,0));
      }

      if (qsum>0) q=1.;
      else if (qsum<0) q=-1.;
      else q=0.;

      // Variables for determining average Bz
      double Bz_avg=0.,Bx,By,Bz;
      unsigned int num_hits=segment->hits.size();

      if (mysegments.size()>1){
	pack1_matched[i]=1;

	DHelicalFit fit;
	if (segment){ 
	  for (unsigned int n=0;n<segment->hits.size();n++){
	    double covxx=segment->hits[n]->covxx;
	    double covyy=segment->hits[n]->covyy;
	    double covxy=segment->hits[n]->covxy;
	    double x=segment->hits[n]->x;
	    double y=segment->hits[n]->y;
	    double z=segment->hits[n]->wire->origin(2);
	    fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	    bfield->GetField(x,y,z,Bx,By,Bz);
	    Bz_avg-=Bz;
	  }
	}
	if (match2){
	  for (unsigned int n=0;n<match2->hits.size();n++){
	    double covxx=match2->hits[n]->covxx;
	    double covyy=match2->hits[n]->covyy;
	    double covxy=match2->hits[n]->covxy;
	    double x=match2->hits[n]->x;
	    double y=match2->hits[n]->y;
	    double z=match2->hits[n]->wire->origin(2);
	    fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	    bfield->GetField(x,y,z,Bx,By,Bz);
	    Bz_avg-=Bz;	    
	  }
	  num_hits+=match2->hits.size();
	}
	if (match3){
	  for (unsigned int n=0;n<match3->hits.size();n++){
	    double covxx=match3->hits[n]->covxx;
	    double covyy=match3->hits[n]->covyy;
	    double covxy=match3->hits[n]->covxy;
	    double x=match3->hits[n]->x;
	    double y=match3->hits[n]->y;
	    double z=match3->hits[n]->wire->origin(2);
	    fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	    bfield->GetField(x,y,z,Bx,By,Bz);
	    Bz_avg-=Bz;
	  }
	  num_hits+=match3->hits.size();
	}
	if (match4){
	  for (unsigned int n=0;n<match4->hits.size();n++){
	    double covxx=match4->hits[n]->covxx;
	    double covyy=match4->hits[n]->covyy;
	    double covxy=match4->hits[n]->covxy;
	    double x=match4->hits[n]->x;
	    double y=match4->hits[n]->y;
	    double z=match4->hits[n]->wire->origin(2);
	    fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	    bfield->GetField(x,y,z,Bx,By,Bz);
	    Bz_avg-=Bz;
	  }
	  num_hits+=match4->hits.size();
	}
	// Fake point at origin
	fit.AddHitXYZ(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	if (fit.FitTrackRiemann(mysegments[0]->rc)==NOERROR){      
	  // Charge
	  //if (q==0) 
	  q=fit.q;
	  // Curvature
	  mysegments[1]->S(0,0)=q/2./fit.r0;
	  // Estimate for azimuthal angle
	  phi0=fit.phi; 
	  mysegments[1]->S(1,0)=phi0;
	  // remaining tracking parameters
	  tanl=fit.tanl;
	  z_vertex=fit.z_vertex;
	
	  mysegments[1]->S(3,0)=tanl;
	  mysegments[1]->S(4,0)=z_vertex;
	  xc=mysegments[1]->xc=fit.x0;
	  yc=mysegments[1]->yc=fit.y0;
	  rc=mysegments[1]->rc=fit.r0;

	  //printf("Match %x %x %x %x\n",segment,match2,match3,match4);
	  //printf("xc %f yc %f rc %f\n",xc,yc,rc);
	  
	  // Try to match to package 2 again.
	  if (match2==NULL && package[1].size()>0 &&
	      (match2=GetTrackMatch(q,zpackage[1],mysegments[1],package[1],
				    match_id))!=NULL){ 
	    // Insert the segment from package 2 into the track
	    mysegments.push_back(match2);

	    // remove the segment from the list 
	    package[1].erase(package[1].begin()+match_id);

	    qsum+=match2->S(0,0)/fabs(match2->S(0,0));

	    // Redo the fit with the additional hits from package 2
	    for (unsigned int n=0;n<match2->hits.size();n++){
	      double covxx=match2->hits[n]->covxx;
	      double covyy=match2->hits[n]->covyy;
	      double covxy=match2->hits[n]->covxy;
	      double x=match2->hits[n]->x;
	      double y=match2->hits[n]->y;
	      double z=match2->hits[n]->wire->origin(2);
	      fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	      bfield->GetField(x,y,z,Bx,By,Bz);
	      Bz_avg-=Bz;
	    }
	    num_hits+=match2->hits.size();
	    
	    if (qsum>0) q=1.;
	    else if (qsum<0) q=-1.;
	    else q=0.;
	    
	    if (fit.FitTrackRiemann(mysegments[1]->rc)==NOERROR){     
	      // Charge
	      //if (q==0) 
	      q=fit.q;
	      // Curvature
	      mysegments[1]->S(0,0)=q/2./fit.r0;
	      // Estimate for azimuthal angle
	      phi0=fit.phi;
	      mysegments[1]->S(1,0)=phi0;
	      // remaining tracking parameters
	      tanl=fit.tanl;
	      z_vertex=fit.z_vertex;
	
	      mysegments[1]->S(3,0)=tanl;
	      mysegments[1]->S(4,0)=z_vertex;
	      xc=mysegments[1]->xc=fit.x0;
	      yc=mysegments[1]->yc=fit.y0;
	      rc=mysegments[1]->rc=fit.r0;
	      
	      //printf("Match %x %x %x %x\n",segment,match2,match3,match4);
	      //printf("xc %f yc %f rc %f\n",xc,yc,rc);
	  
	    }
	  }
	  
	  // Try to match to package 3 again.
	  if (match3==NULL && package[2].size()>0 &&
	      (match3=GetTrackMatch(q,zpackage[2],mysegments[1],package[2],
				    match_id))!=NULL){
	    // Insert the segment from package 3 into the track
	    mysegments.push_back(match3);

	    // remove the segment from the list 
	    package[2].erase(package[2].begin()+match_id);
	    
	    qsum+=match3->S(0,0)/fabs(match3->S(0,0));

	    // Redo the fit with the additional hits from package 3
	    for (unsigned int n=0;n<match3->hits.size();n++){
	      double covxx=match3->hits[n]->covxx;
	      double covyy=match3->hits[n]->covyy;
	      double covxy=match3->hits[n]->covxy;
	      double x=match3->hits[n]->x;
	      double y=match3->hits[n]->y;
	      double z=match3->hits[n]->wire->origin(2);
	      fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	      bfield->GetField(x,y,z,Bx,By,Bz);
	      Bz_avg-=Bz;
	    }
	    num_hits+=match3->hits.size();
	    
	    if (qsum>0) q=1.;
	    else if (qsum<0) q=-1.;
	    else q=0.;

	    if (fit.FitTrackRiemann(mysegments[1]->rc)==NOERROR){     
	      // Charge
	      //if (q==0) 
	      q=fit.q;
	      // Curvature
	      mysegments[2]->S(0,0)=q/2./fit.r0;
	      // Estimate for azimuthal angle
	      phi0=fit.phi;
	      mysegments[2]->S(1,0)=phi0;
	      // remaining tracking parameters
	      tanl=fit.tanl;
	      z_vertex=fit.z_vertex;
	
	      mysegments[2]->S(3,0)=tanl;
	      mysegments[2]->S(4,0)=z_vertex;
	      xc=mysegments[2]->xc=fit.x0;
	      yc=mysegments[2]->yc=fit.y0;
	      rc=mysegments[2]->rc=fit.r0;
	      
	      //printf("Match %x %x %x %x\n",segment,match2,match3,match4);
	      //printf("xc %f yc %f rc %f\n",xc,yc,rc);
	      
	      // If we failed to match to package 4, try again.
	      if (match4==NULL && package[3].size()>0 && 
		  (match4=GetTrackMatch(q,zpackage[3],mysegments[2],package[3],
					match_id))!=NULL){
		// Insert the segment from package 4 into the track
		mysegments.push_back(match4);

		// remove the segment from the list 
		package[3].erase(package[3].begin()+match_id); 
	      }
	    }
	    // Try to match to package 4 again.
	    if (match4==NULL && package[3].size()>0 &&
		(match4=GetTrackMatch(q,zpackage[3],mysegments[1],package[3],
				      match_id))!=NULL){ 
	      
	      // Insert the segment from package 4 into the track
	      mysegments.push_back(match4);

	      // remove the segment from the list
	      package[3].erase(package[3].begin()+match_id);
	    }
	  }
	}
      
	DVector3 mom,pos;
	Bz_avg/=double(num_hits);

	// Try to fix tracks that appear to point upstream
	if (tanl<0){ 
	  double x=segment->hits[0]->x;
	  double y=segment->hits[0]->y;
	  double ratio=sqrt(x*x+y+y)/2./rc;
	  if (ratio<1.){
	    double sperp=2.*rc*asin(ratio);
	    tanl=(segment->hits[0]->wire->origin.z()-Z_VERTEX)/sperp;
	    z_vertex=Z_VERTEX;
	  }
	}	
	double pt=0.003*Bz_avg*rc;
	double theta=M_PI_2-atan(tanl);
	double d=yc-q*rc*cos(phi0);
	
	mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);	  
	pos.SetXYZ(-d*sin(phi0),d*cos(phi0),z_vertex);
	
	if (z_vertex<0)
	  {
	  double sperp=(Z_VERTEX-z_vertex)/tanl;
	  double kappa=q/2./rc;
	  double sin2ks=sin(2.*kappa*sperp);
	  double cos2ks=cos(2.*kappa*sperp);
	  pos(0)=pos.x()+cos(phi0)*sin2ks/2./kappa
	    -sin(phi0)*(1.-cos2ks)/2./kappa;
	  pos(1)=pos.y()+sin(phi0)*sin2ks/2./kappa
	    +cos(phi0)*(1.-cos2ks)/2./kappa;
	  pos(2)=Z_VERTEX;
	  mom(0)=pt*(cos(phi0)*cos2ks-sin(phi0)*sin2ks);
	  mom(1)=pt*(sin(phi0)*cos2ks+cos(phi0)*sin2ks);
	}
		
	// Create new track, starting with the current segment
	DTrackCandidate *track = new DTrackCandidate;
	track->setPosition(pos);
	track->setMomentum(mom);
	track->setCharge(q);
	
	for (unsigned int m=0;m<mysegments.size();m++)
	  track->AddAssociatedObject(mysegments[m]);
	
	_data.push_back(track); 
	
      }
    }
  }

  // Prune segments in package 1 that have been matched to other segments
  vector<DFDCSegment*>pack1_left_over;     
  for (unsigned int i=0;i<package[0].size();i++){
    if (pack1_matched[i]!=1) pack1_left_over.push_back(package[0][i]);
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
      double z_vertex=segment->S(4,0);
      xc=segment->xc;
      yc=segment->yc;
      rc=segment->rc;
      // Sign of the charge
      double q=segment->S(0,0)/fabs(segment->S(0,0));
      
      double qsum=q;
      
      // Start filling vector of segments belonging to current track    
      vector<DFDCSegment*>mysegments; 
      mysegments.push_back(segment);

      // Try matching to package 3
      if (package[2].size()>0 && 
	  (match3=GetTrackMatch(q,zpackage[2],segment,package[2],match_id))
	  !=NULL){
	// Insert the segment from package 3 into the track
	mysegments.push_back(match3);

	// remove the segment from the list 
	package[2].erase(package[2].begin()+match_id);

	qsum+=match3->S(0,0)/fabs(match3->S(0,0));

	// Try matching to package 4
	if (package[3].size()>0 && 
	    (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	     match_id))!=NULL){
	  // Insert the segment from package 4 into the track 
	  mysegments.push_back(match4);

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
	mysegments.push_back(match4);

	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);

	qsum+=match4->S(0,0)/fabs(match4->S(0,0));
      }
      
      if (qsum>0) q=1.;
      else if (qsum<0) q=-1.;
      else q=0.;

      // Variables for determining average Bz
      double Bz_avg=0.,Bx,By,Bz;
      unsigned int num_hits=segment->hits.size();

      if (mysegments.size()>1){
	DHelicalFit fit;
	if (segment){ 
	  for (unsigned int n=0;n<segment->hits.size();n++){
	    double covxx=segment->hits[n]->covxx;
	    double covyy=segment->hits[n]->covyy;
	    double covxy=segment->hits[n]->covxy;
	    double x=segment->hits[n]->x;
	    double y=segment->hits[n]->y;
	    double z=segment->hits[n]->wire->origin(2);
	    fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	    bfield->GetField(x,y,z,Bx,By,Bz);
	    Bz_avg-=Bz;
	  }
	}
	if (match3){
	  for (unsigned int n=0;n<match3->hits.size();n++){
	    double covxx=match3->hits[n]->covxx;
	    double covyy=match3->hits[n]->covyy;
	    double covxy=match3->hits[n]->covxy;
	    double x=match3->hits[n]->x;
	    double y=match3->hits[n]->y;
	    double z=match3->hits[n]->wire->origin(2);
	    fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	    bfield->GetField(x,y,z,Bx,By,Bz);
	    Bz_avg-=Bz;	    
	  }
	  num_hits+=match3->hits.size();
	}
	if (match4){
	  for (unsigned int n=0;n<match4->hits.size();n++){
	    double covxx=2.*match4->hits[n]->covxx;
	    double covyy=2.*match4->hits[n]->covyy;
	    double covxy=2.*match4->hits[n]->covxy;
	    double x=match4->hits[n]->x;
	    double y=match4->hits[n]->y;
	    double z=match4->hits[n]->wire->origin(2);
	    fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);  
	    bfield->GetField(x,y,z,Bx,By,Bz);
	    Bz_avg-=Bz;
	  }
	  num_hits+=match4->hits.size();
	}
	// Fake point at origin
	fit.AddHitXYZ(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	if (fit.FitTrackRiemann(mysegments[0]->rc)==NOERROR){
	  // Charge
	  //if (q==0) 
	  q=fit.q;
	  // Curvature
	  mysegments[1]->S(0,0)=q/2./fit.r0;
	  // Estimate for azimuthal angle
	  phi0=fit.phi;
	  mysegments[1]->S(1,0)=phi0;
	  // remaining tracking parameters
	  tanl=fit.tanl;
	  z_vertex=fit.z_vertex;

	  mysegments[1]->S(3,0)=tanl;
	  mysegments[1]->S(4,0)=z_vertex;
	  xc=mysegments[1]->xc=fit.x0;
	  yc=mysegments[1]->yc=fit.y0;
	  rc=mysegments[1]->rc=fit.r0;
	  
	  // Try to match to package 4 again.
	  if (match4==NULL && package[3].size()>0 &&
	      (match4=GetTrackMatch(q,zpackage[3],mysegments[1],package[3],match_id))
	      !=NULL){
	    // Insert the segment from package 4 into the track
	    mysegments.push_back(match4);

	    // remove the segment from the list 
	    package[3].erase(package[3].begin()+match_id);
	  }

	  // Try to match to package 1 again	  
	  if (pack1_left_over.size()>0){
	    DFDCSegment *match1=NULL;
	    if ((match1=GetTrackMatch(q,zpackage[0],mysegments[1],
				      pack1_left_over,
				      match_id))!=NULL){
	      mysegments.push_back(match1);
	      pack1_left_over.erase(pack1_left_over.begin()+match_id);

	      // Refit with additional hits from package 1
	      for (unsigned int n=0;n<match1->hits.size();n++){
		double covxx=match1->hits[n]->covxx;
		double covyy=match1->hits[n]->covyy;
		double covxy=match1->hits[n]->covxy;
		double x=match1->hits[n]->x;
		double y=match1->hits[n]->y;
		double z=match1->hits[n]->wire->origin(2);
		fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
		bfield->GetField(x,y,z,Bx,By,Bz);
		Bz_avg-=Bz;
	      }
	      num_hits+=match1->hits.size();

	      if (fit.FitTrackRiemann(mysegments[1]->rc)==NOERROR){
		// Charge
		q=fit.q;
		// Curvature
		mysegments[1]->S(0,0)=q/2./fit.r0;
		// Estimate for azimuthal angle
		phi0=fit.phi;
		mysegments[1]->S(1,0)=phi0;
		// remaining tracking parameters
		tanl=fit.tanl;
		z_vertex=fit.z_vertex;
		
		mysegments[1]->S(3,0)=tanl;
		mysegments[1]->S(4,0)=z_vertex;
		xc=mysegments[1]->xc=fit.x0;
		yc=mysegments[1]->yc=fit.y0;
		rc=mysegments[1]->rc=fit.r0;
	      }
	    }
	  }
	}
      } // if mysegments.size()>1
      else{
	for (unsigned int n=0;n<segment->hits.size();n++){
	  double x=segment->hits[n]->x;
	  double y=segment->hits[n]->y;
	  double z=segment->hits[n]->wire->origin(2);
	  bfield->GetField(x,y,z,Bx,By,Bz);
	  Bz_avg-=Bz;
	}
      }

      DVector3 mom,pos;
      Bz_avg/=double(num_hits);
      
      // Try to fix tracks that appear to point upstream
      if (tanl<0){
	double x=segment->hits[0]->x;
	double y=segment->hits[0]->y;
	double ratio=sqrt(x*x+y+y)/2./rc;
	if (ratio<1.){
	  double sperp=2.*rc*asin(ratio);
	  tanl=(segment->hits[0]->wire->origin.z()-Z_VERTEX)/sperp;
	  z_vertex=Z_VERTEX;
	}
      }	
      double pt=0.003*Bz_avg*rc;
      double theta=M_PI_2-atan(tanl);
      double d=yc-q*rc*cos(phi0);

      mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);
      pos.SetXYZ(-d*sin(phi0),d*cos(phi0),z_vertex);

      if (z_vertex<0){
	double sperp=(Z_VERTEX-z_vertex)/tanl;
	double kappa=q/2./rc;
	double sin2ks=sin(2.*kappa*sperp);
	double cos2ks=cos(2.*kappa*sperp);
	pos(0)=pos.x()+cos(phi0)*sin2ks/2./kappa
	  -sin(phi0)*(1.-cos2ks)/2./kappa;
	pos(1)=pos.y()+sin(phi0)*sin2ks/2./kappa
	  +cos(phi0)*(1.-cos2ks)/2./kappa;
	pos(2)=Z_VERTEX;
	mom(0)=pt*(cos(phi0)*cos2ks-sin(phi0)*sin2ks);
	mom(1)=pt*(sin(phi0)*cos2ks+cos(phi0)*sin2ks);
      }

      // Create new track, starting with the current segment
      DTrackCandidate *track = new DTrackCandidate;
      track->setPosition(pos);
      track->setMomentum(mom);
      track->setCharge(q);

      for (unsigned int m=0;m<mysegments.size();m++)
	track->AddAssociatedObject(mysegments[m]);

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
      double z_vertex=segment->S(4,0);  
      xc=segment->xc;
      yc=segment->yc;
      rc=segment->rc;
      // Sign of the charge
      double q=segment->S(0,0)/fabs(segment->S(0,0));

      // Start filling vector of segments belonging to current track    
      vector<DFDCSegment*>mysegments; 
      mysegments.push_back(segment);
      double qsum=q;
      
      // Try matching to package 4
      if (package[3].size()>0 && 
	  (match4=GetTrackMatch(q,zpackage[3],segment,package[3],match_id))
	  !=NULL){
	// Insert the segment from package 4 into the track 
	mysegments.push_back(match4);
	
	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);
	
	qsum+=match4->S(0,0)/fabs(match4->S(0,0));
      }	
      
      if (qsum>0) q=1.;
      else if (qsum<0) q=-1.;
      else q=0.;   

      // Variables for determining average Bz
      double Bz_avg=0.,Bx,By,Bz;
      unsigned int num_hits=segment->hits.size();
          
      if (mysegments.size()>1){
	DHelicalFit fit;
	for (unsigned int m=0;m<mysegments.size();m++){
	  for (unsigned int n=0;n<mysegments[m]->hits.size();n++){
	    DFDCPseudo *hit=mysegments[m]->hits[n];
	    double covxx=hit->covxx;
	    double covyy=hit->covyy;
	    double covxy=hit->covxy;
	    //double hit_chi2=mysegments[m]->track[n].chi2;
	    /*if (hit_chi2>1.){
	      covxx*=hit_chi2;      
	      covyy*=hit_chi2;
	      covxy*=hit_chi2;
	      }*/
	    fit.AddHitXYZ(hit->x,hit->y,hit->wire->origin(2),covxx,covyy,covxy);
	    bfield->GetField(hit->x,hit->y,hit->wire->origin(2),Bx,By,Bz);
	    Bz_avg-=Bz;	    
	  }
	  num_hits+=mysegments[m]->hits.size();
	}
	// Fake point at origin
	fit.AddHitXYZ(0.,0.,Z_VERTEX,BEAM_VAR,BEAM_VAR,0.);
	if (fit.FitTrackRiemann(mysegments[0]->rc)==NOERROR){     	
	  // Charge
	  //if (q==0) 
	  q=fit.q;
	  // Estimate for azimuthal angle
	  phi0=fit.phi;
	  mysegments[0]->S(1,0)=phi0;
	  // remaining tracking parameters
	  tanl=fit.tanl;
	  z_vertex=fit.z_vertex;

	  mysegments[0]->S(3,0)=tanl;
	  mysegments[0]->S(4,0)=z_vertex;
	  xc=mysegments[0]->xc=fit.x0;
	  yc=mysegments[0]->yc=fit.y0;
	  rc=mysegments[0]->rc=fit.r0;
	}
      
	// Try to match to package 2 again.
	if (package[1].size()>0 &&
	    (match2=GetTrackMatch(q,zpackage[1],mysegments[0],package[1],
				  match_id))!=NULL){
	  // Insert the segment from package 2 into the track 
	  mysegments.push_back(match2);
	  
	  // remove the segment from the list 
	  package[1].erase(package[1].begin()+match_id);
	  
	  // Redo the fit with the additional hits from package 2
	  for (unsigned int n=0;n<match2->hits.size();n++){
	    double covxx=match2->hits[n]->covxx;
	    double covyy=match2->hits[n]->covyy;
	    double covxy=match2->hits[n]->covxy;
	    double x=match2->hits[n]->x;
	    double y=match2->hits[n]->y;
	    double z=match2->hits[n]->wire->origin(2);
	    fit.AddHitXYZ(x,y,z,covxx,covyy,covxy);
	    bfield->GetField(x,y,z,Bx,By,Bz);
	    Bz_avg-=Bz;	    
	  }
	  num_hits+=match2->hits.size();

	  if (fit.FitTrackRiemann(mysegments[0]->rc)==NOERROR){     
	    // Charge
	    //if (q==0) 
	    q=fit.q;
	    // Curvature
	    mysegments[0]->S(0,0)=q/2./fit.r0;
	    // Estimate for azimuthal angle
	    phi0=fit.phi;
	    mysegments[0]->S(1,0)=phi0;
	    // remaining tracking parameters
	    tanl=fit.tanl;
	    z_vertex=fit.z_vertex;
	    
	    mysegments[0]->S(3,0)=tanl;
	    mysegments[0]->S(4,0)=z_vertex;
	    xc=mysegments[0]->xc=fit.x0;
	    yc=mysegments[0]->yc=fit.y0;
	    rc=mysegments[0]->rc=fit.r0;
	  }
	}

	// Try to match to package 1 again	  
	if (pack1_left_over.size()>0){
	  DFDCSegment *match1=NULL;
	  if ((match1=GetTrackMatch(q,zpackage[0],mysegments[0],
				    pack1_left_over,
				    match_id))!=NULL){
	    mysegments.push_back(match1);
	    pack1_left_over.erase(pack1_left_over.begin()+match_id);
	  }
	}
      } // if mysegments.size()>1
      else{
	for (unsigned int n=0;n<segment->hits.size();n++){
	  double x=segment->hits[n]->x;
	  double y=segment->hits[n]->y;
	  double z=segment->hits[n]->wire->origin(2);
	  bfield->GetField(x,y,z,Bx,By,Bz);
	  Bz_avg-=Bz;
	}
      }
	
      DVector3 mom,pos;
      Bz_avg/=double(num_hits);

      // Try to fix tracks that appear to point downstream
      if (tanl<0){ 
	double x=segment->hits[0]->x;
	double y=segment->hits[0]->y;
	double ratio=sqrt(x*x+y+y)/2./rc;
	if (ratio<1.){
	  double sperp=2.*rc*asin(ratio);
	  tanl=(segment->hits[0]->wire->origin.z()-Z_VERTEX)/sperp;
	  z_vertex=Z_VERTEX;
	}
      }	
      double pt=0.003*Bz_avg*rc;
      double theta=M_PI_2-atan(tanl);
      double d=yc-q*rc*cos(phi0);

      mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);
      pos.SetXYZ(-d*sin(phi0),d*cos(phi0),z_vertex);

      if (z_vertex<0){
	double sperp=(Z_VERTEX-z_vertex)/tanl;
	double kappa=q/2./rc;
	double sin2ks=sin(2.*kappa*sperp);
	double cos2ks=cos(2.*kappa*sperp);
	pos(0)=pos.x()+cos(phi0)*sin2ks/2./kappa
	  -sin(phi0)*(1.-cos2ks)/2./kappa;
	pos(1)=pos.y()+sin(phi0)*sin2ks/2./kappa
	  +cos(phi0)*(1.-cos2ks)/2./kappa;
	pos(2)=Z_VERTEX;
	mom(0)=pt*(cos(phi0)*cos2ks-sin(phi0)*sin2ks);
	mom(1)=pt*(sin(phi0)*cos2ks+cos(phi0)*sin2ks);
      }
      
      // Create new track, starting with the current segment
      DTrackCandidate *track = new DTrackCandidate;
      track->setPosition(pos);
      track->setMomentum(mom);      
      track->setCharge(q);
      
      for (unsigned int m=0;m<mysegments.size();m++)
	track->AddAssociatedObject(mysegments[m]);
      
      _data.push_back(track); 
    }
  }

  // Now collect stray segments in package 4
  for (unsigned int k=0;k<package[3].size();k++){
    DFDCSegment *segment=package[3][k];
    
    DVector3 pos,mom;
    double tanl=segment->S(3,0);
    double phi0=segment->S(1,0);
     // Circle parameters
    xc=segment->xc;
    yc=segment->yc;
    rc=segment->rc;
    // Sign of the charge
    double q=segment->S(0,0)/fabs(segment->S(0,0));
       
    double Bz_avg=0.;
    // Compute average magnitic field for the segment
    for (unsigned int m=0;m<segment->hits.size();m++){
      DFDCPseudo *hit=segment->hits[m];
      double Bx,By,Bz;
      bfield->GetField(hit->x,hit->y,hit->wire->origin(2),Bx,By,Bz);
      Bz_avg-=Bz;
    }
    Bz_avg/=double(segment->hits.size());
    double pt=0.003*Bz_avg*rc;
    double theta=M_PI_2-atan(tanl);
    double d=yc-q*rc*cos(phi0);
    double z_vertex=segment->S(4,0);
    
    mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);
    pos.SetXYZ(-d*sin(phi0),d*cos(phi0),z_vertex);

    if (z_vertex<0){
      double sperp=(Z_VERTEX-z_vertex)/tanl;
      double kappa=q/2./rc;
      double sin2ks=sin(2.*kappa*sperp);
      double cos2ks=cos(2.*kappa*sperp);
      pos(0)=pos.x()+cos(phi0)*sin2ks/2./kappa
	-sin(phi0)*(1.-cos2ks)/2./kappa;
      pos(1)=pos.y()+sin(phi0)*sin2ks/2./kappa
	+cos(phi0)*(1.-cos2ks)/2./kappa;
      pos(2)=Z_VERTEX;
      mom(0)=pt*(cos(phi0)*cos2ks-sin(phi0)*sin2ks);
      mom(1)=pt*(sin(phi0)*cos2ks+cos(phi0)*sin2ks);
    }

    
    // Create new track, starting with the current segment
    DTrackCandidate *track = new DTrackCandidate;
    track->setPosition(pos);
    track->setMomentum(mom);    
    track->setCharge(q);

    track->AddAssociatedObject(segment);

    _data.push_back(track); 
  }


  // Finally, output stray segments in package 1
  for (unsigned int k=0;k<pack1_left_over.size();k++){
    DFDCSegment *segment=pack1_left_over[k];
    
    DVector3 pos,mom;
    double tanl=segment->S(3,0);
    double phi0=segment->S(1,0);
     // Circle parameters
    xc=segment->xc;
    yc=segment->yc;
    rc=segment->rc;
    // Sign of the charge
    double q=segment->S(0,0)/fabs(segment->S(0,0));
       
    double Bz_avg=0.;
    // Compute average magnitic field for the segment
    for (unsigned int m=0;m<segment->hits.size();m++){
      DFDCPseudo *hit=segment->hits[m];
      double Bx,By,Bz;
      bfield->GetField(hit->x,hit->y,hit->wire->origin(2),Bx,By,Bz);
      Bz_avg-=Bz;
    }
    Bz_avg/=double(segment->hits.size());
    double pt=0.003*Bz_avg*rc;
    double theta=M_PI_2-atan(tanl);
    double d=yc-q*rc*cos(phi0);
    double z_vertex=segment->S(4,0); 

    // Momentum and "vertex" position
    mom.SetMagThetaPhi(pt/sin(theta),theta,phi0);
    pos.SetXYZ(-d*sin(phi0),d*cos(phi0),z_vertex);

    if (z_vertex<0){
      double sperp=(Z_VERTEX-z_vertex)/tanl;
      double kappa=q/2./rc;
      double sin2ks=sin(2.*kappa*sperp);
      double cos2ks=cos(2.*kappa*sperp);
      pos(0)=pos.x()+cos(phi0)*sin2ks/2./kappa
	-sin(phi0)*(1.-cos2ks)/2./kappa;
      pos(1)=pos.y()+sin(phi0)*sin2ks/2./kappa
	+cos(phi0)*(1.-cos2ks)/2./kappa;
      pos(2)=Z_VERTEX;
      mom(0)=pt*(cos(phi0)*cos2ks-sin(phi0)*sin2ks);
      mom(1)=pt*(sin(phi0)*cos2ks+cos(phi0)*sin2ks);
    }
    
    // Create new track, starting with the current segment
    DTrackCandidate *track = new DTrackCandidate;
    track->setPosition(pos);
    track->setMomentum(mom);    
    track->setCharge(q);

    track->AddAssociatedObject(segment);

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
  if (z<pos.z()) mom=-1.0*mom;

  // Match to the next package by swimming the track through the field
  double diff_min=1000.,diff;
  double s;
  if (stepper.SwimToPlane(pos,mom,origin,norm,&s)==false){
    for (unsigned int j=0;j<package.size();j++){
      DFDCSegment *segment2=package[j];
 
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
          double y2=segment->hits[0]->y;
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

  // Since we are assuming that the particle is coming from the target,
  // it is possible that the particle could have undergone more than a full
  // revolution and still end up producing the current segment, in which 
  // case if the vertex position is more or less correct, the dip angle is 
  // too large.  Try matching again with an adjusted tanl.
  if (match==NULL){
    diff_min=1000.;
    double tanl=segment->S(3,0);
    double sperp=(z-segment->S(4,0))/tanl;
    segment->S(3,0)=tanl*sperp/(sperp+2.*segment->rc*M_PI);
    if (GetPositionAndMomentum(segment,pos,mom)!=NOERROR) return NULL;
    if (z<pos.z()) mom=-1.0*mom;
    
    for (unsigned int j=0;j<package.size();j++){
      DFDCSegment *segment2=package[j];
 
      double x2=segment2->hits[segment2->hits.size()-1]->x;
      double y2=segment2->hits[segment2->hits.size()-1]->y;
      diff=sqrt((pos(0)-x2)*(pos(0)-x2)+(pos(1)-y2)*(pos(1)-y2));
      if (diff<diff_min&&diff<MATCH_RADIUS(mom.Mag())){
	diff_min=diff;
	match=segment2;
	match_id=j;
      }
    }
    // Restore old value
    segment->S(3,0)=tanl;
  }

  // If matching in the forward direction did not work, try swimming and
  // matching backwards... with modified tanl values
  if (match==NULL){
    diff_min=1000.;
    for (unsigned int i=0;i<package.size();i++){
      DFDCSegment *segment2=package[i];
      double tanl=segment2->S(3,0);
      double sperp=(z-segment2->S(4,0))/tanl;
      segment2->S(3,0)=tanl*sperp/(sperp+2.*segment2->rc*M_PI);
      
      if (GetPositionAndMomentum(segment2,pos,mom)==NOERROR){
        mom=-1.0*mom;
        origin(2)=segment->hits[0]->wire->origin(2);
        if (stepper.SwimToPlane(pos,mom,origin,norm,&s)==false){
          double x2=segment->hits[0]->x;
          double y2=segment->hits[0]->y;
          diff=sqrt((pos(0)-x2)*(pos(0)-x2)+(pos(1)-y2)*(pos(1)-y2));
          if (diff<diff_min&&diff<MATCH_RADIUS(mom.Mag())){
	    diff_min=diff;
	    match=segment2;
            match_id=i;
          }
        }	
      } 
      // Restore old value
      segment2->S(3,0)=tanl;
    }
  }


  if(DEBUG_HISTS){
    match_dist_fdc->Fill(mom.Mag(),diff_min);
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
  //  if (sqrt(x*x+y*y)>FDC_OUTER_RADIUS) return VALUE_OUT_OF_RANGE;

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
  double Bx,By,Bz;
  bfield->GetField(x,y,z,Bx,By,Bz);
  Bz=fabs(Bz);
  
  // Momentum
  double px=(cosp*cos2ks-sinp*sin2ks)*0.003*Bz/2./kappa;
  double py=(sinp*cos2ks+cosp*sin2ks)*0.003*Bz/2./kappa;
  double pz=0.003*Bz*tanl/2./kappa;

  //if (sqrt(px*px+py*py)>PT_MAX) return VALUE_OUT_OF_RANGE;

  pos.SetXYZ(x,y,z);
  mom.SetXYZ(px,py,pz);

  return NOERROR;
}



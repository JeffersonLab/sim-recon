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
#include "FDC/DFDCPseudo_factory.h"
#include "FDC/DFDCSegment_factory.h"
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
///
/// DTrackCandidate_factory_FDCCathodes::brun():
///
jerror_t DTrackCandidate_factory_FDCCathodes::brun(JEventLoop* eventLoop, 
						   int runnumber) {
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  
  USE_FDC=true;
  if (!dgeom->GetFDCZ(z_wires)){
    _DBG_<< "FDC geometry not available!" <<endl;
    USE_FDC=false;
  }

  // Get the position of the CDC downstream endplate from DGeometry
  double endplate_dz,endplate_rmin,endplate_rmax;
  dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z+=endplate_dz;
  dgeom->GetTargetZ(TARGET_Z);

  DEBUG_HISTS=false;
  gPARMS->SetDefaultParameter("TRKFIND:DEBUG_HISTS", DEBUG_HISTS);

  APPLY_MOMENTUM_CORRECTION=false;
  gPARMS->SetDefaultParameter("TRKFIND:APPLY_MOMENTUM_CORRECTION",APPLY_MOMENTUM_CORRECTION);
  p_factor1=1.61*M_PI/180.;
  p_factor2=-0.0766;

  string description = "If hit wih largest R is less than this, then a ";
  description += "fake point will be added when fitting the parameters ";
  description += "for the track candidate in the 'FDCCathodes' factory. ";
  description += "The point will be on the beamline (x,y) = (0,0) and ";
  description += "at a z location determined from the geometry center of ";
  description += "target (via DGeometry::GetTargetZ()";
  MAX_R_VERTEX_LIMIT = 10.0;
  gPARMS->SetDefaultParameter("TRKFIND:MAX_R_VERTEX_LIMIT", MAX_R_VERTEX_LIMIT, description);

  if(DEBUG_HISTS) {
    dapp->Lock();
    match_dist_fdc=(TH2F*)gROOT->FindObject("match_dist_fdc");
    if (!match_dist_fdc) 
      match_dist_fdc=new TH2F("match_dist_fdc",
		  "Matching distance for connecting FDC segments",
			      50,0.,7,100,0,25.);
    Hcircle_fit_prob=(TH1F*)gROOT->FindObject("Hcircle_fit_prob");
    if (!Hcircle_fit_prob){
      Hcircle_fit_prob=new TH1F("Hcircle_fit_prob","CL for circle fit",
				1000,0,1);
    }
    dapp->Unlock();
  }
    
  // Initialize the stepper
  stepper=new DMagneticFieldStepper(bfield);
  stepper->SetStepSize(1.0);

  return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DTrackCandidate_factory_FDCCathodes::erun(void)
{
  if (stepper) delete stepper;
  return NOERROR;
}
//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_FDCCathodes::fini(void)
{
  
  if (stepper) delete stepper;
  return NOERROR;
}


// local rotine for sorting track candidates
inline bool DTrackCandidate_cmp(const DTrackCandidate *a,const DTrackCandidate *b){
  return (a->momentum().Mag()>b->momentum().Mag());
}

// Local routine for sorting segments by charge and curvature
inline bool DTrackCandidate_segment_cmp(const DFDCSegment *a, const DFDCSegment *b){
  //  double k1=a->S(0,0),k2=b->S(0,0);
  //double q1=k1/fabs(k1),q2=k2/fabs(k2);
  //if (q1!=q2) return q1<q2;
  //return fabs(k1)<fabs(k2); 
  if (a->q!=b->q) return a->q<b->q;
  return a->rc>b->rc;
}

//------------------
// evnt:  main segment linking routine
//------------------
jerror_t DTrackCandidate_factory_FDCCathodes::evnt(JEventLoop *loop, int eventnumber)
{
  if (!USE_FDC) return NOERROR;

  vector<const DFDCSegment*>segments;
  eventLoop->Get(segments);
  // abort if there are no segments
  if (segments.size()==0.) return NOERROR;

  // Clear vector of fit results.
  fit_results.clear();
  
  // Sort segments by charge and curvature
  std::sort(segments.begin(), segments.end(), DTrackCandidate_segment_cmp);

  // Group segments by package
  vector<DFDCSegment*>packages[4];  
  for (unsigned int i=0;i<segments.size();i++){
    const DFDCSegment *segment=segments[i];
    int packno=(segment->hits[0]->wire->layer-1)/6;
    packages[packno].push_back((DFDCSegment*)segment);
  }

  // Loop over all the packages to match to segments in packages downstream
  // of the current package
  for (unsigned int i=0;i<4;i++){
    if (packages[i].size()>0) LinkSegments(i,packages);
  }

  // Now collect stray segments
  for (unsigned int i=0;i<4;i++){
    for (unsigned int k=0;k<packages[i].size();k++){
      DFDCSegment *segment=packages[i][k];

      // Make sure that this segment does not match one of the previously 
      // linked tracks
      if (LinkSegment(segment)==false){   
	DVector3 pos,mom;
	// Circle parameters
	xc=segment->xc;
	yc=segment->yc;
	rc=segment->rc;
	// Sign of the charge
	q=segment->q;
	// Angle parameter
	tanl=segment->tanl;
	
	// Try to fix relatively high momentum tracks in the very forward 
	// direction that look like low momentum tracks due to small pt.
	// Assume that the particle came from the center of the target.
	const DFDCPseudo *first_hit=segment->hits[0];
	double max_r=first_hit->xy.Mod();
	unsigned int index=0;
	for (unsigned int m=1;m<segment->hits.size();m++){
	  const DFDCPseudo *hit=segment->hits[m];
	  double r=hit->xy.Mod();
	  if (r>max_r){
	    max_r=r;
	    index=m;
	  }
	}
	if (rc<0.5*max_r && max_r<6.0){
	  const DFDCPseudo *hit=segment->hits[index];
	  tanl=(hit->wire->origin.z()-TARGET_Z)/max_r;
	  rc=0.5*max_r;
	  xc=0.5*hit->xy.X();
	  yc=0.5*hit->xy.Y();
	}	
	// |Bz| at the first hit in the segment
	double Bz=-bfield->GetBz(first_hit->xy.X(),first_hit->xy.Y(),
				 first_hit->wire->origin.z());
	
	pos.SetXYZ(first_hit->xy.X(),first_hit->xy.Y(),
		   first_hit->wire->origin.z());    
	GetPositionAndMomentum(Bz,pos,mom);
	
	// Empirical correction to the momentum 
	if (APPLY_MOMENTUM_CORRECTION){
	  double p_mag=mom.Mag();
	  mom.SetMag(p_mag*(1.+p_factor1/mom.Theta()+p_factor2));
	}
	
	// Create new track, starting with the current segment
	DTrackCandidate *track = new DTrackCandidate;
	track->setPosition(pos);
	track->setMomentum(mom);    
	track->setCharge(q);

	track->AddAssociatedObject(segment);

	// Add another fit class to the list for use by a later track...
	DHelicalFit fit;
	fit.r0=segment->rc;
	for (unsigned int is=0;is<segment->hits.size();is++){
	  fit.AddHit(segment->hits[is]);
	}
	fit_results.push_back(fit);

	_data.push_back(track);
      }
    }
  }
   
  // Sort the tracks by momentum
  sort(_data.begin(),_data.end(),DTrackCandidate_cmp);

  return NOERROR;
}


// Routine to do a crude match between fdc points and a helical approximation to
// the trajectory
double DTrackCandidate_factory_FDCCathodes::DocaToHelix(const DFDCPseudo *hit){
  double sperp=(hit->wire->origin.z()-zs)*cotl;
  double twoks=twokappa*sperp;
  double sin2ks=sin(twoks);
  double cos2ks=cos(twoks);
  double one_minus_cos2ks=1.-cos2ks;
  double x=xs+(cosphi*sin2ks-sinphi*one_minus_cos2ks)*one_over_twokappa;
  double y=ys+(sinphi*sin2ks+cosphi*one_minus_cos2ks)*one_over_twokappa;
  double dx=x-hit->xy.X();
  double dy=y-hit->xy.Y();

  return sqrt(dx*dx+dy*dy);
}

// Propagate track from one package to the next and look for a match to a 
// segment in the new package
DFDCSegment *DTrackCandidate_factory_FDCCathodes::GetTrackMatch(DFDCSegment *segment,
								vector<DFDCSegment*>package,
								unsigned int &match_id){
  DFDCSegment *match=NULL;
  
  // Get the position and momentum at the exit of the package for the 
  // current segment
  GetPositionAndMomentum(segment);

  // Match to the next package
  double diff_min=1e6,diff;
  for (unsigned int j=0;j<package.size();j++){
    DFDCSegment *segment2=package[j];
    for (unsigned int k=0;k<segment2->hits.size();k++){
      diff=DocaToHelix(segment2->hits[k]);
      if (diff<diff_min&&diff<Match(p)){
	diff_min=diff;
	match=segment2;
	match_id=j;
      }
    }
  }
  
  // If matching in the forward direction did not work, try 
  // matching backwards...
  if (match==NULL){
    diff_min=1e6;
    for (unsigned int i=0;i<package.size();i++){
      DFDCSegment *segment2=package[i];
      GetPositionAndMomentum(segment2);
      for (unsigned int k=0;k<segment->hits.size();k++){
	diff=DocaToHelix(segment->hits[k]);
	if (diff<diff_min&&diff<Match(p)){
	  diff_min=diff;
	  match=segment2;
	  match_id=i;
	}       
      }
    }
  }

  if(DEBUG_HISTS){
    match_dist_fdc->Fill(p,diff_min);
  }
  return match;
}

// Obtain position and momentum at the exit of a given package using the 
// helical track model.
//
jerror_t DTrackCandidate_factory_FDCCathodes::GetPositionAndMomentum(const DFDCSegment *segment){
  // Position of track segment at last hit plane of package
  xs=segment->xc+segment->rc*cos(segment->Phi1);
  ys=segment->yc+segment->rc*sin(segment->Phi1);
  zs=segment->hits[0]->wire->origin.z();

  // Track parameters
  //double kappa=segment->q/(2.*segment->rc);
  double my_phi0=segment->phi0;
  double my_tanl=segment->tanl;
  double z0=segment->z_vertex;

  // Useful intermediate variables
  double cosp=cos(my_phi0);
  double sinp=sin(my_phi0);
  double twoks=segment->q*(zs-z0)/(my_tanl*segment->rc);
  double sin2ks=sin(twoks);
  double cos2ks=cos(twoks); 

  // Get Bfield
  double Bz=fabs(bfield->GetBz(xs,ys,zs));

  // Momentum
  p=0.003*Bz*segment->rc/cos(atan(my_tanl));
  cosphi=cosp*cos2ks-sinp*sin2ks;
  sinphi=sinp*cos2ks+cosp*sin2ks;
  twokappa=segment->q/segment->rc;
  one_over_twokappa=1./twokappa;
  cotl=1./my_tanl;

  return NOERROR;
}



// Obtain position and momentum at the exit of a given package using the 
// helical track model.
//
jerror_t DTrackCandidate_factory_FDCCathodes::GetPositionAndMomentum(DFDCSegment *segment,
					      DVector3 &pos, DVector3 &mom){
  // Position of track segment at last hit plane of package
  double x=segment->xc+segment->rc*cos(segment->Phi1);
  double y=segment->yc+segment->rc*sin(segment->Phi1);
  double z=segment->hits[0]->wire->origin.z();
  pos.SetXYZ(x,y,z);

  zs=z;
  ys=y;
  xs=x;

  // Make sure that the position makes sense!
  //  if (sqrt(x*x+y*y)>FDC_OUTER_RADIUS) return VALUE_OUT_OF_RANGE;

  // Track parameters
  //double kappa=segment->q/(2.*segment->rc);
  double my_phi0=segment->phi0;
  double my_tanl=segment->tanl;
  double z0=segment->z_vertex;

  // Useful intermediate variables
  double cosp=cos(my_phi0);
  double sinp=sin(my_phi0);
  // double twoks=2.*kappa*(z-z0)/my_tanl;
  double twoks=segment->q*(z-z0)/(my_tanl*segment->rc);
  double sin2ks=sin(twoks);
  double cos2ks=cos(twoks); 

  // Get Bfield
  double Bz=fabs(bfield->GetBz(x,y,z));

  // Momentum
  double pt=0.003*Bz*segment->rc;
  cosphi=cosp*cos2ks-sinp*sin2ks;
  sinphi=sinp*cos2ks+cosp*sin2ks;
  twokappa=segment->q/segment->rc;
  one_over_twokappa=1./twokappa;
  cotl=1./my_tanl;
  mom.SetXYZ(pt*(cosp*cos2ks-sinp*sin2ks),pt*(sinp*cos2ks+cosp*sin2ks),
	     pt*my_tanl);

  return NOERROR;
}

// Routine to return momentum and position given the helical parameters and the
// z-component of the magnetic field
jerror_t DTrackCandidate_factory_FDCCathodes::GetPositionAndMomentum(const double Bz_avg,DVector3 &pos,DVector3 &mom){
  double pt=0.003*Bz_avg*rc;
  double phi1=atan2(pos.y()-yc,pos.x()-xc);
  double dphi_s=q*(pos.z()-endplate_z)/(rc*tanl);
  double dphi1=phi1-dphi_s;// was -
  double x=xc+rc*cos(dphi1);
  double y=yc+rc*sin(dphi1);
  pos.SetXYZ(x,y,endplate_z);

  dphi1*=-1.;
  if (q<0) dphi1+=M_PI;

  double px=pt*sin(dphi1);
  double py=pt*cos(dphi1);
  double pz=pt*tanl;
  mom.SetXYZ(px,py,pz);

  return NOERROR;
}

double DTrackCandidate_factory_FDCCathodes::GetCharge(const DVector3 &pos,
						 const DFDCSegment *segment
						      ){
  // Magnitude of change in phi to go from the reference plane to another point
  // in the segment
  double dphi=(segment->hits[0]->wire->origin.z()-pos.z())/(rc*tanl);

  // Phi at reference plane
  double Phi1=atan2(pos.y()-yc,pos.x()-xc);

  // Positive and negative changes in phi
  double phiplus=Phi1+dphi;
  double phiminus=Phi1-dphi;
  DVector2 plus(xc+rc*cos(phiplus),yc+rc*sin(phiplus));
  DVector2 minus(xc+rc*cos(phiminus),yc+rc*sin(phiminus));

  // Compute differences 
  double d2plus=(plus-segment->hits[0]->xy).Mod2();
  double d2minus=(minus-segment->hits[0]->xy).Mod2();

  // Look for smallest difference to determine q
  if (d2minus<d2plus){
    return -1.;
  }

  return 1.;
}

// Routine to loop over segments in one of the packages, linking them with 
// segments in the packages downstream of this package
void DTrackCandidate_factory_FDCCathodes::LinkSegments(unsigned int pack1,
					vector<DFDCSegment *>packages[4]){ 
  unsigned int match_id=0;
  unsigned int pack2=pack1+1;
  unsigned int pack3=pack2+1;
  unsigned int pack4=pack3+1;
  bool pack2_exists=(pack2<4);
  bool pack3_exists=(pack3<4);
  bool pack4_exists=(pack4<4);

  // Keep track of the segments in the first package that have been matched to
  // other segments
  unsigned int num_pack1=packages[pack1].size();
  vector<int>matched(num_pack1);

  // Loop over the segments in the first package
  for (unsigned int i=0;i<num_pack1;i++){
    DFDCSegment *segment=packages[pack1][i];
    DFDCSegment *match2=NULL;
    DFDCSegment *match3=NULL;
    DFDCSegment *match4=NULL;
 
    // Tracking parameters from first segment
    xc=segment->xc;
    yc=segment->yc;
    rc=segment->rc;
    tanl=segment->tanl;

    // Sign of the charge
    q=segment->q;
    
    // Start filling vector of segments belonging to current track    
    vector<DFDCSegment*>mysegments; 
    mysegments.push_back(segment);
    
    // Try matching to package 2
    if (pack2_exists && packages[pack2].size()>0 
	&& (match2=GetTrackMatch(segment,packages[pack2],match_id))!=NULL){

      // Insert the segment from package 2 into the track 
      mysegments.push_back(match2);
         
      // remove the segment from the list 
      packages[pack2].erase(packages[pack2].begin()+match_id);
      
      // Try matching to package 3 starting from package 2
      if (pack3_exists && packages[pack3].size()>0
	  && (match3=GetTrackMatch(match2,packages[pack3],match_id))!=NULL){

	// Insert the segment from package 3 into the track 
	mysegments.push_back(match3);
	
	// remove the segment from the list 
	packages[pack3].erase(packages[pack3].begin()+match_id);
	
	// Try matching to package 4 starting from package 3
	if (pack4_exists && packages[pack4].size()>0 &&  
	    (match4=GetTrackMatch(match3,packages[pack4],match_id))!=NULL){
	  
	  // Insert the segment from package 4 into the track 
	  mysegments.push_back(match4);
	  
	  // remove the segment from the list 
	  packages[pack4].erase(packages[pack4].begin()+match_id);
	}
      }
      // No match in package 3, try for 4
      else if(pack4_exists && packages[pack4].size()>0 && 
	      (match4=GetTrackMatch(match2,packages[pack4],match_id))!=NULL){

	// Insert the segment from package 4 into the track 
	mysegments.push_back(match4);
	
	// remove the segment from the list 
	packages[pack4].erase(packages[pack4].begin()+match_id);
      }
    }
    // No match in package 2, try for 3
    else if (pack3_exists && packages[pack3].size()>0 && 
	     (match3=GetTrackMatch(segment,packages[pack3],match_id))!=NULL){

      // Insert the segment from package 3 into the track
      mysegments.push_back(match3);
	
      // remove the segment from the list 
      packages[pack3].erase(packages[pack3].begin()+match_id);
      
      // Try matching to package 4
      if (pack4_exists && packages[pack4].size()>0 && 
	  (match4=GetTrackMatch(match3,packages[pack4],match_id))!=NULL){

	// Insert the segment from package 4 into the track 
	mysegments.push_back(match4);
	
	// remove the segment from the list 
	packages[pack4].erase(packages[pack4].begin()+match_id);
      }
    }    
    // No match to package 2 or 3, try 4
    else if (pack4_exists && packages[pack4].size()>0 && 
	     (match4=GetTrackMatch(segment,packages[pack4], match_id))!=NULL){

      // Insert the segment from package 4 into the track 
      mysegments.push_back(match4);
	
      // remove the segment from the list 
      packages[pack4].erase(packages[pack4].begin()+match_id);
    }
       
    // If we found a match, redo the helical fit with all the hits
    if (mysegments.size()>1){
      // Mark the segment in the first package as matched
      matched[i]=1;
    
      // Variables for determining average Bz
      double Bz_avg=0.;
      unsigned int num_hits=segment->hits.size();
      
      DHelicalFit fit;
      double max_r=0.,my_x=0.,my_y=0.,my_z=0.;
      double min_r=1000.;
      // First add the hits to the fit class
      for (unsigned int m=0;m<mysegments.size();m++){
	for (unsigned int n=0;n<mysegments[m]->hits.size();n++){
	  const DFDCPseudo *hit=mysegments[m]->hits[n];
	  fit.AddHit(hit);
	  Bz_avg-=bfield->GetBz(hit->xy.X(),hit->xy.Y(),
				hit->wire->origin.z());
	  double r=hit->xy.Mod();
	  if (r<min_r) min_r=r;
	  if (r>max_r){
	    max_r=r;
	    my_x=hit->xy.X();
	    my_y=hit->xy.Y();
	    my_z=hit->wire->origin.z();
	  }
	}
	num_hits+=mysegments[m]->hits.size();
      }
      // Fake point at origin
      if (max_r<MAX_R_VERTEX_LIMIT) fit.AddHitXYZ(0.,0.,TARGET_Z,BEAM_VAR,BEAM_VAR,0.);
      // Do the fit
      if (fit.FitCircleAndLineRiemann(rc)==NOERROR){      
	if (DEBUG_HISTS){
	  Hcircle_fit_prob->Fill(TMath::Prob(fit.chisq,fit.ndof));
	}
	if (fit.r0>0.5*(max_r-min_r)){
	  // New track parameters
	  tanl=fit.tanl;
	  xc=fit.x0;
	  yc=fit.y0;
	  rc=fit.r0;
	}
	else{

	  // At least use two segments to make a crude approximation for tanl..
	  tanl=(mysegments[1]->hits[0]->wire->origin.z()
		-mysegments[0]->hits[0]->wire->origin.z())
	    /(mysegments[1]->hits[0]->xy.Mod()
	      -mysegments[0]->hits[0]->xy.Mod());
	}
      }
   
      // Try to fix relatively high momentum tracks in the very forward 
      // direction that look like low momentum tracks due to small pt.
      // Assume that the particle came from the center of the target.
      if (max_r<6.0 && rc<0.5*max_r){
	tanl=(my_z-TARGET_Z)/max_r;
	rc=0.5*max_r;
	xc=0.5*my_x;
	yc=0.5*my_y;	 
	fit.r0=rc;
      }		

      // Position of a hit in the segment
      const DFDCPseudo *fdchit=segment->hits[0];
      DVector3 pos(fdchit->xy.X(),fdchit->xy.Y(),fdchit->wire->origin.z());
      
      // Create new track, starting with the current segment
      DTrackCandidate *track = new DTrackCandidate;

      if (match2){
	q=GetCharge(pos,match2);
      }
      else if (match3){
	q=GetCharge(pos,match3);
      }
      else if (match4){
	q=GetCharge(pos,match4);
      }

      DVector3 mom;
      Bz_avg/=double(num_hits);
      GetPositionAndMomentum(Bz_avg,pos,mom);
	
      // Empirical correction to the momentum
      if (APPLY_MOMENTUM_CORRECTION){
	double p_mag=mom.Mag();
	mom.SetMag(p_mag*(1.+p_factor1/mom.Theta()+p_factor2));
      }
      
      track->setCharge(q);
      track->setPosition(pos);
      track->setMomentum(mom);
      
      for (unsigned int m=0;m<mysegments.size();m++)
	track->AddAssociatedObject(mysegments[m]);
      
      _data.push_back(track); 
      fit_results.push_back(fit);
    }
  }
  
  // Cull the segments that have been matched already from the list of segments
  // in the first package
  vector<DFDCSegment*>left_over;
  for (unsigned int i=0;i<num_pack1;i++){
    if (matched[i]!=1) left_over.push_back(packages[pack1][i]);
  }
  packages[pack1].assign(left_over.begin(),left_over.end());

}

// Routine that tries to link a stray segment with an already existing track
// candidate
bool DTrackCandidate_factory_FDCCathodes::LinkSegment(const DFDCSegment *segment){
  int new_packno=(segment->hits[0]->wire->layer-1)/6;
  for (unsigned int i=0;i<_data.size();i++){
    vector<const DFDCSegment *>segments;
    _data[i]->GetT(segments);

    // if a track already has a matched segment in this package, skip...
    for (unsigned int j=0;j<segments.size();j++){
      int packno=(segments[j]->hits[0]->wire->layer-1)/6;
      if (packno==new_packno) return false;
    }

    DVector3 pos=_data[i]->position();
    DVector3 mom=_data[i]->momentum();
    if (GetTrackMatch(_data[i]->charge(),pos,mom,segment)){
      const DFDCPseudo *first_hit=segment->hits[0];
      
      // Add the segment to the track as an associated object
      _data[i]->AddAssociatedObject(segment);
      
      // Add hits to fit_results for this track and compute the average Bz
      double Bz_avg=-bfield->GetBz(first_hit->xy.X(),first_hit->xy.Y(),
				   first_hit->wire->origin.z());
      fit_results[i].AddHit(first_hit);
      double min_r=1000.,max_r=0.;
      for (unsigned int k=1;k<segment->hits.size();k++){
	const DFDCPseudo *hit=segment->hits[k];
	fit_results[i].AddHit(hit);
	Bz_avg-=bfield->GetBz(hit->xy.X(),hit->xy.Y(),hit->wire->origin.z());
	double r=hit->xy.Mod();
	if (r>max_r) max_r=r;
	if (r<min_r) min_r=r;
      }
      Bz_avg/=double(segment->hits.size());
      
      // Redo the helical fit with the additional hits
      double old_rc=fit_results[i].r0;
      if (fit_results[i].FitCircleAndLineRiemann(old_rc)==NOERROR){      	
	// New circle radius
	rc=fit_results[i].r0;
	if (rc>0.5*(max_r-min_r)){
	  // Other new track parameters
	  tanl=fit_results[i].tanl;
	  xc=fit_results[i].x0;
	  yc=fit_results[i].y0;
	  q=_data[i]->charge();
	  
	  pos.SetXYZ(first_hit->xy.X(),first_hit->xy.Y(),
		     first_hit->wire->origin.z());  
	  GetPositionAndMomentum(Bz_avg,pos,mom);
	  
	  _data[i]->setPosition(pos);
	  _data[i]->setMomentum(mom);
	}
      }

      return true;
    }
  }
  return false;
}

// Routine for matching to a segment using the stepper
bool DTrackCandidate_factory_FDCCathodes::GetTrackMatch(double q,
							DVector3 &pos,
							DVector3 &mom,
							const DFDCSegment *segment){
  const DVector3 norm(0,0,1.);
  stepper->SetCharge(q);
  
  const DFDCPseudo *hit=segment->hits[0];
  if (stepper->SwimToPlane(pos,mom,hit->wire->origin,norm,NULL)==false){
    double dx=hit->xy.X()-pos.x();
    double dy=hit->xy.Y()-pos.y();
    double d=sqrt(dx*dx+dy*dy);
    if (d<Match(mom.Mag())) return true;
  }
  return false;
}


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
#include <JANA/JCalibration.h>
#include "FDC/DFDCPseudo_factory.h"
#include "FDC/DFDCSegment_factory.h"
#include "DHelicalFit.h"
#include "DHoughFind.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>

#define MAX_SEGMENTS 20
#define HALF_PACKAGE 6.0
#define FDC_OUTER_RADIUS 50.0 
#define BEAM_VAR 0.0208 // cm^2
#define HIT_CHI2_CUT 10.0
///
/// DTrackCandidate_factory_FDCCathodes::brun():
///
jerror_t DTrackCandidate_factory_FDCCathodes::brun(JEventLoop* eventLoop, 
						   int32_t runnumber) {
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield(runnumber);
  FactorForSenseOfRotation=(bfield->GetBz(0.,0.,65.)>0.)?-1.:1.;

  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  
  USE_FDC=true;
  if (!dgeom->GetFDCZ(z_wires)){
    _DBG_<< "FDC geometry not available!" <<endl;
    USE_FDC=false;
  }
  /*
  zpack[0]=z_wires[0]-1.;
  zpack[1]=z_wires[6]-1.;
  zpack[2]=z_wires[12]-1.;
  zpack[3]=z_wires[18]-1.;
  */

  // Get the position of the CDC downstream endplate from DGeometry
  double endplate_dz,endplate_rmin,endplate_rmax;
  dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z+=endplate_dz;


  JCalibration *jcalib = dapp->GetJCalibration(runnumber);
  map<string, double> targetparms;
  if (jcalib->Get("TARGET/target_parms",targetparms)==false){
    TARGET_Z = targetparms["TARGET_Z_POSITION"];
  }
  else{
    dgeom->GetTargetZ(TARGET_Z);
  }

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
  MAX_R_VERTEX_LIMIT = 50.0;
  gPARMS->SetDefaultParameter("TRKFIND:MAX_R_VERTEX_LIMIT", MAX_R_VERTEX_LIMIT, description);

  MATCHING_PHI_CUT=10.0;
  gPARMS->SetDefaultParameter("TRKFIND:MATCHING_PHI_CUT", MATCHING_PHI_CUT);

  if(DEBUG_HISTS) {
    dapp->Lock();
    match_dist_fdc=(TH2F*)gROOT->FindObject("match_dist_fdc");
    if (!match_dist_fdc){ 
      match_dist_fdc=new TH2F("match_dist_fdc",
			      "Matching distance for connecting FDC segments",
			      50,0.,7,500,0,100.);
    }
    match_center_dist2=(TH2F*)gROOT->FindObject("match_center_dist2");
    if (!match_center_dist2){
      match_center_dist2=new TH2F("match_center_dist2","matching distance squared between two circle centers vs p",50,0,1.5,100,0,100);
      match_center_dist2->SetXTitle("p [GeV/c]");
      match_center_dist2->SetYTitle("(#Deltad)^{2} [cm^{2}]");
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
  if (stepper) {
    delete stepper;
    stepper = nullptr;
  }

  return NOERROR;
}
//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_FDCCathodes::fini(void)
{
  if (stepper) {
    delete stepper;
    stepper = nullptr;
  }  

  return NOERROR;
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


inline bool DTrackCandidate_segment_cmp_by_z(const DFDCSegment *a, 
					     const DFDCSegment *b){  
  return (a->hits[0]->wire->origin.z()<b->hits[0]->wire->origin.z());
}


//------------------
// evnt:  main segment linking routine
//------------------
jerror_t DTrackCandidate_factory_FDCCathodes::evnt(JEventLoop *loop, uint64_t eventnumber)
{
  if (!USE_FDC) return NOERROR;

  vector<const DFDCSegment*>segments;
  eventLoop->Get(segments);

  // abort if there are no segments
  if (segments.size()==0.) return NOERROR;

  std::sort(segments.begin(), segments.end(), DTrackCandidate_segment_cmp);

  // Group segments by package
  vector<DFDCSegment*>packages[4];
  for (unsigned int i=0;i<segments.size();i++){
    const DFDCSegment *segment=segments[i];
    packages[segment->package].push_back((DFDCSegment*)segment);
  }
  // Keep track of the segments in each package that have been paired with
  // other segments
  vector<vector<int> >is_paired;
  for (unsigned int i=0;i<4;i++){
    vector<int>temp(packages[i].size());
    is_paired.push_back(temp);
  } 
  // Loop over all the packages to match to segments in packages downstream
  // of the current package 
  vector<pair<const DFDCSegment*,const DFDCSegment*> >paired_segments;
  for (unsigned int i=0;i<3;i++){
    if (packages[i].size()>0) LinkSegments(i,packages,paired_segments,is_paired);
  }
  
  // Link pairs of segments into groups of three linked segments
  vector<vector<const DFDCSegment *> >triplets;
  vector<int>is_tripled(paired_segments.size());
  for (unsigned int i=0;i<paired_segments.size();i++){
    for (unsigned int j=0;j<paired_segments.size();j++){
      if (i==j) continue;
      if (paired_segments[i].second==paired_segments[j].first){
	is_tripled[i]=1;
	is_tripled[j]=1;

	vector<const DFDCSegment *>triplet;
	triplet.push_back(paired_segments[i].first);
	triplet.push_back(paired_segments[i].second);
	triplet.push_back(paired_segments[j].second);
	triplets.push_back(triplet);
      }  
    }
  } 
  // Link triplets with pairs to form groups of four linked segments
  vector<int>is_quadrupled(triplets.size());
  vector<vector<const DFDCSegment *> >quadruplets;
  for (unsigned int i=0;i<triplets.size();i++){
    for (unsigned int j=0;j<paired_segments.size();j++){
      if (triplets[i][2]==paired_segments[j].first){
	is_quadrupled[i]=1;
	is_tripled[j]=1;

	vector<const DFDCSegment*>quadruplet=triplets[i];
	quadruplet.push_back(paired_segments[j].second);
	quadruplets.push_back(quadruplet);
      } 
    }
  }
  // Mark all triplets that are subsets of quadruplets that have not been marked
  // previously.
  for (unsigned int i=0;i<quadruplets.size();i++){
    for (unsigned int j=0;j<triplets.size();j++){
      if (is_quadrupled[j]==0){
	unsigned int num=0;
	for (unsigned int k=0;k<4;k++){
	  for (unsigned int n=0;n<3;n++){
	    if (quadruplets[i][k]==triplets[j][n]) num++;
	  }
	}
	if (num==3) is_quadrupled[j]=1;
      }
    }
  }

  // Start gathering groups into a list of linked segments to elevate to track
  // candidates
  vector<vector<const DFDCSegment *> >mytracks;
  if (quadruplets.size()==1){
    mytracks.push_back(quadruplets[0]);
  }
  else if (quadruplets.size()>1){    
    // Because segments could have been added to the triplets on either end,
    // we need to check for clones
    vector<int>is_clone(quadruplets.size());
    for (unsigned int i=0;i<is_clone.size()-1;i++){
      for (unsigned int j=i+1;j<is_clone.size();j++){
	unsigned int num=0;
	for (unsigned int k=0;k<4;k++){
	  if (quadruplets[i][k]==quadruplets[j][k]) num++;
	}
	if (num==4){
	  is_clone[j]=1;
	  printf("Got clone!\n");
	}
      }
    }
    for (unsigned int i=0;i<quadruplets.size();i++){
      if (is_clone[i]==0){
	mytracks.push_back(quadruplets[i]);
      }
    }
  }
  // If we could not link some of the pairs together, create two-segment 
  // "tracks"
  for (unsigned int i=0;i<is_tripled.size();i++){
    if (is_tripled[i]==0){
      vector<const DFDCSegment *>mytrack;
      mytrack.push_back(paired_segments[i].first);
      mytrack.push_back(paired_segments[i].second);
      mytracks.push_back(mytrack);
    }
  }
  // if we could not link some of the triplets to other segments, create 
  // three-segment "tracks"
  for (unsigned int i=0;i<triplets.size();i++){
    if (is_quadrupled[i]==0){
      mytracks.push_back(triplets[i]);
    }
  }

  // For each set of matched segments, redo the helical fit with all the hits 
  // and create a new track candidate
  for (unsigned int i=0;i<mytracks.size();i++){  
    // Create the fit object and add the hits
    DHelicalFit fit;
    double max_r=0.;
    rc=0.; // create a guess for rc
    for (unsigned int m=0;m<mytracks[i].size();m++){
      rc+=mytracks[i][m]->rc;
      for (unsigned int n=0;n<mytracks[i][m]->hits.size();n++){
	const DFDCPseudo *hit=mytracks[i][m]->hits[n];
	fit.AddHit(hit);
	
	double r=hit->xy.Mod();
	if (r>max_r){
	  max_r=r;
	}
      }
    }
    rc/=double(mytracks[i].size());
    // Fake point at origin
    bool use_fake_point=false;
    if (max_r<MAX_R_VERTEX_LIMIT){
      fit.AddHitXYZ(0.,0.,TARGET_Z,BEAM_VAR,BEAM_VAR,0.);
      use_fake_point=true;
    }
    // Do the fit
    if (fit.FitTrackRiemann(rc)==NOERROR){    
      // New track parameters
      tanl=fit.tanl;
      xc=fit.x0;
      yc=fit.y0;
      rc=fit.r0;
      q=FactorForSenseOfRotation*fit.h;
    }
    
    // Try to fix relatively high momentum tracks in the very forward 
    // direction that look like low momentum tracks due to small pt.
    // Assume that the particle came from the center of the target.
    
    if (rc<0.5*max_r && max_r<10.0){
      if (use_fake_point) fit.PruneHit(0); 
      fit.FitTrack_FixedZvertex(TARGET_Z);
      tanl=fit.tanl;
      rc=fit.r0;
      xc=fit.x0;
      yc=fit.y0;
      fit.FindSenseOfRotation();
      q=FactorForSenseOfRotation*fit.h;
    }      
    
    // Create new track, starting with the most upstream segment
    DTrackCandidate *track = new DTrackCandidate;
    //circle fit parameters
    track->rc=rc;
    track->xc=xc;
    track->yc=yc;

    // Get the momentum and position just upstream of first hit
    DVector3 mom,pos;
    GetPositionAndMomentum(mytracks[i],pos,mom);
    
    // Empirical correction to the momentum
    if (APPLY_MOMENTUM_CORRECTION){
      double p_mag=mom.Mag();
      mom.SetMag(p_mag*(1.+p_factor1/mom.Theta()+p_factor2));
    }
    
    track->chisq=fit.chisq;
    track->Ndof=fit.ndof;
    track->setPID((q > 0.0) ? PiPlus : PiMinus);
    track->setPosition(pos);
    track->setMomentum(mom);
    
    for (unsigned int m=0;m<mytracks[i].size();m++){
      track->AddAssociatedObject(mytracks[i][m]);
    }
    
    _data.push_back(track); 
 
  }

  // Now try to attach stray segments to existing tracks
  for (unsigned int i=0;i<4;i++){
    for (unsigned int k=0;k<packages[i].size();k++){
      DFDCSegment *segment=packages[i][k];
      if (is_paired[i][k]==0 && LinkStraySegment(segment)) is_paired[i][k]=1;
    }
  }

  // Create track stubs for unused segments
  for (unsigned int j=0;j<4;j++){
    for (unsigned int i=0;i<packages[j].size();i++){
      if (is_paired[j][i]==0){
	const DFDCSegment* segment=packages[j][i];
	
	// Get the momentum and position at a specific z position
	DVector3 mom, pos;
	GetPositionAndMomentum(segment,pos,mom);

	// Empirical correction to the momentum 
	if (APPLY_MOMENTUM_CORRECTION){
	  double p_mag=mom.Mag();
	  mom.SetMag(p_mag*(1.+p_factor1/mom.Theta()+p_factor2));
	}
	
	// Create new track, starting with the current segment
	DTrackCandidate *track = new DTrackCandidate;
	track->rc=segment->rc;
	track->xc=segment->xc;
	track->yc=segment->yc;
	
	track->setPosition(pos);
	track->setMomentum(mom);
	track->setPID((segment->q > 0.0) ? PiPlus : PiMinus);
	track->Ndof=segment->Ndof;
	track->chisq=segment->chisq;
      
	track->AddAssociatedObject(segment);
	
	_data.push_back(track);
      }
    }
  }

  return NOERROR;
}


// Routine to do a crude match between fdc points and a helical approximation to
// the trajectory
double DTrackCandidate_factory_FDCCathodes::DocaSqToHelix(const DFDCPseudo *hit){
  double sperp=(hit->wire->origin.z()-zs)*cotl;
  double twoks=twokappa*sperp;
  double sin2ks=sin(twoks);
  double cos2ks=cos(twoks);
  double one_minus_cos2ks=1.-cos2ks;
 
  double x=xs+(cosphi*sin2ks-sinphi*one_minus_cos2ks)*one_over_twokappa;
  double y=ys+(sinphi*sin2ks+cosphi*one_minus_cos2ks)*one_over_twokappa;
  double dx=x-hit->xy.X();
  double dy=y-hit->xy.Y();

  return (dx*dx+dy*dy);
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
  double my_p=p;
  
  // Match to the next package
  double doca2_min=1e6,doca2;
  for (unsigned int j=0;j<package.size();j++){
    DFDCSegment *segment2=package[j];
    doca2=DocaSqToHelix(segment2->hits[0]);

    if (doca2<doca2_min){
      doca2_min=doca2;
      if(doca2<Match(p)){
	match=segment2;
	match_id=j;
      }
    }
  }

  if(DEBUG_HISTS){
    match_dist_fdc->Fill(p,doca2_min);
  }
  if (match!=NULL) return match;

  // If matching in the forward direction did not work, try 
  // matching backwards...
  doca2_min=1e6;
  for (unsigned int i=0;i<package.size();i++){
    DFDCSegment *segment2=package[i];
    GetPositionAndMomentum(segment2);
    doca2=DocaSqToHelix(segment->hits[segment->hits.size()-1]);
    if (doca2<doca2_min){
      doca2_min=doca2;
      if (doca2<Match(p)){
	match=segment2;
	match_id=i;
      }
    }       
  }
  if (match!=NULL) return match;

  // Match by centers of circles
  double circle_center_diff2_min=1e6;
  for (unsigned int j=0;j<package.size();j++){
    DFDCSegment *segment2=package[j];
    
    double dx=segment->xc-segment2->xc;
    double dy=segment->yc-segment2->yc;
    double circle_center_diff2=dx*dx+dy*dy;
    
    if (circle_center_diff2<circle_center_diff2_min){
      circle_center_diff2_min=circle_center_diff2;
      if (circle_center_diff2_min<4.0){
	match=segment2;
	match_id=j;
      }
    }
  }
  if (DEBUG_HISTS){
    match_center_dist2->Fill(my_p,circle_center_diff2_min);
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
  twokappa=FactorForSenseOfRotation*segment->q/segment->rc;

  one_over_twokappa=1./twokappa;
  cotl=1./my_tanl;

  // Useful intermediate variables
  double cosp=cos(my_phi0);
  double sinp=sin(my_phi0);
  double twoks=twokappa*(zs-z0)*cotl;
  double sin2ks=sin(twoks);
  double cos2ks=cos(twoks); 
  
  // Get Bfield
  double Bz=fabs(bfield->GetBz(xs,ys,zs));

  // Momentum
  p=0.003*Bz*segment->rc/cos(atan(my_tanl));

  cosphi=cosp*cos2ks-sinp*sin2ks;
  sinphi=sinp*cos2ks+cosp*sin2ks;

  return NOERROR;
}

// Routine to return momentum and position given the helical parameters and the
// z-component of the magnetic field
jerror_t 
DTrackCandidate_factory_FDCCathodes::GetPositionAndMomentum(
				        vector<const DFDCSegment *>segments,
							    DVector3 &pos,
							    DVector3 &mom){
  // Hit in the most upstream package
  const DFDCPseudo *hit=segments[0]->hits[segments[0]->hits.size()-1];
  double zhit=hit->wire->origin.z();
  double xhit=hit->xy.X();
  double yhit=hit->xy.Y();
  
  // Position
  double dz=1.;
  double zmin=zhit-dz;
  double phi1=atan2(yhit-yc,xhit-xc);
  double q_over_rc_tanl=q/(rc*tanl);
  double dphi_s=dz*q_over_rc_tanl;
  double dphi1=phi1-dphi_s;// was -
  double x=xc+rc*cos(dphi1);
  double y=yc+rc*sin(dphi1);
  pos.SetXYZ(x,y,zmin);

  dphi1*=-1.;
  if (FactorForSenseOfRotation*q<0) dphi1+=M_PI;

  // Find the average Bz
  double Bz=0.;
  double z=zmin;
  unsigned int num_segments=segments.size();
  double zmax=segments[num_segments-1]->hits[0]->wire->origin.z();
  unsigned int num_samples=20*num_segments;
  double one_over_denom=1./double(num_samples);
  dz=(zmax-zmin)*one_over_denom;
  for (unsigned int i=0;i<num_samples;i++){
    double my_dphi=phi1+(z-zmin)*q_over_rc_tanl;
    x=xc+rc*cos(my_dphi);
    y=yc+rc*sin(my_dphi);
    Bz+=bfield->GetBz(x,y,z);

    z+=dz;
  }
  Bz=fabs(Bz)*one_over_denom;
  
  
  // Momentum 
  double pt=0.003*Bz*rc; 
  double px=pt*sin(dphi1);
  double py=pt*cos(dphi1);
  double pz=pt*tanl;
  mom.SetXYZ(px,py,pz);

  return NOERROR;
}

// Routine to return momentum and position given the helical parameters and the
// z-component of the magnetic field
jerror_t 
DTrackCandidate_factory_FDCCathodes::GetPositionAndMomentum(
							    const DFDCSegment *segment,
							    DVector3 &pos,
							    DVector3 &mom){
  // Hit in the most upstream package
  const DFDCPseudo *hit=segment->hits[segment->hits.size()-1];
  double zhit=hit->wire->origin.z();
  double xhit=hit->xy.X();
  double yhit=hit->xy.Y();
  
  // Position
  double dz=1.;
  double zmin=zhit-dz;
  double phi1=atan2(yhit-segment->yc,xhit-segment->xc);
  double q_over_rc_tanl=segment->q/(segment->rc*segment->tanl);
  double dphi_s=dz*q_over_rc_tanl;
  double dphi1=phi1-dphi_s;// was -
  double x=segment->xc+segment->rc*cos(dphi1);
  double y=segment->yc+segment->rc*sin(dphi1);
  pos.SetXYZ(x,y,zmin);

  dphi1*=-1.;
  if (FactorForSenseOfRotation*segment->q<0) dphi1+=M_PI;

  // Find Bz at x,y,zmin
  double Bz=bfield->GetBz(x,y,zmin);
  
  // Momentum 
  double pt=0.003*Bz*segment->rc; 
  double px=pt*sin(dphi1);
  double py=pt*cos(dphi1);
  double pz=pt*segment->tanl;
  mom.SetXYZ(px,py,pz);

  return NOERROR;
}


// Routine to loop over segments in one of the packages, linking them with 
// segments in the package downstream of this package
void DTrackCandidate_factory_FDCCathodes::LinkSegments(unsigned int pack1,
						       vector<DFDCSegment *>packages[4], vector<pair<const DFDCSegment*,const DFDCSegment*> >&paired_segments,
						       vector<vector<int> >&is_paired){ 
  unsigned int match_id=0;
  unsigned int pack2=pack1+1;

  // Loop over the segments in package "pack1"
  for (unsigned int i=0;i<packages[pack1].size();i++){
    DFDCSegment *segment=packages[pack1][i];
    DFDCSegment *match2=NULL;

    // Try matching to the next package
    if (packages[pack2].size()>0 
	&& (match2=GetTrackMatch(segment,packages[pack2],match_id))!=NULL){
      
      pair<const DFDCSegment*,const DFDCSegment*> mypair;
      mypair.first=segment;
      mypair.second=match2;
      paired_segments.push_back(mypair);
      is_paired[pack2][match_id]=1;
      is_paired[pack1][i]=1;
    }
  
  }
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
    double d2=dx*dx+dy*dy;

    if (d2<Match(mom.Mag())) return true;
  }
  return false;
}

// Routine that tries to link a stray segment with an already existing track
// candidate
bool DTrackCandidate_factory_FDCCathodes::LinkStraySegment(const DFDCSegment *segment){
  // Loop over existing candidates looking for potential holes
  for (unsigned int i=0;i<_data.size();i++){
    bool got_segment_in_package=false;

    // Get the segments already associated with this track 
    vector<const DFDCSegment*>segments;
    _data[i]->GetT(segments);
    // Flag if segment is in a package that has already been used for this 
    // candidate
    for (unsigned int j=0;j<segments.size();j++){
      if (segments[j]->package==segment->package){
	got_segment_in_package=true;
	break;
      }
    }
    if (got_segment_in_package==false){
      // Try to link this segment to an existing candidate
      DVector3 pos=_data[i]->position();
      DVector3 mom=_data[i]->momentum();

      // Switch the direction of the momentum if we would need to backtrack to 
      // get to the segment
      if (segment->hits[0]->wire->origin.z()<pos.z()){
	mom=-1.0*mom;
      }

      if (GetTrackMatch(_data[i]->charge(),pos,mom,segment)){
	// Add the segment as an associated object to _data[i]
	_data[i]->AddAssociatedObject(segment);
   
	// Add the new segment and sort by z
	segments.push_back(segment);
	sort(segments.begin(),segments.end(),DTrackCandidate_segment_cmp_by_z);
	
	// Create fit object and add hits
	DHelicalFit fit;
	for (unsigned int m=0;m<segments.size();m++){
	  for (unsigned int k=0;k<segments[m]->hits.size();k++){
	    const DFDCPseudo *hit=segments[m]->hits[k];
	    fit.AddHit(hit);
	  }
	}
	
	// Redo the helical fit with the additional hits
	if (fit.FitTrackRiemann(segment->rc)==NOERROR){      	
	  rc=fit.r0;
	  tanl=fit.tanl;
	  xc=fit.x0;
	  yc=fit.y0;
	  q=FactorForSenseOfRotation*fit.h;
	  // Get position and momentum just upstream of first hit
	  GetPositionAndMomentum(segments,pos,mom);

	  // Empirical correction to the momentum
	  if (APPLY_MOMENTUM_CORRECTION){
	    double p_mag=mom.Mag();
	    mom.SetMag(p_mag*(1.+p_factor1/mom.Theta()+p_factor2));
	  }
	  
	  _data[i]->chisq=fit.chisq;
	  _data[i]->Ndof=fit.ndof;
	  _data[i]->setPID((q > 0.0) ? PiPlus : PiMinus);
	  _data[i]->setPosition(pos);
	  _data[i]->setMomentum(mom); 
	}

	return true;
      }
    }
  }
  return false;
}

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
#include "DHelicalFit.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>

//#define MATCH_RADIUS(p) (2.79+2.88/(p)/(p))
#define MATCH_RADIUS(p) (2.75+2./(p+0.25)/(p+0.25))
#define MAX_SEGMENTS 20
#define HALF_PACKAGE 6.0
#define FDC_OUTER_RADIUS 50.0 
#define BEAM_VAR 0.0208 // cm^2
#define HIT_CHI2_CUT 10.0
///
/// DTrackCandidate_factory_FDCCathodes::brun():
///
jerror_t DTrackCandidate_factory_FDCCathodes::brun(JEventLoop* eventLoop, 
						   int runnumber) {
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();
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
  MAX_R_VERTEX_LIMIT = 50.0;
  gPARMS->SetDefaultParameter("TRKFIND:MAX_R_VERTEX_LIMIT", MAX_R_VERTEX_LIMIT, description);

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
jerror_t DTrackCandidate_factory_FDCCathodes::evnt(JEventLoop *loop, int eventnumber)
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

  // Loop over all the packages to match to segments in packages downstream
  // of the current package
  for (unsigned int i=0;i<4;i++){
    if (packages[i].size()>0) LinkSegments(i,packages);
  }

  // Now collect stray segments
  // if (false)
  for (unsigned int i=0;i<4;i++){
    for (unsigned int k=0;k<packages[i].size();k++){
      DFDCSegment *segment=packages[i][k];

      // Make sure that this segment does not match one of the previously 
      // linked tracks
      if (LinkStraySegment(segment)==false){   
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
	const DFDCPseudo *segment_hit=segment->hits[segment->hits.size()-1];
	if (rc<10.0 && segment_hit->xy.Mod()<10.0 ){
	  DHelicalFit fit;
	  for (unsigned int n=0;n<segment->hits.size();n++){
	    const DFDCPseudo *hit=segment->hits[n];
	    fit.AddHit(hit);
	  }
	  fit.FitTrack_FixedZvertex(TARGET_Z);
	  tanl=fit.tanl;
	  rc=fit.r0;
	  xc=fit.x0;
	  yc=fit.y0;
	  fit.FindSenseOfRotation();
	  q=FactorForSenseOfRotation*fit.h;
	}	

	// Get the momentum and position at a specific z position
	DVector3 mom, pos(segment_hit->xy.X(),segment_hit->xy.Y(),
			  segment_hit->wire->origin.z()); 
	//Get position and momentum just upstream of this hit
	GetPositionAndMomentum(pos,mom);
	
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

  // If matching in the forward direction did not work, try 
  // matching backwards...
  if (match==NULL){
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
  }

  
  if (match==NULL){
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



// Obtain position and momentum at the exit of a given package using the 
// helical track model.
//
jerror_t DTrackCandidate_factory_FDCCathodes::GetPositionAndMomentum(DFDCSegment *segment,
					      DVector3 &pos, DVector3 &mom){
  // Position of track segment at last hit plane of package
  double x=segment->xc+segment->rc*cos(segment->Phi1);
  double y=segment->yc+segment->rc*sin(segment->Phi1);
  double z=segment->hits[0]->wire->origin.z();

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
  twokappa=FactorForSenseOfRotation*segment->q/segment->rc; 
  one_over_twokappa=1./twokappa;
  cotl=1./my_tanl;

  // Useful intermediate variables
  double cosp=cos(my_phi0);
  double sinp=sin(my_phi0);
  double twoks=twokappa*(z-z0)*cotl;
  double sin2ks=sin(twoks);
  double cos2ks=cos(twoks); 

  // Get Bfield
  double Bz=fabs(bfield->GetBz(x,y,z));

  // Momentum
  double pt=0.003*Bz*segment->rc;
  cosphi=cosp*cos2ks-sinp*sin2ks;
  sinphi=sinp*cos2ks+cosp*sin2ks;
  mom.SetXYZ(pt*cosphi,pt*sinphi,pt*my_tanl);

  return NOERROR;
}

// Routine to return momentum and position given the helical parameters and the
// z-component of the magnetic field
jerror_t 
DTrackCandidate_factory_FDCCathodes::GetPositionAndMomentum(DVector3 &pos,
							    DVector3 &mom){
  // Position
  double phi1=atan2(pos.y()-yc,pos.x()-xc);
  double q_over_rc_tanl=q/(rc*tanl);
  double dz=1.;
  double dphi_s=dz*q_over_rc_tanl;
  double dphi1=phi1-dphi_s;// was -
  double x=xc+rc*cos(dphi1);
  double y=yc+rc*sin(dphi1);
  double z=pos.z()-dz;
  pos.SetXYZ(x,y,z);

  dphi1*=-1.;
  if (FactorForSenseOfRotation*q<0) dphi1+=M_PI;

  // Find Bz
  double Bz=fabs(bfield->GetBz(x,y,z));

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
    vector<const DFDCSegment*>mysegments; 
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
    }
       
    // If we found a match, redo the helical fit with all the hits
    if (mysegments.size()>1){
      // Mark the segment in the first package as matched
      matched[i]=1;

      // Create the fit object and add the hits
      DHelicalFit fit;
      double max_r=0.,z_at_max_r=0.,x_at_max_r=0.,y_at_max_r=0.;
      for (unsigned int m=0;m<mysegments.size();m++){
	for (unsigned int n=0;n<mysegments[m]->hits.size();n++){
	  const DFDCPseudo *hit=mysegments[m]->hits[n];
	  fit.AddHit(hit);

	  double r=hit->xy.Mod();
	  if (r>max_r){
	    max_r=r;
	    z_at_max_r=hit->wire->origin.z();
	    x_at_max_r=hit->xy.X();
	    y_at_max_r=hit->xy.Y();
	  }
	}
      }
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

      // Create new track, starting with the current segment
      DTrackCandidate *track = new DTrackCandidate;

      // Get the momentum and position just upstream of first hit
      DVector3 mom,pos;
      GetPositionAndMomentum(mysegments,pos,mom);
    	
      // Empirical correction to the momentum
      if (APPLY_MOMENTUM_CORRECTION){
	double p_mag=mom.Mag();
	mom.SetMag(p_mag*(1.+p_factor1/mom.Theta()+p_factor2));
      }
      
      track->chisq=fit.chisq;
      track->Ndof=fit.ndof;
      track->setCharge(q);
      track->setPosition(pos);
      track->setMomentum(mom);
      
      for (unsigned int m=0;m<mysegments.size();m++){
	track->AddAssociatedObject(mysegments[m]);
      }
	
      _data.push_back(track); 
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
	  _data[i]->setCharge(q);
	  _data[i]->setPosition(pos);
	_data[i]->setMomentum(mom); 
	}

	return true;
      }
    }
  }
  return false;
}

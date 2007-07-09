//*****************************************************************************
// DTrackLinker_factory.cc - factory linking segments into tracks
//*****************************************************************************

#include "DANA/DApplication.h"
#include "DTrackLinker_factory.h"
#include "DMagneticFieldStepper.h"
#include "FDC/DFDCPseudo_factory.h"
#include "FDC/DFDCSegment_factory.h"

#define MATCH_RADIUS 5.0
#define MAX_SEGMENTS 20

bool DTrackLinker_cmp(const DFDCPseudo *a,const DFDCPseudo *b){
  return a->wire->layer > b->wire->layer;
}

///
/// DTrackLinker_factory::DTrackLinker_factory():
/// default constructor -- opens log file
///
DTrackLinker_factory::DTrackLinker_factory() {
        logFile = new ofstream("DTrackLinker_factory.log");
        _log = new JStreamLog(*logFile, "TRACKLINKER");
        *_log << "File initialized." << endMsg;
}

///
/// DTrackLinker_factory::~DTrackLinker_factory():
/// default destructor -- closes log file
///
DTrackLinker_factory::~DTrackLinker_factory() {
        logFile->close();
        delete logFile;
        delete _log;
}
///
/// DTrackLinker_factory::brun():
///
jerror_t DTrackLinker_factory::brun(JEventLoop* eventLoop, int eventNo) {
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();

  return NOERROR;
}
///
/// DTrackLinker_factory::evnt():
/// Routine where segments are linked into track candidates
///
jerror_t DTrackLinker_factory::evnt(JEventLoop* eventLoop, int eventNo) {
  vector<const DFDCSegment*>segments;
  eventLoop->Get(segments);
  
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
      // Create new track, starting with the current segment
      DTrackLinker *track = new DTrackLinker;
      DFDCSegment *segment=package[0][i];
      
      // Sign of the charge
      double kappa=segment->S(0,0);
      double q=kappa/fabs(kappa);
      
      // Start filling vector of points belonging to current track    
      vector<DFDCPseudo*>points; 
      points.assign(segment->hits.begin(),segment->hits.end());

      // Check that the tangent of the dip angle makes sense for FDC hits
      if (segment->S(3,0)<0) continue;
      
      // Try matching to package 2
      if (package[1].size()>0 && 
	  (match2=GetTrackMatch(q,zpackage[1],segment,package[1],match_id))
	  !=NULL){
	// Insert the points in the segment from package 2 into the track 
	points.insert(points.begin(),match2->hits.begin(),match2->hits.end());
	
	// remove the segment from the list 
	package[1].erase(package[1].begin()+match_id);

	// Try matching to package 3
	if (package[2].size()>0 && 
	    (match3=GetTrackMatch(q,zpackage[2],match2,package[2],match_id))
	    !=NULL){
	  // Insert the points in the segment from package 3 into the track 
	  points.insert(points.begin(),match3->hits.begin(),
			match3->hits.end());

	  // remove the segment from the list 
	  package[2].erase(package[2].begin()+match_id);

	  // Try matching to package 4
	  if (package[3].size()>0 && 
	      (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	       match_id))!=NULL){
	    // Insert the points in the segment from package 4 into the track 
	    points.insert(points.begin(),match4->hits.begin(),
			  match4->hits.end());

	    // remove the segment from the list 
	    package[3].erase(package[3].begin()+match_id);
	  }
	}
	// No match in package 3, try for 4
	else if(package[3].size()>0 && 
		(match4=GetTrackMatch(q,zpackage[3],match2,package[3],
		 match_id))!=NULL){
	  // Insert the points in the segment from package 4 into the track 
	  points.insert(points.begin(),match4->hits.begin(),
			match4->hits.end());

	  // remove the segment from the list 
	  package[3].erase(package[3].begin()+match_id);
	}
      }
      // No match in package 2, try for 3
      else if (package[2].size()>0 && 
	       (match3=GetTrackMatch(q,zpackage[2],segment,package[2],
		match_id))!=NULL){
	// Insert the points in the segment from package 3 into the track 
	points.insert(points.begin(),match3->hits.begin(),match3->hits.end());

	// remove the segment from the list 
	package[2].erase(package[2].begin()+match_id);
	
	// Try matching to package 4
	if (package[3].size()>0 && 
	    (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	     match_id))!=NULL){
	  // Insert the points in the segment from package 4 into the track 
	  points.insert(points.begin(),match4->hits.begin(),
			match4->hits.end());

	  // remove the segment from the list 
	  package[3].erase(package[3].begin()+match_id);
	}
      }    
      // No match to package 2 or 3, try 4
      else if (package[3].size()>0 && 
	       (match4=GetTrackMatch(q,zpackage[3],segment,package[3],
		match_id))!=NULL){
	// Insert the points in the segment from package 4 into the track 
	points.insert(points.begin(),match4->hits.begin(),
		      match4->hits.end());
	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);
      }
      
      // Sort hits by layer, starting at the most downstream plane
      std::sort(points.begin(),points.end(),DTrackLinker_cmp);
      track->points=points;
      
      track->S.ResizeTo(segment->S);
      track->S=segment->S;
      track->cov.ResizeTo(segment->cov);
      track->cov=segment->cov;
      
      _data.push_back(track); 
    }
  }
 
  // Next try to link segments starting at package 2
  if (package[1].size()>0 ){
    // Loop over segments in the 2nd package, matching them to segments in 
    // the third and fourth (most downstream) packages.
    for (unsigned int i=0;i<package[1].size();i++){
      // Create new track, starting with the current segment
      DTrackLinker *track = new DTrackLinker;
      DFDCSegment *segment=package[1][i];
      
      // Sign of the charge
      double kappa=segment->S(0,0);
      double q=kappa/fabs(kappa);
      
      // Start filling vector of points belonging to current track    
      vector<DFDCPseudo*>points; 
      points.assign(segment->hits.begin(),segment->hits.end());

      // Check that the tangent of the dip angle makes sense for FDC hits
      if (segment->S(3,0)<0) continue;

      // Try matching to package 3
      if (package[2].size()>0 && 
	  (match3=GetTrackMatch(q,zpackage[2],segment,package[2],match_id))
	  !=NULL){
	// Insert the points in the segment from package 3 into the track 
	points.insert(points.begin(),match3->hits.begin(),match3->hits.end());

	// remove the segment from the list 
	package[2].erase(package[2].begin()+match_id);

	// Try matching to package 4
	if (package[3].size()>0 && 
	    (match4=GetTrackMatch(q,zpackage[3],match3,package[3],
	     match_id))!=NULL){
	  // Insert the points in the segment from package 4 into the track 
	  points.insert(points.begin(),match4->hits.begin(),
			match4->hits.end());

	  // remove the segment from the list 
	  package[3].erase(package[3].begin()+match_id);
	}	
      }
      // No match in 3, try for 4
      else if (package[3].size()>0 && 
	       (match4=GetTrackMatch(q,zpackage[3],segment,package[3],
		match_id))!=NULL){
	// Insert the points in the segment from package 4 into the track 
	points.insert(points.begin(),match4->hits.begin(),match4->hits.end());

	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);
      }
    
      // Sort hits by layer, starting at the most downstream plane
      std::sort(points.begin(),points.end(),DTrackLinker_cmp);
      track->points=points;
      
      track->S.ResizeTo(segment->S);
      track->S=segment->S;
      track->cov.ResizeTo(segment->cov);
      track->cov=segment->cov;
      
      _data.push_back(track); 
    }
  }
  
  // Next try to link segments starting at package 3
  if(package[2].size()>0){
    // Loop over segments in the 3rd package, matching them to segments in 
    // the fourth (most downstream) packages.
    for (unsigned int i=0;i<package[2].size();i++){
      // Create new track, starting with the current segment
      DTrackLinker *track = new DTrackLinker;
      DFDCSegment *segment=package[2][i];
      
      // Sign of the charge
      double kappa=segment->S(0,0);
      double q=kappa/fabs(kappa);
      
      // Start filling vector of points belonging to current track    
      vector<DFDCPseudo*>points; 
      points.assign(segment->hits.begin(),segment->hits.end());

      // Check that the tangent of the dip angle makes sense for FDC hits
      if (segment->S(3,0)<0) continue;
      
      // Try matching to package 4
      if (package[3].size()>0 && 
	  (match4=GetTrackMatch(q,zpackage[3],segment,package[3],match_id))
	  !=NULL){
	// Insert the points in the segment from package 4 into the track 
	points.insert(points.begin(),match4->hits.begin(),
		      match4->hits.end());
	
	// remove the segment from the list 
	package[3].erase(package[3].begin()+match_id);
      }	
      
      // Sort hits by layer, starting at the most downstream plane
      std::sort(points.begin(),points.end(),DTrackLinker_cmp);
      track->points=points;
      
      track->S.ResizeTo(segment->S);
      track->S=segment->S;
      track->cov.ResizeTo(segment->cov);
      track->cov=segment->cov;
      
      _data.push_back(track); 
    }
  }

  // Now collect stray segments in package 4
  for (unsigned int k=0;k<package[3].size();k++){
    DTrackLinker *track = new DTrackLinker;
    DFDCSegment *segment=package[3][k];
    
    track->points=segment->hits;
    track->S.ResizeTo(segment->S);
    track->S=segment->S;
    track->cov.ResizeTo(segment->cov);
    track->cov=segment->cov;
    
    _data.push_back(track); 
  }

  return NOERROR;
}

// Swim track from one package to the next and look for a match to a segment
// in the new package
DFDCSegment *DTrackLinker_factory::GetTrackMatch(double q,double z, 
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
  GetPositionAndMomentum(segment,pos,mom);

  // Match to the next package by swimming the track through the field
  double diff_min=1000.,diff;
  if (stepper.SwimToPlane(pos,mom,origin,norm)==false){
    for (unsigned int j=0;j<package.size();j++){
      DFDCSegment *segment2=package[j];
	  
      // Skip to next segment if the sign of the charge is wrong
      double kappa=segment2->S(0,0);
      double q2=kappa/fabs(kappa);

      if (q2!=q) continue;
     
      // Check that the tangent of the dip angle makes sense for FDC hits
      if (segment2->S(3,0)<0) continue;
 
      double x2=segment2->hits[segment2->hits.size()-1]->x;
      double y2=segment2->hits[segment2->hits.size()-1]->y;
      diff=sqrt((pos(0)-x2)*(pos(0)-x2)+(pos(1)-y2)*(pos(1)-y2));

      if (diff<diff_min&&diff<MATCH_RADIUS){
	diff_min=diff;
	match=segment2;
	match_id=j;
      }
    }
  }
 
  return match;
}


// Obtain position and momentum at the exit of a given package using the 
// helical track model.
//
jerror_t DTrackLinker_factory::GetPositionAndMomentum(DFDCSegment *segment,
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

  pos.SetXYZ(x,y,z);
  mom.SetXYZ(px,py,pz);

  return NOERROR;
}

//------------------
// toString
//------------------
const string DTrackLinker_factory::toString(void)
{
  // Ensure our Get method has been called so _data is up to date
  Get();
  if(_data.size()<=0)return string(); // don't print anything if we have no data!
  
  return string();

}

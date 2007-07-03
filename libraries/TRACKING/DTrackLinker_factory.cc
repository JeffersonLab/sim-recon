//*****************************************************************************
// DTrackLinker_factory.cc - factory linking segments into tracks
//*****************************************************************************

#include "DANA/DApplication.h"
#include "DTrackLinker_factory.h"
#include "DMagneticFieldStepper.h"
#include "FDC/DFDCPseudo_factory.h"
#include "FDC/DFDCSegment_factory.h"

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
    
  DVector3 norm;  // normal to FDC planes
  norm.SetXYZ(0.,0.,1.);
  DFDCSegment *match=NULL;

  // Loop over segments in the first package, matching them to segments in the
  // second, third, and fourth (most downstream) packages.
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
 
    // Get the position and momentum at the exit of the package for the current
    // segment
    DVector3 pos,mom;
    GetPositionAndMomentum(segment,pos,mom);

    // Match to the next package by swimming the track through the field
    double diff_min=1000.,diff;
    for (unsigned int j=0;j<package[1].size();j++){
      DFDCSegment *segment2=package[1][j];
      DVector3 origin=segment2->hits[segment2->hits.size()-1]->wire->origin;

      // Skip to next segment if the sign of the charge is wrong
      double kappa=segment2->S(0,0);
      double q2=kappa/fabs(kappa);
      if (q2!=q) continue;

      DMagneticFieldStepper stepper(bfield, q); 
      if (stepper.SwimToPlane(pos,mom,origin,norm)==false){
	double x2=segment2->hits[segment2->hits.size()-1]->x;
	double y2=segment2->hits[segment2->hits.size()-1]->y;
	diff=sqrt((pos(0)-x2)*(pos(0)-x2)+(pos(1)-y2)*(pos(1)-y2));
	if (diff<diff_min){
	  diff_min=diff;
	  match=segment2;
	}
      }
    }
 
    // if we found a match, move on to the 3rd package
    if (match){
      // Insert the points in the segment from package 2 into the track 
      points.insert(points.begin(),match->hits.begin(),match->hits.end());

      // Get position and momentum at last plane in package 2
      GetPositionAndMomentum(match,pos,mom);
    
      // Reinitialize match pointer
      match=NULL; 
            
      // Match to the next package by swimming the track through the field
      diff_min=1000;
      for (unsigned int j=0;j<package[2].size();j++){
	DFDCSegment *segment2=package[2][j];
	DVector3 origin=segment2->hits[segment2->hits.size()-1]->wire->origin;

	// Skip to next segment if the sign of the charge is wrong
	double kappa=segment2->S(0,0);
	double q2=kappa/fabs(kappa);
        if (q2!=q) continue;

	DMagneticFieldStepper stepper(bfield, q); 
	if (stepper.SwimToPlane(pos,mom,origin,norm)==false){
	  double x2=segment2->hits[segment2->hits.size()-1]->x;
	  double y2=segment2->hits[segment2->hits.size()-1]->y;
	  diff=sqrt((pos(0)-x2)*(pos(0)-x2)+(pos(1)-y2)*(pos(1)-y2));
	  if (diff<diff_min){
	    diff_min=diff;
	    match=segment2;
	  }
	}
      }

      // If we have a match, move on to the 4th package
      if (match){
	// Insert the points in the segment from package 3 into the track 
	points.insert(points.begin(),match->hits.begin(),match->hits.end());

	// Get position and momentum at last plane in package 3
	GetPositionAndMomentum(match,pos,mom);
    
	// Reinitialize match pointer
	match=NULL; 

	// Match to the next package by swimming the track through the field
	diff_min=1000;
	for (unsigned int j=0;j<package[3].size();j++){
	  DFDCSegment *segment2=package[3][j];
	  DVector3 origin
	    =segment2->hits[segment2->hits.size()-1]->wire->origin;

	  // Skip to next segment if the sign of the charge is wrong
	  double kappa=segment2->S(0,0);
	  double q2=kappa/fabs(kappa);
          if (q2!=q) continue;

	  DMagneticFieldStepper stepper(bfield, q); 
	  if (stepper.SwimToPlane(pos,mom,origin,norm)==false){
	    double x2=segment2->hits[segment2->hits.size()-1]->x;
	    double y2=segment2->hits[segment2->hits.size()-1]->y;
	    diff=sqrt((pos(0)-x2)*(pos(0)-x2)+(pos(1)-y2)*(pos(1)-y2));
	    if (diff<diff_min){
	      diff_min=diff;
	      match=segment2;
	    }
	  }
	}
      
	if (match){
	  // Insert the points in the segment from package 4 into the track 
	  points.insert(points.begin(),match->hits.begin(),match->hits.end());
	}
      }
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
   
  return NOERROR;
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

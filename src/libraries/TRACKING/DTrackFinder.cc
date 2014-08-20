// $Id$
//
//    File: DTrackFinder.cc
// Created: Fri Aug 15 09:43:08 EDT 2014
// Creator: staylor (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "DTrackFinder.h"

#define CDC_MATCH_RADIUS 5.0

bool DTrackFinder_cdc_hit_cmp(const DCDCTrackHit *a,const DCDCTrackHit *b){
  return(a->wire->origin.Y()>b->wire->origin.Y());
}


//---------------------------------
// DTrackFinder    (Constructor)
//---------------------------------
DTrackFinder::DTrackFinder(JEventLoop *loop)
{
 
}

//---------------------------------
// ~DTrackFinder    (Destructor)
//---------------------------------
DTrackFinder::~DTrackFinder()
{

}

void DTrackFinder::Reset(void){
   
  axial_hits.clear();
  stereo_hits.clear();
  axial_segments.clear();
  cdc_tracks.clear();

}




void DTrackFinder::AddHit(const DCDCTrackHit *hit){  
  int ring=hit->wire->ring;
  if (ring<=4) axial_hits.push_back(cdc_hit_t(hit));
  else if (ring<=12) stereo_hits.push_back(cdc_hit_t(hit));
  else if (ring<=16) axial_hits.push_back(cdc_hit_t(hit));
  else if (ring<=24) stereo_hits.push_back(cdc_hit_t(hit));
  else axial_hits.push_back(cdc_hit_t(hit)); 
}


// Find segments in cdc axial layers
bool DTrackFinder::FindAxialSegments(void){
  if (axial_hits.size()==0) return false;

  // Group adjacent axial_hits into pairs
  vector<pair<unsigned int,unsigned int> > pairs;
  for (unsigned int i=0;i<axial_hits.size()-1;i++){
    for (unsigned int j=i+1;j<axial_hits.size();j++){ 
      const DCDCWire *first_wire=axial_hits[i].hit->wire;
      const DCDCWire *second_wire=axial_hits[j].hit->wire;
      int r1=first_wire->ring;
      int r2=second_wire->ring;
      int s1=first_wire->straw;
      int s2=second_wire->straw;
      double d=(first_wire->origin-second_wire->origin).Perp();

      if ((abs(r1-r2)<=2 && d<CDC_MATCH_RADIUS) 
	  || (abs(r1-r2)==0 && abs(s1-s2)==1)){
	pair <unsigned int,unsigned int> mypair(i,j);
	pairs.push_back(mypair);
      }
    }
  }
  // Link pairs of axial_hits together into segments
  for (unsigned int i=0;i<pairs.size();i++){
    if (axial_hits[pairs[i].first].used==false 
	&& axial_hits[pairs[i].second].used==false){
      vector<const DCDCTrackHit *>neighbors;
      unsigned int old=i;
      unsigned int old_first=pairs[old].first;
      unsigned int old_second=pairs[old].second;
      axial_hits[old_first].used=true;
      axial_hits[old_second].used=true;
      neighbors.push_back(axial_hits[old_first].hit);
      neighbors.push_back(axial_hits[old_second].hit);
      for (unsigned int j=i+1;j<pairs.size();j++){
	unsigned int first=pairs[j].first;
	unsigned int second=pairs[j].second;
	old_first=pairs[old].first;
	old_second=pairs[old].second;
	if ((axial_hits[old_first].used || axial_hits[old_second].used)
	    && (first==old_first || first==old_second || second==old_second
		|| second==old_first)){
	  if (axial_hits[first].used==false){
	    axial_hits[first].used=true;
	    neighbors.push_back(axial_hits[first].hit);
	  }
	  if (axial_hits[second].used==false){
	    axial_hits[second].used=true;
	    neighbors.push_back(axial_hits[second].hit);
	  }  
	  if (axial_hits[old_first].used==false){
	    axial_hits[old_first].used=true;
	    neighbors.push_back(axial_hits[old_first].hit);
	  }
	  if (axial_hits[old_second].used==false){
	    axial_hits[old_second].used=true;
	    neighbors.push_back(axial_hits[old_second].hit);
	  }
	}
	old=j;
      }
 
      sort(neighbors.begin(),neighbors.end(),DTrackFinder_cdc_hit_cmp);

      DVector3 dir=neighbors[neighbors.size()-1]->wire->origin
	-neighbors[0]->wire->origin;
      dir.SetMag(1.);
      
      dir.Print();

      axial_segments.push_back(cdc_segment_t(neighbors,dir));
    }

  }
  
  return true;
}



// Link axial segments together to form track candidates and match to stereo 
// hits
void DTrackFinder::LinkCDCSegments(void){
  unsigned int num_axial=axial_segments.size();
  for (unsigned int i=0;i<num_axial-1;i++){
    if (axial_segments[i].matched==false){
      DTrackFinder::cdc_track_t mytrack(axial_segments[i].hits);

      DVector3 pos0=axial_segments[i].hits[0]->wire->origin;
      DVector3 vhat=axial_segments[i].dir;

      for (unsigned int j=i+1;j<num_axial;j++){
	if (axial_segments[j].matched==false){
	  DVector3 pos1=axial_segments[j].hits[0]->wire->origin;
	  DVector3 dir1=axial_segments[j].hits[0]->wire->udir;
	  DVector3 diff=pos1-pos0;
	  double s=diff.Dot(vhat);
	  double d=(diff-s*vhat).Mag();

	  if (d<CDC_MATCH_RADIUS){
	    axial_segments[j].matched=true;	   
	    mytrack.axial_hits.insert(mytrack.axial_hits.end(),
				  axial_segments[j].hits.begin(),
				  axial_segments[j].hits.end());
	    sort(mytrack.axial_hits.begin(),mytrack.axial_hits.end(),
		 DTrackFinder_cdc_hit_cmp);
	    
	    vhat=mytrack.axial_hits[mytrack.axial_hits.size()-1]->wire->origin
	      -mytrack.axial_hits[0]->wire->origin;
	    vhat.SetMag(1.);
	  }
	}
      }
      //  Position of the first axial wire in the track  
      pos0=mytrack.axial_hits[0]->wire->origin;

      // Grab axial hits not associated with segments
      bool got_match=false;
      for (unsigned int j=0;j<axial_hits.size();j++){
	if (axial_hits[j].used==false){
	  if (MatchCDCHit(vhat,pos0,axial_hits[j].hit)){
	    axial_hits[j].used=true;
	    mytrack.axial_hits.push_back(axial_hits[j].hit);
	    got_match=true;
	  }
	}
      }
      // Resort if we added axial hits and recompute direction vector
      if (got_match){
	sort(mytrack.axial_hits.begin(),mytrack.axial_hits.end(),
	     DTrackFinder_cdc_hit_cmp);
	
	vhat=mytrack.axial_hits[mytrack.axial_hits.size()-1]->wire->origin
	  -mytrack.axial_hits[0]->wire->origin;
	vhat.SetMag(1.);   
	pos0=mytrack.axial_hits[0]->wire->origin;
      }
      
    
      // Now try to associate stereo hits with this track
      for (unsigned int j=0;j<stereo_hits.size();j++){
	if (stereo_hits[j].used==false){
	  if (MatchCDCHit(vhat,pos0,stereo_hits[j].hit)){
	    stereo_hits[j].used=true;
	    mytrack.stereo_hits.push_back(stereo_hits[j].hit);
	  }
	}
      }
      size_t num_stereo=mytrack.stereo_hits.size();
      size_t num_axial=mytrack.axial_hits.size();
      if (num_stereo>0 && num_stereo+num_axial>4){
	mytrack.dir=vhat;
	printf("-----------\n");
	vhat.Print();
	if (mytrack.FindStateVector()==NOERROR){
	  mytrack.S.Print();
	  cdc_tracks.push_back(mytrack);
	}
      }
    }
  }

}


// Match a CDC hit with a line starting at pos0 going in the vhat direction
bool DTrackFinder::MatchCDCHit(const DVector3 &vhat,const DVector3 &pos0,
			       const DCDCTrackHit *hit){  
  DVector3 pos1=hit->wire->origin;
  DVector3 uhat=hit->wire->udir;
  DVector3 diff=pos1-pos0;
  double vhat_dot_uhat=vhat.Dot(uhat);
  double scale=1./(1.-vhat_dot_uhat*vhat_dot_uhat);
  double s=scale*(vhat_dot_uhat*diff.Dot(vhat)-diff.Dot(uhat));
  double t=scale*(diff.Dot(vhat)-vhat_dot_uhat*diff.Dot(uhat));
  double d=(diff+s*uhat-t*vhat).Mag();

  if (d<CDC_MATCH_RADIUS) return true;

  return false;
}

// Compute initial guess for state vector (x,y,tx,ty) for a track in the CDC
// by fitting a line to the intersections between the line in the xy plane and 
// the stereo wires.
jerror_t DTrackFinder::cdc_track_t::FindStateVector(void){  
  // Parameters for line in x-y plane
  double vx=this->dir.x();
  double vy=this->dir.y();
  DVector3 pos0=this->axial_hits[0]->wire->origin;
  double xa=pos0.x();
  double ya=pos0.y();

  double sumv=0,sumx=0,sumy=0,sumz=0,sumxx=0,sumyy=0,sumxz=0,sumyz=0;
  for (unsigned int i=0;i<this->stereo_hits.size();i++){
    // Intersection of line in xy-plane with this stereo straw
    DVector3 origin_s=this->stereo_hits[i]->wire->origin;
    DVector3 dir_s=this->stereo_hits[i]->wire->udir;
    double ux_s=dir_s.x();
    double uy_s=dir_s.y();
    double dx=xa-origin_s.x();
    double dy=ya-origin_s.y();
    double s=(dx*vy-dy*vx)/(ux_s*vy-uy_s*vx);
    DVector3 pos1=origin_s+s*dir_s;
    double x=pos1.x(),y=pos1.y(),z=pos1.z();
    
    if (z>17.0 && z<167.0){ // Check for CDC dimensions
      sumv+=1.;
      sumx+=x;
      sumxx+=x*x;
      sumy+=y;
      sumyy+=y*y;
      sumz+=z;
      sumxz+=x*z;
      sumyz+=y*z;
    }
  }
  const double EPS=1e-3;
  double xdenom=sumv*sumxz-sumx*sumz;
  if (fabs(xdenom)<EPS) return VALUE_OUT_OF_RANGE;
 
  double ydenom=sumv*sumyz-sumy*sumz;
  if (fabs(ydenom)<EPS) return VALUE_OUT_OF_RANGE;
 
  double xtemp=sumv*sumxx-sumx*sumx;
  double xslope=xtemp/xdenom;
  double ytemp=sumv*sumyy-sumy*sumy;
  double yslope=ytemp/ydenom;

  //  double z0x=(sumxx*sumz-sumx*sumxz)/xtemp;
  double z0y=(sumyy*sumz-sumy*sumyz)/ytemp;
  
  // Increment just beyond point largest in y
  double delta_z=(yslope>0)?0.5:-0.5;

  //Starting z position
  this->z=z0y+ya/yslope+delta_z;

  this->S(state_x)=xa+xslope*delta_z;
  this->S(state_y)=ya+yslope*delta_z;
  this->S(state_tx)=xslope;
  this->S(state_ty)=yslope;

  return NOERROR;

}

// Given state vector S, find doca to wire given by origin and wdir
double DTrackFinder::FindDoca(double z,const DMatrix4x1 &S,const DVector3 &wdir,
			      const DVector3 &origin) const{
  DVector3 pos(S(state_x),S(state_y),z);
  DVector3 diff=pos-origin;
  
  DVector3 uhat(S(state_tx),S(state_ty),1.);
  uhat.SetMag(1.); 
  DVector3 vhat=wdir;
  //  vhat.SetMag(1.);

  double vhat_dot_diff=diff.Dot(vhat);
  double uhat_dot_diff=diff.Dot(uhat);
  double uhat_dot_vhat=uhat.Dot(vhat);
  double D=1.-uhat_dot_vhat*uhat_dot_vhat;
  double N=uhat_dot_vhat*vhat_dot_diff-uhat_dot_diff;
  double N1=vhat_dot_diff-uhat_dot_vhat*uhat_dot_diff;
  double scale=1./D;
  double s=scale*N;
  double t=scale*N1;
  
  diff+=s*uhat-t*vhat;
  return diff.Mag();
}




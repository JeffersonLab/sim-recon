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

bool DTrackFinder_fdc_hit_cmp(const DTrackFinder::fdc_hit_t &a,
			      const DTrackFinder::fdc_hit_t &b){
  return (a.hit->wire->origin.z()>b.hit->wire->origin.z());
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

  fdc_hits.clear();
  for (unsigned int i=0;i<4;i++) fdc_segments[i].clear();
  fdc_tracks.clear();

}

void DTrackFinder::AddHit(const DFDCPseudo *hit){
  fdc_hits.push_back(fdc_hit_t(hit));
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

      axial_segments.push_back(cdc_segment_t(neighbors,dir));
    }

  }
  
  return true;
}



// Link axial segments together to form track candidates and match to stereo 
// hits
bool DTrackFinder::LinkCDCSegments(void){
  unsigned int num_axial=axial_segments.size();
  if (num_axial<1) return false;
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
	if (mytrack.FindStateVector()==NOERROR){
	  cdc_tracks.push_back(mytrack);
	}
      }
    }
  }
  return true;
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
			      const DVector3 &origin,DVector3 *poca) const{
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
  
  if (poca!=NULL) *poca=pos+s*uhat;

  diff+=s*uhat-t*vhat;
  return diff.Mag();
}


// Find segments by associating adjacent hits within a package together.
bool DTrackFinder::FindFDCSegments(void){
  if (fdc_hits.size()==0) return false;
  unsigned int num_hits=fdc_hits.size();
  const double MATCH_RADIUS=2.0;
  const double ADJACENT_MATCH_RADIUS=1.0;

  // Order points by z
  sort(fdc_hits.begin(),fdc_hits.end(),DTrackFinder_fdc_hit_cmp);

  // Put indices for the first point in each plane before the most downstream
  // plane in the vector x_list.
  double old_z=fdc_hits[0].hit->wire->origin.z();
  vector<unsigned int>x_list;
  x_list.push_back(0);
  for (unsigned int i=0;i<num_hits;i++){
    if (fdc_hits[i].hit->wire->origin.z()!=old_z){
      x_list.push_back(i);
    }
    old_z=fdc_hits[i].hit->wire->origin.z();
  }
  x_list.push_back(num_hits); 

  unsigned int start=0;
  // loop over the start indices, starting with the first plane
  while (start<x_list.size()-1){
    // Now loop over the list of track segment start fdc_hits
    for (unsigned int i=x_list[start];i<x_list[start+1];i++){
      if (fdc_hits[i].used==false){
	fdc_hits[i].used=true;
	
	// Point in the current plane in the package 
	DVector2 XY=fdc_hits[i].hit->xy;
	double z=fdc_hits[i].hit->wire->origin.z();
	
	// Create list of nearest neighbors
	vector<const DFDCPseudo*>neighbors;
	neighbors.push_back(fdc_hits[i].hit);
	unsigned int match=0;
	double delta,delta_min=1000.;
	for (unsigned int k=0;k<x_list.size()-1;k++){
	  delta_min=1000.;
	  match=0;
	  for (unsigned int m=x_list[k];m<x_list[k+1];m++){
	    delta=(XY-fdc_hits[m].hit->xy).Mod();
	    double delta_z=fabs(z-fdc_hits[m].hit->wire->origin.z());
	    if (delta<delta_min){
	      delta_min=delta;
	      if (delta<MATCH_RADIUS && delta_z<10.0) match=m;
	    }
	  }
	  if (fdc_hits[match].used==false){
	    XY=fdc_hits[match].hit->xy;
	    fdc_hits[match].used=true;
	    neighbors.push_back(fdc_hits[match].hit);	  
	  }
	}
	unsigned int num_neighbors=neighbors.size();

	// Look for hits adjacent to the ones we have in our segment candidate
	for (unsigned int k=0;k<num_hits;k++){
	  if (fdc_hits[k].used==false){
	    for (unsigned int j=0;j<num_neighbors;j++){
	      delta=(fdc_hits[k].hit->xy-neighbors[j]->xy).Mod();

	      if (delta<ADJACENT_MATCH_RADIUS && 
		  abs(neighbors[j]->wire->wire-fdc_hits[k].hit->wire->wire)<=1
		  && neighbors[j]->wire->origin.z()==fdc_hits[k].hit->wire->origin.z()){
	        fdc_hits[k].used=true;
		neighbors.push_back(fdc_hits[k].hit);
	      }      
	    }
	  }
	} // loop looking for hits adjacent to hits on segment

	if (neighbors.size()>3){
	  unsigned int packNum=(neighbors[0]->wire->layer-1)/6;
	  fdc_segments[packNum].push_back(fdc_segment_t(neighbors));
	}
      }
    }// loop over start points within a plane
    
    // Look for a new plane to start looking for a segment
    while (start<x_list.size()-1){
      if (fdc_hits[x_list[start]].used==false) break;
      start++;
    }

  }

  return true;
}

// Link segments from package to package by doing straight-line projections
bool DTrackFinder::LinkFDCSegments(void){
  // Cuts for matching
  const double MATCH_RADIUS=2.0;
  const double LINK_MATCH_RADIUS=7.0;

  // Vector to store hits for the linked segments
  vector<const DFDCPseudo *>myhits;

  // loop over packages
  for (unsigned int i=0;i<4;i++){
    for (unsigned int j=0;j<fdc_segments[i].size();j++){
      if (fdc_segments[i][j].matched==false){
	unsigned i_plus_1=i+1; 
	if (i_plus_1<4){
	  double tx=fdc_segments[i][j].S(state_tx);
	  double ty=fdc_segments[i][j].S(state_ty);
	  double x0=fdc_segments[i][j].S(state_x);
	  double y0=fdc_segments[i][j].S(state_y);
	  
	  for (unsigned int k=0;k<fdc_segments[i_plus_1].size();k++){
	    if (fdc_segments[i_plus_1][k].matched==false){
	      double z=fdc_segments[i_plus_1][k].hits[0]->wire->origin.z();
	      DVector2 proj(x0+tx*z,y0+ty*z);
	    
	      if ((proj-fdc_segments[i_plus_1][k].hits[0]->xy).Mod()<LINK_MATCH_RADIUS){
		fdc_segments[i_plus_1][k].matched=true;
		myhits.insert(myhits.end(),fdc_segments[i_plus_1][k].hits.begin(),
				fdc_segments[i_plus_1][k].hits.end());
		
		unsigned int i_plus_2=i_plus_1+1;
		if (i_plus_2<4){
		  tx=fdc_segments[i_plus_1][k].S(state_tx);
		  ty=fdc_segments[i_plus_1][k].S(state_ty);
		  x0=fdc_segments[i_plus_1][k].S(state_x);
		  y0=fdc_segments[i_plus_1][k].S(state_y);
		  
		  for (unsigned int m=0;m<fdc_segments[i_plus_2].size();m++){
		    if (fdc_segments[i_plus_2][m].matched==false){
		      z=fdc_segments[i_plus_2][m].hits[0]->wire->origin.z();
		      proj.Set(x0+tx*z,y0+ty*z);
		      
		      if ((proj-fdc_segments[i_plus_2][m].hits[0]->xy).Mod()<LINK_MATCH_RADIUS){
			fdc_segments[i_plus_2][m].matched=true;
			myhits.insert(myhits.end(),fdc_segments[i_plus_2][m].hits.begin(),
				      fdc_segments[i_plus_2][m].hits.end());
			
			unsigned int i_plus_3=i_plus_2+1;
			if (i_plus_3<4){
			  tx=fdc_segments[i_plus_2][m].S(state_tx);
			  ty=fdc_segments[i_plus_2][m].S(state_ty);
			  x0=fdc_segments[i_plus_2][m].S(state_x);
			  y0=fdc_segments[i_plus_2][m].S(state_y);
			  
			  for (unsigned int n=0;n<fdc_segments[i_plus_3].size();n++){
			    if (fdc_segments[i_plus_3][n].matched==false){
			      z=fdc_segments[i_plus_3][n].hits[0]->wire->origin.z();
			      proj.Set(x0+tx*z,y0+ty*z);
			      
			      if ((proj-fdc_segments[i_plus_3][n].hits[0]->xy).Mod()<LINK_MATCH_RADIUS){
				fdc_segments[i_plus_3][n].matched=true;
				myhits.insert(myhits.end(),fdc_segments[i_plus_3][n].hits.begin(),
					      fdc_segments[i_plus_3][n].hits.end());
				
				break;
			      } // matched a segment
			    }
			  }  // loop over last set of segments
			} // if we have another package to loop over
			break;
		      } // matched a segment
		    }
		  } // loop over second-to-last set of segments
		}
		break;
	      } // matched a segment
	    }
	  } // loop over third-to-last set of segments
	}	
	if (myhits.size()>0){ 
	  myhits.insert(myhits.begin(),fdc_segments[i][j].hits.begin(),fdc_segments[i][j].hits.end());	
	  fdc_tracks.push_back(fdc_segment_t(myhits));
	}	  
	myhits.clear();
      } // check if we have already used this segment
    } // loop over first set of segments
  } // loop over packages

  // Try to link tracks together
  if (fdc_tracks.size()>1){
    for (unsigned int i=0;i<fdc_tracks.size()-1;i++){
      DMatrix4x1 S=fdc_tracks[i].S;
      size_t last_index_1=fdc_tracks[i].hits.size()-1;
      int first_pack_1=(fdc_tracks[i].hits[0]->wire->layer-1)/6;
      int last_pack_1=(fdc_tracks[i].hits[last_index_1]->wire->layer-1)/6;
      for (unsigned int j=i+1;j<fdc_tracks.size();j++){
	size_t last_index_2=fdc_tracks[j].hits.size()-1;
	int first_pack_2=(fdc_tracks[j].hits[0]->wire->layer-1)/6;
	int last_pack_2=(fdc_tracks[j].hits[last_index_2]->wire->layer-1)/6;

	if (last_pack_1<first_pack_2 || first_pack_1 > last_pack_2){
	  double z=fdc_tracks[j].hits[0]->wire->origin.z();
	  DVector2 proj(S(state_x)+z*S(state_tx),S(state_y)+z*S(state_ty));
	  double diff=(fdc_tracks[j].hits[0]->xy-proj).Mod();

	  if (diff<MATCH_RADIUS){
	    // Combine the hits from the two tracks and recompute the state 
	    // vector S
	    if (last_pack_1<first_pack_2){
	      fdc_tracks[i].hits.insert(fdc_tracks[i].hits.end(),
					fdc_tracks[j].hits.begin(),
					fdc_tracks[j].hits.end());
	    }
	    else{
	      fdc_tracks[i].hits.insert(fdc_tracks[i].hits.begin(),
					fdc_tracks[j].hits.begin(),
					fdc_tracks[j].hits.end());
	    }
	    fdc_tracks[i].FindStateVector();
	    
	    // Drop the second track from the list 
	    fdc_tracks.erase(fdc_tracks.begin()+j);
	    break;
	  }
	}
      }
    } // loop over tracks 
  } // check if we have more than one track

  // Try to attach unmatched segments to tracks
  for (unsigned int i=0;i<fdc_tracks.size();i++){
    DMatrix4x1 S=fdc_tracks[i].S;
    size_t last_index_1=fdc_tracks[i].hits.size()-1;
    int first_pack_1=(fdc_tracks[i].hits[0]->wire->layer-1)/6;
    int last_pack_1=(fdc_tracks[i].hits[last_index_1]->wire->layer-1)/6; 
    for (unsigned int j=0;j<4;j++){
      for (unsigned int k=0;k<fdc_segments[j].size();k++){
	if (fdc_segments[j][k].matched==false){
	  int pack_2=(fdc_segments[j][k].hits[0]->wire->layer-1)/6;
	  if (pack_2<first_pack_1 || pack_2>last_pack_1){
	    double z=fdc_segments[j][k].hits[0]->wire->origin.z();
	    DVector2 proj(S(state_x)+z*S(state_tx),S(state_y)+z*S(state_ty));
	    double diff=(fdc_segments[j][k].hits[0]->xy-proj).Mod();

	    if (diff<MATCH_RADIUS){
	      fdc_segments[j][k].matched=true;
	      
	      // Add hits and recompute S
	      if (pack_2<first_pack_1){
		fdc_tracks[i].hits.insert(fdc_tracks[i].hits.begin(),
					  fdc_segments[j][k].hits.begin(),
					  fdc_segments[j][k].hits.end());
	      }
	      else {
		fdc_tracks[i].hits.insert(fdc_tracks[i].hits.end(),
					  fdc_segments[j][k].hits.begin(),
					  fdc_segments[j][k].hits.end());
		
	      }
	      fdc_tracks[i].FindStateVector();
	    }
	  }
	} // check if already matched to other segments
      } // loop over segments in package
    } // loop over packages
  } //loop over existing tracks
  

  return true;
}



// Use linear regression on the hits to obtain a first guess for the state
// vector.  Method taken from Numerical Recipes in C.
DMatrix4x1 
DTrackFinder::fdc_segment_t::FindStateVector(void) const {
  double S1=0.;
  double S1z=0.;
  double S1y=0.;
  double S1zz=0.;
  double S1zy=0.;  
  double S2=0.;
  double S2z=0.;
  double S2x=0.;
  double S2zz=0.;
  double S2zx=0.;

  double sig2v=0.25; // rough guess;

  for (unsigned int i=0;i<hits.size();i++){
    double cosa=hits[i]->wire->udir.y();
    double sina=hits[i]->wire->udir.x();
    double x=hits[i]->xy.X();
    double y=hits[i]->xy.Y();
    double z=hits[i]->wire->origin.z();
    double sig2x=cosa*cosa/12+sina*sina*sig2v;
    double sig2y=sina*sina/12+cosa*cosa*sig2v;
    double one_over_var1=1/sig2y;
    double one_over_var2=1/sig2x;

    S1+=one_over_var1;
    S1z+=z*one_over_var1;
    S1y+=y*one_over_var1;
    S1zz+=z*z*one_over_var1;
    S1zy+=z*y*one_over_var1;    
    
    S2+=one_over_var2;
    S2z+=z*one_over_var2;
    S2x+=x*one_over_var2;
    S2zz+=z*z*one_over_var2;
    S2zx+=z*x*one_over_var2;
  }
  double D1=S1*S1zz-S1z*S1z;
  double y_intercept=(S1zz*S1y-S1z*S1zy)/D1;
  double y_slope=(S1*S1zy-S1z*S1y)/D1;
  double D2=S2*S2zz-S2z*S2z;
  double x_intercept=(S2zz*S2x-S2z*S2zx)/D2;
  double x_slope=(S2*S2zx-S2z*S2x)/D2;

  return DMatrix4x1(x_intercept,y_intercept,x_slope,y_slope);
}

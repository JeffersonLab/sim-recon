// $Id$
//
//    File: DVertex_factory.cc
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <PID/DChargedTrack.h>
#include <TMath.h>
#include "DVertex_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DVertex_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DVertex_factory::brun(jana::JEventLoop *loop, int runnumber)
{
  // Initialize to reasonable defaults
  target_length = 0.0;
  target_z = 0.0;
  
  // Get Target parameters from XML
  DApplication *dapp = dynamic_cast<DApplication*> (loop->GetJApplication());
  DGeometry *geom = dapp ? dapp->GetDGeometry(runnumber):NULL;
  if(geom){
    geom->GetTargetLength(target_length);
    geom->GetTargetZ(target_z);
  }
  
  // Get the particle ID algorithms
  vector<const DParticleID *> pid_algorithms;
  loop->Get(pid_algorithms);
  if(pid_algorithms.size()<1){
    _DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
    return RESOURCE_UNAVAILABLE;
  }
  // Drop the const qualifier from the DParticleID pointer (I'm surely going to hell for this!)
  pid_algorithm = const_cast<DParticleID*>(pid_algorithms[0]);
  
  // Warn user if something happened that caused us NOT to get a pid_algorithm object pointer
  if(!pid_algorithm){
    _DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
    return RESOURCE_UNAVAILABLE;
  }
  

  // Specify the size of the vertexInfo pool.
  // This is not a strict limit. Rather, it is the size the pool will be reduced
  // to if it has grown larger on the previous event and the current event does
  // not require this many. This prevents memory-leak-like behavior when running
  // many threads where each thread keeps allocating bigger and bigger pools as
  // it comes across slightly busier and busier events.
  MAX_VERTEXINFOS = 10;
  
  Nbinst = 600;
  tmin = -100.0;
  tmax = 500.0;
  Nbinsz = 50;
  zmin = 0.0;
  zmax = 100.0;
  
       
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DVertex_factory::evnt(JEventLoop *loop, int eventnumber)
{
  this->eventnumber=eventnumber;

  // Get list of all charged tracks
  vector<const DTrackTimeBased*> tracks;
  loop->Get(tracks);
  
  // group according to candidate id
  vector<vector<const DTrackTimeBased*> >tracks_by_candidate;
  pid_algorithm->GroupTracks(tracks,tracks_by_candidate);
  
  // Make list of vertices
  MakeVertices(tracks_by_candidate);

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DVertex_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DVertex_factory::fini(void)
{
	return NOERROR;
}


// Form vertices from grouping charged particle tracks together according to
// proximity in time and position of the closest approach to the beam line
jerror_t DVertex_factory::MakeVertices(vector<vector<const DTrackTimeBased*> >&tracks_by_candidate){
  // To minimize memory usage and time in allocation, we maintain a
  // pool of vertexInfo_t objects. Make sure the pool is large enough to hold
  // all of the charged particles we have this event. 
  unsigned int Nparticles_total =tracks_by_candidate.size();
  for(unsigned int i=vertexInfos_pool.size(); i<Nparticles_total; i++){
    vertexInfo_t *pi = new vertexInfo_t();
    
    pi->SetLimits(tmin, tmax, zmin, zmax, Nbinst, Nbinsz);
    vertexInfos_pool.push_back(pi);
  }
  // Periodically delete some vertexInfo_t objects if the pool gets too large.
  // This prevents memory-leakage-like behavor.
  if((Nparticles_total < MAX_VERTEXINFOS) && (vertexInfos_pool.size()>MAX_VERTEXINFOS)){
    for(unsigned int i=MAX_VERTEXINFOS; i<vertexInfos_pool.size(); i++)delete vertexInfos_pool[i];
    vertexInfos_pool.resize(MAX_VERTEXINFOS);
  }
  
  // Vector to hold list of vertexInfo_t objects for all charged tracks
  vector<vertexInfo_t*> vertices;

  // Assign and fill vertexInfo_t objects for each charged track
  for(unsigned int i=0; i<tracks_by_candidate.size(); i++){
    vertexInfo_t *pi = vertexInfos_pool[vertices.size()];
    // Use the fit result with the highest figure of merit
    FillVertexInfoChargedTrack(pi, &tracks_by_candidate[i]);
    vertices.push_back(pi);
  }

  // Each charged particle has a histogram of t vs.z
  // values filled using approriate uncertainties (no covariance). We
  // can now use this list to identify resonances in the t/z plane which
  // indicate a vertex location. Particles within 3sigma in both t and
  // z of the resonance will be grouped together as belonging to the 
  // same vertex. We loop until all particles have been assigned
  // to a group, even if that means assigning particles to their own
  // "group of one".
  
  // Loop until all particles have been assigned to a group.
  vector< vector<vertexInfo_t *> > groups;
  while(!AllInGroups(vertices)){
    // Make a list of all particles that have not been assigned
    // to a group
    vector<const DHoughFind*> unassigned;
    for(unsigned int i=0; i<vertices.size(); i++){
      if(!vertices[i]->is_in_group)unassigned.push_back(vertices[i]);
    }
    // Find the maximum t,z coordinate by adding all unassigned
    // particle's histos together
    DVector2 maxloc = DHoughFind::GetMaxBinLocation(unassigned);
    
    if(debug_level>0)_DBG_<<"Location of maximum: t="<<maxloc.X()<<"  z="<<maxloc.Y()<<endl;		

    // Loop over all unassigned particles, assigning any within
    // 3 sigma in both t and z to the new group. We loop over
    // the vertices vector just because it saves a dynamic_cast
    // if we were to use the unassigned vector.
    vector<vertexInfo_t *> new_group;
    for(unsigned int i=0; i<vertices.size(); i++){
      vertexInfo_t *pi = vertices[i];
      if(pi->is_in_group)continue;
			
      double delta_t = fabs(maxloc.X() - pi->t);
      if(delta_t/pi->sigmat > 3.0)continue;
      
      double delta_z = fabs(maxloc.Y() - pi->z);
      if(delta_z/pi->sigmaz > 3.0)continue;
      
      // Assign this particle to the group
      pi->is_in_group=true;
      new_group.push_back(pi);
    }
    
    // At this point it's possible (but hopefully unlikely) that the
    // maximum in the t,z sum histo was generated at an in-between place
    // with no single particle nearby. In that case, the new_group is
    // empty, even though there are unassigned particles. The best we
    // can do here is to assign one particle to the new_group and hope
    // that the next iteration groups the remaining ones appropriately.
    // To try and minimize the chances of placing a particle from the
    // L1 trigger event in its own group, we choose the particle with a time
    // the furthest away from t=0.
    if(new_group.size()==0){
      vertexInfo_t *pi_with_max_t = NULL;
      double delta_t_max=0.0;
      for(unsigned int i=0; i<vertices.size(); i++){
	vertexInfo_t *pi = vertices[i];
	if(pi->is_in_group)continue;
			
	double delta_t = fabs(maxloc.X() - pi->t);
	if(delta_t>delta_t_max || pi_with_max_t==NULL){
	  delta_t_max = delta_t;
	  pi_with_max_t = pi;
	}
      }
			
      if(pi_with_max_t==NULL){
	_DBG_<<"pi_with_max_t==NULL. This should never happen! Complain to davidl@jlab.org"<<endl;
	_DBG_<<"event number: "<<eventnumber<<endl;
	break;
			}
			
      new_group.push_back(pi_with_max_t);
    }
		
    // Set the is_in_group flags for all of the members of the new group
    for(unsigned int i=0; i<new_group.size(); i++)new_group[i]->is_in_group = true;
    
    // Add the new group to the list of groups
    groups.push_back(new_group);  
  }

  // OK, we've now grouped the particles together into groups. Create a new
  // DVertex for each group	
  for (unsigned int i=0;i<groups.size();i++){
    vector<vertexInfo_t *> &group = groups[i];

    DVertex *my_vertex = new DVertex;
    my_vertex->x.SetXYZT(0.0, 0.0, target_z,0.0);
    my_vertex->cov.ResizeTo(3,3);
    my_vertex->beamline_used = false;

    // Simply average POCA to beamline for all tracks	
    DVector3 temp;
    double t0=0.;
    double sum_invar=0, invar=0;
    for(unsigned int j=0; j<group.size(); j++){
      const DTrackTimeBased *trk=(*group[j]->hypotheses)[0];   
      invar=1./(group[j]->sigmat*group[j]->sigmat);
      sum_invar+=invar;
      t0+=group[j]->t*invar;
      temp+=trk->position();
    }
    t0*=1./sum_invar;
    temp*=1./double(group.size());
    my_vertex->x.SetVect(temp);
    my_vertex->x.SetT(t0);

    // Add list of tracks used to create this vertex
    for(unsigned int j=0; j<group.size(); j++){
      vector<DVertex::track_info_t>track_infos;
      vector<const DTrackTimeBased*>tracks=*(group[j]->hypotheses);
      for (unsigned int m=0;m<tracks.size();m++){
	DVertex::track_info_t my_track_info;
	my_track_info.track=tracks[m];
	my_track_info.FOM=my_track_info.tprojected=0.;
	track_infos.push_back(my_track_info);
      }
      my_vertex->hypotheses.push_back(track_infos);
    } 
    
    _data.push_back(my_vertex);
  }

  
  return NOERROR;
}


//------------------
// FillVertexInfoChargedTrack
//------------------
void DVertex_factory::FillVertexInfoChargedTrack(DVertex_factory::vertexInfo_t *pi, vector<const DTrackTimeBased *>*hypotheses)
{
	pi->Reset();
	pi->hypotheses = hypotheses;

	pi->t = (*hypotheses)[0]->t0();
	pi->sigmat = (*hypotheses)[0]->t0_err();
	pi->z = (*hypotheses)[0]->z();
	pi->sigmaz = 0.8/sin((*hypotheses)[0]->momentum().Theta()); // in cm.  For now, use 3mm wide angle track resolution scaled by sin(theta)

	pi->Fill(pi->t, pi->sigmat, pi->z, pi->sigmaz);
}

//------------------
// AllInGroups
//------------------
bool DVertex_factory::AllInGroups(vector<vertexInfo_t*> &vertices)
{
  for(unsigned int i=0; i<vertices.size(); i++){
    if(!vertices[i]->is_in_group)return false;
  }
  return true;
}


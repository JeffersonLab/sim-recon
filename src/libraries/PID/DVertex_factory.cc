// $Id$
//
//    File: DVertex_factory.cc
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <TROOT.h>
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
  
  //DEBUG_HISTS=true;
  DEBUG_HISTS=false;
  if (DEBUG_HISTS){
    dapp->Lock();
    
    fcal_match= (TH2F *)gROOT->FindObject("fcal_match");
    if (!fcal_match) fcal_match=new TH2F("fcal_match","#delta r vs p for FCAL/track matching",100,0,10,100,0,20);
    
    fcal_dt=(TH2F *)gROOT->FindObject("fcal_dt");
    if (!fcal_dt) fcal_dt=new TH2F("fcal_dt","projected time at target for FCAL vs momentum",100,0,10,100,-10,10.);

    
    dapp->Unlock();
  }

       
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

  // Get ToF points
  vector<const DTOFPoint *>tof_points;
  loop->Get(tof_points);

  // Get BCAL and FCAL clusters
  vector<const DBCALShower*>bcal_clusters;
  eventLoop->Get(bcal_clusters);
  vector<const DFCALCluster*>fcal_clusters;
  eventLoop->Get(fcal_clusters);


  // Find the time at the vertex by locking to the RF clock and match the
  // tracks with the outer detectors.  If a match is found, compute a FOM for
  // quality of the PID from the time projected from the outer detector back
  // toward the vertex.
  for (unsigned int i=0;i<_data.size();i++){
    // Match vertex time to RF bucket
    double corrected_rf_time=(_data[i]->x.Z()-target_z)/SPEED_OF_LIGHT;
    double tdiff=corrected_rf_time-_data[i]->x.T();
    if (tdiff<-1.){
      _data[i]->x.SetT(corrected_rf_time-2.);
    }
    else if (tdiff>1.){
      _data[i]->x.SetT(corrected_rf_time+2.);
    }
    else{
      _data[i]->x.SetT(corrected_rf_time);
    }

    // Here we loop over the mass hypotheses for each track
    for (unsigned int j=0;j<_data[i]->hypotheses.size();j++){
      vector<DVertex::track_info_t>tracks=_data[i]->hypotheses[j];
      for (unsigned int k=0;k<tracks.size();k++){
	bool matched_outer_detector=false;
	double tproj=0.,dmin=1000.;
	unsigned int bcal_id=0,tof_id=0,fcal_id=0;
	// Try matching the track with hits in the outer detectors
	if (pid_algorithm->MatchToBCAL(tracks[k].track,bcal_clusters,tproj,
				       bcal_id)
	    ==NOERROR){
	  matched_outer_detector=true;
	  tracks[k].tprojected=tproj;
	  tdiff=tproj-_data[i]->x.T();
	  double p=tracks[k].track->momentum().Mag();
	  double bcal_sigma=0.00255*pow(p,-2.52)+0.220;
	  double bcal_chi2=tdiff*tdiff/(bcal_sigma*bcal_sigma);
	  tracks[k].FOM=TMath::Prob(tracks[k].track->chi2_dedx+bcal_chi2,2);
	}
	else if (pid_algorithm->MatchToTOF(tracks[k].track,tof_points,tproj,
					   tof_id)
		 ==NOERROR){
	  matched_outer_detector=true;
	  tracks[k].tprojected=tproj;
	  tdiff=tproj-_data[i]->x.T();
	  double tof_sigma=0.08;
	  double tof_chi2=tdiff*tdiff/(tof_sigma*tof_sigma);
	  tracks[k].FOM=TMath::Prob(tracks[k].track->chi2_dedx+tof_chi2,2);
	}
	if (pid_algorithm->MatchToFCAL(tracks[k].track,fcal_clusters,tproj,
				       fcal_id,dmin)
	    ==NOERROR){
	  if (matched_outer_detector==false){
	    matched_outer_detector=true;
	    tracks[k].tprojected=tproj-2.218; // correction determine from fit to simulated data
	    tdiff=tproj-2.218-_data[i]->x.T();
	    double fcal_sigma=0.6; // straight-line fit to high momentum data
	    double fcal_chi2=tdiff*tdiff/(fcal_sigma*fcal_sigma);
	    tracks[k].FOM=TMath::Prob(tracks[k].track->chi2_dedx+fcal_chi2,2);
	  }
	  if (DEBUG_HISTS){
	    TH2F *fcal_dt=(TH2F*)gROOT->FindObject("fcal_dt");
	    if (fcal_dt) fcal_dt->Fill(tracks[k].track->momentum().Mag(),tproj);
	  }

	}
	if (DEBUG_HISTS){
	  TH2F *fcal_match=(TH2F *) gROOT->FindObject("fcal_match");
	  if (fcal_match) fcal_match->Fill(tracks[k].track->momentum().Mag(),
					   dmin);
	}
	  

	// If we did not succeed in matching to an outer detector, we are are
	// left with only the track info (dedx, chi2) for the FOM...
	if (matched_outer_detector==false){
	  tracks[k].FOM=tracks[k].track->FOM;
	}
      }
    }

  }



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


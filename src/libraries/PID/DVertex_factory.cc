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
#include "TRACKING/DTrackFitter.h"
using namespace jana;

bool static DVertex_hypothesis_cmp(DVertex::track_info_t a,
				   DVertex::track_info_t b)
{
  return (a.FOM>b.FOM);
}


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

  // Configuration parameters
  GROUP_NUM_SIGMAS_TIME = 100.0; // originally was 3, but changed as temporary measure
  GROUP_NUM_SIGMAS_Z    = 100.0; // originally was 3, but changed as temporary measure
  DEBUG_HISTS           = false; 
  gPARMS->SetDefaultParameter("PID:GROUP_NUM_SIGMAS_TIME", GROUP_NUM_SIGMAS_TIME, "Number of sigmas particles (tracks and or photons) can be apart in time and still be associated with same vertex");
  gPARMS->SetDefaultParameter("PID:GROUP_NUM_SIGMAS_Z"   , GROUP_NUM_SIGMAS_Z   , "Number of sigmas particles (tracks and or photons) can be apart in z and still be associated with same vertex");
  gPARMS->SetDefaultParameter("PID:DEBUG_HISTS"   , DEBUG_HISTS   , "Enable generation and filling if PID related debugging histos");

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
  
  if (DEBUG_HISTS){
    dapp->Lock();
    
    fcal_match= (TH2F *)gROOT->FindObject("fcal_match");
    if (!fcal_match) fcal_match=new TH2F("fcal_match","#delta r vs p for FCAL/track matching",100,0,10,100,0,20);
    
    fcal_dt=(TH2F *)gROOT->FindObject("fcal_dt");
    if (!fcal_dt) fcal_dt=new TH2F("fcal_dt","projected time at target for FCAL vs momentum",100,0,10,100,-10,10.);

    Nsigmas_t_photons=(TH1F *)gROOT->FindObject("Nsigmas_t_photons");
    if (!Nsigmas_t_photons) Nsigmas_t_photons=new TH1F("Nsigmas_t_photons","Number of sigmas in t (evnt)",5000, 0.0 , 500.0);

    Nsigmas_t_particles=(TH1F *)gROOT->FindObject("Nsigmas_t_particles");
    if (!Nsigmas_t_particles) Nsigmas_t_particles=new TH1F("Nsigmas_t_particles","Number of sigmas in t (AssignParticlesToGroups)",5000, 0.0 , 500.0);

    Nsigmas_z_particles=(TH1F *)gROOT->FindObject("Nsigmas_z_particles");
    if (!Nsigmas_z_particles) Nsigmas_z_particles=new TH1F("Nsigmas_z_particles","Number of sigmas in z (AssignParticlesToGroups)",5000, 0.0 , 500.0);

    
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

  // clear the list of orphan showers
  orphan_showers.clear();

  // Get list of all charged tracks
  vector<const DTrackTimeBased*> tracks;
  loop->Get(tracks);
  
  // Get ToF points
  vector<const DTOFPoint *>tof_points;
  loop->Get(tof_points);
  
  // Get BCAL and FCAL showers
  vector<const DBCALShower*>bcal_showers;
  eventLoop->Get(bcal_showers, "KLOE" );
  vector<const DFCALCluster*>fcal_showers;
  eventLoop->Get(fcal_showers);

  // group tracks according to candidate id
  vector<vector<const DTrackTimeBased*> >tracks_by_candidate;
  pid_algorithm->GroupTracks(tracks,tracks_by_candidate);

  // To minimize memory usage and time in allocation, we maintain a
  // pool of vertexInfo_t objects. Make sure the pool is large enough to hold
  // all of the particles we have in this event. 
  unsigned int Nparticles_total =tracks_by_candidate.size()+bcal_showers.size()
    +fcal_showers.size();
  for(unsigned int i=vertexInfos_pool.size(); i<Nparticles_total; i++){
    vertexInfo_t *vi = new vertexInfo_t();
    
    vi->SetLimits(tmin, tmax, zmin, zmax, Nbinst, Nbinsz);
    vertexInfos_pool.push_back(vi);
  }
  // Periodically delete some vertexInfo_t objects if the pool gets too large.
  // This prevents memory-leakage-like behavor.
  if((Nparticles_total < MAX_VERTEXINFOS) && (vertexInfos_pool.size()>MAX_VERTEXINFOS)){
    for(unsigned int i=MAX_VERTEXINFOS; i<vertexInfos_pool.size(); i++)delete vertexInfos_pool[i];
    vertexInfos_pool.resize(MAX_VERTEXINFOS);
  }
  
  // Make list of vertices
  MakeVertices(tracks_by_candidate);

  // Creat a vector to keep track of the BCAL showers that have been matched 
  // to tracks and keep track of how many are left
  unsigned int remaining_bcal_showers=bcal_showers.size();
  vector<int>bcal_matches(remaining_bcal_showers);

  // Creat a vector to keep track of the FCAL showers that have been matched 
  // to tracks and keep track of how many are left
  unsigned int remaining_fcal_showers=fcal_showers.size();
  vector<int>fcal_matches(remaining_fcal_showers);

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
      vector<DVertex::track_info_t>&tracks=_data[i]->hypotheses[j];
      for (unsigned int k=0;k<tracks.size();k++){
	bool matched_outer_detector=false;
	double tproj=0.,dmin=1000.;
	unsigned int bcal_id=0,tof_id=0,fcal_id=0;
	// Try matching the track with hits in the outer detectors
	if (pid_algorithm->MatchToBCAL(tracks[k].track->rt,
				       DTrackFitter::kTimeBased,
				       bcal_showers,tproj,
				       bcal_id)
	    ==NOERROR){
	  // Store this shower in the showers vector
	  if (bcal_matches[bcal_id]!=1){
	    _data[i]->showers.push_back(DVertex::shower_info_t(bcal_showers[bcal_id],
							       tracks[k].track));
	    // Decrement the shower counter
	    remaining_bcal_showers--;
	  }
	  bcal_matches[bcal_id]=1;
	  matched_outer_detector=true;
	  tracks[k].tprojected=tproj;
	  tdiff=tproj-_data[i]->x.T();
	  double p=tracks[k].track->momentum().Mag();
	  double bcal_sigma=0.00255*pow(p,-2.52)+0.220;
	  double bcal_chi2=tdiff*tdiff/(bcal_sigma*bcal_sigma);
	  tracks[k].FOM=TMath::Prob(tracks[k].track->chi2_dedx+bcal_chi2,2);


	}
	else if (pid_algorithm->MatchToTOF(tracks[k].track->rt,
					   DTrackFitter::kTimeBased,
					   tof_points,tproj,
					   tof_id)
		 ==NOERROR){
	  matched_outer_detector=true;
	  tracks[k].tprojected=tproj;
	  tdiff=tproj-_data[i]->x.T();
	  double tof_sigma=0.08;
	  double tof_chi2=tdiff*tdiff/(tof_sigma*tof_sigma);
	  tracks[k].FOM=TMath::Prob(tracks[k].track->chi2_dedx+tof_chi2,2);
	}
	if (pid_algorithm->MatchToFCAL(tracks[k].track->rt,
				       DTrackFitter::kTimeBased,
				       fcal_showers,tproj,
				       fcal_id,dmin)
	    ==NOERROR){
	  if (matched_outer_detector==false){
	    // Store this shower in the showers vector
	    if (fcal_matches[fcal_id]!=1){
	      _data[i]->showers.push_back(DVertex::shower_info_t(fcal_showers[fcal_id],
								 tracks[k].track));
	      // Decrement the shower counter
	      remaining_fcal_showers--;
	    }
	    fcal_matches[fcal_id]=1;

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
      // Sort hypotheses according to figure of merit
      sort(tracks.begin(),tracks.end(),DVertex_hypothesis_cmp);

    }// loop over tracks 
  }

  // Vector to hold list of vertexInfo_t objects for all unmatched showers
  vector<vertexInfo_t*> vertices;  
  if (remaining_bcal_showers>0){
    for (unsigned int i=0;i<bcal_showers.size();i++){
      if (bcal_matches[i]==0){
	vertexInfo_t *vi = vertexInfos_pool[vertices.size()];
	FillVertexInfoBCAL(vi,bcal_showers[i],i);
	vertices.push_back(vi);	
      }
    }
  }

  if (remaining_fcal_showers>0){
    for (unsigned int i=0;i<fcal_showers.size();i++){
      if (fcal_matches[i]==0){
	vertexInfo_t *vi = vertexInfos_pool[vertices.size()];
	FillVertexInfoFCAL(vi,fcal_showers[i],i);
	vertices.push_back(vi);	
      }
    }
  }

  // Now try to associate BCAL and FCAL showers with the vertex
  for (unsigned int i=0;i<_data.size();i++){
    // Correct for flight time to this particular vertex
    for (unsigned int k=0;k<vertices.size();k++){
      vertexInfo_t *vi=vertices[k];
      if (vi->bcal!=NULL && vi->is_matched_to_vertex==false){
	double tflight=(DVector3(vi->bcal->x,vi->bcal->y,vi->bcal->z)
			-_data[i]->x.Vect()).Mag()/SPEED_OF_LIGHT;
	vi->t-=tflight;
	vi->Fill(vi->t,vi->sigmat,vi->z,vi->sigmaz);
      }
      if (vi->fcal!=NULL && vi->is_matched_to_vertex==false){
	double tflight=(vi->fcal->getCentroid()-_data[i]->x.Vect()).Mag()
	  /SPEED_OF_LIGHT;
	vi->t-=tflight;
	vi->Fill(vi->t,vi->sigmat,vi->z,vi->sigmaz);
      }
    }
    // Group photons together by time
    vector< vector<vertexInfo_t *> > groups;
    AssignParticlesToGroups(vertices,groups);

    // Try to associate a group of photons with this vertex
    double t0=0.;
    double sum_invar=0, invar=0;
    for (unsigned int k=0;k<groups.size();k++){
      vector<vertexInfo_t *> &group = groups[k];
      for (unsigned int m=0;m<group.size();m++){
	invar=1./(group[m]->sigmat*group[m]->sigmat);
	sum_invar+=invar;
	t0+=group[m]->t*invar;
      }  
      t0*=1./sum_invar;
      double t_sigma_photons=1./sqrt(sum_invar);
      double t_sigma_tot=sqrt(t_sigma_photons*t_sigma_photons
			      +_data[i]->t_sigma*_data[i]->t_sigma);

      double Nsigmas_t = fabs(t0-_data[i]->x.T())/t_sigma_tot;
      if(DEBUG_HISTS)Nsigmas_t_photons->Fill(Nsigmas_t);
      if (Nsigmas_t < GROUP_NUM_SIGMAS_TIME){
	for (unsigned int m=0;m<group.size();m++){
	  _data[i]->showers.push_back(DVertex::shower_info_t(group[m]->bcal,group[m]->fcal));

	  // Flag that we have associated the showers in this group with a 
	  // vertex
	  group[m]->is_matched_to_vertex=true;
	  // Decrement the bcal and fcal shower counters
	  if (group[m]->fcal!=NULL){
	    fcal_matches[group[m]->fcal_index]=1;
	    remaining_fcal_showers--;
	  }
	  if (group[m]->bcal!=NULL){
	    bcal_matches[group[m]->bcal_index]=1;
	    remaining_bcal_showers--;
	  }
	}
      }
    }

    // Clear the groups vector so that new groups can be formed
    groups.clear();
    // Allow those photons that have not been matched to a vertex already 
    // to form a new group with a different vertex
    for (unsigned int k=0;k<vertices.size();k++){
      vertexInfo_t *vi=vertices[k];
      if (vi->is_matched_to_vertex==false) vi->is_in_group=false;
    }
  }

  // Deal with bcal and fcal showers that have not already been matched to 
  // tracks or associated with a vertex.  Since we have no vertex information 
  // for these showers, assume they came from the center of the target.  These 
  // photons are put in a separate container called "orphaned_showers".
  if (remaining_bcal_showers>0 || remaining_fcal_showers>0){
    DVector3 target_center(0.,0.,65.);
    vertices.clear();
    if (remaining_bcal_showers>0){
      for (unsigned int i=0;i<bcal_showers.size();i++){
	if (bcal_matches[i]!=1){
	  vertexInfo_t *vi = vertexInfos_pool[vertices.size()];
	  FillVertexInfoBCAL(vi,bcal_showers[i],i);
	  double tflight=(DVector3(vi->bcal->x,vi->bcal->y,vi->bcal->z)
			  -target_center).Mag()/SPEED_OF_LIGHT;
	  vi->t-=tflight;
	  vi->Fill(vi->t,vi->sigmat,vi->z,vi->sigmaz);
	  vertices.push_back(vi);	
	}
      }
      
    }
    if (remaining_fcal_showers>0){
      for (unsigned int i=0;i<fcal_showers.size();i++){
	if (fcal_matches[i]!=1){
	  vertexInfo_t *vi = vertexInfos_pool[vertices.size()];
	  FillVertexInfoFCAL(vi,fcal_showers[i],i);
	  double tflight=(vi->fcal->getCentroid()-target_center).Mag()
	  /SPEED_OF_LIGHT;
	  vi->t-=tflight;
	  vi->Fill(vi->t,vi->sigmat,vi->z,vi->sigmaz);
	  vertices.push_back(vi);		  
	}
      }
    }
    // group the photons together
    vector< vector<vertexInfo_t *> > groups;
    AssignParticlesToGroups(vertices,groups);
    for (unsigned int k=0;k<groups.size();k++){
      vector<DVertex::shower_info_t >shower_group;
      vector<vertexInfo_t *> &group = groups[k];
      for (unsigned int m=0;m<group.size();m++){
	shower_group.push_back(DVertex::shower_info_t(group[m]->bcal,group[m]->fcal));
      }
      orphan_showers.push_back(shower_group);
      shower_group.clear();
    }
  }
  
  // At this point, any orphan showers have been grouped together
  // into lists with each list containing the showers believed to 
  // have come from the same vertex. A DVertex object needs to be
  // created for each of these lists. The DVertex will be located
  // at the center of the target for lack of any better information.
  for(unsigned int i=0; i<orphan_showers.size(); i++)
  {
	DVertex *my_vertex = new DVertex;
	my_vertex->x.SetXYZT(0.0, 0.0, target_z,0.0);
	my_vertex->cov.ResizeTo(3,3);
    	my_vertex->beamline_used = true;
        
        // Set errors for photon-only vertex. The error in "t"
        // is given by the length of the target and the amount
        // of time it would take a photon to traverse the target.
        // This may be very unrealistic since the time is set to
        // "0" which for simulated data is too perfect.
        my_vertex->t_sigma = 30.0/sqrt(12.0) * (1.0/30.0); // 30cm target, 1/(30 cm/ns)
        
        // For the covariance of the position, use beam dimensions
        // at the target (rough). The radial dimension of the beam
        // is ~3mm in diameter. For both x and y, assume 0.3/2.0
        // and for z assume 30cm/sqrt(12). We could also divide the
        // radial error by sqrt(12), but this is already too accurate
        // compared to vertices obtained through tracking!
        double sigma_r2 = pow(0.3/2.0,2.0);
        double sigma_z2 = pow(30.0/sqrt(12.0),2.0);
        DMatrix &cov = my_vertex->cov;
        cov[0][0] = sigma_r2;  cov[0][1] =   0.0;     cov[0][2] =   0.0;
        cov[1][0] =   0.0;     cov[1][1] = sigma_r2;  cov[1][2] =   0.0;
        cov[2][0] =   0.0;     cov[2][1] =   0.0;     cov[2][2] = sigma_z2;
        
        my_vertex->showers = orphan_showers[i]; // copy entire list of showers for this vertex
        
        _data.push_back(my_vertex);
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
  // Vector to hold list of vertexInfo_t objects for all charged tracks
  vector<vertexInfo_t*> vertices;

  // Assign and fill vertexInfo_t objects for each charged track
  for(unsigned int i=0; i<tracks_by_candidate.size(); i++){
    vertexInfo_t *vi = vertexInfos_pool[vertices.size()];
    // Use the fit result with the highest figure of merit
    FillVertexInfoChargedTrack(vi, &tracks_by_candidate[i]);
    vertices.push_back(vi);
  }
  // Group tracks together
  vector< vector<vertexInfo_t *> > groups;
  AssignParticlesToGroups(vertices,groups);

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
    my_vertex->t_sigma=0.; // <------ this needs to be fixed

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
// FillVertexInfoBCAL
//------------------
#define EPS 1.e-8
void DVertex_factory::FillVertexInfoBCAL(DVertex_factory::vertexInfo_t *vi,
					 const DBCALShower *bcal,
					 unsigned int index){
   vi->Reset();
   vi->bcal=bcal;
   vi->bcal_index=index;
   vi->fcal_index=0;
   vi->fcal=NULL;
   vi->z=0;  //will be filled in later
   vi->sigmaz=30.0/sqrt(12); // in cm. Use length of target for z-resolution of photons
   vi->t=bcal->t; // will be propagated to vertex later
   vi->sigmat = bcal->tErr;
   // somtimes KLOE algorithm returns a zero RMS..
   if (vi->sigmat<EPS) vi->sigmat=0.5; // guess!
   
}

//------------------
// FillVertexInfoFCAL
//------------------
void DVertex_factory::FillVertexInfoFCAL(DVertex_factory::vertexInfo_t *vi,
					 const DFCALCluster *fcal,
					 unsigned int index){
   vi->Reset();
   vi->fcal=fcal;
   vi->fcal_index=index;
   vi->bcal_index=0;
   vi->bcal=NULL;
   vi->z=0;  //will be filled in later
   vi->sigmaz=30.0/sqrt(12); // in cm. Use length of target for z-resolution of photons
   vi->t=fcal->getTime(); // will be propagated to vertex later
   vi->sigmat=0.5;		   
}

//------------------
// FillVertexInfoChargedTrack
//------------------
void DVertex_factory::FillVertexInfoChargedTrack(DVertex_factory::vertexInfo_t *vi, vector<const DTrackTimeBased *>*hypotheses)
{
	vi->Reset();
	vi->hypotheses = hypotheses;

	vi->t = (*hypotheses)[0]->t0();
	vi->sigmat = (*hypotheses)[0]->t0_err();
	vi->z = (*hypotheses)[0]->z();
	vi->sigmaz = 0.8/sin((*hypotheses)[0]->momentum().Theta()); // in cm.  For now, use 3mm wide angle track resolution scaled by sin(theta)

	vi->Fill(vi->t, vi->sigmat, vi->z, vi->sigmaz);
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

// Group "particles" (either tracks or photon candidates) together by time 
// and z position
void DVertex_factory::AssignParticlesToGroups(vector<vertexInfo_t*> &vertices,
					      vector< vector<vertexInfo_t *> > &groups
					      ){
  
  // Each particle has a histogram of t vs.z
  // values filled using approriate uncertainties (no covariance). We
  // can now use this list to identify resonances in the t/z plane which
  // indicate a vertex location. Particles within 3sigma in both t and
  // z of the resonance will be grouped together as belonging to the 
  // same vertex. We loop until all particles have been assigned
  // to a group, even if that means assigning particles to their own
  // "group of one".
  
  // Loop until all particles have been assigned to a group.
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
      vertexInfo_t *vi = vertices[i];
      if(vi->is_in_group)continue;
      if (vi->is_matched_to_vertex) continue;
			
      double delta_t = fabs(maxloc.X() - vi->t);
      double Nsigmas_t = delta_t/vi->sigmat;
      if(DEBUG_HISTS)Nsigmas_t_particles->Fill(Nsigmas_t);
      if(Nsigmas_t > GROUP_NUM_SIGMAS_TIME)continue;
      
      double delta_z = fabs(maxloc.Y() - vi->z);
      double Nsigmas_z = delta_z/vi->sigmaz;
      if(DEBUG_HISTS)Nsigmas_z_particles->Fill(Nsigmas_z);
      if(Nsigmas_z > GROUP_NUM_SIGMAS_Z)continue;
      
      // Assign this particle to the group
      vi->is_in_group=true;
      new_group.push_back(vi);
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
      vertexInfo_t *vi_with_max_t = NULL;
      double delta_t_max=0.0;
      for(unsigned int i=0; i<vertices.size(); i++){
	vertexInfo_t *vi = vertices[i];
	if(vi->is_in_group)continue;
			
	double delta_t = fabs(maxloc.X() - vi->t);
	if(delta_t>delta_t_max || vi_with_max_t==NULL){
	  delta_t_max = delta_t;
	  vi_with_max_t = vi;
	}
      }
			
      if(vi_with_max_t==NULL){
	_DBG_<<"vi_with_max_t==NULL. This should never happen! Complain to davidl@jlab.org"<<endl;
	_DBG_<<"event number: "<<eventnumber<<endl;
	break;
			}
			
      new_group.push_back(vi_with_max_t);
    }
		
    // Set the is_in_group flags for all of the members of the new group
    for(unsigned int i=0; i<new_group.size(); i++)new_group[i]->is_in_group = true;
    
    // Add the new group to the list of groups
    groups.push_back(new_group);  
  }
}

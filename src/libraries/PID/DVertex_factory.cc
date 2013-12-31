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
#include <TRACKING/DReferenceTrajectory.h>
using namespace jana;


inline bool DVertexSortTracks(const DChargedTrack *a,const DChargedTrack *b){
  const DChargedTrackHypothesis *hyp_a=a->dChargedTrackHypotheses[0];
  const DChargedTrackHypothesis *hyp_b=b->dChargedTrackHypotheses[0];
  
  return (hyp_a->momentum().Mag()>hyp_b->momentum().Mag());
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
  loop->GetSingle(dAnalysisUtilities);

  // Get Target parameters from XML
  DApplication *locApplication = dynamic_cast<DApplication*> (loop->GetJApplication());
  DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(runnumber):NULL;
  dTargetCenter_Z = 0.0;
  if(locGeometry)
    locGeometry->GetTargetZ(dTargetCenter_Z);

  // Configuration parameters
  GROUP_NUM_SIGMAS_TIME = 100.0; // originally was 3, but changed as temporary measure
  GROUP_NUM_SIGMAS_Z    = 100.0; // originally was 3, but changed as temporary measure
  DEBUG_HISTS           = false; 
  gPARMS->SetDefaultParameter("PID:GROUP_NUM_SIGMAS_TIME", GROUP_NUM_SIGMAS_TIME, "Number of sigmas particles (tracks and or photons) can be apart in time and still be associated with same vertex");
  gPARMS->SetDefaultParameter("PID:GROUP_NUM_SIGMAS_Z"   , GROUP_NUM_SIGMAS_Z   , "Number of sigmas particles (tracks and or photons) can be apart in z and still be associated with same vertex");
  gPARMS->SetDefaultParameter("PID:DEBUG_HISTS"   , DEBUG_HISTS   , "Enable generation and filling if PID related debugging histos");

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
    locApplication->Lock();

    Nsigmas_t_particles=(TH1F *)gROOT->FindObject("Nsigmas_t_particles");
    if (!Nsigmas_t_particles) Nsigmas_t_particles=new TH1F("Nsigmas_t_particles","Number of sigmas in t (AssignParticlesToGroups)",5000, 0.0 , 500.0);

    Nsigmas_z_particles=(TH1F *)gROOT->FindObject("Nsigmas_z_particles");
    if (!Nsigmas_z_particles) Nsigmas_z_particles=new TH1F("Nsigmas_z_particles","Number of sigmas in z (AssignParticlesToGroups)",5000, 0.0 , 500.0);
    
    locApplication->Unlock();
  }

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DVertex_factory::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DChargedTrack *> locChargedTracks;
	loop->Get(locChargedTracks);

	// Sort the tracks such that the one with the highest momentum is first
	sort(locChargedTracks.begin(),locChargedTracks.end(),DVertexSortTracks);
	
#if 1
	vector<vertex_t>groups;  // vector to keep track of combinations of tracks
	unsigned int numTracks=locChargedTracks.size();
	if (numTracks>0){
	  unsigned int last_index=numTracks-1;
	  for (unsigned int j=0;j<last_index;j++){
	    const DChargedTrackHypothesis *track1=locChargedTracks[j]->dChargedTrackHypotheses[0];
	    DKinematicData track1_kd=*track1;
	    
	    // Try to match this track with the other tracks
	    unsigned int trackbits=0x1u<<j;   // store track indices in a bit pattern
	    DVector3 avg_pos;
	    // double weight=0;
	    DVector3 weight;
	    double t0=0.;
	    double t0_err=0.,one_over_t0_var=0.;
	    double t0_weight=0.;
	    // boolean to make sure we use the current track only once when computing averages
	    double use_track_in_average=true; 
	    // loop over the remaining tracks
	    for (unsigned int i=j+1;i<numTracks;i++){
	      const DChargedTrackHypothesis *track2=locChargedTracks[i]->dChargedTrackHypotheses[0];
	      DKinematicData track2_kd=*track2;
	      
	      // Find the intersection between the two tracks
	      DVector3 intersection;
	      double doca=1000.,var_doca=1000.;
			doca = dAnalysisUtilities->Calc_DOCAVertex(&track1_kd, &track2_kd, intersection);

	      // Add the track to the group if the combination of the current two tracks meets a minimum 
	      // probability criterion.
	      double prob=TMath::Prob(doca*doca/var_doca,1);
	      if (prob>0.01){
		trackbits|=0x1u<<i;
		
		// Accumulate values for computing vertex time and position
		if (use_track_in_average){
		  t0_err=track1_kd.t0_err();
		  one_over_t0_var=1./(t0_err*t0_err);
		  t0+=one_over_t0_var*track1_kd.t0();
		  t0_weight+=one_over_t0_var;

		  DVector3 pos1=track1_kd.position();
		  DMatrixDSym cov1=track1_kd.errorMatrix();
		  DVector3 dpos(pos1.x()/cov1(3,3),pos1.y()/cov1(4,4),pos1.z()/cov1(5,5));
		  DVector3 dweight(1./cov1(3,3),1/cov1(4,4),1/cov1(5,5));
		  
		  avg_pos+=dpos;
		  weight+=dweight;
		  
		  use_track_in_average=false;
		}
		
		t0_err=track2_kd.t0_err();
		one_over_t0_var=1./(t0_err*t0_err);
		t0+=one_over_t0_var*track2_kd.t0();
		t0_weight+=one_over_t0_var;

		DVector3 pos2=track2_kd.position();
		DMatrixDSym cov2=track2_kd.errorMatrix();
		DVector3 dpos(pos2.x()/cov2(3,3),pos2.y()/cov2(4,4),pos2.z()/cov2(5,5));
		DVector3 dweight(1./cov2(3,3),1/cov2(4,4),1/cov2(5,5));

		avg_pos+=dpos;
		weight+=dweight;
	      }
	    }
	    if (t0_weight==0.){
	      DVector3 pos1=track1_kd.position();
	      DMatrixDSym cov1=track1_kd.errorMatrix();
	      avg_pos.SetXYZ(pos1.x()/cov1(3,3),pos1.y()/cov1(4,4),pos1.z()/cov1(5,5));
	      weight.SetXYZ(1./cov1(3,3),1/cov1(4,4),1/cov1(5,5));
	      t0=track1_kd.t0();
	      double t0_err=track1_kd.t0_err();
	      t0_weight=1./(t0_err*t0_err);
	    }
	    groups.push_back(vertex_t(trackbits,avg_pos,weight,t0,t0_weight));
	  }
	  // Cut out groups that contain tracks that are members of other groups
	  vector<unsigned int>groups_to_cut;
	  if (groups.size()>1){
	    for (unsigned int i=0;i<groups.size()-1;i++){
	      for (unsigned int j=i+1;j<groups.size();j++){ 
		if (groups[i].trackbits&groups[j].trackbits){
		  groups[i].trackbits|=groups[j].trackbits;
		  groups_to_cut.push_back(j);
		}
	      }
	    }	 
	  }
	  vector<vertex_t> groups_to_keep;
	  for (unsigned int i=0;i<groups.size();i++){
	    bool use_group=true;
	    for (unsigned int j=0;j<groups_to_cut.size();j++){
	      if (i==groups_to_cut[j]) use_group=false;
	    }
	    if (use_group) groups_to_keep.push_back(groups[i]);
	  }
	  unsigned int all_used_tracks=0; // Keep track of all grouped tracks
	  for (unsigned int i=0;i<groups_to_keep.size();i++){
	    all_used_tracks|=groups_to_keep[i].trackbits;

	    DVertex *locVertex = new DVertex;

	    double weight_x=groups_to_keep[i].weight.x();
	    double weight_y=groups_to_keep[i].weight.y();
	    double weight_z=groups_to_keep[i].weight.z();
	    
	    // Vertex position
	    locVertex->dSpacetimeVertex.SetX(groups_to_keep[i].pos.x()/weight_x);
	    locVertex->dSpacetimeVertex.SetY(groups_to_keep[i].pos.y()/weight_y);
	    locVertex->dSpacetimeVertex.SetZ(groups_to_keep[i].pos.z()/weight_z);

	    // Covariance matrix for vertex position
	    locVertex->locCovarianceMatrix.ResizeTo(3,3);
	    locVertex->locCovarianceMatrix(0,0)=1./weight_x;
	    locVertex->locCovarianceMatrix(1,1)=1./weight_y;
	    locVertex->locCovarianceMatrix(2,2)=1./weight_z;
	    
	    // Vertex time
	    // T=0 corresponds to the center of the target, so correct the average t0 accordingly
	    double t0=groups_to_keep[i].t0/groups_to_keep[i].t0_weight
	      +(dTargetCenter_Z-locVertex->dSpacetimeVertex.Z())/SPEED_OF_LIGHT;
	    locVertex->dSpacetimeVertex.SetT(t0);
	    double var_t0=1./groups_to_keep[i].t0_weight+1./(weight_z*SPEED_OF_LIGHT*SPEED_OF_LIGHT);
	    locVertex->dTimeUncertainty=sqrt(var_t0);

	    for (unsigned int k=0;k<numTracks;k++){
	      if ((0x1u<<k)&groups_to_keep[i].trackbits){
		locVertex->dChargedTracks.push_back(locChargedTracks[k]);
	      }
	    }

	    _data.push_back(locVertex);
	  }
	  // Add single tracks unmatched to other tracks to the vertex list
	  if ((all_used_tracks & (0x1u<<last_index)) == 0){
	    const DChargedTrackHypothesis *locTrack=locChargedTracks[last_index]->dChargedTrackHypotheses[0];

	    DVertex *locVertex = new DVertex;
	    locVertex->locCovarianceMatrix.ResizeTo(3,3);
	    locVertex->dSpacetimeVertex.SetVect(locTrack->position());
	    locVertex->dSpacetimeVertex.SetT(locTrack->t0());
	    locVertex->dTimeUncertainty=locTrack->t0_err();
	    locVertex->dChargedTracks.push_back(locChargedTracks[0]);
	    
	    _data.push_back(locVertex);
	    

	  }

	}
	  
#else
	// To minimize memory usage and time in allocation, we maintain a
	// pool of vertexInfo_t objects. Make sure the pool is large enough to hold
	// all of the particles we have in this event. 
	unsigned int loc_i;
	for(loc_i = dVertexInfoPool.size(); loc_i < locChargedTracks.size(); loc_i++){
		vertexInfo_t *locVertexInfo = new vertexInfo_t();
		locVertexInfo->SetLimits(tmin, tmax, zmin, zmax, Nbinst, Nbinsz);
		dVertexInfoPool.push_back(locVertexInfo);
	}

	// Periodically delete some vertexInfo_t objects if the pool gets too large.
	// This prevents memory-leakage-like behavor.
	if((locChargedTracks.size() < MAX_VERTEXINFOS) && (dVertexInfoPool.size() > MAX_VERTEXINFOS)){
		for(loc_i = MAX_VERTEXINFOS; loc_i < dVertexInfoPool.size(); loc_i++) delete dVertexInfoPool[loc_i];
		dVertexInfoPool.resize(MAX_VERTEXINFOS);
	}

	// Make list of vertices
	MakeVertices(locChargedTracks);
#endif


	// If no vertices but neutral showers present, create vertex at center of target
	if(_data.size() == 0){
		vector<const DNeutralShower *> locNeutralShowers;
		loop->Get(locNeutralShowers);
		if(locNeutralShowers.size() > 0){ //use center of target
			DVertex *locVertex = new DVertex;
			locVertex->dTimeUncertainty = 0.0;
			DLorentzVector locSpacetimeVertex(0.0, 0.0, dTargetCenter_Z, 0.0);
			locVertex->dSpacetimeVertex = locSpacetimeVertex;
			_data.push_back(locVertex);
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
jerror_t DVertex_factory::MakeVertices(vector<const DChargedTrack*> &locChargedTracks){
	unsigned int loc_i, loc_j;
	// Vector to hold list of vertexInfo_t objects for all charged tracks
	vector<vertexInfo_t*> locVertexInfos;

	// Assign and fill vertexInfo_t objects for each charged track
	for(loc_i = 0; loc_i < locChargedTracks.size(); loc_i++){
		vertexInfo_t *locVertexInfo = dVertexInfoPool[locVertexInfos.size()];
		// Use the fit result with the highest figure of merit
		FillVertexInfoChargedTrack(locVertexInfo, locChargedTracks[loc_i]);
		locVertexInfos.push_back(locVertexInfo);
	}
	// Group tracks together
	vector< vector<vertexInfo_t *> > locVertexInfoGroups;
	AssignParticlesToGroups(locVertexInfos, locVertexInfoGroups);

	// OK, we've now grouped the particles together into groups. Create a new
	// DVertex for each group	
	for (loc_i = 0; loc_i < locVertexInfoGroups.size(); loc_i++){
		vector<vertexInfo_t *> &locVertexInfoGroup = locVertexInfoGroups[loc_i];

		DVertex *locVertex = new DVertex;
		locVertex->locCovarianceMatrix.ResizeTo(3,3);

		// Simply average POCA to beamline for all tracks	
		DVector3 temp;
		double t0 = 0.;
		double sum_invar = 0.0, invar = 0.0;
		for(loc_j = 0; loc_j < locVertexInfoGroup.size(); loc_j++){
			invar = 1.0/(locVertexInfoGroup[loc_j]->sigmat*locVertexInfoGroup[loc_j]->sigmat);
			sum_invar += invar;
			t0 += locVertexInfoGroup[loc_j]->t*invar;
			temp += locVertexInfoGroup[loc_j]->dChargedTrack->dChargedTrackHypotheses[0]->position();
		}
		t0 *= 1.0/sum_invar;
		temp *= 1.0/double(locVertexInfoGroup.size());
		locVertex->dSpacetimeVertex.SetVect(temp);
		locVertex->dSpacetimeVertex.SetT(t0);
		locVertex->dTimeUncertainty = 0.; // <------ this needs to be fixed

		// Add list of tracks used to create this vertex
		for(loc_j = 0; loc_j < locVertexInfoGroup.size(); loc_j++)
			locVertex->dChargedTracks.push_back(locVertexInfoGroup[loc_j]->dChargedTrack);
		_data.push_back(locVertex);
	}
	
	return NOERROR;
}

//------------------
// FillVertexInfoChargedTrack
//------------------
void DVertex_factory::FillVertexInfoChargedTrack(DVertex_factory::vertexInfo_t *locVertexInfo, const DChargedTrack *locChargedTrack)
{
	locVertexInfo->Reset();
	locVertexInfo->dChargedTrack = locChargedTrack;
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->dChargedTrackHypotheses[0];

	locVertexInfo->t = locChargedTrackHypothesis->t0();
	locVertexInfo->sigmat = locChargedTrackHypothesis->t0_err();
	locVertexInfo->z = locChargedTrackHypothesis->z();
	locVertexInfo->sigmaz = 0.8/sin(locChargedTrackHypothesis->momentum().Theta()); // in cm.  For now, use 3mm wide angle track resolution scaled by sin(theta)

	locVertexInfo->Fill(locVertexInfo->t, locVertexInfo->sigmat, locVertexInfo->z, locVertexInfo->sigmaz);
}

// Group particles together by time and z position
void DVertex_factory::AssignParticlesToGroups(vector<vertexInfo_t*> &locVertexInfos, vector< vector<vertexInfo_t *> > &locVertexInfoGroups){
	unsigned int loc_i;

	// Each particle has a histogram of t vs.z
	// values filled using approriate uncertainties (no covariance). We
	// can now use this list to identify resonances in the t/z plane which
	// indicate a vertex location. Particles within 3sigma in both t and
	// z of the resonance will be grouped together as belonging to the 
	// same vertex. We loop until all particles have been assigned
	// to a group, even if that means assigning particles to their own
	// "group of one".
	
	// Loop until all particles have been assigned to a group.
	while(!AllInGroups(locVertexInfos)){
		// Make a list of all particles that have not been assigned to a group
		vector<const DHoughFind*> locUnassignedHoughs;
		for(loc_i = 0; loc_i < locVertexInfos.size(); loc_i++){
			if(!locVertexInfos[loc_i]->is_in_group) locUnassignedHoughs.push_back(locVertexInfos[loc_i]);
		}
		// Find the maximum t,z coordinate by adding all unassigned
		// particle's histos together
		DVector2 maxloc = DHoughFind::GetMaxBinLocation(locUnassignedHoughs);
		
		if(debug_level>0)_DBG_<<"Location of maximum: t="<<maxloc.X()<<"  z="<<maxloc.Y()<<endl;		

		// Loop over all unassigned particles, assigning any within
		// 3 sigma in both t and z to the new group. We loop over
		// the locVertexInfos vector just because it saves a dynamic_cast
		// if we were to use the unassigned vector.
		vector<vertexInfo_t *> locVertexInfoGroup;
		for(loc_i = 0; loc_i < locVertexInfos.size(); loc_i++){
			vertexInfo_t *locVertexInfo = locVertexInfos[loc_i];
			if (locVertexInfo->is_in_group) continue;
			if (locVertexInfo->is_matched_to_vertex) continue;
			
			double delta_t = fabs(maxloc.X() - locVertexInfo->t);
			double Nsigmas_t = delta_t/locVertexInfo->sigmat;
			if(DEBUG_HISTS) Nsigmas_t_particles->Fill(Nsigmas_t);
			if(Nsigmas_t > GROUP_NUM_SIGMAS_TIME) continue;
			
			double delta_z = fabs(maxloc.Y() - locVertexInfo->z);
			double Nsigmas_z = delta_z/locVertexInfo->sigmaz;
			if(DEBUG_HISTS) Nsigmas_z_particles->Fill(Nsigmas_z);
			if(Nsigmas_z > GROUP_NUM_SIGMAS_Z) continue;
			
			// Assign this particle to the group
			locVertexInfo->is_in_group = true;
			locVertexInfoGroup.push_back(locVertexInfo);
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
		if(locVertexInfoGroup.size() == 0){
			vertexInfo_t *vi_with_max_t = NULL;
			double delta_t_max = 0.0;
			for(loc_i = 0; loc_i < locVertexInfos.size(); loc_i++){
				vertexInfo_t *locVertexInfo = locVertexInfos[loc_i];
				if(locVertexInfo->is_in_group) continue;
				double delta_t = fabs(maxloc.X() - locVertexInfo->t);
				if(delta_t>delta_t_max || vi_with_max_t==NULL){
					delta_t_max = delta_t;
					vi_with_max_t = locVertexInfo;
				}
			}
			if(vi_with_max_t == NULL){
				_DBG_<<"vi_with_max_t==NULL. This should never happen! Complain to davidl@jlab.org"<<endl;
				break;
			}
			locVertexInfoGroup.push_back(vi_with_max_t);
		}
		
		// Set the is_in_group flags for all of the members of the new group
		for(loc_i = 0; loc_i < locVertexInfoGroup.size(); loc_i++) locVertexInfoGroup[loc_i]->is_in_group = true;
		
		// Add the new group to the list of groups
		locVertexInfoGroups.push_back(locVertexInfoGroup);
	}
}

//------------------
// AllInGroups
//------------------
bool DVertex_factory::AllInGroups(vector<vertexInfo_t*> &locVertexInfos)
{
	for(unsigned int loc_i = 0; loc_i < locVertexInfos.size(); loc_i++){
		if(!locVertexInfos[loc_i]->is_in_group) return false;
	}
	return true;
}


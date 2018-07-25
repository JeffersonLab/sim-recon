// $Id$
//
//    File: DTrackTimeBased_factory_THROWN.cc
// Created: Wed Nov 18 06:25:19 EST 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#include <cmath>
using namespace std;

#include <DANA/DApplication.h>

#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <TRACKING/DTrackWireBased.h>

#include "DTrackTimeBased_factory_THROWN.h"
#include "DMCThrown.h"
#include "DReferenceTrajectory.h"
#include "DRandom.h"
#include "DMatrix.h"
#include "DTrackHitSelector.h"


//------------------
// DTrackTimeBased_factory_THROWN
//------------------
DTrackTimeBased_factory_THROWN::DTrackTimeBased_factory_THROWN()
{
	fitter = NULL;
	hitselector=NULL;
	
	RootGeom = NULL;
	geom = NULL;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory_THROWN::brun(jana::JEventLoop *loop, int32_t runnumber)
{
	// Get pointer to DTrackFitter object that actually fits a track
	vector<const DTrackFitter *> fitters;
	loop->Get(fitters);
	if(fitters.size()<1){
		_DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	
	// Drop the const qualifier from the DTrackFitter pointer (I'm surely going to hell for this!)
	fitter = const_cast<DTrackFitter*>(fitters[0]);

	// Warn user if something happened that caused us NOT to get a fitter object pointer
	if(!fitter){
		_DBG_<<"ERROR: Unable to get a DTrackFitter object! Chisq for DTrackTimeBased:THROWN will NOT be calculated!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	// Get pointer to DTrackHitSelector object
	vector<const DTrackHitSelector *> hitselectors;
	loop->Get(hitselectors);
	if(hitselectors.size()<1){
		_DBG_<<"ERROR: Unable to get a DTrackHitSelector object! NO DTrackTimeBased:THROWN objects will be created!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	hitselector = hitselectors[0];

	// Set DGeometry pointers so it can be used by the DReferenceTrajectory class
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	geom = dapp->GetDGeometry(runnumber);
	
	// Set magnetic field pointer
	bfield = dapp->GetBfield();

	// load DC geometry
    geom->GetCDCWires(cdcwires);
    //   geom->GetCDCRmid(cdc_rmid); // THIS ISN'T IMPLEMENTED!!
    // extract the "mean" radius of each ring from the wire data
    for(int ring=0; ring<cdcwires.size(); ring++)
  		cdc_rmid.push_back( cdcwires[ring][0]->origin.Perp() );
  	double loc_endplate_dz, loc_endplate_rmin, loc_endplate_rmax;
  	geom->GetCDCEndplate(cdc_endplate_z, loc_endplate_dz, loc_endplate_rmin, loc_endplate_rmax);

   	// Get z positions of fdc wire planes
   	geom->GetFDCZ(fdc_z_wires);
   	// for now, assume the z extent of a package is the difference between the positions
   	// of the two wire planes.  save half of this distance
   	fdc_package_size = (fdc_z_wires[1]-fdc_z_wires[0]) / 2.;
   	geom->GetFDCRmin(fdc_rmin_packages);
   	geom->GetFDCRmax(fdc_rmax);

	// Get the particle ID algorithms
	loop->GetSingle(dParticleID);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackTimeBased_factory_THROWN::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCPseudo*> fdcpseudos;
	vector<const DTrackWireBased*> wbtracks;
	loop->Get(mcthrowns);
	loop->Get(cdctrackhits);
	loop->Get(fdcpseudos);
	loop->Get(wbtracks, "THROWN");

	for(unsigned int i=0; i< mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];

		if(fabs(thrown->charge())<1)continue;

		// First, copy over the DKinematicData part
		DTrackTimeBased *track = new DTrackTimeBased();
	  *static_cast<DKinematicData*>(track) = *static_cast<const DKinematicData*>(thrown);

		// Set PID		
      track->setPID(thrown->PID());

		// Add DMCThrown as associated object
		track->AddAssociatedObject(thrown);

		// We need to swim a reference trajectory here. To avoid the overhead
		// of allocating/deallocating them every event, we keep a pool and
		// re-use them. If the pool is not big enough, then add one to the
		// pool.
      unsigned int locNumInitialReferenceTrajectories = rt_pool.size();
		if(rt_pool.size()<=_data.size()){
			// This is a little ugly, but only gets called a few times throughout the life of the process
			// Note: these never get deleted, even at the end of process.			
			rt_pool.push_back(new DReferenceTrajectory(bfield));
		}

		DReferenceTrajectory *rt = rt_pool[_data.size()];
      if(locNumInitialReferenceTrajectories == rt_pool.size()) //didn't create a new one
        rt->Reset();
      rt->q = track->charge();
      //		track->rt = rt;
		DVector3 pos = track->position();
		DVector3 mom = track->momentum();
		rt->SetMass(thrown->mass());
		rt->SetDGeometry(geom);
		rt->SetDRootGeom(RootGeom);
		rt->Swim(pos, mom, track->charge());

		// Find hits that should be on this track and add them as associated objects
		vector<const DCDCTrackHit*> cdchits;
		vector<const DFDCPseudo*> fdchits;
		if(hitselector && cdctrackhits.size()>0)hitselector->GetCDCHits(DTrackHitSelector::kHelical, rt, cdctrackhits, cdchits);
		if(hitselector && fdcpseudos.size()>0)hitselector->GetFDCHits(DTrackHitSelector::kHelical, rt, fdcpseudos, fdchits);
		for(unsigned int j=0; j<cdchits.size(); ++j)track->AddAssociatedObject(cdchits[j]);
		for(unsigned int j=0; j<fdchits.size(); ++j)track->AddAssociatedObject(fdchits[j]);

		track->measured_cdc_hits_on_track = cdchits.size();
 	    track->measured_fdc_hits_on_track = fdchits.size();

		// Since we have swum a DReferenceTrajectory, use the same algorithm that's used 
		// for the internal reference trajectories in DTrackFitterKalmanSIMD
		// It's a little nasty, but should work for now

		set<const DCDCWire *> expected_hit_straws;
		set<int> expected_hit_fdc_planes;

		for(int i=0; i<rt->Nswim_steps; i++) {
			double z = rt->swim_steps[i].origin.Z();
			double r = rt->swim_steps[i].origin.Perp();

			// CHECK HITS IN CDC
			if(z <= cdc_endplate_z) {
			// figure out the radial position of the point to see which ring it's in
			double r = rt->swim_steps[i].origin.Perp();
			int ring=0;
			for(; ring<cdc_rmid.size(); ring++) {
				if( (r<cdc_rmid[ring]-0.78) || (fabs(r-cdc_rmid[ring])<0.78) )
					break;
			}
			if(ring == cdc_rmid.size()) ring--;
			//_DBG_ << "ring = " << ring << endl;
			//_DBG_ << "ring = " << ring << "  stereo = " << cdcwires[ring][0]->stereo << endl;
			int best_straw=0;
			double best_dist_diff=fabs((rt->swim_steps[i].origin 
				- cdcwires[ring][0]->origin).Mag());		
	    	// match based on straw center
	    	for(int straw=1; straw<cdcwires[ring].size(); straw++) {
	    		DVector3 wire_position = cdcwires[ring][straw]->origin;  // start with the nominal wire center
	    		// now take into account the z dependence due to the stereo angle
	    		double dz = rt->swim_steps[i].origin.Z() - cdcwires[ring][straw]->origin.Z();
	    		double ds = dz*tan(cdcwires[ring][straw]->stereo);
	    		wire_position += DVector3(-ds*sin(cdcwires[ring][straw]->origin.Phi()), ds*cos(cdcwires[ring][straw]->origin.Phi()), dz);
	    		double diff = fabs((rt->swim_steps[i].origin
					- wire_position).Mag());
				if( diff < best_dist_diff )
					best_straw = straw;
	    	}
	    
	    	expected_hit_straws.insert(cdcwires[ring][best_straw]);
			}
			
			// CHECK HITS IN FDC
			if( z>=fdc_z_wires[0] && z<=fdc_z_wires[fdc_z_wires.size()-1]) {
			// check to make sure that the track goes through the sensitive region of the FDC
			// assume one hit per plane

			// see if we're in the "sensitive area" of a package
			for(int plane=0; plane<fdc_z_wires.size(); plane++) {
				int package = plane/6;
				if(fabs(z-fdc_z_wires[plane]) < fdc_package_size) {
					if( r<fdc_rmax && r>fdc_rmin_packages[package]) {
						expected_hit_fdc_planes.insert(plane);
					}
					break; // found the right plane
				}
 			}
			}
		}
	
		track->potential_cdc_hits_on_track = expected_hit_straws.size();
		track->potential_fdc_hits_on_track = expected_hit_fdc_planes.size();

   		//_DBG_ << " CDC hits/potential hits " << cdchits.size() << "/" << track->potential_cdc_hits_on_track 
        //	 << "  FDC hits/potential hits " << fdchits.size() << "/" << track->potential_fdc_hits_on_track  << endl;

		// We want to get chisq and Ndof values for this track using the hits from above.
		// We do this using the DTrackFitter object. This more or less guarantees that the
		// chisq calculation is done in the same way as it is for track fitting. Note
		// that no fitting is actually done here so this should be reasonably fast
		if(fitter){
			fitter->Reset();
			fitter->AddHits(cdchits);
			fitter->AddHits(fdchits);
			double chisq;
			int Ndof;
			vector<DTrackFitter::pull_t> pulls;
			fitter->ChiSq(DTrackFitter::kTimeBased, rt, &chisq, &Ndof, &pulls);
			track->chisq = chisq;
			track->Ndof = Ndof;
			track->pulls = pulls;
		}else{
			track->chisq = 0.0;
			track->Ndof = 0;
		}
		
		// For this to work properly with DChargedTrack, we need to put something
		// in for the candidateid and the FOM (figure of merit) used to decide if
		// this is the right mass hypothesis for the candidate. Since there is no
		// candidate and we *know* it's the right hypothesis, we set the candidateid
		// to the thrown object's oid and set the FOM to 1.
		track->candidateid = thrown->id;
		track->trackid = thrown->id;
		track->FOM = 1.0;

		// Add wire-based track as associated object. Even though they should
		// be in the same order, we verify it is the correct one by checking
		// that the candidateid is the same as ours (i.e. the same thrown track)
		for(unsigned int j=0; j<wbtracks.size(); j++){
			if(wbtracks[j]->candidateid == track->candidateid){
				track->AddAssociatedObject(wbtracks[j]);
				break;
			}
		}

      // Set MC Hit-matching information
      track->dMCThrownMatchMyID = thrown->myid;
      track->dNumHitsMatchedToThrown = track->Ndof + 5;

		_data.push_back(track);
	}

	// Set CDC ring & FDC plane hit patterns
	for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i)
	{
		vector<const DCDCTrackHit*> locCDCTrackHits;
		_data[loc_i]->Get(locCDCTrackHits);

		vector<const DFDCPseudo*> locFDCPseudos;
		_data[loc_i]->Get(locFDCPseudos);

		_data[loc_i]->dCDCRings = dParticleID->Get_CDCRingBitPattern(locCDCTrackHits);
		_data[loc_i]->dFDCPlanes = dParticleID->Get_FDCPlaneBitPattern(locFDCPseudos);
		
	}

	return NOERROR;
}


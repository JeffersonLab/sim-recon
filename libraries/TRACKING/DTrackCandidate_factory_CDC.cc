// $Id$
//
//    File: DTrackCandidate_factory_CDC.cc
// Created: Thu Sep  6 14:47:48 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#include <algorithm>
using namespace std;

#include <DVector2.h>
#include <CDC/DCDCTrackHit.h>
#include "DQuickFit.h"
#include "DRiemannFit.h"

#include "DTrackCandidate_factory_CDC.h"

bool CDCSortByRdecreasing(DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit1, DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit2) {
	// use the ring number to sort by R(decreasing) and then straw(increasing)
	if(hit1->hit->wire->ring == hit2->hit->wire->ring){
		return hit1->hit->wire->straw < hit2->hit->wire->straw;
	}
	return hit1->hit->wire->ring > hit2->hit->wire->ring;
}
bool CDCSortByRincreasing(DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit1, DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit2) {
	// use the ring number to sort by R
	return hit1->hit->wire->ring < hit2->hit->wire->ring;
}
bool CDCSortByZincreasing(DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit1, DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit2) {
	// use the z_stereo position to sort
	return hit1->z_stereo < hit2->z_stereo;
}
bool CDCSortByStereoPhiincreasing(DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit1, DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit2) {
	// use the z_stereo position to sort
	return hit1->phi_stereo < hit2->phi_stereo;
}

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_CDC::init(void)
{
	MAX_SUBSEED_STRAW_DIFF = 1;
	MIN_SUBSEED_HITS = 3;
	MIN_SEED_HITS  = 3;
	MAX_SUBSEED_LINKED_HITS = 12;
	MIN_SEED_DIST = 4.0; // cm
	MAX_HIT_DIST = 3.0; // cm
	MAX_HIT_CIRCLE_DIST = 0.8; // cm
	MAX_SEED_TIME_DIFF = 450.0; // ns
	MAX_CIRCLE_CLONE_FILTER_FRAC = 0.2;
	MAX_STEREO_PHI_DELTA = 0.35; // rad
	TARGET_Z_MIN = 50.0;
	TARGET_Z_MAX = 80.0;
	
	MAX_HIT_DIST2 = MAX_HIT_DIST*MAX_HIT_DIST;

	// Initialize cdchits_by_superlayer with empty vectors for each superlayer
	vector<DCDCTrkHit*> mt;
	for(int i=0; i<5; i++)cdchits_by_superlayer.push_back(mt);
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_CDC::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTrackCandidate_factory_CDC::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_CDC::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_CDC::evnt(JEventLoop *loop, int eventnumber)
{
	// Get CDC hits
	GetCDCHits(loop);
	
	// Find all seeds in all 3 axial superlayers
	vector<DCDCSeed> seeds_sl5;
	vector<DCDCSeed> seeds_sl3;
	vector<DCDCSeed> seeds_sl1;
	if(debug_level>3)_DBG_<<"========= SL5 =========="<<endl;
	FindSeeds(cdchits_by_superlayer[5-1], seeds_sl5);
	if(debug_level>3)_DBG_<<"========= SL3 =========="<<endl;
	FindSeeds(cdchits_by_superlayer[3-1], seeds_sl3);
	if(debug_level>3)_DBG_<<"========= SL1 =========="<<endl;
	FindSeeds(cdchits_by_superlayer[1-1], seeds_sl1);
	
	// Link seeds using average phi of seed hits
	vector<DCDCSeed> seeds_tmp, seeds;
	LinkSeeds(seeds_sl5, seeds_sl3, seeds_tmp, MAX_SUBSEED_LINKED_HITS);
	LinkSeeds(seeds_tmp, seeds_sl1, seeds, 2*MAX_SUBSEED_LINKED_HITS);
	
	// Fit linked seeds to circles
	for(unsigned int i=0; i<seeds.size(); i++)seeds[i].valid = FitCircle(seeds[i]);
	
	// Filter out duplicates of seeds by clearing their "valid" flags
	FilterCloneSeeds(seeds);

	// Extend seeds into stereo layers and do fit to find z and theta
	for(unsigned int i=0; i<seeds.size(); i++){
		DCDCSeed &seed = seeds[i];
		if(debug_level>3)_DBG_<<"----- Seed "<<i<<" ------"<<endl;
		if(debug_level>3)_DBG_<<"seed.fit.phi="<<seed.fit.phi<<endl;
		if(!seed.valid)continue;

		// Add stereo hits to seed
		AddStereoHits(cdchits_by_superlayer[2-1], seed);
		AddStereoHits(cdchits_by_superlayer[4-1], seed);
		
		// If no stereo hits were found for this seed, then
		// we can't fit it.
		if(seed.stereo_hits.size()==0)continue;
		
		// Fit stereo hits to get theta and vertex z position
		FindThetaZ(seed);
		if(!seed.valid)continue;

		// The following is from a fit of ratio of thrown to reconstructed
		// transverse momentum vs. theta for the 1400A field
		//double par[] = {0.984463, 0.150759, -0.414933, 0.257472, -0.055801};
		//double theta = seed.theta;
		//double ff = par[0]+theta*(par[1]+theta*(par[2]+theta*(par[3]+theta*par[4])));
		double p_trans = seed.fit.p_trans*0.95;
		double phi = seed.fit.phi;
		double q = seed.fit.q;
		double theta = seed.theta;

		//Make a track candidate from results
		DTrackCandidate *can = new DTrackCandidate;
		DVector3 pos, mom;
		pos.SetXYZ(0.0, 0.0, seed.z_vertex);
		mom.SetMagThetaPhi(p_trans/sin(theta), theta, phi);
		//pos.SetXYZ(0.0, 0.0, 65.0);
		//mom.SetMagThetaPhi(p_trans, M_PI/2.0, phi);
		can->setPosition(pos);
		can->setMomentum(mom);
		//can->setCharge(seed.q);
		can->setCharge(q);

		//can->q = can->charge();
		//can->phi = mom.Phi();
		//can->theta = mom.Theta();
		//can->p_trans = mom.Pt();
		//can->p = mom.Mag();
		//can->z_vertex = pos.Z();
		//if(debug_level>3)_DBG_<<"can->phi="<<can->phi<<"  phi="<<phi<<" can->theta="<<can->theta<<" theta="<<theta<<endl;
		
		_data.push_back(can);
	}

	return NOERROR;
}

//------------------
// GetCDCHits
//------------------
void DTrackCandidate_factory_CDC::GetCDCHits(JEventLoop *loop)
{
	// Delete old hits
	for(unsigned int i=0; i<cdctrkhits.size(); i++)delete cdctrkhits[i];
	cdctrkhits.clear();
	for(unsigned int i=0; i<cdchits_by_superlayer.size(); i++)cdchits_by_superlayer[i].clear();

	// Get the "raw" hits. These already have the wire associated with them.
	vector<const DCDCTrackHit*> cdctrackhits;
	loop->Get(cdctrackhits);
	
	// If there are no hits, then bail now
	if(cdctrackhits.size()==0)return;
	
	// Create DCDCTrkHit objects out of these.
	for(unsigned int i=0; i<cdctrackhits.size(); i++){
		DCDCTrkHit *cdctrkhit = new DCDCTrkHit;
		cdctrkhit->hit = cdctrackhits[i];
		cdctrkhit->flags = cdctrkhit->hit->wire->stereo==0.0 ? NONE:IS_STEREO;
		cdctrkhit->flags |= NOISE; // (see below)

		// Add to "master" list
		cdctrkhits.push_back(cdctrkhit);
		
		// Sort into list of hits by superlayer
		if(cdctrkhit->hit->wire->ring<=3)cdchits_by_superlayer[0].push_back(cdctrkhit);
		else if(cdctrkhit->hit->wire->ring<= 7)cdchits_by_superlayer[1].push_back(cdctrkhit);
		else if(cdctrkhit->hit->wire->ring<=12)cdchits_by_superlayer[2].push_back(cdctrkhit);
		else if(cdctrkhit->hit->wire->ring<=16)cdchits_by_superlayer[3].push_back(cdctrkhit);
		else if(cdctrkhit->hit->wire->ring<=25)cdchits_by_superlayer[4].push_back(cdctrkhit);
	}
	
	// Sort the individual superlayer lists by decreasing values of R
	for(unsigned int i=0; i<cdchits_by_superlayer.size(); i++){
		sort(cdchits_by_superlayer[i].begin(), cdchits_by_superlayer[i].end(), CDCSortByRdecreasing);
	}
	
	// Filter out noise hits. All hits are initially flagged as "noise".
	// Hits with a neighbor within MAX_HIT_DIST have their noise flags cleared.
	for(unsigned int i=0; i<cdctrkhits.size(); i++){ // cut above should ensure cdctrkhits.size() is at least 1
		DCDCTrkHit *trkhit1 = cdctrkhits[i];
		if(!(trkhit1->flags & NOISE))continue; // this hit already not marked for noise
		for(unsigned int j=0; j<cdctrkhits.size(); j++){
			if(j==i)continue;
			double d2 = trkhit1->Dist2(cdctrkhits[j]);
			if(d2<9.0*MAX_HIT_DIST2){
				trkhit1->flags &= ~NOISE;
				cdctrkhits[j]->flags &= ~NOISE;
				break;
			}// if
		}// j
	}// i
}

//------------------
// FindSeeds
//------------------
void DTrackCandidate_factory_CDC::FindSeeds(vector<DCDCTrkHit*> &hits, vector<DCDCSeed> &seeds)
{
	// Sort through hits ring by ring to find sub-seeds from neighboring wires
	// in the same ring. What we want is a list of sub-seeds for each ring.
	// Each sub-seed is a list of adjacent (or nearly adjacent) hits. Note
	// that hits are ordered by ring, then straw. We ignore the boundary in 
	// straw number for now which will cause some sub-seeds to be listed as
	// 2 separate ones. These will be recombined later.
	DCDCSeed subseed;
	vector<DCDCSeed > subseeds;
	vector<vector<DCDCSeed > > rings;
	int last_ring = -1;
	for(unsigned int i=0; i<hits.size(); i++){
		DCDCTrkHit *trkhit = hits[i];
		
		// Clear USED, IN_SEED, IN_LINE, IN_TRACK flags for all these hits
		trkhit->flags &= ~(USED|IN_SEED|IN_LINE|IN_TRACK);
		
		//if(trkhit->flags & NOISE)continue;

		// Check if ring number has changed. 
		if(trkhit->hit->wire->ring != last_ring){
			if(subseed.hits.size()!=0)subseeds.push_back(subseed);
			if(subseeds.size()!=0)rings.push_back(subseeds);
			if(debug_level>3)_DBG_<<"  subseed hits:"<<subseed.hits.size()<<"  subseeds:"<<subseeds.size()<<endl;
			subseeds.clear();
			subseed.hits.clear();
			subseed.hits.push_back(trkhit);
			subseed.linked=false;
			last_ring = trkhit->hit->wire->ring;
			continue;
		}
		
		// Check if this hit is a neighbor of the last hit added to the subseed
		if(abs(subseed.hits[subseed.hits.size()-1]->hit->wire->straw - trkhit->hit->wire->straw)>MAX_SUBSEED_STRAW_DIFF){
			if(subseed.hits.size()!=0)subseeds.push_back(subseed);
			if(debug_level>3)_DBG_<<"  subseed hits:"<<subseed.hits.size()<<endl;
			subseed.hits.clear();
			subseed.linked=false;
		}
		subseed.hits.push_back(trkhit);
	}
	if(subseed.hits.size()!=0)subseeds.push_back(subseed);
	if(subseeds.size()!=0)rings.push_back(subseeds);
	if(debug_level>3)_DBG_<<"  subseed hits:"<<subseed.hits.size()<<"  subseeds:"<<subseeds.size()<<endl;
	if(debug_level>3)_DBG_<<"rings: "<<rings.size()<<endl;

	// If we have no "rings", then there must be no seeds. Bail now
	if(rings.size()==0)return;
	
	// Loop over rings
	for(ringiter ring =rings.begin(); ring!=rings.end(); ring++){
		vector<DCDCSeed> &subseeds = *ring;
		ringiter next_ring = ring;
		next_ring++;
		
		// Loop over subseeds of this ring
		for(unsigned int j=0; j<subseeds.size(); j++){
			if(subseeds[j].linked)continue;
			
			// This subseed was never used. Start a new seed
			vector<DCDCSeed*> parent;
			parent.push_back(&subseeds[j]);
			LinkSubSeeds(parent, next_ring, rings.end(), seeds);
		}
	}
}

//------------------
// LinkSubSeeds
//------------------
void DTrackCandidate_factory_CDC::LinkSubSeeds(vector<DCDCSeed*> &parent, ringiter ring, ringiter ringend, vector<DCDCSeed> &seeds)
{
	/// Combine subseeds from layers in a single superlayer into seeds.
	///
	/// This a a re-entrant routine (i.e. it calls itself recursively). Upon
	/// entry, <i>parent</i> contains a list of pointers to all of the subseeds
	/// from the rings outside of <i>ring</i> that are to be combined into
	/// a seed. This will search through all subseeds of <i>ring</i> and if
	/// any are found that can extend the parent, a copy of parent is made,
	/// the current subseed of this ring is added to it, and then it is
	/// passed on to another call to this routine. If no matches are found
	/// (which will necessarily be the case for the inner most layer), then
	/// the subseeds in <i>parent</i> will be combined into a single seed and
	/// that added to <i>seeds</i>.
	
	// Make sure parent has at least one subseed
	if(parent.size()==0){
		_DBG_<<"parent has no subseeds!!"<<endl;
		return;
	}
	
	// Set flag to keep track of whether this is the end of the seed or not
	bool seed_extended = false;
	if(ring!=ringend){

		// Last subseed in parent list is the one we need to compare to
		DCDCSeed *parent_subseed = parent[parent.size()-1];
		double r_parent = parent_subseed->hits[0]->hit->wire->origin.Perp();
	
		// Loop over subseeds in this ring
		vector<DCDCSeed> &subseeds = (*ring);
		ring++; // increment ring iterator to point to next level down in case we recall ouself below
		
		for(unsigned int i=0; i<subseeds.size(); i++){
			// Find the difference between this subseed's distance from the beamline
			// and the parent's.
			double dr = r_parent - subseeds[i].hits[0]->hit->wire->origin.Perp();
		
			// Check if this subseed is close enough to the parent's to extend it
			if(parent_subseed->MinDist2(subseeds[i])<(dr*dr+MAX_HIT_DIST2)){
				vector<DCDCSeed*> myparent = parent;
				myparent.push_back(&subseeds[i]);
				subseeds[i].linked = true;
				LinkSubSeeds(myparent, ring, ringend, seeds);
				seed_extended = true;
			}
		}
	}
	
	// Check if this is the end of the line.
	if(!seed_extended){
		// This is the end of this seed. Combine the subseeds into a single
		// seed and add it to the list.
		DCDCSeed seed;
		seed.hits.clear();
		for(unsigned int i=0; i<parent.size(); i++){
			seed.Merge(*parent[i]);
		}
		if(seed.hits.size()>=MIN_SEED_HITS)seeds.push_back(seed);
	}
}

//------------------
// LinkSeeds
//------------------
void DTrackCandidate_factory_CDC::LinkSeeds(vector<DCDCSeed> &in_seeds1, vector<DCDCSeed> &in_seeds2, vector<DCDCSeed> &seeds, unsigned int max_linked_hits)
{
	/// Loop over the two input sets of seeds and compare the phi angles of
	/// their first and last hits. The closest combination is checked to see
	/// if we should combine them into a new a new seed. Any seeds from either
	/// list that are not linked, but have enough hits to be a full seed in and
	/// of themselves will be added to the <i>seeds</i> list.
	if(debug_level>3)_DBG_<<"Linking seeds: Num seeds in list 1:"<<in_seeds1.size()<<" Num seeds in list 2:"<<in_seeds2.size()<<endl;

	// Clear all "linked" flags
	for(unsigned int i=0; i< in_seeds1.size(); i++)in_seeds1[i].linked=false;
	for(unsigned int i=0; i< in_seeds2.size(); i++)in_seeds2[i].linked=false;

	for(unsigned int i=0; i< in_seeds1.size(); i++){
		vector<DCDCTrkHit*> &hits1 = in_seeds1[i].hits;
		if(hits1.size()<1)continue;
		
		for(unsigned int j=0; j< in_seeds2.size(); j++){
			vector<DCDCTrkHit*> &hits2 = in_seeds2[j].hits;
			if(hits2.size()<1)continue;
			
			// Find the hits from the seeds that are closest in R
			const DCDCWire *wire1a = hits1[0]->hit->wire;
			const DCDCWire *wire1b = hits1[hits1.size()-1]->hit->wire;
			const DCDCWire *wire2a = hits2[0]->hit->wire;
			const DCDCWire *wire2b = hits2[hits2.size()-1]->hit->wire;
			const DVector3 *pos1, *pos2;
			if(wire1a->ring > wire2a->ring){
				// seed1 is outside seed2
				pos1 = &(wire1a->ring > wire1b->ring ? wire1b:wire1a)->origin;
				pos2 = &(wire2a->ring < wire2b->ring ? wire2b:wire2a)->origin;
			}else{
				// seed2 is outside seed1
				pos1 = &(wire1a->ring < wire1b->ring ? wire1b:wire1a)->origin;
				pos2 = &(wire2a->ring > wire2b->ring ? wire2b:wire2a)->origin;
			}

			// Compare the phi values of the first and last points of the two seeds
			// to see if they are close enough to link.
			double dphi = fabs(pos1->Phi() - pos2->Phi());
			while(dphi>M_PI)dphi-=2.0*M_PI;
			if(fabs(dphi)<M_PI/6.0){

				DCDCSeed seed = in_seeds1[i];
				seed.Merge(in_seeds2[j]);
				seeds.push_back(seed);
				
				// OK, at this point we have linked the seeds and we *may* want to
				// to set thier "linked" flags. The linked flags are really used
				// to decide later if the partial seed should be promoted to
				// a full seed. The idea here is to only flag the seed as "linked"
				// if the link is strong, otherwise leaved it marked as unlinked
				// so it has the option of becoming another track candidate at
				// the next level. Note that even if we set a linked flag here,
				// the partial seed can still be linked with other partials in
				// this loop.
				//
				// The criteria for a strong link are:
				// 1. the number of hits in the partner's partial seed is not
				//    greater than max_linked_hits (which is passed to us as an
				//    argument)
				// 2. The minimum distance between the partial seeds' ends is less
				//    than 3*MAX_HIT_DIST in the phi direction.
				
				// Find distance between seeds' ends
				double dr = fabs(pos1->Perp() - pos2->Perp());
				DVector3 d = *pos1 - *pos2;
				double dist_phi2 = d.Perp2() - dr*dr;
				
				if(dist_phi2 < 9.0*MAX_HIT_DIST2){
				
					// Mark the seeds has having been linked. Here we do something
					// a little funny. We only mark a seed as linked if the other
					// seed doesn't have too many hits. This is because we can get
					// "super seeds" from a single superlayer made from many hits,
					// often due to crossing tracks.
					if(hits1.size()<=max_linked_hits)in_seeds2[j].linked = true;
					if(hits2.size()<=max_linked_hits)in_seeds1[i].linked = true;
				}
				if(debug_level>3)_DBG_<<"  linked seeds "<<i<<" ("<<hits1.size()<<" hits) and "<<j<<" ("<<hits2.size()<<" hits)  dist_phi="<<sqrt(dist_phi2)<<endl;
			}else{
				if(debug_level>3)_DBG_<<"  rejected link between "<<i<<" and "<<j<<" due to avg. phi (fabs("<<in_seeds1[i].phi_avg<<" - "<<in_seeds2[j].phi_avg<<")>="<<M_PI/6.0<<")"<<endl;
			}
		}
	}
	
	// Any seeds that weren't linked from either list should be promoted
	// to full seeds.
	for(unsigned int i=0; i< in_seeds1.size(); i++){
		if(!in_seeds1[i].linked){
			seeds.push_back(in_seeds1[i]);
			if(debug_level>3)_DBG_<<"Promoting seed "<<i<<" from list 1"<<endl;
		}
	}
	for(unsigned int i=0; i< in_seeds2.size(); i++){
		if(!in_seeds2[i].linked){
			seeds.push_back(in_seeds2[i]);
			if(debug_level>3)_DBG_<<"Promoting seed "<<i<<" from list 2"<<endl;
		}
	}
	
	// Calculate average time for all axial hits in each seed and set
	// all seeds to a charge of +1
	for(unsigned int i=0; i< seeds.size(); i++){
		DCDCSeed &seed = seeds[i];
		seed.tdrift_avg = 0.0;
		unsigned int Nin_time=0;
		for(unsigned int j=0; j<seed.hits.size(); j++){
			double tdrift = seed.hits[j]->hit->tdrift;
			if(tdrift>500.0){
				seed.hits[j]->flags |= OUT_OF_TIME;
				continue;
			}
			seed.tdrift_avg += tdrift;
			Nin_time++;
		}
		seed.tdrift_avg /= (double)Nin_time;
	}
}

//------------------
// FitCircle
//------------------
bool DTrackCandidate_factory_CDC::FitCircle(DCDCSeed &seed)
{
	/// Do a quick fit of the 2-D projection of the axial hits
	/// in seed (seed should contain *only* axial hits) to a 
	/// circle. Determine the sign of the charge (and correspondingly
	/// the initial phi angle) assuming that
	/// the majority of hits come from the outgoing part of the
	/// track. If the resulting circle passes within
	/// MAX_HIT_DIST of more the majority of hits, then true
	/// is returned indicating success. Otherwise, false is
	/// returned indicating failure and that the seed should be
	/// discarded.
	
	// Loop over hits in seed and add them to the seed's DQuickFit object
	seed.fit.Clear();
	for(unsigned int j=0; j<seed.hits.size(); j++){
		if(seed.hits[j]->flags&OUT_OF_TIME)continue;
		const DVector3 &pos = seed.hits[j]->hit->wire->origin;
		seed.fit.AddHitXYZ(pos.X(), pos.Y(), pos.Z());
	}
	seed.fit.FitCircle();
	seed.fit.GuessChargeFromCircleFit();

	// Check if majority of hits are close to circle.
	double x0 = seed.fit.x0;
	double y0 = seed.fit.y0;
	double r0 = sqrt(x0*x0 + y0*y0);
	unsigned int N=0;
	for(unsigned int i=0; i<seed.hits.size(); i++){
		if(seed.hits[i]->flags&OUT_OF_TIME){
			if(debug_level>6)_DBG_<<"discarding out of time hit"<<endl;
			continue;
		}
		double dx = seed.hits[i]->hit->wire->origin.X() - x0;
		double dy = seed.hits[i]->hit->wire->origin.Y() - y0;
		double d = sqrt(dx*dx + dy*dy);
		if(debug_level>10)_DBG_<<"dist="<<d-r0<<endl;
		if(fabs(d-r0)<=MAX_HIT_DIST)N++;
	}
	
	if(debug_level>3)_DBG_<<"Circle fit: Nhits="<<seed.fit.GetNhits()<<"  p_trans="<<seed.fit.p_trans<<" q="<<seed.fit.q<<" N="<<N<<" phi="<<seed.fit.phi<<endl;
	if(N<MIN_SEED_HITS){
	if(debug_level>3)_DBG_<<"Rejected circle fit due to too few hits on track (N="<<N<<" MIN_SEED_HITS="<<MIN_SEED_HITS<<")"<<endl;
		return false;
	}
	
	if(N<seed.hits.size()/2){
	if(debug_level>3)_DBG_<<"Rejected circle fit due to minority of hits on track (N="<<N<<" seed.hits.size()/2="<<seed.hits.size()/2<<")"<<endl;
		return false;
	}
	return true;
}

//------------------
// FilterCloneSeeds
//------------------
void DTrackCandidate_factory_CDC::FilterCloneSeeds(vector<DCDCSeed> &seeds)
{
	/// Look for clones of seeds and flag all but one as not valid.
	/// Two tracks are considered clones if their trajectories are
	/// really close.
	
	/// The "closeness" of the two trajectories is determined using
	/// two methods, either of which can flag the two seeds as clones.
	///
	/// For the first method, we use the "asymmetry" of the circle
	/// centers. The asymmetry is defined as the distance between
	/// the two circle's centers divided by the sum of the magnitudes
	/// of the circle's centers relative to the origin. If this is below
	/// some threshold (say 0.05) then the tracks are seeds are
	/// considered clones.
	///
	/// For the second method, we check the distance between the
	/// trajectories at a point near the last hit and a point near
	/// one of the middle hits of one of the seeds. Since all
	/// seeds pass through the origin, this effectively gives 3 
	/// points on the circle(s). If the total sum of the distances
	/// between these points and the other seed's track are less than
	/// some small value, the two are considered clones. Note that
	/// currently, the choice of "middle" hits is done simply by
	/// looking at the N/2th hit. This may not actually be near the
	/// middle of the R values so a more sophisticated algorithm
	/// may need to be implemented later.

	for(unsigned int i=0; i<seeds.size(); i++){
		DCDCSeed &seed1 = seeds[i];
		if(!seed1.valid)continue;
		
		// Calculate point on seed1 closest to last(furthest out) hit on seed1
		DVector2 r1(seed1.fit.x0, seed1.fit.y0);
		const DVector3 &wire_pos1 = seed1.hits[seed1.hits.size()-1]->hit->wire->origin;
		DVector2 alpha1(wire_pos1.X(), wire_pos1.Y());
		alpha1 -= r1;
		alpha1 *= r1.Mod()/alpha1.Mod();
		
		// Calculate point on seed1 at about the midpoint hit on seed1
		const DVector3 &wire_pos2 = seed1.hits[seed1.hits.size()/2]->hit->wire->origin;
		DVector2 alpha2(wire_pos2.X(), wire_pos2.Y());
		alpha2 -= r1;
		alpha2 *= r1.Mod()/alpha2.Mod();

		for(unsigned int j=i+1; j<seeds.size(); j++){
			DCDCSeed &seed2 = seeds[j];
			if(!seed2.valid)continue;
			
			// Flag to indicate these are clones
			bool are_clones = false;

			// If the two circle fits are significantly different (i.e.
			// their x0, y0 coordinates aren't close), then skip assume
			// these seeds aren't duplicates.
			DVector2 r2(seed2.fit.x0, seed2.fit.y0);
			DVector2 dr = r1-r2;
			double asym = dr.Mod()/(r1.Mod()+r2.Mod());
			if(debug_level>3)_DBG_<<"asym="<<asym<<endl;
			if(asym<0.05) are_clones = true;

			if(!are_clones){
				// Find total absolute distance from alpha1 and seed2 circle and
				// alpha2 and seed2 circle.
				DVector2 d1 = dr+alpha1;
				DVector2 d2 = dr+alpha2;
				double d = fabs(d1.Mod()-r2.Mod()) + fabs(d2.Mod() - r2.Mod());
				if(d<1.5) are_clones = true;
				if(debug_level>3)_DBG_<<"d="<<d<<endl;
			}
			
			// Remove a clone in necessary
			if(are_clones){
				// These seeds appear to be clones of one another. Mark one as not valid.
				if(seed1.fit.chisq<=seed2.fit.chisq){
					seed2.valid = false;
				}else{
					seed1.valid = false;
				}
				if(debug_level>3)_DBG_<<"Filtering clone circle seed (seed1.fit.chisq="<<seed1.fit.chisq<<" seed2.fit.chisq="<<seed2.fit.chisq<<endl;
			}
		}
	}
}

//------------------
// AddStereoHits
//------------------
void DTrackCandidate_factory_CDC::AddStereoHits(vector<DCDCTrkHit*> &stereo_hits, DCDCSeed &seed)
{
	// To find the z coordinate, we look at the 2D projection of the
	// stereo wire and find the intersection point of that with the
	// circle found in FitSeed().

	// Loop over stereo hits to find z-values for any that cross this circle
	for(unsigned int i=0; i<stereo_hits.size(); i++){
		DCDCTrkHit *trkhit = stereo_hits[i];
		
		// Don't consider hits that are out of time with the seed
		if(fabs(trkhit->hit->tdrift-seed.tdrift_avg)>MAX_SEED_TIME_DIFF){
			if(debug_level>3)_DBG_<<"Out of time stereo hit: seed.tdrift_avg="<<seed.tdrift_avg<<" tdrift="<<trkhit->hit->tdrift<<endl;
			continue;
		}
		
		// Calculate intersection points between circle and stereo wire
		const DCDCWire *wire = trkhit->hit->wire;
		DVector2 r1(wire->origin.X(), wire->origin.Y());
		DVector2 r2(wire->udir.X(), wire->udir.Y());
		DVector2 R(seed.fit.x0, seed.fit.y0);
		double r2_mod2 = r2.Mod2();
		if(r2_mod2<1.0E-10){
			_DBG_<<"wire in CDC stereo hit list is not stereo!"<<endl;
			continue; // this must not be a stereo wire!
		}

		double a = r2_mod2;
		double b = 2.0*r2*(r1-R);
		double c = r1.Mod2()-2.0*r1*R;
		double A = b*b - 4.0*a*c;

		// Check that this wire intersects this circle
		if(A<0.0)continue; // line along wire does not intersect circle, ever.

		// Calculate intersection points for two roots of alpha
		double B = sqrt(A);
		double alpha1 = (-b - B)/(2.0*a);
		double alpha2 = (-b + B)/(2.0*a);
		if(debug_level>3)_DBG_<<"alpha1="<<alpha1<<" alpha2="<<alpha2<<endl;
		
		// At this point we must decide which value of alpha to use. The
		// proper way would likely involve either trying both possibilities
		// and taking the one that gave the better chi-sq for a line fit
		// of phi vs. z or looking at the surrounding axial layers
		// and using the value which puts the hit closest to those.
		// For now, we just use the value closest to zero (i.e. closest to
		// the center of the wire).
		double alpha = sqrt(r2_mod2)*(fabs(alpha1)<fabs(alpha2) ? alpha1:alpha2);
		
		// Now we must convert the alpha value into a z-value. To do this,
		// we use the known theta angle of the wire. The distance along the
		// wire from the center in 3D is given by:
		//
		//   s = alpha/sin(theta_wire)
		//
		// The z coordinate of the hit is given by:
		//
		//    z = z1 + s*z2
		//
		// where z1 is the z coordinate of the center of the wire and z2
		// is the z component of the direction of the wire. i.e.
		//
		//    z2 = cos(theta_wire)
		//
		// This means  sin(theta_wire) = sqrt(1 - (z2)^2)
		//double z2 = wire->udir.Z();
		double s = alpha/fabs(sin(wire->stereo));
		if(debug_level>3)_DBG_<<"alpha="<<alpha<<" s="<<s<<" ring="<<wire->ring<<" straw="<<wire->straw<<" stereo="<<wire->stereo<<endl;
		if(fabs(s) > wire->L/2.0)continue; // if wire doesn't cross circle, skip hit
		
		// Add this hit to the stereo hit list
		DVector3 pos = wire->origin + s*wire->udir;
		trkhit->x_stereo = pos.X();
		trkhit->y_stereo = pos.Y();
		trkhit->z_stereo = pos.Z();
		trkhit->phi_stereo = atan2(trkhit->y_stereo-R.Y(), trkhit->x_stereo-R.X());
		R*=-1.0; // make R point from center of circle to beamline instead of other way around
		if(debug_level>3){
			_DBG_<<" --- wire->udir X, Y, Z = "<<wire->udir.X()<<", "<<wire->udir.Y()<<", "<<wire->udir.Z()<<endl;
			_DBG_<<" -- ring="<<wire->ring<<" trkhit->z_stereo="<<trkhit->z_stereo<<" trkhit->y_stereo="<<trkhit->y_stereo<<" trkhit->x_stereo="<<trkhit->x_stereo<<endl;
			_DBG_<<" -- phi_stereo="<<trkhit->phi_stereo<<"  R.Phi()="<<R.Phi()<<"  (X,Y)=("<<R.X()<<", "<<R.Y()<<")"<<endl;
			_DBG__;
		}
		trkhit->phi_stereo -= R.Phi(); // make angle relative to beamline
		
		// We want this to go either from 0 to +2pi for positive charge, or
		// 0 to -2pi for negative.
		double phi_hi = seed.fit.q>0.0 ? +2.0*M_PI:0.0;
		double phi_lo = seed.fit.q>0.0 ? 0.0:-2.0*M_PI;
		while(trkhit->phi_stereo<phi_lo)trkhit->phi_stereo+=2.0*M_PI;
		while(trkhit->phi_stereo>phi_hi)trkhit->phi_stereo-=2.0*M_PI;
		trkhit->flags |= VALID_STEREO;
		seed.stereo_hits.push_back(trkhit);
	}
}

//------------------
// FindThetaZ
//------------------
void DTrackCandidate_factory_CDC::FindThetaZ(DCDCSeed &seed)
{
	FindTheta(seed, TARGET_Z_MIN, TARGET_Z_MAX);
	FindZ(seed, seed.theta_min, seed.theta_max);
	
	// If z_vertex is not inside the target limits, then flag this
	// seed as invalid.
	if(seed.z_vertex<TARGET_Z_MIN || seed.z_vertex>TARGET_Z_MAX){
		if(debug_level>3)_DBG_<<"Seed z-vertex outside of target range (z="<<seed.z_vertex<<" TARGET_Z_MIN="<<TARGET_Z_MIN<<" TARGET_Z_MAX="<<TARGET_Z_MAX<<endl;
		seed.valid=false;
	}
	
	return;
}

//------------------
// FindTheta
//------------------
void DTrackCandidate_factory_CDC::FindTheta(DCDCSeed &seed, double target_z_min, double target_z_max)
{
	/// Find the theta value using the valid stereo hits from <i>seed</i>. The values
	/// for phi_stereo and z_stereo are assumed to be valid as is the status of the
	/// VALID_STEREO bit in flags. The value of seed.fit.r0 is also used to calculate
	/// theta.
	///
	/// This uses a histogramming technique that looks at the overlaps of the
	/// angle ranges subtended by each hit between the given target limits.
	/// The overlaps usually lead to a range of values for theta. The limits
	/// of these are stored in the theta_min and theta_max fields of the seed.
	/// The centroid of the range is stored in the theta field.
	
	// We use a simple array to store our histogram here. We don't want to use
	// ROOT histograms because they are not thread safe.
	unsigned int Nbins = 1000;
	unsigned int hist[Nbins];
	for(unsigned int i=0; i<Nbins; i++)hist[i] = 0; // clear histogram
	double bin_width = 2.0*M_PI/(double)Nbins;
	double hist_low_limit = -M_PI; // lower edge of histogram limits
	
	// Loop over valid hits, filling the histogram
	double r0 = seed.fit.r0;
	for(unsigned int i=0; i<seed.stereo_hits.size(); i++){
		DCDCTrkHit *trkhit = seed.stereo_hits[i];
		if(!trkhit->flags&VALID_STEREO)continue;
		
		// Calculate upper and lower limits in theta
		double alpha = r0*trkhit->phi_stereo;
		if(seed.fit.q<0.0)alpha = -alpha;
		double tmin = atan2(alpha, trkhit->z_stereo - target_z_min);
		double tmax = atan2(alpha, trkhit->z_stereo - target_z_max);
		if(tmin>tmax){
			double tmp = tmin;
			tmin=tmax;
			tmax=tmp;
		}
		if(debug_level>3)_DBG_<<" -- phi_stereo="<<trkhit->phi_stereo<<" z_stereo="<<trkhit->z_stereo<<endl;
		if(debug_level>3)_DBG_<<" -- tmin="<<tmin<<"  tmax="<<tmax<<endl;
		
		// Find index of bins corresponding to tmin and tmax
		unsigned int imin = (unsigned int)floor((tmin-hist_low_limit)/bin_width);
		unsigned int imax = (unsigned int)floor((tmax-hist_low_limit)/bin_width);
		
		// If entire range of this hit is outside of the histogram limit
		// then discard this hit.
		if(imax<0 || imin>=Nbins)continue;
		
		// Clip limits of bin range to our histogram limits
		if(imin<0)imin=0;
		if(imin>=Nbins)imin=Nbins-1;
		if(imax<0)imax=0;
		if(imax>=Nbins)imax=Nbins-1;
		
		// Increment all bins between imin and imax
		for(unsigned int j=imin; j<=imax; j++)hist[j]++;
	}
	
	// Look for the indexes of the plateau
	unsigned int istart=0;
	unsigned int iend=0;
	for(unsigned int i=1; i<Nbins; i++){
		if(hist[i]>hist[istart]){
			istart = i;
			if(debug_level>3)_DBG_<<" -- istart="<<istart<<" (theta="<<hist_low_limit + bin_width*(0.5+(double)istart)<<" , N="<<hist[i]<<")"<<endl;
		}
		if(hist[i] == hist[istart])iend = i;
	}
	
	// If there are no entries in the histogram, then flag this seed as invalid
	if(hist[istart]==0.0)seed.valid=false;
	
	// Calculate theta limits
	seed.theta_min = hist_low_limit + bin_width*(0.5+(double)istart);
	seed.theta_max = hist_low_limit + bin_width*(0.5+(double)iend);
	seed.theta = (seed.theta_max + seed.theta_min)/2.0;
	if(debug_level>3)_DBG_<<"istart="<<istart<<" iend="<<iend<<" theta_min="<<seed.theta_min<<" theta_max="<<seed.theta_max<<endl;
}

//------------------
// FindZ
//------------------
void DTrackCandidate_factory_CDC::FindZ(DCDCSeed &seed, double theta_min, double theta_max)
{
	/// Find the z value of the vertex using the valid stereo hits from <i>seed</i>. The values
	/// for phi_stereo and z_stereo are assumed to be valid as is the status of the
	/// VALID_STEREO bit in flags.
	///
	/// This uses a histogramming technique that looks at the overlaps of the
	/// z ranges subtended by each hit between the given theta limits.
	/// The overlaps usually lead to a range of values for z_vertex. The limits
	/// of these are stored in the z_min and z_max fields of the seed.
	/// The centroid of the range is stored in the z_vertex field.
	
	// We use a simple array to store our histogram here. We don't want to use
	// ROOT histograms because they are not thread safe.
	unsigned int Nbins = 300;
	unsigned int hist[Nbins];
	for(unsigned int i=0; i<Nbins; i++)hist[i] = 0; // clear histogram
	double bin_width = 0.5; // bins are 5mm
	double hist_low_limit = 0.0; // lower edge of histogram limits
	
	// Loop over valid hits, filling the histogram
	double r0 = seed.fit.r0;
	double tan_alpha_min = tan(theta_min)/r0;
	double tan_alpha_max = tan(theta_max)/r0;
	for(unsigned int i=0; i<seed.stereo_hits.size(); i++){
		DCDCTrkHit *trkhit = seed.stereo_hits[i];
		if(!trkhit->flags&VALID_STEREO)continue;
		
		// Calculate upper and lower limits in z
		double q_sign = seed.fit.q>0.0 ? +1.0:-1.0;
		double zmin = trkhit->z_stereo - q_sign*trkhit->phi_stereo/tan_alpha_min;
		double zmax = trkhit->z_stereo - q_sign*trkhit->phi_stereo/tan_alpha_max;
		if(zmin>zmax){
			double tmp = zmin;
			zmin=zmax;
			zmax=tmp;
		}
		if(debug_level>3)_DBG_<<" -- phi_stereo="<<trkhit->phi_stereo<<" z_stereo="<<trkhit->z_stereo<<endl;
		if(debug_level>3)_DBG_<<" -- zmin="<<zmin<<"  zmax="<<zmax<<endl;
		
		// Find index of bins corresponding to tmin and tmax
		unsigned int imin = (unsigned int)floor((zmin-hist_low_limit)/bin_width);
		unsigned int imax = (unsigned int)floor((zmax-hist_low_limit)/bin_width);
		
		// If entire range of this hit is outside of the histogram limit
		// then discard this hit.
		if(imax<=0 || imin>=Nbins)continue;
		
		// Clip limits of bin range to our histogram limits
		if(imin<0)imin=0;
		if(imin>=Nbins)imin=Nbins-1;
		if(imax<0)imax=0;
		if(imax>=Nbins)imax=Nbins-1;
		
		// Increment all bins between imin and imax
		for(unsigned int j=imin; j<=imax; j++)hist[j]++;
	}
	
	// Look for the indexes of the plateau
	unsigned int istart=0;
	unsigned int iend=0;
	for(unsigned int i=1; i<Nbins; i++){
		if(hist[i]>hist[istart]){
			istart = i;
			if(debug_level>3)_DBG_<<" -- istart="<<istart<<" (z="<<hist_low_limit + bin_width*(0.5+(double)istart)<<" , N="<<hist[i]<<")"<<endl;
		}
		if(hist[i] == hist[istart])iend = i;
	}

	// If there are no entries in the histogram, then flag this seed as invalid
	if(hist[istart]==0.0)seed.valid=false;
	
	// Calculate z limits
	seed.z_min = hist_low_limit + bin_width*(0.5+(double)istart);
	seed.z_max = hist_low_limit + bin_width*(0.5+(double)iend);
	seed.z_vertex = (seed.z_max + seed.z_min)/2.0;
	if(debug_level>3)_DBG_<<"istart="<<istart<<" iend="<<iend<<" z_min="<<seed.z_min<<" z_max="<<seed.z_max<<" hits[istart]="<<hist[istart]<<endl;
}

//------------------
// NumEligibleSeedHits
//------------------
int DTrackCandidate_factory_CDC::NumEligibleSeedHits(vector<DCDCTrkHit*> &hits)
{
	int N=0;
	for(unsigned int i=0; i<hits.size(); i++){
		if((hits[i]->flags & (USED | CANT_BE_IN_SEED | OUT_OF_TIME)) == 0)N++;
	}
	
	return N;
}

//------------------
// DCDCSeed::Merge
//------------------
void DTrackCandidate_factory_CDC::DCDCSeed::Merge(DCDCSeed& seed)
{
	phi_avg *= (double)hits.size();
	for(unsigned int i=0; i<seed.hits.size(); i++){
		hits.push_back(seed.hits[i]);
		phi_avg += seed.hits[i]->hit->wire->phi;
	}
	phi_avg/=(double)hits.size();
}

//------------------
// DCDCSeed::MinDist2
//------------------
double DTrackCandidate_factory_CDC::DCDCSeed::MinDist2(DCDCSeed& seed)
{
	/// Returns the minimum distance squared between this and the given seed.
	/// Only the first and last hits of each seed's hit list are used
	/// to calculate a maximum of 4 distances (minimum of 1) of which, the
	/// smallest is returned.
	if(hits.size()==0 || seed.hits.size()==0){
		_DBG_<<"Number of seed hits 0! (Nthis="<<hits.size()<<" ,Nseed="<<seed.hits.size()<<")"<<endl;
		return 1.0E10;
	}
	
	// This is kind of messy looking, but it avoids extra calls to
	// DCDCTrkHit::Dist2 that aren't necessary.
	double d2min = hits[0]->Dist2(seed.hits[0]);
	if(hits.size()>0){
		double d2 = hits[hits.size()-1]->Dist2(seed.hits[0]);
		if(d2<d2min)d2min = d2;
		
		if(seed.hits.size()>0){
			double d2 = hits[hits.size()-1]->Dist2(seed.hits[seed.hits.size()-1]);
			if(d2<d2min)d2min = d2;
		}
	}
	if(seed.hits.size()>0){
		double d2 = hits[0]->Dist2(seed.hits[seed.hits.size()-1]);
		if(d2<d2min)d2min = d2;
	}
	
	return d2min;
}

//------------------
// toString
//------------------
const string DTrackCandidate_factory_CDC::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("      id: Nhits: q:     p:       theta:   phi:  p_trans:   x:     y:     z:    dz/dphi:");

	for(unsigned int i=0; i<_data.size(); i++){
		DTrackCandidate *trackcandidate = _data[i];
		printnewrow();
		
		printcol("%lx",    trackcandidate->id);
#if 0
		printcol("%d",    trackcandidate->hitid.size());
		printcol("%+d", (int)trackcandidate->q);
		printcol("%3.3f", trackcandidate->p);
		printcol("%1.3f", trackcandidate->theta);
		printcol("%1.3f", trackcandidate->phi);
		printcol("%3.2f", trackcandidate->p_trans);
		printcol("%2.2f", trackcandidate->x0);
		printcol("%2.2f", trackcandidate->y0);
		printcol("%2.2f", trackcandidate->z_vertex);
		printcol("%1.3f", trackcandidate->dzdphi);
#endif
		printrow();
	}
	
	return _table;
}

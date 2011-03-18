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
#include "DHelicalFit.h"
#include <HDGEOMETRY/DGeometry.h>
#include <HDGEOMETRY/DMagneticFieldMap.h>
#include <FDC/DFDCPseudo.h>
#define BeamRMS 0.1
#define Z_TARGET 65.0

#include "DTrackCandidate_factory_CDC.h"

class DVector3_with_perp:public DVector3
{
	public:
		DVector3_with_perp():DVector3(){CalcPerp();}
		DVector3_with_perp(double x, double y, double z):DVector3(x,y,z){CalcPerp();}
		double CalcPerp(void){
			perp = Perp();
			return perp;
		}
		double perp;
};

bool SortIntersections(const DVector3_with_perp &a,const DVector3_with_perp &b){
  if (a.perp<b.perp) return true;
  return false;
}

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
	MAX_ALLOWED_CDC_HITS = 1000;
	MAX_SUBSEED_STRAW_DIFF = 1;
	MIN_SEED_HITS  = 2;
	MAX_SUBSEED_LINKED_HITS = 12;
	MAX_RING_SUBSEED_HITS = 4;
	MAX_HIT_DIST = 4.0; // cm
	MAX_SEED_TIME_DIFF = 450.0; // ns
	MAX_CDC_MATCH_ANGLE = 20.0; // degrees
	MAX_FDC_MATCH_ANGLE = 40.0; // degrees
	MAX_SEED_LINK_ANGLE = M_PI/6.0*57.3; // degrees
	TARGET_Z_MIN = 50.0;
	TARGET_Z_MAX = 80.0;
	DEBUG_LEVEL = 0;
	
	// Initialize cdchits_by_superlayer with empty vectors for each superlayer
	vector<DCDCTrkHit*> mt;
	for(int i=0; i<5; i++)cdchits_by_superlayer.push_back(mt);
	
	// Set the layer numbers defining the superlayer boundaries
	superlayer_boundaries.push_back( 4);
	superlayer_boundaries.push_back(12);
	superlayer_boundaries.push_back(16);
	superlayer_boundaries.push_back(24);
	superlayer_boundaries.push_back(28);
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_CDC::brun(JEventLoop *eventLoop, int runnumber)
{
	gPARMS->SetDefaultParameter("TRKFIND:MAX_ALLOWED_CDC_HITS", MAX_ALLOWED_CDC_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SUBSEED_STRAW_DIFF", MAX_SUBSEED_STRAW_DIFF);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_SEED_HITS", MIN_SEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SUBSEED_LINKED_HITS", MAX_SUBSEED_LINKED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_RING_SUBSEED_HITS", MAX_RING_SUBSEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_HIT_DIST", MAX_HIT_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_TIME_DIFF", MAX_SEED_TIME_DIFF);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_CDC_MATCH_ANGLE", MAX_CDC_MATCH_ANGLE);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_FDC_MATCH_ANGLE", MAX_FDC_MATCH_ANGLE);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_LINK_ANGLE", MAX_SEED_LINK_ANGLE);
	gPARMS->SetDefaultParameter("TRKFIND:DEBUG_LEVEL", DEBUG_LEVEL);
	
	MAX_HIT_DIST2 = MAX_HIT_DIST*MAX_HIT_DIST;

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
	if(DEBUG_LEVEL>3)_DBG_<<"========= SL5 =========="<<endl;
	FindSeeds(cdchits_by_superlayer[5-1], seeds_sl5);
	if(DEBUG_LEVEL>3)_DBG_<<"========= SL3 =========="<<endl;
	FindSeeds(cdchits_by_superlayer[3-1], seeds_sl3);
	if(DEBUG_LEVEL>3)_DBG_<<"========= SL1 =========="<<endl;
	FindSeeds(cdchits_by_superlayer[1-1], seeds_sl1);
	
	// Link seeds using average phi of seed hits
	vector<DCDCSeed> seeds_tmp, seeds;
	LinkSeeds(seeds_sl5, seeds_sl3, seeds_tmp, MAX_SUBSEED_LINKED_HITS);
	LinkSeeds(seeds_tmp, seeds_sl1, seeds, 2*MAX_SUBSEED_LINKED_HITS);
	
	// Check to add lone hits in SL3 to seeds with only SL1 hits
	PickupUnmatched(seeds);
	
	// Drop seeds containing hits from only a single axial superlayer besides SL1
	DropIncompleteSeeds(seeds);
	
	// Fit linked seeds to circles
	for(unsigned int i=0; i<seeds.size(); i++){
		if(DEBUG_LEVEL>5)_DBG_<<"-- Starting fit for seed "<<i<<" Nhits="<<seeds[i].hits.size()<<" phi_avg="<<seeds[i].phi_avg<<endl;
		seeds[i].valid = FitCircle(seeds[i]);
	}
	
	// Filter out duplicates of seeds by clearing their "valid" flags
	FilterCloneSeeds(seeds);

	// Extend seeds into stereo layers and do fit to find z and theta
	for(unsigned int i=0; i<seeds.size(); i++){
		DCDCSeed &seed = seeds[i];
		if(DEBUG_LEVEL>3)_DBG_<<"----- Seed "<<i<<" ------"<<endl;
		if(DEBUG_LEVEL>3)_DBG_<<"seed.fit.phi="<<seed.fit.phi<<endl;
		if(!seed.valid)continue;

		// Add stereo hits to seed
		AddStereoHits(cdchits_by_superlayer[2-1], seed);
		AddStereoHits(cdchits_by_superlayer[4-1], seed);
		
		// If no stereo or FDC hits were found for this seed, then
		// we can't fit it.
		if(seed.stereo_hits.size()==0 && seed.fdchits.size()==0)continue;

		if (FindThetaZRegression(seed)!=NOERROR){
		  // If the linear regression doesn't work try the histogramming method
		  // Fit stereo hits to get theta and vertex z position
		  FindThetaZ(seed);
		  if(!seed.valid){
		    continue;

		    // Assume that the track came from one end or the other 
		    // of the target and use a point in one of the stereo 
		    // layers to estimate tanl
		    if (seed.z_vertex>TARGET_Z_MAX)
		      seed.z_vertex=TARGET_Z_MAX;
		    else
		      seed.z_vertex=TARGET_Z_MIN;
		    double x=seed.stereo_hits[0].x_stereo;
		    double y=seed.stereo_hits[0].y_stereo;
		    double ratio=sqrt(x*x+y*y)/2./seed.fit.r0;
		    if (ratio<1.){
		      double tanl=(seed.stereo_hits[0].z_stereo-seed.z_vertex)/
			(2.*seed.fit.r0*asin(ratio));
		      seed.theta=M_PI_2-atan(tanl);
		    }
		    
		  }
		}
	
		// The following is from a fit of ratio of thrown to reconstructed
		// transverse momentum vs. theta for the 1400A field
		//double par[] = {0.984463, 0.150759, -0.414933, 0.257472, -0.055801};
		//double theta = seed.theta;
		//double ff = par[0]+theta*(par[1]+theta*(par[2]+theta*(par[3]+theta*par[4])));
		double p_trans = seed.fit.p_trans*seed.FindAverageBz(loop)/2.0;
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
	
		for (unsigned int n=0;n<seed.hits.size();n++){
		  const DCDCTrackHit *cdchit=(seed.hits[n])->hit;
		  can->AddAssociatedObject(cdchit);
		}
		for (unsigned int n=0;n<seed.stereo_hits.size();n++){
		  const DCDCTrackHit *cdchit=(seed.stereo_hits[n]).hit;
		  can->AddAssociatedObject(cdchit);
		}

	
		//can->q = can->charge();
		//can->phi = mom.Phi();
		//can->theta = mom.Theta();
		//can->p_trans = mom.Pt();
		//can->p = mom.Mag();
		//can->z_vertex = pos.Z();
		if(DEBUG_LEVEL>3)_DBG_<<"Final Candidate parameters: p="<<mom.Mag()<<" theta="<<mom.Theta()<<"  phi="<<mom.Phi()<<" z="<<pos.Z()<<endl;
		
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
	
	// If there are too many hits, bail with a warning message
	if(cdctrackhits.size()>MAX_ALLOWED_CDC_HITS){
		_DBG_<<"Too many hits in CDC ("<<cdctrackhits.size()<<", max="<<MAX_ALLOWED_CDC_HITS<<")! Track finding in CDC bypassed for event "<<loop->GetJEvent().GetEventNumber()<<endl;
		cdctrackhits.clear();
		return;
	}
	
	// Create DCDCTrkHit objects out of these.
	for(unsigned int i=0; i<cdctrackhits.size(); i++){
	  DCDCTrkHit *cdctrkhit = new DCDCTrkHit;
	  cdctrkhit->hit = cdctrackhits[i];
	  cdctrkhit->flags = cdctrkhit->hit->wire->stereo==0.0 ? NONE:IS_STEREO;
	  cdctrkhit->flags |= NOISE; // (see below)
	  
	  // Add to "master" list
	  cdctrkhits.push_back(cdctrkhit);
	  
		// Sort into list of hits by superlayer
#if 1
		for(unsigned int j=0; j<superlayer_boundaries.size(); j++){
			if(cdctrkhit->hit->wire->ring<=superlayer_boundaries[j]){
				cdchits_by_superlayer[j].push_back(cdctrkhit);
				break;
			}
		}
#else
	  if(cdctrkhit->hit->wire->ring<=4)cdchits_by_superlayer[0].push_back(cdctrkhit);
	  else if(cdctrkhit->hit->wire->ring<= 12)cdchits_by_superlayer[1].push_back(cdctrkhit);
	  else if(cdctrkhit->hit->wire->ring<=16)cdchits_by_superlayer[2].push_back(cdctrkhit);
	  else if(cdctrkhit->hit->wire->ring<=24)cdchits_by_superlayer[3].push_back(cdctrkhit);
	  else if(cdctrkhit->hit->wire->ring<=28)cdchits_by_superlayer[4].push_back(cdctrkhit);
#endif
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
			if(DEBUG_LEVEL>3)_DBG_<<"  subseed hits:"<<subseed.hits.size()<<"  subseeds:"<<subseeds.size()<<endl;
			subseeds.clear();
			subseed.hits.clear();
			subseed.hits.push_back(trkhit);
			subseed.linked=false;
			last_ring = trkhit->hit->wire->ring;
			continue;
		}
		
		// Check if this hit is a neighbor of the last hit added to the subseed
		if((unsigned int)abs(subseed.hits[subseed.hits.size()-1]->hit->wire->straw - trkhit->hit->wire->straw)>MAX_SUBSEED_STRAW_DIFF){
			if(subseed.hits.size()!=0)subseeds.push_back(subseed);
			if(DEBUG_LEVEL>3)_DBG_<<"  subseed hits:"<<subseed.hits.size()<<endl;
			subseed.hits.clear();
			subseed.linked=false;
		}
		subseed.hits.push_back(trkhit);
	}
	if(subseed.hits.size()!=0)subseeds.push_back(subseed);
	if(subseeds.size()!=0)rings.push_back(subseeds);
	if(DEBUG_LEVEL>3)_DBG_<<"  subseed hits:"<<subseed.hits.size()<<"  subseeds:"<<subseeds.size()<<endl;
	if(DEBUG_LEVEL>3)_DBG_<<"rings: "<<rings.size()<<endl;

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
			if(subseeds[j].hits.size()>MAX_RING_SUBSEED_HITS){
				if(DEBUG_LEVEL>4)_DBG_<<"rejecting subseed for having too many hits in a single layer ("<<subseeds[j].hits.size()<<" max="<<MAX_RING_SUBSEED_HITS<<")"<<endl;
				continue;
			}

			// This subseed was never used. Start a new seed
			vector<DCDCSeed*> parent;
			parent.push_back(&subseeds[j]);
			LinkSubSeeds(parent, next_ring, rings.end(), seeds);
		}
	}
	
	// Mark all seeds as valid
	for(unsigned int i=0; i<seeds.size(); i++){
		seeds[i].valid = true;
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
		if(seed.hits.size()>=MIN_SEED_HITS){
			seeds.push_back(seed);
		}else{
			if(DEBUG_LEVEL>1)_DBG_<<"rejecting seed due to too few hits (have "<<seed.hits.size()<<" need "<<MIN_SEED_HITS<<")"<<endl;
		}
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
	if(DEBUG_LEVEL>3)_DBG_<<"Linking seeds: Num seeds in list 1:"<<in_seeds1.size()<<" Num seeds in list 2:"<<in_seeds2.size()<<endl;

	// Clear all "linked" flags
	for(unsigned int i=0; i< in_seeds1.size(); i++)in_seeds1[i].linked=false;
	for(unsigned int i=0; i< in_seeds2.size(); i++)in_seeds2[i].linked=false;

	for(unsigned int i=0; i< in_seeds1.size(); i++){
		vector<DCDCTrkHit*> &hits1 = in_seeds1[i].hits;
		if(DEBUG_LEVEL>5)_DBG_<<"hits1.size()="<<hits1.size()<<endl;
		if(hits1.size()<1)continue;
		
		for(unsigned int j=0; j< in_seeds2.size(); j++){
			vector<DCDCTrkHit*> &hits2 = in_seeds2[j].hits;
			if(DEBUG_LEVEL>5)_DBG_<<" hits2.size()="<<hits2.size()<<endl;
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
			if(fabs(dphi*57.3)<MAX_SEED_LINK_ANGLE){

				DCDCSeed seed = in_seeds1[i];
				seed.Merge(in_seeds2[j]);
				seeds.push_back(seed);
				
				// Mark all of the hits in the new, merged seed as USED. Some will
				// already have the flag set by having been previously merged. The
				// USED flag is used later when matching SL1 only seeds with stray
				// hits in either SL3 or the FDC.
				for(unsigned int k=0; k<seed.hits.size(); k++)seed.hits[k]->flags |= USED;
				
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
				}else{
					if(DEBUG_LEVEL>3)_DBG_<<"   not marking seeds as linked (dist_phi="<<sqrt(dist_phi2)<<" cm > "<<sqrt(9.0*MAX_HIT_DIST2)<<" cm)"<<endl;
				}
				if(DEBUG_LEVEL>3)_DBG_<<"  linked seeds "<<i<<" ("<<hits1.size()<<" hits) and "<<j<<" ("<<hits2.size()<<" hits)  dist_phi="<<sqrt(dist_phi2)<<endl;
			}else{
				if(DEBUG_LEVEL>3)_DBG_<<"  rejected link between "<<i<<" and "<<j<<" due to avg. phi (fabs("<<in_seeds1[i].phi_avg*57.3<<" - "<<in_seeds2[j].phi_avg*57.3<<")>="<<MAX_SEED_LINK_ANGLE<<")"<<endl;
			}
		}
	}
	
	// Any seeds that weren't linked from either list should be promoted
	// to full seeds.
	for(unsigned int i=0; i< in_seeds1.size(); i++){
		if(!in_seeds1[i].linked){
			seeds.push_back(in_seeds1[i]);
			if(DEBUG_LEVEL>3)_DBG_<<"Promoting seed "<<i<<" from list 1"<<endl;
		}
	}
	for(unsigned int i=0; i< in_seeds2.size(); i++){
		if(!in_seeds2[i].linked){
			seeds.push_back(in_seeds2[i]);
			if(DEBUG_LEVEL>3)_DBG_<<"Promoting seed "<<i<<" from list 2"<<endl;
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
	
	// Loop over hits in seed and add them to the seed's DHelicalFit object
	seed.fit.Clear();
	for(unsigned int j=0; j<seed.hits.size(); j++){
		if(seed.hits[j]->flags&OUT_OF_TIME)continue;
		const DVector3 &pos = seed.hits[j]->hit->wire->origin;
		seed.fit.AddHitXYZ(pos.X(), pos.Y(),pos.Z());
	}
	
	// Add in any FDC hits
	for(unsigned int j=0; j<seed.fdchits.size(); j++){
		const DFDCPseudo *hit = seed.fdchits[j];
		seed.fit.AddHitXYZ(hit->xy.X(), hit->xy.Y(), hit->wire->origin.Z());
	}

	// Try and fit the circle using a Riemann fit. If 
	// it fails, try a basic fit with QuickFit.
	if(seed.fit.FitCircleRiemann(BeamRMS)!=NOERROR || seed.fit.chisq>20){
	  if(DEBUG_LEVEL>3)_DBG_<<"Riemann fit failed. Attempting regression fit..."<<endl;
	  if(seed.fit.FitCircle()!=NOERROR || seed.fit.chisq>20){
	    if(DEBUG_LEVEL>3)_DBG_<<"Regression circle fit failed. Trying straight line."<<endl;
	    if(seed.fit.FitCircleStraightTrack()!=NOERROR){
	      if(DEBUG_LEVEL>3)_DBG_<<"Trying straight line fit failed. Giving up."<<endl;
	      return false;
	    }
	  }else{
	    seed.fit.GuessChargeFromCircleFit(); // for Riemann fit
		}
	}else{
		seed.fit.GuessChargeFromCircleFit(); // for regression fit
	}

	// Check if majority of hits are close to circle.
	double x0 = seed.fit.x0;
	double y0 = seed.fit.y0;
	double r0 = seed.fit.r0;
	unsigned int N=0;
	for(unsigned int i=0; i<seed.hits.size(); i++){
		if(seed.hits[i]->flags&OUT_OF_TIME){
			if(DEBUG_LEVEL>6)_DBG_<<"discarding out of time hit"<<endl;
			continue;
		}
		double dx = seed.hits[i]->hit->wire->origin.X() - x0;
		double dy = seed.hits[i]->hit->wire->origin.Y() - y0;
		double d = sqrt(dx*dx + dy*dy);
		if(DEBUG_LEVEL>15)_DBG_<<"dist="<<d-r0<<endl;
		if(fabs(d-r0)<=MAX_HIT_DIST)N++;
	}
	for(unsigned int i=0; i<seed.fdchits.size(); i++){
		const DFDCPseudo *hit = seed.fdchits[i];
		double dx = hit->xy.X() - x0;
		double dy = hit->xy.Y() - y0;
		double d = sqrt(dx*dx + dy*dy);
		if(DEBUG_LEVEL>10)_DBG_<<"dist="<<d-r0<<endl;
		if(fabs(d-r0)<=MAX_HIT_DIST)N++;
	}
	
	if(DEBUG_LEVEL>3)_DBG_<<"Circle fit: Nhits="<<seed.fit.GetNhits()<<"  p_trans="<<seed.fit.p_trans<<" q="<<seed.fit.q<<" N="<<N<<" phi="<<seed.fit.phi<<endl;
	if(N<MIN_SEED_HITS){
	if(DEBUG_LEVEL>3)_DBG_<<"Rejected circle fit due to too few hits on track (N="<<N<<" MIN_SEED_HITS="<<MIN_SEED_HITS<<")"<<endl;
		return false;
	}
	
	if(N<seed.hits.size()/2){
	if(DEBUG_LEVEL>3)_DBG_<<"Rejected circle fit due to minority of hits on track (N="<<N<<" seed.hits.size()/2="<<seed.hits.size()/2<<")"<<endl;
		return false;
	}
	return true;
}

//------------------
// PickupUnmatched
//------------------
void DTrackCandidate_factory_CDC::PickupUnmatched(vector<DCDCSeed> &seeds)
{
	/// Look for single hits in superlayers that did not have enough
	/// neighbors to make seeds of their own, but could be part of
	/// a seed made in another axial superlayer. This is mainly here
	/// to address tracks that have a single hit in SL3 that is
	/// otherwise discarded. Seeds with only hits in SL1 often
	/// don't enough information to form a good candidate so we try
	/// and pick up these single hits in order to get a good handle
	/// on the candidate.

	// Loop through seeds, looking for ones with only SL1 hits
	for(unsigned int i=0; i<seeds.size(); i++){
		if(!seeds[i].valid)continue;
		vector<DCDCTrkHit*> &hits = seeds[i].hits;
		if(hits.size()<1)continue;
		bool has_non_SL1_hit = false;
		for(unsigned int j=0; j<hits.size(); j++){
			if(hits[j]->hit->wire->ring>superlayer_boundaries[0]){
				has_non_SL1_hit = true;
				break;
			}
		}
		if(has_non_SL1_hit)continue;
		if(DEBUG_LEVEL>1)_DBG_<<"seed "<<i<<" has only SL1 hits. Looking for stray hits to match it with ..."<<endl;
		
		// This seed must have only SL1 hits, look for unused SL3 hits
		// at a phi angle consistent with the outermost hit.
		bool added_hit = false;
		const DCDCWire *sl1_wire = hits[hits.size()-1]->hit->wire;
		for(unsigned int j=0; j<cdctrkhits.size(); j++){
			DCDCTrkHit *sl3_hit = cdctrkhits[j];
			if(sl3_hit->flags&USED)continue;
			if(sl3_hit->hit->wire->ring<=superlayer_boundaries[1])continue;
			if(sl3_hit->hit->wire->ring>superlayer_boundaries[2])continue;
			double a = sl1_wire->origin.Angle(sl3_hit->hit->wire->origin);
			if(fabs(a)<MAX_CDC_MATCH_ANGLE/57.3){
				sl3_hit->flags |= USED;
				hits.push_back(sl3_hit);
				added_hit = true;
				if(DEBUG_LEVEL>1)_DBG_<<"Adding stray SL3 hit to SL1 seed a="<<a*57.3<<"  flags="<<sl3_hit->flags<<"  ring="<<sl3_hit->hit->wire->ring<<endl;
			}
		}

		if(added_hit){
			if(DEBUG_LEVEL>1)_DBG_<<"  found SL3 CDC hit to add to seed "<<i<<endl;
			continue; // We found at least one SL3 hit that was added to the seed
		}
		
		if(DEBUG_LEVEL>1)_DBG_<<"  looking for DFDCPseudo hit to add to seed "<<i<<endl;

		// Tracks at around 17 degrees can miss SL3 (or at least
		// inefficiency can cause the 1 or 2 hits we would get there
		// to be lost). In principle, the first package of the FDC
		// can get enough hits to form a candidate that is then matched
		// either matched with the CDC candidate in DTrackCandidate_factory.cc
		// or extends close enough to the real track in the CDC to be
		// used outright.
		//
		// There are a few tracks though that will not have enough hits in
		// package 1 of the FDC to form a candidate. For these, we need to
		// pick up any lone hits in the FDC and try and add them to 
		// our SL1 seed in order to get a better candiate.
		//
		// Note that this is being done while looking at single track
		// events with no noise. In the presence of noice hits, using
		// single hits like this can be a problem.
		vector<const DFDCPseudo*> fdchits;
		eventLoop->Get(fdchits);
		vector<const DFDCPseudo*> fdc1_hits_matched;
		for(unsigned int j=0; j<fdchits.size(); j++){
			const DFDCPseudo *fdchit = fdchits[j];
			if(fdchit->wire->layer>6)continue; // only interested in 1st package
			DVector2 sl1(sl1_wire->origin.X(), sl1_wire->origin.Y());
			DVector2 fdc1=fdchit->xy;
			double a = sl1.DeltaPhi(fdc1);

			if(fabs(a)<MAX_FDC_MATCH_ANGLE/57.3){
				fdc1_hits_matched.push_back(fdchit);
				if(DEBUG_LEVEL>3)_DBG_<<"Potential match of FDC to CDC candidate delta phi="<<a*57.3<<endl;
			}
		}
		
		// If any matches were found, add them to the seed.
		// Note: We may want to put a limit here so this is only done if
		// too few hits are found for the FDC segment factory
		// to form a segment. That will be left as a future optimization.
		if(fdc1_hits_matched.size()>0){
			if(DEBUG_LEVEL>1)_DBG_<<"Matched "<<fdc1_hits_matched.size()<<" hits from FDC package 1"<<endl;
			seeds[i].fdchits = fdc1_hits_matched;
		}else{
			if(DEBUG_LEVEL>1)_DBG_<<"  no matching stray hits found to add to seed "<<i<<endl;
		}
	}
}

//------------------
// DropIncompleteSeeds
//------------------
void DTrackCandidate_factory_CDC::DropIncompleteSeeds(vector<DCDCSeed> &seeds)
{
	/// Look for seeds that have hits only in SL3 or only in SL5 and drop them.

	int iSL1_lo = 1;
	int iSL1_hi = superlayer_boundaries[0];
	int iSL3_lo = superlayer_boundaries[1]+1;
	int iSL3_hi = superlayer_boundaries[2];
	int iSL5_lo = superlayer_boundaries[3]+1;
	int iSL5_hi = superlayer_boundaries[3];

	// Loop through seeds, looking for ones with only SL1 hits
	for(unsigned int i=0; i<seeds.size(); i++){
		if(!seeds[i].valid)continue;
		vector<DCDCTrkHit*> &hits = seeds[i].hits;
		if(hits.size()<1)continue;
		bool has_SL1_hit = false;
		bool has_SL3_hit = false;
		bool has_SL5_hit = false;
		for(unsigned int j=0; j<hits.size(); j++){
			int ring = hits[j]->hit->wire->ring;
			if(ring<iSL1_lo || ring>iSL1_hi)has_SL1_hit = true;
			if(ring<iSL3_lo || ring>iSL3_hi)has_SL3_hit = true;
			if(ring<iSL5_lo || ring>iSL5_hi)has_SL5_hit = true;
			if(has_SL1_hit)break;
		}
		if(has_SL1_hit)continue;
		if(has_SL3_hit && has_SL5_hit)continue;
		
		// Seed contains hits only from superlayer 3 or only from superlayer 5
		// Flag the seed as invalid so it is ignored later
		seeds[i].valid = false;
		if(DEBUG_LEVEL>1)_DBG_<<"Dropping seed "<<i<<": hits in oly SL3 or only SL5 (has_SL3_hit="<<has_SL3_hit<<", has_SL5_hit="<<has_SL5_hit<<")"<<endl;
	}
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
			if(DEBUG_LEVEL>3)_DBG_<<"asym="<<asym<<endl;
			if(asym<0.05) are_clones = true;

			if(!are_clones){
				// Find total absolute distance from alpha1 and seed2 circle and
				// alpha2 and seed2 circle.
				DVector2 d1 = dr+alpha1;
				DVector2 d2 = dr+alpha2;
				double d = fabs(d1.Mod()-r2.Mod()) + fabs(d2.Mod() - r2.Mod());
				if(d<1.5) are_clones = true;
				if(DEBUG_LEVEL>3)_DBG_<<"d="<<d<<endl;
			}
			// Check that the charges are the same
			if (are_clones){
			  if (seed1.fit.q!=seed2.fit.q) are_clones=false;
			  
			}
			// Remove a clone in necessary
			if(are_clones){
				// These seeds appear to be clones of one another. Mark one as not valid.
			  
				if(seed1.fit.chisq<=seed2.fit.chisq){
					seed2.valid = false;
				}else{
					seed1.valid = false;
				}
			  
				if(DEBUG_LEVEL>3)_DBG_<<"Filtering clone circle seed (seed1.fit.chisq="<<seed1.fit.chisq<<" seed2.fit.chisq="<<seed2.fit.chisq<<endl;
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
			if(DEBUG_LEVEL>3)_DBG_<<"Out of time stereo hit: seed.tdrift_avg="<<seed.tdrift_avg<<" tdrift="<<trkhit->hit->tdrift<<endl;
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

		if(DEBUG_LEVEL>15)_DBG_<<"alpha1="<<alpha1<<" alpha2="<<alpha2<<endl;
		
		// At this point we must decide which value of alpha to use. The
		// proper way would likely involve either trying both possibilities
		// and taking the one that gave the better chi-sq for a line fit
		// of phi vs. z or looking at the surrounding axial layers
		// and using the value which puts the hit closest to those.
		// For now, we just use the value closest to zero (i.e. closest to
		// the center of the wire).
		//double alpha = sqrt(r2_mod2)*(fabs(alpha1)<fabs(alpha2) ? alpha1:alpha2);
		
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
		//double s = alpha/fabs(sin(wire->stereo));
		// The above code unnecessarily multiples and divides by sin(theta_wire)
	        double s=(fabs(alpha1)<fabs(alpha2) ? alpha1:alpha2);

		//if(DEBUG_LEVEL>3)_DBG_<<"alpha="<<alpha<<" s="<<s<<" ring="<<wire->ring<<" straw="<<wire->straw<<" stereo="<<wire->stereo<<endl;
		if(DEBUG_LEVEL>15)_DBG_<<"s="<<s<<" ring="<<wire->ring<<" straw="<<wire->straw<<" stereo="<<wire->stereo<<endl;
		if(fabs(s) > wire->L/2.0)continue; // if wire doesn't cross circle, skip hit
		
	
		// Determine 3d space point for this hit
		DCDCTrkHit mytrkhit = *trkhit; // we need to make a copy since x_stereo,... are unique based on the cadidate
		DVector3 pos = wire->origin + s*wire->udir;
		mytrkhit.x_stereo = pos.X();
		mytrkhit.y_stereo = pos.Y();
		mytrkhit.z_stereo = pos.Z();
		
		// Verify that this hit is reasonably close to the axial hits in X and Y
		// If not, then drop this hit for this seed.
		double phi_diff = fabs(pos.Phi() - seed.phi_avg);
		if(phi_diff>M_PI)phi_diff = 2.0*M_PI - phi_diff;
		phi_diff*=57.3; // convert to degrees
		if(fabs(phi_diff)>MAX_SEED_LINK_ANGLE){ // yes, the fabs is needed ...
			if(DEBUG_LEVEL>4)_DBG_<<"rejecting stereo hit at phi="<<pos.Phi()<<" for being too far away from axial hits(phi_diff="<<phi_diff<<")"<<endl;
			continue;
		}
		mytrkhit.phi_stereo = atan2(mytrkhit.y_stereo-R.Y(), mytrkhit.x_stereo-R.X());
		R*=-1.0; // make R point from center of circle to beamline instead of other way around
		if(DEBUG_LEVEL>15){
			_DBG_<<" --- wire->udir X, Y, Z = "<<wire->udir.X()<<", "<<wire->udir.Y()<<", "<<wire->udir.Z()<<endl;
			_DBG_<<" -- ring="<<wire->ring<<" trkhit->z_stereo="<<trkhit->z_stereo<<" trkhit->y_stereo="<<trkhit->y_stereo<<" trkhit->x_stereo="<<trkhit->x_stereo<<endl;
			_DBG_<<" -- phi_stereo="<<trkhit->phi_stereo<<"  R.Phi()="<<R.Phi()<<"  (X,Y)=("<<R.X()<<", "<<R.Y()<<")"<<endl;
			_DBG__;
		}
		mytrkhit.phi_stereo -= R.Phi(); // make angle relative to beamline
		
		// We want this to go either from 0 to +2pi for positive charge, or
		// 0 to -2pi for negative.
		double phi_hi = seed.fit.q>0.0 ? +2.0*M_PI:0.0;
		double phi_lo = seed.fit.q>0.0 ? 0.0:-2.0*M_PI;
		while(mytrkhit.phi_stereo<phi_lo){
		  mytrkhit.phi_stereo+=2.0*M_PI;
		}
		while(mytrkhit.phi_stereo>phi_hi){
		  mytrkhit.phi_stereo-=2.0*M_PI;
		}
		mytrkhit.flags |= VALID_STEREO;
		seed.stereo_hits.push_back(mytrkhit);
		if(DEBUG_LEVEL>10)_DBG_<<"Adding CDC stereo hit: ring="<<mytrkhit.hit->wire->ring<<" straw="<<mytrkhit.hit->wire->straw<<endl;
	}
	if(DEBUG_LEVEL>5)_DBG_<<"Num stereo hits: "<<seed.stereo_hits.size()<<endl;
}

//------------------
// FindThetaZ
//------------------
void DTrackCandidate_factory_CDC::FindThetaZ(DCDCSeed &seed)
{
	// Decide whether to use the histogramming method or the
	// straight track method base on the transverse momentum
if(DEBUG_LEVEL>2)_DBG_<<"p_trans:"<<seed.fit.p_trans<<endl;
//if(seed.fit.p_trans<3.0)
 if (true)
	  {
		FindTheta(seed, TARGET_Z_MIN, TARGET_Z_MAX);
		FindZ(seed, seed.theta_min, seed.theta_max);
	}else{
		FindThetaZStraightTrack(seed);
		
		// The value of p_trans that was originally determined
		// by DQuickFit::FitCircleStraightTrack() can often lead
		// to a total momentum larger than the beam momentum.
		// Re-estimate p_trans limiting the search with the
		// theta we just found.
		seed.fit.SearchPtrans(9.0*sin(seed.theta), 0.1);
		//seed.fit.QuickPtrans();
		if(DEBUG_LEVEL>2)_DBG_<<"p_trans from search:"<<seed.fit.p_trans<<"  (ptot="<<seed.fit.p_trans/sin(seed.theta)<<")"<<endl;
	}
	
	// If z_vertex is not inside the target limits, then flag this
	// seed as invalid.
	if(seed.z_vertex<TARGET_Z_MIN || seed.z_vertex>TARGET_Z_MAX){
		if(DEBUG_LEVEL>3)_DBG_<<"Seed z-vertex outside of target range (z="<<seed.z_vertex<<" TARGET_Z_MIN="<<TARGET_Z_MIN<<" TARGET_Z_MAX="<<TARGET_Z_MAX<<endl;
		seed.valid=false;
	}
	
	return;
}

//------------------
// FindThetaZStraightTrack
//------------------
void DTrackCandidate_factory_CDC::FindThetaZStraightTrack(DCDCSeed &seed)
{
	/// In the case of high momentum tracks, the values of phi_stereo will
	/// all be close together (and near 0 or +-2pi). For these cases,
	/// we calculate theta from r and z points which should be pretty linear.
	
	if(DEBUG_LEVEL>3)_DBG_<<"Fitting theta and z for high transverse momentum track (p_trans="<<seed.fit.p_trans<<")"<<endl;

	// Do a linear regression of R vs. Z for stereo hits. First, we make
	// a list of the hits we'll use.
	vector<double> r;
	vector<double> z;
	for(unsigned int i=0; i<seed.stereo_hits.size(); i++){
		DCDCTrkHit *trkhit = &seed.stereo_hits[i];
		if(!trkhit->flags&VALID_STEREO)continue;
		
		//double R = sqrt(pow(trkhit->x_stereo,2.0) + pow(trkhit->y_stereo,2.0));
		double R=sqrt(trkhit->x_stereo*trkhit->x_stereo
			      +trkhit->y_stereo*trkhit->y_stereo);
		r.push_back(R);
		z.push_back(trkhit->z_stereo);
	}

	for(unsigned int i=0; i<seed.fdchits.size(); i++){
		const DFDCPseudo *fdchit = seed.fdchits[i];
		//double R = sqrt(pow((double)fdchit->x,2.0) + pow((double)fdchit->y,2.0));
		//double R=sqrt(fdchit->x*fdchit->x+fdchit->y*fdchit->y);
		double R=fdchit->xy.Mod();
		r.push_back(R);
		z.push_back(fdchit->wire->origin.Z());
	}

	// Add center of target as a point to constrain the fit a little more
	double z_target = (TARGET_Z_MIN+TARGET_Z_MAX)/2.0;
	r.push_back(0.0);
	z.push_back(z_target);
	
	// Calculate average z and r
	double Ravg=0.0, Zavg=0.0;
	for(unsigned int i=0; i<r.size(); i++){
		Ravg += r[i];
		Zavg += z[i];
	}
	Ravg /= (double)r.size();
	Zavg /= (double)z.size();

	// Now, do the regression
	double Szr=0.0, Szz=0.0, Srr=0.0;;
	for(unsigned int i=0; i<r.size(); i++){
		if(DEBUG_LEVEL>4)_DBG_<<"r="<<r[i]<<"\t z="<<z[i]<<"  Szz_i="<<pow((z[i] - Zavg), 2.0)<<" Srr_i="<<pow((r[i] - Ravg),2.0)<<endl;
		double dz=z[i] - Zavg;
		double dr=r[i] - Ravg;
		/*
		Szr += (z[i] - Zavg)*(r[i] - Ravg);
		Szz += pow((z[i] - Zavg), 2);
		Srr += pow((r[i] - Ravg), 2);
		*/
		Szr += dz*dr;
		Szz += dz*dz;
		Srr += dr*dr;
	}

	if(Szz>Srr){
		// track is more horizontal
		seed.z_vertex = Zavg - Ravg*Szz/Szr;
		seed.theta = atan2(Szr, Szz);
	}else{
		// track is more vertical
		seed.z_vertex = Zavg - Ravg*Szr/Srr;
		seed.theta = M_PI_2 - atan2(Szr, Srr);
	}
	if(seed.theta<0.0){
		if(DEBUG_LEVEL>3)_DBG_<<"theta less than 0 ("<<seed.theta<<") this simply shoudn't happen!"<<endl;
		seed.theta+=M_PI;
	}
	if(DEBUG_LEVEL>3)_DBG_<<"theta="<<seed.theta<<"  z_vertex="<<seed.z_vertex<<endl;
	if(DEBUG_LEVEL>4){
		_DBG_<<"Szz="<<Szz<<"  Srr="<<Srr<<endl;
		double z_vertex = Zavg - Ravg*Szz/Szr;
		double theta = atan2(Szr, Szz);
		_DBG_<<"for Szz>Srr : theta="<<theta<<"  z_vertex="<<z_vertex<<endl;
		z_vertex = Zavg - Ravg*Szr/Srr;
		theta = M_PI_2 - atan2(Szr, Srr);
		_DBG_<<"for Szz<=Srr : theta="<<theta<<"  z_vertex="<<z_vertex<<endl;
	}
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
	
	// Loop over CDC hits, filling the histogram
	double r0 = seed.fit.r0;
	for(unsigned int i=0; i<seed.stereo_hits.size(); i++){
		DCDCTrkHit *trkhit = &seed.stereo_hits[i];
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
		if(DEBUG_LEVEL>3)_DBG_<<" -- phi_stereo="<<trkhit->phi_stereo<<" z_stereo="<<trkhit->z_stereo<<"  alpha="<<alpha<<endl;
		if(DEBUG_LEVEL>3)_DBG_<<" -- tmin="<<tmin<<"  tmax="<<tmax<<endl;
		
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
	
	// Loop over FDC hits, filling the histogram
	for(unsigned int i=0; i<seed.fdchits.size(); i++){
		const DFDCPseudo *fdchit = seed.fdchits[i];
		
		// Calculate upper and lower limits in theta
		DVector2 R(seed.fit.x0, seed.fit.y0);
		double phi_stereo = atan2(fdchit->xy.Y()-R.Y(), 
					  fdchit->xy.X()-R.X());
		double alpha = -r0*phi_stereo;
		if(seed.fit.q<0.0)alpha = -alpha;
		double tmin = atan2(alpha, fdchit->wire->origin.Z() - target_z_min);
		double tmax = atan2(alpha, fdchit->wire->origin.Z() - target_z_max);
		if(tmin>tmax){
			double tmp = tmin;
			tmin=tmax;
			tmax=tmp;
		}
		if(DEBUG_LEVEL>3)_DBG_<<" -- phi_stereo="<<phi_stereo<<" z_stereo="<<fdchit->wire->origin.Z()<<"  alpha="<<alpha<<endl;
		if(DEBUG_LEVEL>3)_DBG_<<" -- tmin="<<tmin<<"  tmax="<<tmax<<endl;
		
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
			if(DEBUG_LEVEL>3)_DBG_<<" -- istart="<<istart<<" (theta="<<hist_low_limit + bin_width*(0.5+(double)istart)<<" , N="<<hist[i]<<")"<<endl;
		}
		if(hist[i] == hist[istart])iend = i;
	}
	
	// If there are no entries in the histogram, then flag this seed as invalid
	if(hist[istart]==0.0)seed.valid=false;
	
	// Calculate theta limits
	seed.theta_min = hist_low_limit + bin_width*(0.5+(double)istart);
	seed.theta_max = hist_low_limit + bin_width*(0.5+(double)iend);
	seed.theta = (seed.theta_max + seed.theta_min)/2.0;
	if(DEBUG_LEVEL>3)_DBG_<<"istart="<<istart<<" iend="<<iend<<" theta_min="<<seed.theta_min<<" theta_max="<<seed.theta_max<<endl;
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
	
	// Loop over CDC hits, filling the histogram
	double r0 = seed.fit.r0;
	double tan_alpha_min = tan(theta_min)/r0;
	double tan_alpha_max = tan(theta_max)/r0;
	for(unsigned int i=0; i<seed.stereo_hits.size(); i++){
		DCDCTrkHit *trkhit = &seed.stereo_hits[i];
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
		if(DEBUG_LEVEL>3)_DBG_<<" -- phi_stereo="<<trkhit->phi_stereo<<" z_stereo="<<trkhit->z_stereo<<endl;
		if(DEBUG_LEVEL>3)_DBG_<<" -- zmin="<<zmin<<"  zmax="<<zmax<<endl;
		
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

	// Loop over FDC hits, filling the histogram
	for(unsigned int i=0; i<seed.fdchits.size(); i++){
		const DFDCPseudo *fdchit = seed.fdchits[i];
		
		// Calculate upper and lower limits in z
		DVector2 R(seed.fit.x0, seed.fit.y0);
		double phi_stereo = atan2(fdchit->xy.Y()-R.Y(), 
					  fdchit->xy.X()-R.X());
		double z = fdchit->wire->origin.Z();
		double q_sign = seed.fit.q>0.0 ? +1.0:-1.0;
		double zmin = z - q_sign*phi_stereo/tan_alpha_min;
		double zmax = z - q_sign*phi_stereo/tan_alpha_max;
		if(zmin>zmax){
			double tmp = zmin;
			zmin=zmax;
			zmax=tmp;
		}
		if(DEBUG_LEVEL>3)_DBG_<<" -- phi_stereo="<<phi_stereo<<" z_stereo="<<z<<endl;
		if(DEBUG_LEVEL>3)_DBG_<<" -- zmin="<<zmin<<"  zmax="<<zmax<<endl;
		
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
			if(DEBUG_LEVEL>3)_DBG_<<" -- istart="<<istart<<" (z="<<hist_low_limit + bin_width*(0.5+(double)istart)<<" , N="<<hist[i]<<")"<<endl;
		}
		if(hist[i] == hist[istart])iend = i;
	}

	// If there are no entries in the histogram, then flag this seed as invalid
	if(hist[istart]==0.0)seed.valid=false;
	
	// Calculate z limits
	seed.z_min = hist_low_limit + bin_width*(0.5+(double)istart);
	seed.z_max = hist_low_limit + bin_width*(0.5+(double)iend);
	seed.z_vertex = (seed.z_max + seed.z_min)/2.0;
	if(DEBUG_LEVEL>3)_DBG_<<"istart="<<istart<<" iend="<<iend<<" z_min="<<seed.z_min<<" z_max="<<seed.z_max<<" hits[istart]="<<hist[istart]<<endl;
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
// DCDCSeed::DCDCSeed
//------------------
DTrackCandidate_factory_CDC::DCDCSeed::DCDCSeed()
{	
	phi_avg = 0.0;
	tdrift_avg = 0.0;
	linked = false;
	valid = false;
	theta = 0.0;
	z_vertex = 0.0;
	q = 0.0;
	theta_min = 0.0;
	theta_max = 0.0;
	z_min = 0.0;
	z_max = 0.0;
}

//------------------
// DCDCSeed::Merge
//------------------
void DTrackCandidate_factory_CDC::DCDCSeed::Merge(DCDCSeed& seed)
{	
	// Need avg. phi based on hits from both seeds
	double x = (double)hits.size()*cos(phi_avg);
	double y = (double)hits.size()*sin(phi_avg);

	for(unsigned int i=0; i<seed.hits.size(); i++){
		hits.push_back(seed.hits[i]);
		
		x += cos(seed.hits[i]->hit->wire->phi);
		y += sin(seed.hits[i]->hit->wire->phi);
	}

	phi_avg = atan2(y,x);
	if(phi_avg<0.0)phi_avg += 2.0*M_PI;
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
// DCDCSeed::FindAverageBz
//------------------
double DTrackCandidate_factory_CDC::DCDCSeed::FindAverageBz(JEventLoop *loop)
{
  //return 2.0;
	if(!loop)return 0.0;
	DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp)return 0.0;

	DMagneticFieldMap *bfield = dapp->GetBfield();
	if(!bfield)return 0.0;

	double Bz_sum=0.0;
	for(unsigned int i=0; i<stereo_hits.size(); i++){
	  DCDCTrkHit *hit = &stereo_hits[i];
	  Bz_sum += bfield->GetBz(hit->x_stereo,hit->y_stereo,
				  hit->z_stereo);
	}
	// If there are no stereo hits, fall back on available FDC hits
	for (unsigned int i=0; i<fdchits.size();i++){
	  const DFDCPseudo *hit= fdchits[i];
	  Bz_sum += bfield->GetBz(hit->xy.X(), hit->xy.Y(), 
				  hit->wire->origin.z());
	}
	
	return Bz_sum/(double)(stereo_hits.size()+fdchits.size());
}


//------------------
// FindThetaZRegression
//------------------
// Linear regression to find tan(lambda) and z_vertex
jerror_t DTrackCandidate_factory_CDC::FindThetaZRegression(DCDCSeed &seed){

	if(DEBUG_LEVEL>3)_DBG_<<"Finding theta and z via linear regression method."<<endl;

  if (seed.fit.normal.Mag()==0.) return VALUE_OUT_OF_RANGE;
  // Vector of intersections between the circles of the measurements and the plane intersecting the Riemann surface
  vector<DVector3_with_perp>intersections;
  DVector3_with_perp beam(0,0,65.);
  intersections.push_back(beam);
  
  // CDC stereo hits
  for (unsigned int m=0;m<seed.stereo_hits.size();m++){
    DCDCTrkHit *trkhit=&seed.stereo_hits[m];
    double R2=trkhit->x_stereo*trkhit->x_stereo+trkhit->y_stereo*trkhit->y_stereo;

    DVector3_with_perp intersection;
    DVector3 N=seed.fit.normal;
    //double c0=seed.fit.c_origin;
    double A=seed.fit.c_origin+R2*N.Z();
    double B=N.Perp();
    double C=B*R2-A*A;
    
    if (C>=0) {
      double sqrtC=sqrt(C);
      double x1=(-N.X()*A+N.Y()*sqrtC)/B;
      double y1=(-N.Y()*A-N.X()*sqrtC)/B;   
      double x2=(-N.X()*A-N.Y()*sqrtC)/B;
      double y2=(-N.Y()*A+N.X()*sqrtC)/B;
      
      if (fabs(trkhit->y_stereo-y1)<fabs(trkhit->y_stereo-y2)){
	intersection.SetX(x1);
	intersection.SetY(y1);
      }
      else{
	intersection.SetX(x2);
	intersection.SetY(y2);
      }
      intersection.SetZ(trkhit->z_stereo);
      
      intersections.push_back(intersection);
		if(DEBUG_LEVEL>5)_DBG_<<"Adding CDC hit "<<m<<" z="<<intersection.z()<<endl;      
      
    }

  }

  // FDC pseudo hits
  for (unsigned int m=0;m<seed.fdchits.size();m++){
    const DFDCPseudo *trkhit=seed.fdchits[m];
    //double R2=trkhit->x*trkhit->x+trkhit->y*trkhit->y;
    double R2=trkhit->xy.Mod2();
    DVector3_with_perp intersection;
    DVector3 N=seed.fit.normal;
    //double c0=seed.fit.c_origin;
    double A=seed.fit.c_origin+R2*N.Z();
    double B=N.Perp();
    double C=B*R2-A*A;
    
    if (C>=0) {
      double sqrtC=sqrt(C);
      double x1=(-N.X()*A+N.Y()*sqrtC)/B;
      double y1=(-N.Y()*A-N.X()*sqrtC)/B;   
      double x2=(-N.X()*A-N.Y()*sqrtC)/B;
      double y2=(-N.Y()*A+N.X()*sqrtC)/B;
      
      if (fabs(trkhit->xy.Y()-y1)<fabs(trkhit->xy.Y()-y2)){
	intersection.SetX(x1);
	intersection.SetY(y1);
      }
      else{
	intersection.SetX(x2);
	intersection.SetY(y2);
      }
      intersection.SetZ(trkhit->wire->origin.Z());
      
      intersections.push_back(intersection);
		if(DEBUG_LEVEL>5)_DBG_<<"Adding FDC hit "<<m<<" z="<<intersection.z()<<endl;      
    }

  }
  
  // The SortIntersections method requires the perp field to be valid
  // for all objects. Go through and fill in those fields now.
  for(unsigned int i=0; i<intersections.size(); i++)intersections[i].CalcPerp();
  
  // Now, sort the entries
  sort(intersections.begin(),intersections.end(),SortIntersections);

  // Compute the arc lengths between (0,0) and (xi,yi)
  vector<double>arclengths(intersections.size());
  vector<double>var_z(intersections.size());
  arclengths[0]=0.;
  //double temp=1.6/sin(6.*M_PI/180.);
  //double varz=temp*temp/12.;
  for (unsigned int m=1;m<arclengths.size();m++){
    double diffx=intersections[m].x()-intersections[0].x();
    double diffy=intersections[m].y()-intersections[0].y();
    double chord=sqrt(diffx*diffx+diffy*diffy);
    double s=2*seed.fit.r0*asin(chord/2./seed.fit.r0);
    
    arclengths[m]=s;
    //var_z[m]=varz;
    var_z[m]=19.52;
  }

  //Linear regression to find z0, tanl
  double tanl=0.,z0=0.;
  if (arclengths.size()>2){ // Do fit only if have more than one measurement
    //var_z[0]=30.*30./12.;  // Use length of the target
    var_z[0]=75.;
    double sumv=0.,sumx=0.;
    double sumy=0.,sumxx=0.,sumxy=0.;
    for (unsigned int m=0;m<intersections.size();m++){
      double temp=1./var_z[m];
      sumv+=temp;
      sumx+=arclengths[m]*temp;
      sumy+=intersections[m].z()*temp;;
      sumxx+=arclengths[m]*arclengths[m]*temp;
      sumxy+=arclengths[m]*intersections[m].z()*temp;
    }
    double Delta=sumv*sumxx-sumx*sumx;
    if (Delta==0.) return VALUE_OUT_OF_RANGE;
    
    tanl=(sumv*sumxy-sumx*sumy)/Delta;
    z0=(sumxx*sumy-sumx*sumxy)/Delta;
  }
  else if(arclengths.size()==2){
    z0=Z_TARGET;
    tanl=(intersections[1].z()-z0)/arclengths[1];
  }else{
	if(DEBUG_LEVEL>5)_DBG_<<"Fit failed for theta-z via regressionz due to too few hits with z-info"<<endl;
	return VALUE_OUT_OF_RANGE;
  }
  
  if (z0>TARGET_Z_MAX || z0<TARGET_Z_MIN){
	if(DEBUG_LEVEL>5)_DBG_<<"Fit failed for theta-z via regressionz value out of target range (z="<<z0<<")"<<endl;
	return VALUE_OUT_OF_RANGE;
  }

  seed.theta=M_PI/2-atan(tanl);
  seed.z_vertex=z0;

  return NOERROR;
}

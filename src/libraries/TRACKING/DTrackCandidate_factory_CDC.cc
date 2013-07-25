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

#define BeamRMS 1.0
#define EPS 1e-3
#ifndef M_TWO_PI
#define M_TWO_PI 6.28318530717958647692
#endif

#include "DTrackCandidate_factory_CDC.h"

#define TWO(c)     (0x1u << (c))
#define MASK(c) ((unsigned int)(-1)) / (TWO(TWO(c)) + 1u)
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (TWO(c))) & MASK(c))
 
inline int bitcount (unsigned int n)  {
  n = COUNT(n, 0) ;
  n = COUNT(n, 1) ;
  n = COUNT(n, 2) ;
  n = COUNT(n, 3) ;
  n = COUNT(n, 4) ;
  /*n = COUNT(n, 5) ;    for 64-bit integers */
 return n ;
}


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

inline bool SortIntersections(const intersection_t &a,const intersection_t &b){
  if (a.perp2<b.perp2) return true;
  return false;
}
inline bool CDCSortByRdecreasing2(const DTrackCandidate_factory_CDC::DCDCTrkHit &hit1,
				 const DTrackCandidate_factory_CDC::DCDCTrkHit &hit2){
  // use the ring number to sort by R(decreasing) and then straw(increasing)
  if(hit1.hit->wire->ring == hit2.hit->wire->ring){
    return hit1.hit->wire->straw < hit2.hit->wire->straw;
  }
  return hit1.hit->wire->ring > hit2.hit->wire->ring;
}


inline bool CDCSortByRdecreasing(DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit1, DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit2) {
	// use the ring number to sort by R(decreasing) and then straw(increasing)
	if(hit1->hit->wire->ring == hit2->hit->wire->ring){
		return hit1->hit->wire->straw < hit2->hit->wire->straw;
	}
	return hit1->hit->wire->ring > hit2->hit->wire->ring;
}
inline bool CDCSortByRincreasing(DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit1, DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit2) {
	// use the ring number to sort by R
	return hit1->hit->wire->ring < hit2->hit->wire->ring;
}
inline bool CDCSortByZincreasing(DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit1, DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit2) {
	// use the z_stereo position to sort
	return hit1->z_stereo < hit2->z_stereo;
}
inline bool CDCSortByStereoPhiincreasing(DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit1, DTrackCandidate_factory_CDC::DCDCTrkHit* const &hit2) {
	// use the z_stereo position to sort
	return hit1->phi_stereo < hit2->phi_stereo;
}

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_CDC::init(void)
{
	MAX_ALLOWED_CDC_HITS = 10000;
	MAX_SUBSEED_STRAW_DIFF = 1;
	MIN_SEED_HITS  = 2;
	MAX_SUBSEED_LINKED_HITS = 12;
	MAX_RING_SUBSEED_HITS = 4;
	MAX_HIT_DIST = 4.0; // cm
	MAX_SEED_TIME_DIFF = 1000.0; // ns
	MAX_CDC_MATCH_ANGLE = 10.0; // degrees

	MAX_SEED_LINK_ANGLE = 10.0; // degrees
	TARGET_Z_MIN = 50.0;
	TARGET_Z_MAX = 80.0;
	TARGET_Z=65.0;
	DEBUG_LEVEL = 0;

	VERTEX_Z_MIN=-100.0;
	VERTEX_Z_MAX=+200.0;

	FILTER_SEEDS=true;
	
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
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();

  // get the geometry
  const DGeometry *geom = dapp->GetDGeometry(runnumber);

  geom->GetTargetZ(TARGET_Z);
  double zrange;
  geom->GetTargetLength(zrange);
  TARGET_Z_MIN=TARGET_Z-0.5*zrange;
  TARGET_Z_MAX=TARGET_Z+0.5*zrange;
  

	gPARMS->SetDefaultParameter("TRKFIND:MAX_ALLOWED_CDC_HITS", MAX_ALLOWED_CDC_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SUBSEED_STRAW_DIFF", MAX_SUBSEED_STRAW_DIFF);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_SEED_HITS", MIN_SEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SUBSEED_LINKED_HITS", MAX_SUBSEED_LINKED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_RING_SUBSEED_HITS", MAX_RING_SUBSEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_HIT_DIST", MAX_HIT_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_TIME_DIFF", MAX_SEED_TIME_DIFF);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_CDC_MATCH_ANGLE", MAX_CDC_MATCH_ANGLE);

	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_LINK_ANGLE", MAX_SEED_LINK_ANGLE);
	gPARMS->SetDefaultParameter("TRKFIND:DEBUG_LEVEL", DEBUG_LEVEL);
	gPARMS->SetDefaultParameter("TRKFIND:FILTER_SEEDS",FILTER_SEEDS);
	
	MAX_HIT_DIST2 = MAX_HIT_DIST*MAX_HIT_DIST;
	
	DEBUG_HISTS=false;
	gPARMS->SetDefaultParameter("TRKFIND:DEBUG_HISTS", DEBUG_HISTS);

	if(DEBUG_HISTS) {
	  dapp->Lock();
	  Hdphi_s=(TH2F*)gROOT->FindObject("Hdphi_s");
	  if (!Hdphi_s) 
	    Hdphi_s=new TH2F("Hdphi_s",
			     "Matching #delta#phi for attaching stereo hits",
			      200,0.,1000,100,0,50.);
	  Hdphi_s=(TH2F*)gROOT->FindObject("Hdphi_s");
	  if (!Hdphi) 
	    Hdphi=new TH2F("Hdphi",
			     "Matching #delta#phi for linking seeds",
			      200,0.,1000,100,0,50.);

	  dapp->Unlock();
	}
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
  	for(unsigned int i=0; i<cdctrkhits.size(); i++)delete cdctrkhits[i];
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_CDC::evnt(JEventLoop *loop, int eventnumber)
{
  // Get CDC hits
  unsigned int numcdchits=0;
  if (GetCDCHits(loop,numcdchits)!=NOERROR){
    return RESOURCE_UNAVAILABLE;
  }
	
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
  vector<DCDCSeed> seeds_tmp, seeds_tmpSL5,seeds_tmpSL3_SL5,seeds;
  LinkSeeds(seeds_sl5, seeds_sl3, seeds_tmp, MAX_SUBSEED_LINKED_HITS);

  // Put aside the seeds that only have hits in SL5
  for (unsigned int i=0;i<seeds_tmp.size();i++){
    unsigned int max_ind=seeds_tmp[i].hits.size()-1;
    if (seeds_tmp[i].hits[max_ind]->hit->wire->ring>24){
      seeds_tmpSL5.push_back(seeds_tmp[i]);
    }
    else{
      seeds_tmpSL3_SL5.push_back(seeds_tmp[i]);
    }
  }
  
  // Link seeds between SL1 and SL3+SL5
  LinkSeeds(seeds_tmpSL3_SL5, seeds_sl1, seeds, 2*MAX_SUBSEED_LINKED_HITS);

  // Insert the SL5 seeds into the full list of seeds
  seeds.insert(seeds.begin(),seeds_tmpSL5.begin(),seeds_tmpSL5.end());

  // Check to add lone hits in SL3 to seeds with only SL1 hits
  PickupUnmatched(seeds);
  
  // Drop seeds containing hits from only a single axial superlayer besides SL1
  DropIncompleteSeeds(seeds);
  
  if (seeds.size()>0){
    // Create vectors to store bit pattern of hits used by the seeds
    for (unsigned int i=0;i<seeds.size();i++){
      seeds[i].HitBitPattern.clear();
      seeds[i].HitBitPattern.resize(numcdchits/(8*sizeof(unsigned int))+1);
    } 
    
    // Fit linked seeds to circles
    for(unsigned int i=0; i<seeds.size(); i++){
      if(DEBUG_LEVEL>5)_DBG_<<"-- Starting fit for seed "<<i<<" Nhits="<<seeds[i].hits.size()<<" phi_avg="<<seeds[i].phi_avg<<endl;
      seeds[i].valid = FitCircle(seeds[i]);
    }
      
    // Extend seeds into stereo layers and do fit to find z and theta
    for(unsigned int i=0; i<seeds.size(); i++){
      DCDCSeed &seed = seeds[i];
      if(DEBUG_LEVEL>3)_DBG_<<"----- Seed "<<i<<" ------"<<endl;
      if(DEBUG_LEVEL>3)_DBG_<<"seed.fit.phi="<<seed.fit.phi<<endl;
      if(!seed.valid)continue;
      
      // Add stereo hits to seed
      AddStereoHits(cdchits_by_superlayer[2-1], seed);

      // If the track would pass through the first stereo layers but there are
      // no intersections with the circle, mark seed as invalid.
      unsigned int first_stereo_count=seed.stereo_hits.size();
      if (first_stereo_count==0){
	seed.valid=false;
	continue;
      }
      
      // Find places where the stereo straws transition from + to - (or - to +)
      // angles, add the corresponding intersection points as points on the 
      // circle, and refit...
      if (first_stereo_count>0){
	RefitCircleWithStereoIntersections(seed);
      }
      
      // Go on to the next set of stereo layers
      AddStereoHits(cdchits_by_superlayer[4-1], seed);
  
      // If no stereo hits were found for this seed, then
      // we can't fit it.
      if(seed.stereo_hits.size()==0){
	seed.valid=false;
	continue;
      }
      
      // Find places where the stereo straws transition from + to - (or - to +)
      // angles, add the corresponding intersection points as points on the 
      // circle, and refit...
      if (seed.stereo_hits.size()>first_stereo_count){
	RefitCircleWithStereoIntersections(seed);
      }

      // Next do the line fit to determine tan(lambda)
      DoLineFit(seed);
   
      // Check for charge consistency
      seed.CheckCharge();
    }
  
    // Filter out duplicates of seeds by clearing their "valid" flags
    FilterCloneSeeds(seeds);

    // Try to gather stray stereo hits to add to existing seeds   
    //if (false)
    for(unsigned int i=0; i<seeds.size(); i++){
      DCDCSeed seed = seeds[i];	
      if (seed.valid){
	bool got_hits=false;
	AddStrayStereoHits(cdchits_by_superlayer[2-1],seed,got_hits);
	AddStrayStereoHits(cdchits_by_superlayer[4-1],seed,got_hits);

	if (got_hits){	 
	  // Sort stereo hits by R(decreasing)
	  sort(seed.stereo_hits.begin(),seed.stereo_hits.end(),
	       CDCSortByRdecreasing2);

	  // Find places where the stereo straws transition from + to - 
	  // (or - to +) angles, add the corresponding intersection points 
	  // as points on the circle, and refit...
	  if (RefitCircleWithStereoIntersections(seed)==true){
	    // Try to grab more stereo hits using the revised circle fit
	    AddStrayStereoHits(cdchits_by_superlayer[2-1],seed,got_hits);
	    AddStrayStereoHits(cdchits_by_superlayer[4-1],seed,got_hits);
	
	    // Next do the line fit to determine tan(lambda) 
	    DoLineFit(seed);

	    // Check for charge consistency
	    seed.CheckCharge();

	    // Replace the existing seed with the revised version
	    seeds[i]=seed;
	  } // circe refit succeeded?
	} // did we add any more hits to the seed?
      }
    }
  
    // Put the good seeds in the list of cdc track candidates
    for(unsigned int i=0; i<seeds.size(); i++){
      DCDCSeed seed = seeds[i];
      if (seed.valid){
	// The following is from a fit of ratio of thrown to reconstructed
	// transverse momentum vs. theta for the 1400A field
	//double par[] = {0.984463, 0.150759, -0.414933, 0.257472, -0.055801};
	//double theta = seed.theta;
	//double ff = par[0]+theta*(par[1]+theta*(par[2]+theta*(par[3]+theta*par[4])));
	
	DVector3 pos, mom;		
	GetPositionAndMomentum(seed,pos,mom);
	
	//Make a track candidate from results
	DTrackCandidate *can = new DTrackCandidate;
	
	can->setPosition(pos);
	can->setMomentum(mom);
	//can->setCharge(seed.q);
	can->setCharge(seed.fit.q);
	
	for (unsigned int n=0;n<seed.hits.size();n++){
	  const DCDCTrackHit *cdchit=(seed.hits[n])->hit;
	  can->used_cdc_indexes.push_back(seed.hits[n]->index);
	  can->AddAssociatedObject(cdchit);
	}
	for (unsigned int n=0;n<seed.stereo_hits.size();n++){
	  const DCDCTrackHit *cdchit=(seed.stereo_hits[n]).hit;
	  can->used_cdc_indexes.push_back((seed.stereo_hits[n]).index);
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
    }
  }
  
  return NOERROR;
}

//------------------
// GetCDCHits
//------------------
jerror_t DTrackCandidate_factory_CDC::GetCDCHits(JEventLoop *loop,
						 unsigned int &numhits)
{
	// Delete old hits
	for(unsigned int i=0; i<cdctrkhits.size(); i++)delete cdctrkhits[i];
	cdctrkhits.clear();
	for(unsigned int i=0; i<cdchits_by_superlayer.size(); i++)cdchits_by_superlayer[i].clear();

       
	// Get the "raw" hits. These already have the wire associated with them.
	vector<const DCDCTrackHit*> cdctrackhits;
	loop->Get(cdctrackhits);
	numhits=cdctrackhits.size();
	
	// If there are no hits, then bail now
	if(cdctrackhits.size()==0)return RESOURCE_UNAVAILABLE;
	
	// If there are too many hits, bail with a warning message
	if(cdctrackhits.size()>MAX_ALLOWED_CDC_HITS){
		_DBG_<<"Too many hits in CDC ("<<cdctrackhits.size()<<", max="<<MAX_ALLOWED_CDC_HITS<<")! Track finding in CDC bypassed for event "<<loop->GetJEvent().GetEventNumber()<<endl;
		cdctrackhits.clear();
		return UNRECOVERABLE_ERROR;
	}
	
	// Create DCDCTrkHit objects out of these.	
	int oldwire = -1;
	for(unsigned int i=0; i<cdctrackhits.size(); i++){	  
	  // Add to "master" list
	  // ONLY FIRST HIT OF A WIRE
	  int newwire = cdctrackhits[i]->wire->ring*1000 + cdctrackhits[i]->wire->straw;
	  if (newwire!=oldwire){

	    oldwire = newwire;
	  
	    // Add to "master" list
	    DCDCTrkHit *cdctrkhit = new DCDCTrkHit;
	    cdctrkhit->index=i;
	    cdctrkhit->hit = cdctrackhits[i];
	    cdctrkhit->flags = cdctrkhit->hit->is_stereo==false ? NONE:IS_STEREO;
	    cdctrkhit->flags |= NOISE; // (see below)

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
	return NOERROR;
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

		// Perform a preliminary circle fit for seed1
		DHelicalFit fit;
		for (unsigned int k=0;k<hits1.size();k++){
		  DVector3 pos=hits1[k]->hit->wire->origin;
		  fit.AddHitXYZ(pos.x(),pos.y(),pos.z());
		}
		fit.FitCircleRiemann(TARGET_Z,BeamRMS);
		
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
			// Skip if we would have to cross over into the opposite
			// quadrant to match
			if (pos2->x()*pos1->x()<0 && pos2->y()*pos1->y()<0.) continue;

			// Find the position on the circle from the first seed that
			// is closest to pos2
			double dx=pos2->x()-fit.x0;
			double dy=pos2->y()-fit.y0;
			double one_over_denom=1./sqrt(dx*dx+dy*dy);
			double x=fit.x0+fit.r0*dx*one_over_denom;
			double y=fit.y0+fit.r0*dy*one_over_denom;
			DVector2 mypos(x,y);
			
			// Compare phi values to see if the seeds are close enough
			// to link
			double my_dphi=fabs(mypos.Phi()-pos2->Phi());
			while(my_dphi>M_PI)my_dphi-=M_TWO_PI;
			my_dphi*=57.3;
			if (DEBUG_HISTS){
			  Hdphi->Fill(fit.r0,my_dphi);
			}
			double dphi_hits=pos2->Phi()-pos1->Phi();
			while(dphi_hits>M_PI)dphi_hits-=M_TWO_PI;
			dphi_hits=fabs(57.3*dphi_hits);

			if(fabs(my_dphi)<MAX_SEED_LINK_ANGLE && dphi_hits<90.0
			   ){
				DCDCSeed seed = in_seeds1[i];
				seed.Merge(in_seeds2[j]);
				seeds.push_back(seed);
				
				// Mark all of the hits in the new, merged seed as USED. Some will
				// already have the flag set by having been previously merged. The
				// USED flag is used later when matching SL1 only seeds with stray
				// hits in SL3.
				for(unsigned int k=0; k<seed.hits.size(); k++){
				  seed.hits[k]->flags |= USED;
				}
				
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
			  if(DEBUG_LEVEL>3)_DBG_<<"  rejected link between "<<i<<" and "<<j<<" due to (dphi="<< fabs(my_dphi) << ")>="<<MAX_SEED_LINK_ANGLE<<")"<<endl;
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
			if(tdrift>1000.0){
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
	  unsigned int numbits=8*sizeof(unsigned int);
	  seed.HitBitPattern[seed.hits[j]->index/numbits]
	    |=0x1u<<seed.hits[j]->index%numbits;
	
	  if(seed.hits[j]->flags&OUT_OF_TIME)continue;
	  
		const DVector3 &pos = seed.hits[j]->hit->wire->origin;
		seed.fit.AddHitXYZ(pos.X(), pos.Y(),pos.Z());
	}

	// Try and fit the circle using a Riemann fit. If 
	// it fails, try a basic fit with QuickFit.
	double my_BeamRMS=BeamRMS;
	// If there are a reasonable number of axial hits, de-weight the 
	// influence of the fake target point on the fit
	if (seed.hits.size()>4) my_BeamRMS=2.0*BeamRMS;
	if(seed.fit.FitCircleRiemann(TARGET_Z,my_BeamRMS)!=NOERROR
	   /* || seed.fit.chisq>20*/
	   ){
	  if(DEBUG_LEVEL>3)_DBG_<<"Riemann fit failed. Attempting regression fit..."<<endl;
	  if(seed.fit.FitCircle()!=NOERROR 
	     /*|| seed.fit.chisq>20*/
	     ){
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
			//double a = sl1_wire->origin.Angle(sl3_hit->hit->wire->origin);
			double a=sl1_wire->sdir.Angle(sl3_hit->hit->wire->sdir);

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
  /// several methods, any of which can flag the two seeds as clones.
  ///
  /// For the first method we match the hit pattern of the two seeds 
  /// looking for a fraction of shared hits above some minimum fraction.
  ///
  /// For the second method, we use the "asymmetry" of the circle
  /// centers. The asymmetry is defined as the distance between
  /// the two circle's centers divided by the sum of the magnitudes
  /// of the circle's centers relative to the origin. If this is below
  /// some threshold (say 0.05) then the tracks are seeds are
  /// considered clones.
  ///
  /// For the third method, we check the distance between the
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
    
    // Number of words in the hit pattern for seed 1
    unsigned int numwords=seed1.HitBitPattern.size();
    
    for(unsigned int j=0; j<seeds.size(); j++){
      if (j==i) continue;
      DCDCSeed &seed2 = seeds[j];
      if(!seed2.valid)continue;

      // Flag to indicate these are clones
      bool are_clones = false;
      
      //Look for seeds that share hits and mark a seed as invalid if it shares 
      //more than a certain fraction of hits with another seed and the
      //chisq/df from the circle fit is worse than that of the other fit by a 
      //certain amount. 
      unsigned int numcommon=0;
      unsigned int nhits1=seed1.hits.size()+seed1.stereo_hits.size();
      for (unsigned int k=0;k<numwords;k++){
	numcommon+=bitcount(seed1.HitBitPattern[k]&seed2.HitBitPattern[k]);
      }
      double hitratio=double(numcommon)/double(nhits1);
      if (DEBUG_LEVEL>3)_DBG_ << "Hit ratio = " << hitratio << endl;
      if (hitratio>0.85) are_clones=true;
      else if (hitratio>0.5){
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
      }
      // Remove a clone if necessary
      if(are_clones){
	// These seeds appear to be clones of one another. Mark one as not valid.
	unsigned int nhits2=seed2.hits.size()+seed2.stereo_hits.size();
	if (nhits1>nhits2 &&
	    (seed1.fit.chisq/double(seed1.fit.ndof)
	     <seed2.fit.chisq/double(seed2.fit.ndof))){
	  seed2.valid = false;
	}else{
	  seed1.valid = false;
	}
	
	if(DEBUG_LEVEL>3)
	  _DBG_<<"Filtering clone seed (seed1.fit.chisq="<<seed1.fit.chisq<<" seed2.fit.chisq="<<seed2.fit.chisq<<endl;
      }
    }
  }
}
	
// Calculate intersection point between circle and stereo wire	
jerror_t DTrackCandidate_factory_CDC::GetStereoPosition(const DCDCWire *wire,
							const DHelicalFit &fit,
							DVector3 &pos, 
							double &var_z,
							double d){
  DVector3 origin = wire->origin;
  DVector3 dir = (1./wire->udir.z())*wire->udir;
  double dx=origin.x()-fit.x0;  
  double dy=origin.y()-fit.y0;
  double ux=dir.x();
  double uy=dir.y();
  double temp1=ux*ux+uy*uy;
  double temp2=ux*dy-uy*dx;
  double b=-ux*dx-uy*dy;
  double dr=fit.r0-d;
  double r0_sq=dr*dr;
  double A=r0_sq*temp1-temp2*temp2;

  // Check that this wire intersects this circle
  if(A<0.0) return UNRECOVERABLE_ERROR; // line along wire does not intersect circle, ever.

  // Guess for variance for z: assume straw cell size??
  double temp=1.6/sin(wire->stereo);
  var_z=temp*temp/12.;
  
  // Calculate intersection points for the two roots 
  double B = sqrt(A);
  double dz1 = (b-B)/temp1;
  double dz2 = (b+B)/temp1;

  if(DEBUG_LEVEL>15) _DBG_<<"dz1="<<dz1<<" dz2="<<dz2<<endl;
  
  // At this point we must decide which value of alpha to use. 
  // For now, we just use the value closest to zero (i.e. closest to
  // the center of the wire).
  double dz=dz1;
  if (fabs(dz2)<fabs(dz1)){
    dz=dz2;
  }
		
  // Compute the position for this hit
  pos = origin + dz*dir;
 
  // distance along wire relative to origin
  double s=dz/cos(wire->stereo);
  
  if(DEBUG_LEVEL>15)
    _DBG_<<"s="<<s<<" ring="<<wire->ring<<" straw="<<wire->straw<<" stereo="<<wire->stereo<<endl;
  if(fabs(s) > 0.5*wire->L){
    return VALUE_OUT_OF_RANGE; // intersects beyond range of CDC
  }
 
  return NOERROR;
}



//------------------
// AddStereoHits
//------------------
void DTrackCandidate_factory_CDC::AddStereoHits(vector<DCDCTrkHit*> &stereo_hits, DCDCSeed &seed)
{
  if (stereo_hits.size()==0) return;

	// To find the z coordinate, we look at the 2D projection of the
	// stereo wire and find the intersection point of that with the
	// circle found in FitSeed().

  // Find phi angles of hits that are closest to the stereo layer 
  double phi_inner=0.,phi_outer=0.;
  int stereo_ring=stereo_hits[0]->hit->wire->ring;
  unsigned int last_index=seed.hits.size()-1;
  if (stereo_ring>seed.hits[0]->hit->wire->ring){
    phi_inner=seed.hits[0]->hit->wire->phi;
    phi_outer=phi_inner;
    //    return;
  }
  else if (stereo_ring<seed.hits[last_index]->hit->wire->ring){
    phi_inner=seed.hits[last_index]->hit->wire->phi;
    phi_outer=phi_inner;
    //return;
  }
  else{
    for (unsigned int i=0;i<last_index;i++){
      if (seed.hits[i]->hit->wire->ring>stereo_ring 
	  && seed.hits[i+1]->hit->wire->ring<stereo_ring){ 
	phi_inner=seed.hits[i+1]->hit->wire->phi;
	phi_outer=seed.hits[i]->hit->wire->phi;
      }
    }
  }

  

	// Loop over stereo hits to find z-values for any that cross this circle
	for(unsigned int i=0; i<stereo_hits.size(); i++){
		DCDCTrkHit *trkhit = stereo_hits[i];
		
		// Don't consider hits that are out of time with the seed
		if(fabs(trkhit->hit->tdrift-seed.tdrift_avg)>MAX_SEED_TIME_DIFF){
		  if(DEBUG_LEVEL>3)
		    _DBG_<<"Out of time stereo hit: seed.tdrift_avg="<<seed.tdrift_avg<<" tdrift="<<trkhit->hit->tdrift<<endl;
			continue;
		}
		
		// Calculate intersection points between circle and stereo wire
		const DCDCWire *wire = trkhit->hit->wire;
		DVector3 pos;
		double var_z=0.;
		if (GetStereoPosition(wire,seed.fit,pos,var_z)
		    !=NOERROR) continue;
		
		DCDCTrkHit mytrkhit = *trkhit; // we need to make a copy since x_stereo,... are unique based on the candidate
		mytrkhit.x_stereo = pos.X();
		mytrkhit.y_stereo = pos.Y();
		mytrkhit.z_stereo = pos.Z();
		
		// Verify that this hit is reasonably close to the axial hits in X and Y
		// If not, then drop this hit for this seed.
		double my_phi=pos.Phi();
		double phi_diff1=phi_inner-my_phi;
		if(phi_diff1>M_PI)phi_diff1 = M_TWO_PI - phi_diff1;
		double phi_diff2=phi_outer-my_phi;
		if(phi_diff2>M_PI)phi_diff2 = M_TWO_PI - phi_diff2;
		// convert to degrees and take the absolute value
		phi_diff1=fabs(57.3*phi_diff1);
		phi_diff2=fabs(57.3*phi_diff1);

		if (DEBUG_HISTS){
		  if (phi_diff1>phi_diff2) Hdphi_s->Fill(seed.fit.r0,phi_diff2);
		  else Hdphi_s->Fill(seed.fit.r0,phi_diff1);
		}


		if (phi_diff1>MAX_SEED_LINK_ANGLE && phi_diff2>MAX_SEED_LINK_ANGLE){
		  if(DEBUG_LEVEL>4)_DBG_<<"rejecting stereo hit at phi="<<pos.Phi()<<" for being too far away from axial hits(phi_diff="<< ((phi_diff1<phi_diff2)?phi_diff1:phi_diff2) <<")"<<endl;
			continue;
		} 
		// put the z-variance into mytrkhit
		mytrkhit.var_z=var_z;
		
		// Compute phi for the stereo wire
		DVector2 R(seed.fit.x0, seed.fit.y0);
		mytrkhit.phi_stereo = atan2(mytrkhit.y_stereo-R.Y(), mytrkhit.x_stereo-R.X());
		R*=-1.0; // make R point from center of circle to beamline instead of other way around
	
		mytrkhit.phi_stereo -= R.Phi(); // make angle relative to beamline
		
		// We want this to go either from 0 to +2pi for positive charge, or
		// 0 to -2pi for negative.
		double phi_hi = seed.fit.q>0.0 ? +M_TWO_PI:0.0;
		double phi_lo = seed.fit.q>0.0 ? 0.0:-M_TWO_PI;
		while(mytrkhit.phi_stereo<phi_lo){
		  mytrkhit.phi_stereo+=M_TWO_PI;
		}
		while(mytrkhit.phi_stereo>phi_hi){
		  mytrkhit.phi_stereo-=M_TWO_PI;
		}
		trkhit->flags |= USED;
		mytrkhit.flags |= VALID_STEREO;
		seed.stereo_hits.push_back(mytrkhit);
		if(DEBUG_LEVEL>10)
		  _DBG_<<"Adding CDC stereo hit: ring="<<mytrkhit.hit->wire->ring<<" straw="<<mytrkhit.hit->wire->straw<<endl;

		  if(DEBUG_LEVEL>15)
		  {
		  _DBG_<<" --- wire->udir X, Y, Z = "<<wire->udir.X()<<", "<<wire->udir.Y()<<", "<<wire->udir.Z()<<endl;
		  _DBG_<<" -- ring="<<wire->ring<<" trkhit->z_stereo="
		       <<mytrkhit.z_stereo<<" trkhit->y_stereo="
		       <<mytrkhit.y_stereo<<" trkhit->x_stereo="
		       <<mytrkhit.x_stereo<<endl;
		  _DBG_<<" -- phi_stereo="<<mytrkhit.phi_stereo
		       <<"  R.Phi()="<<R.Phi()<<"  (X,Y)=("<<R.X()<<", "
		       <<R.Y()<<")"<<endl;
		  _DBG__;
		}

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
		FindTheta(seed, VERTEX_Z_MIN, VERTEX_Z_MAX);
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
 /*
	if(seed.z_vertex<TARGET_Z_MIN || seed.z_vertex>TARGET_Z_MAX){
		if(DEBUG_LEVEL>3)_DBG_<<"Seed z-vertex outside of target range (z="<<seed.z_vertex<<" TARGET_Z_MIN="<<TARGET_Z_MIN<<" TARGET_Z_MAX="<<TARGET_Z_MAX<<endl;
		seed.valid=false;
	}
	
 */
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

	// Add center of target as a point to constrain the fit a little more
	r.push_back(0.0);
	z.push_back(TARGET_Z);
	
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
	double bin_width = M_TWO_PI/(double)Nbins;
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
	if(phi_avg<0.0)phi_avg += M_TWO_PI;
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

// Use a stereo hit on the track and the helical model for the trajectory to 
// guess the charge
void DTrackCandidate_factory_CDC::DCDCSeed::CheckCharge(){  
   // Get circle parameters
  double rc=fit.r0;
  double xc=fit.x0;
  double yc=fit.y0;
  // Stereo hit position
  DVector2 xy2(stereo_hits[0].x_stereo,stereo_hits[0].y_stereo);
  // Compute phi rotation from "vertex" to cdc hit
  double Phi1=atan2(-yc,-xc);
  double chord=xy2.Mod();
  double ratio=chord/(2.*rc);
  if (ratio>1.) return;
  double dphi=2.*asin(ratio);
  
  // Positive and negative changes in phi
  double phiplus=Phi1+dphi;
  double phiminus=Phi1-dphi;
  DVector2 plus(xc+rc*cos(phiplus),yc+rc*sin(phiplus));
  DVector2 minus(xc+rc*cos(phiminus),yc+rc*sin(phiminus));

  // Compute differences 
  double d2plus=(plus-xy2).Mod2();
  double d2minus=(minus-xy2).Mod2();
  // Guess for charge based on differences
  double my_q=1.;
  if (d2minus<d2plus){
    my_q=-1.;
  }
  // Fix initial guess if there is no agreement
  if (my_q!=fit.q){
    fit.q=my_q;
    fit.phi+=M_PI;
    if (fit.phi>M_PI) fit.phi-=M_TWO_PI;
  }
  
}

//------------------
// DCDCSeed::FindAverageBz
//------------------
double DTrackCandidate_factory_CDC::DCDCSeed::FindAverageBz( const DMagneticFieldMap *bfield
)
{
   //return 2.0;
	if(!bfield)return 0.0;

	double Bz_sum=0.0;
	for(unsigned int i=0; i<stereo_hits.size(); i++){
	  DCDCTrkHit *hit = &stereo_hits[i];
	  Bz_sum += bfield->GetBz(hit->x_stereo,hit->y_stereo,
				  hit->z_stereo);
	}
	
	return Bz_sum/(double)(stereo_hits.size());
}


//------------------
// FindThetaZRegression
//------------------
// Linear regression to find tan(lambda) and z_vertex.
// This method assumes that there are errors in both the z positions and 
// the arc lengths.
// Algorithm from Numerical Recipes in C (2nd. ed.), pp. 668-669.
jerror_t DTrackCandidate_factory_CDC::FindThetaZRegression(DCDCSeed &seed){

  if(DEBUG_LEVEL>3)_DBG_<<"Finding theta and z via linear regression method."<<endl;
  
  DHelicalFit *myfit=&seed.fit;
  if (myfit->normal.Mag()==0.) return VALUE_OUT_OF_RANGE;
  // Vector of intersections between the circle and the stereo wires
  vector<intersection_t>intersections;
			       
  // CDC stereo hits
  int old_ring=-1;
  for (unsigned int m=0;m<seed.stereo_hits.size();m++){
    unsigned int numbits=8*sizeof(unsigned int);
    seed.HitBitPattern[seed.stereo_hits[m].index/numbits]
      |=1<<seed.stereo_hits[m].index%numbits;

    DCDCTrkHit *trkhit=&seed.stereo_hits[m];

    if (trkhit->hit->wire->ring!=old_ring){
      //DVector3_with_perp intersection;
      intersection_t intersection;
      intersection.x=trkhit->x_stereo;
      intersection.y=trkhit->y_stereo;
      intersection.perp2=intersection.x*intersection.x+intersection.y*intersection.y;
      intersection.z=trkhit->z_stereo;
      intersection.var_z=trkhit->var_z;
      
      intersections.push_back(intersection);  
    }
    old_ring=trkhit->hit->wire->ring;
  }
   
  // Now, sort the entries
  sort(intersections.begin(),intersections.end(),SortIntersections);
 
  // Compute the arc lengths between the origin in x and y and (xi,yi)
  vector<double>arclengths(intersections.size()); 
  vector<double>ratios(intersections.size());
  double xc=myfit->x0;
  double yc=myfit->y0;
  double rc=myfit->r0;
  double two_rc=2.*rc;
  // Find POCA to beam line
  double myphi=atan2(yc,xc);
  double y0=yc-rc*sin(myphi);
  double x0=xc-rc*cos(myphi);

  // Arc length to first measurement
  double diffx=intersections[0].x-x0;
  double diffy=intersections[0].y-y0;
  double chord=sqrt(diffx*diffx+diffy*diffy);  
  double ratio=chord/two_rc;
  double s=(ratio<1.)?two_rc*asin(ratio):M_PI_2*two_rc;
  arclengths[0]=s;
  ratios[0]=ratio;

  // Find arc lengths for the rest of the stereo hits
  for (unsigned int m=1;m<arclengths.size();m++){
    unsigned int m_minus_1=m-1;
    diffx=intersections[m].x-intersections[m_minus_1].x;
    diffy=intersections[m].y-intersections[m_minus_1].y;
    chord=sqrt(diffx*diffx+diffy*diffy);  
    ratio=chord/two_rc;
    if (ratio>0.999) return VALUE_OUT_OF_RANGE;
    double ds=two_rc*asin(ratio);
    s+=ds;
    arclengths[m]=s;
    ratios[m]=ratio;
  }

  //Linear regression to find z0, tanl
  double tanl=0.,z0=0.;
  if (arclengths.size()>1){ // Do fit only if have more than one measurement
    DCDCLineFit fit;
    unsigned int n=fit.n=intersections.size();
    fit.s.resize(n);
    fit.var_s.resize(n);
    fit.z.resize(n);
    fit.var_z.resize(n);
    fit.w.resize(n);

    // Find average variances for z and s
    double avg_var_s=0.,avg_var_z=0.;
    double var_r=1.6*1.6/12.;  // assume cell size
    for (unsigned int m=0;m<n;m++){
      fit.s[m]=arclengths[m];     
      fit.var_s[m]=var_r/(1.-ratios[m]*ratios[m]);

      avg_var_s+=fit.var_s[m];
      avg_var_z+=intersections[m].var_z;

      if(DEBUG_LEVEL>5)
	_DBG_<<"Using CDC hit "<<m<<" z="<<intersections[m].z << " s=" 
	     << arclengths[m] <<endl;      
      
    }
    // Scale z errors according to the ratio of the average variances
    double scale2=avg_var_s/avg_var_z;
    double scale=sqrt(scale2);
    vector<double>weight(n);
    for (unsigned int m=0;m<n;m++){
      fit.z[m]=scale*intersections[m].z;
      fit.var_z[m]=scale2*intersections[m].var_z;
      weight[m]=fit.var_s[m]+fit.var_z[m];
    }
    // Perform preliminary fit to find the (scaled) slope tanl
    double sumv=0.,sumx=0.;
    double sumy=0.,sumxx=0.,sumxy=0.;
    for (unsigned int m=0;m<n;m++){
      //double temp=1./var_z[m];
      double temp=1./weight[m];
      sumv+=temp;
      sumx+=arclengths[m]*temp;
      sumy+=fit.z[m]*temp;
      sumxx+=arclengths[m]*arclengths[m]*temp;
      sumxy+=arclengths[m]*fit.z[m]*temp;
    }
    double Delta=sumv*sumxx-sumx*sumx;
    if (Delta==0.) return VALUE_OUT_OF_RANGE;   
    tanl=(sumv*sumxy-sumx*sumy)/Delta;
    fit.z0=(sumxx*sumy-sumx*sumxy)/Delta;

    // Convert tanl to an angle and create two other reference angles
    double angle[3];
    angle[0]=0.;
    angle[1]=atan(tanl);
    angle[2]=1.571;
    // Compute chi^2 values for line fits with these three angles
    double ch[3];
    for (unsigned int m=0;m<3;m++){
      ch[m]=fit.ChiXY(angle[m]);
    }
    // Bracket the minimum chi^2
    fit.BracketMinimumChisq(angle[0],angle[1],angle[2],ch[0],ch[1],ch[2]);
    // Find the minimum chi^2 using Brent's method and compute the best value
    // for lambda
    double lambda=0.;
    fit.FindMinimumChisq(angle[0],angle[1],angle[2],lambda);
    // Undo the scaling 
    z0=fit.z0/scale;
    tanl=tan(lambda)/scale;
  }
  else{
    z0=TARGET_Z;
    tanl=(intersections[0].z-z0)/arclengths[0];
  }
 
  if (z0>VERTEX_Z_MAX || z0<VERTEX_Z_MIN){
	if(DEBUG_LEVEL>5)_DBG_<<"Fit failed for theta-z via regression value out of target range (z="<<z0<<")"<<endl;
	return VALUE_OUT_OF_RANGE;
  }
  
  seed.fit.tanl=tanl;
  seed.theta=M_PI_2-atan(tanl);
  seed.z_vertex=z0;

  return NOERROR;
}

// Get the position and momentum at a fixed radius from the beam line
jerror_t DTrackCandidate_factory_CDC::GetPositionAndMomentum(DCDCSeed &seed,
							     DVector3 &pos,
							     DVector3 &mom){
  // Direction tangent
  double tanl=tan(M_PI_2-seed.theta);

  // Squared radius of cylinder outside start counter but inside CDC inner 
  // radius
  double r2=90.0;
 
  // Circle parameters
  double xc=seed.fit.x0;
  double yc=seed.fit.y0;
  double rc=seed.fit.r0;
  double rc2=rc*rc;
  double xc2=xc*xc;
  double yc2=yc*yc;
  double xc2_plus_yc2=xc2+yc2;

  // variables needed for intersection of circles
  double a=(r2-xc2_plus_yc2-rc2)/(2.*rc);
  double temp1=yc*sqrt(xc2_plus_yc2-a*a);
  double temp2=xc*a;
  double cosphi_plus=(temp2+temp1)/xc2_plus_yc2;
  double cosphi_minus=(temp2-temp1)/xc2_plus_yc2;

  // Check for intersections
  if(!isfinite(temp1) || !isfinite(temp2)){
    // We did not find an intersection between the two circles, so return 
    // sensible defaults for pos and mom
    double my_seed_phi=seed.fit.phi;
    double my_center_phi=atan2(yc,xc);
    double xv=xc-rc*cos(my_center_phi);
    double yv=yc-rc*sin(my_center_phi);
    pos.SetXYZ(xv,yv,seed.z_vertex);

    double pt=-0.5*seed.fit.p_trans*bfield->GetBz(pos.x(),pos.y(),pos.z());
    mom.SetXYZ(pt*cos(my_seed_phi),pt*sin(my_seed_phi),pt*tanl);

    return NOERROR;
  }

  // if we have intersections, find both solutions
  double phi_plus=-acos(cosphi_plus);
  double phi_minus=-acos(cosphi_minus);
  double x_plus=xc+rc*cosphi_plus;
  double x_minus=xc+rc*cosphi_minus;
  double y_plus=yc+rc*sin(phi_plus);
  double y_minus=yc+rc*sin(phi_minus);

  // if the resulting radial position on the circle from the fit does not agree
  // with the radius to which we are matching, we have the wrong sign for phi+ 
  // or phi-
  double r2_plus=x_plus*x_plus+y_plus*y_plus;
  double r2_minus=x_minus*x_minus+y_minus*y_minus;  
  if (fabs(r2-r2_plus)>EPS){
    phi_plus*=-1.;
    y_plus=yc+rc*sin(phi_plus);
  }
  if (fabs(r2-r2_minus)>EPS){
    phi_minus*=-1.;
    y_minus=yc+rc*sin(phi_minus);
  }

  // Choose phi- or phi+ depending on proximity to one of the cdc hits
  double xwire=seed.hits[0]->hit->wire->origin.x();
  double ywire=seed.hits[0]->hit->wire->origin.y();
  double dx=x_minus-xwire;
  double dy=y_minus-ywire;
  double d2_minus=dx*dx+dy*dy;
  dx=x_plus-xwire;
  dy=y_plus-ywire;
  double d2_plus=dx*dx+dy*dy;
  if (d2_plus>d2_minus){
    phi_minus*=-1.;
    if (seed.fit.q<0) phi_minus+=M_PI;  
    double myphi=atan2(yc,xc);
    double xv=xc-rc*cos(myphi);
    double yv=yc-rc*sin(myphi);
    double dx=x_minus-xv;
    double dy=y_minus-yv;
    double chord=sqrt(dx*dx+dy*dy);
    double two_rc=2.*rc;
    double ratio=chord/two_rc;
    double ds=(ratio<1.)?(two_rc*asin(ratio)):(two_rc*M_PI_2);
    pos.SetXYZ(x_minus,y_minus,seed.z_vertex+ds*tanl);

    double pt=-0.5*seed.fit.p_trans*bfield->GetBz(pos.x(),pos.y(),pos.z());
    mom.SetXYZ(pt*sin(phi_minus),pt*cos(phi_minus),pt*tanl);

  }
  else{
    phi_plus*=-1.;
    if (seed.fit.q<0) phi_plus+=M_PI;
    double myphi=atan2(yc,xc);
    double xv=xc-rc*cos(myphi);
    double yv=yc-rc*sin(myphi);
    double dx=x_plus-xv;
    double dy=y_plus-yv;
    double chord=sqrt(dx*dx+dy*dy);
    double two_rc=2.*rc;
    double ratio=chord/two_rc;
    double ds=(ratio<1.)?(two_rc*asin(ratio)):(two_rc*M_PI_2);
    pos.SetXYZ(x_plus,y_plus,seed.z_vertex+ds*tanl); 

    double pt=-0.5*seed.fit.p_trans*bfield->GetBz(pos.x(),pos.y(),pos.z());
    mom.SetXYZ(pt*sin(phi_plus),pt*cos(phi_plus),pt*tanl);

  }

  return NOERROR;
}


//-------------------------------------------------------------------------
// Routines for fitting a line to the stereo data
//-------------------------------------------------------------------------
// Compute the chi^2 for a line fit given errors in both s and z.  Also 
// computes current best guess for the scaled intercept z0. 
// Taken from Numerical Recipes in C (2nd ed.), p. 670.
double DTrackCandidate_factory_CDC::DCDCLineFit::ChiXY(double lambda){
  double tanl=tan(lambda);
  double sumw=0.,avg_s=0.,avg_z=0.,my_chi2=0.;
  for (unsigned i=0;i<n;i++){
    w[i]=1./(tanl*tanl*var_s[i]+var_z[i]);
    sumw+=w[i];
    avg_s+=w[i]*s[i];
    avg_z+=w[i]*z[i];
  }
  avg_s/=sumw;
  avg_z/=sumw;
  z0=avg_z-tanl*avg_s;
  for (unsigned int i=0;i<n;i++){
    double temp=z[i]-z0-tanl*s[i];
    my_chi2+=w[i]*temp*temp;
  }
  return my_chi2;
}

// Routine to bracket the minimum chi^2, from Numerical Recipes in C (2nd ed.),
// pp. 400-401.
#define SHFT(w,x,y,z) (w)=(x);(x)=(y);(y)=(z)
#define SIGN(x,y) ((y)>=0.0 ? fabs(x):-fabs(x))
void DTrackCandidate_factory_CDC
::DCDCLineFit::BracketMinimumChisq(double &a,double &b,double &c,
				   double &chi2a,double &chi2b,
				   double &chi2c){
  const double GOLD=1.618034;
  const double GLIMIT=100.0;

  chi2a=ChiXY(a);
  chi2b=ChiXY(b);
  double chi2u=0.;
  if (chi2b>chi2a){
    double dummy=0.;
    SHFT(dummy,a,b,dummy);
    SHFT(dummy,chi2b,chi2a,dummy);
  }
  c=b+GOLD*(b-a);
  chi2c=ChiXY(c);
  while (chi2b>chi2c){
    double r=(b-a)*(chi2b-chi2c);
    double q=(b-c)*(chi2b-chi2a);
    double q_minus_r=q-r;
    double fabs_q_minus_r=fabs(q_minus_r);
    double max=(fabs_q_minus_r>1.e-20)?fabs_q_minus_r:1.e-20;
    double u=b-((b-c)*q-(b-a)*r)/(2.*SIGN(max,q_minus_r));
    double ulim=b+GLIMIT*(c-b);
    if ((b-u)*(u-c)>0.0){
      chi2u=ChiXY(u);
      if (chi2u<chi2c){
	a=b;
	b=u;
	chi2a=chi2b;
	chi2b=chi2u;
	return;
      }
      else if (chi2u>chi2b){
	c=u;
	chi2c=chi2u;
	return;
      }
      u=c+GOLD*(c-b);
      chi2u=ChiXY(u);     
    }
    else if ((c-u)*(u-ulim)>0.0){
      chi2u=ChiXY(u);
      if (chi2u<chi2c){
	SHFT(b,c,u,c+GOLD*(c-b));
	SHFT(chi2b,chi2c,chi2u,ChiXY(u));
      }
    }
    else if ((u-ulim)*(ulim-c)>=0.0){
      u=ulim;
      chi2u=ChiXY(u);
    }
    else{
      u=c+GOLD*(c-b);
      chi2u=ChiXY(u);
    }
    SHFT(a,b,c,u);
    SHFT(chi2a,chi2b,chi2c,chi2u);
  }
}

// Use Brent's algorithm to find the "true" (within tolerance) minimum chi^2
// after bracketting. Taken from Numerical Recipes in C (2nd. Ed.), pp. 404-405.
double DTrackCandidate_factory_CDC
::DCDCLineFit::FindMinimumChisq(double ax,double bx, double cx, double &xmin){
  const double CGOLD=0.3819660;
  double a=(ax<cx)?ax:cx;
  double b=(ax>cx)?ax:cx;
  double x=bx,w=bx,v=bx;
  double fx=ChiXY(x),fv=fx,fw=fx,fu=0.;
  double tol=1e-3,err=0.0,d=0.,u=0.;
  for (int iter=1;iter<=100;iter++){
    double xm=0.5*(a+b);
    double tol1=tol*fabs(x)+1e-10;
    double tol2=2.0*tol1;
    if (fabs(x-xm)<=(tol2-0.5*(b-a))){
      xmin=x;
      return fx;
    }
    if (fabs(err)>tol1){
      double r=(x-w)*(fx-fv);
      double q=(x-v)*(fx-fw);
      double p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q>0.0) p=-p;
      q=fabs(q);
      double etemp=err;
      err=d;
      if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x)|| p>=q*(b-x)){
	d=CGOLD*(err=(x>=xm ? a-x : b-x));
      }
      else{
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
      }      
    } else {
      d=CGOLD*(err=(x>=xm ? a-x : b-x ));
    }
    u=(fabs(d)>=tol1 ? x+d : x+SIGN(tol1,d));
    fu=ChiXY(u);
    if (fu<=fx){
      if (u>=x) a=x; 
      else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    }
    else{
      if (u<x) a=u;
      else b=u;
      if (fu<=fw || w==x){
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      }
      else if (fu<=fv || v==x || v==w){
	v=u;
	fv=fu;
      }
    }
  }
  xmin=x;
  return fx;
}

// Find places where the stereo straws transition from + to - (or - to +)
// angles, add the corresponding intersection points as points on the 
// circle, and refit...
bool DTrackCandidate_factory_CDC::RefitCircleWithStereoIntersections(DCDCSeed &seed){
  vector<DVector3>intersections;
  for (unsigned int is=0;is<seed.stereo_hits.size()-1;is++){
    unsigned int inext=is+1;
    const DCDCWire *first_wire=seed.stereo_hits[is].hit->wire;
    const DCDCWire *second_wire=seed.stereo_hits[inext].hit->wire;
    if ((first_wire->ring>8 && second_wire->ring<=8)
	|| (first_wire->ring>20 && second_wire->ring<=20)){
      DVector3 u0=first_wire->origin;
      DVector3 udir=first_wire->udir;  
      DVector3 v0=second_wire->origin;
      DVector3 vdir=second_wire->udir;
      DVector3 diff=u0-v0;
      double u_dot_v=udir.Dot(vdir);
      double u_dot_diff=udir.Dot(diff);
      double v_dot_diff=vdir.Dot(diff);
      double scale=1./(1.-u_dot_v*u_dot_v);
      double ul=scale*(u_dot_v*v_dot_diff-u_dot_diff);
      double vl=scale*(v_dot_diff-u_dot_v*u_dot_diff);
      DVector3 pos=0.5*(u0+ul*udir+v0+vl*vdir);
      intersections.push_back(pos);
    }
  }
  if (intersections.size()>0){
    DHelicalFit myfit(seed.fit);
    for (unsigned int k=0;k<intersections.size();k++){
      myfit.AddHitXYZ(intersections[k].X(),intersections[k].Y(),
		      intersections[k].Z());
    }
    if (myfit.FitCircleRiemann(seed.fit.r0)==NOERROR){
      seed.fit.x0=myfit.x0;
      seed.fit.y0=myfit.y0;
      seed.fit.r0=myfit.r0;
      //seed.fit.p_trans=myfit.p_trans;
      seed.fit.p_trans=2.0*0.003*myfit.r0; // with |Bz|=2.0, will be fixed later
	      
      double myphi=atan2(myfit.y0,myfit.x0)-M_PI_2;
      if(myphi<0)myphi+=M_TWO_PI;
      if(myphi>=M_TWO_PI)myphi-=M_TWO_PI;
      seed.fit.phi=myphi;
      
      for (unsigned int js=0;js<seed.stereo_hits.size();js++){
	DCDCTrkHit *stereo=&seed.stereo_hits[js];
	const DCDCWire *wire=stereo->hit->wire;
	double var_z=0.;
	DVector3 pos;
	if (GetStereoPosition(wire,seed.fit,pos,var_z)==NOERROR){
	  stereo->x_stereo=pos.X();
	  stereo->y_stereo=pos.Y();
	  stereo->z_stereo=pos.Z();
	  stereo->var_z=var_z;	
	  // Compute phi for the stereo wire
	  DVector2 R(seed.fit.x0, seed.fit.y0);
	  stereo->phi_stereo = atan2(stereo->y_stereo-R.Y(),
				     stereo->x_stereo-R.X());
	  R*=-1.0; // make R point from center of circle to beamline instead of other way around
	  if(DEBUG_LEVEL>15){
	    _DBG_<<" -- ring="<<wire->ring
		 <<" trkhit->z_stereo="<<stereo->z_stereo
		 <<" trkhit->y_stereo="<<stereo->y_stereo
		 <<" trkhit->x_stereo="<<stereo->x_stereo<<endl;
	    _DBG_<<" -- phi_stereo="<<stereo->phi_stereo<<"  R.Phi()="<<R.Phi()<<"  (X,Y)=("<<R.X()<<", "<<R.Y()<<")"<<endl;
	    _DBG__;
	  }
	  stereo->phi_stereo -= R.Phi(); // make angle relative to beamline
	  
	  // We want this to go either from 0 to +2pi for positive charge, or
	  // 0 to -2pi for negative.
	  double phi_hi = seed.fit.q>0.0 ? +M_TWO_PI:0.0;
	  double phi_lo = seed.fit.q>0.0 ? 0.0:-M_TWO_PI;
	  while(stereo->phi_stereo<phi_lo){
	    stereo->phi_stereo+=M_TWO_PI;
	  }
	  while(stereo->phi_stereo>phi_hi){
	    stereo->phi_stereo-=M_TWO_PI;
	  }
	}
      } // loop over stereo hits
      return true;
    } // circle fit
  } // got +/- intersection?
  return false;
}


// Add stray stereo hits to the track.  This code is intended to capture those 
// hits for which the apparent intersection point is beyond the extent of the CDC
// but should have been matched to the track if the radius of the circle fit were
// closer to the true value.  The proximity of the hits to the hits already on the
// track is more restrictive than the other AddStereoHits routine to minimize false
// matches.
void DTrackCandidate_factory_CDC::AddStrayStereoHits(vector<DCDCTrkHit*> &stereo_hits,
						     DCDCSeed &seed,
						     bool &got_hits){
  for (unsigned int i=0;i<stereo_hits.size();i++){
    DCDCTrkHit *trkhit=stereo_hits[i]; 
    const DCDCWire *new_wire=trkhit->hit->wire;
    const DCDCWire *wire_in_seed=seed.stereo_hits[0].hit->wire;
    if (new_wire->ring<=wire_in_seed->ring) break;

    if (!(trkhit->flags & USED) && fabs(new_wire->straw-wire_in_seed->straw)<=5){
      DVector3 pos;
      double var_z;
      if (GetStereoPosition(new_wire,seed.fit,pos,var_z)!=UNRECOVERABLE_ERROR){
	DCDCTrkHit mytrkhit = *trkhit; // we need to make a copy since x_stereo,... are unique based on the candidate
	mytrkhit.x_stereo = pos.X();
	mytrkhit.y_stereo = pos.Y();
	mytrkhit.z_stereo = pos.Z();

	// put the z-variance into mytrkhit
	mytrkhit.var_z=var_z;
		
	// Compute phi for the stereo wire
	DVector2 R(seed.fit.x0, seed.fit.y0);
	mytrkhit.phi_stereo = atan2(mytrkhit.y_stereo-R.Y(), mytrkhit.x_stereo-R.X());
	R*=-1.0; // make R point from center of circle to beamline instead of other way around
	
	mytrkhit.phi_stereo -= R.Phi(); // make angle relative to beamline
		
	// We want this to go either from 0 to +2pi for positive charge, or
	// 0 to -2pi for negative.
	double phi_hi = seed.fit.q>0.0 ? +M_TWO_PI:0.0;
	double phi_lo = seed.fit.q>0.0 ? 0.0:-M_TWO_PI;
	while(mytrkhit.phi_stereo<phi_lo){
	  mytrkhit.phi_stereo+=M_TWO_PI;
	}
	while(mytrkhit.phi_stereo>phi_hi){
	  mytrkhit.phi_stereo-=M_TWO_PI;
	}
	trkhit->flags |= USED;
	mytrkhit.flags |= VALID_STEREO;
	seed.stereo_hits.push_back(mytrkhit);

	got_hits=true;

	if(DEBUG_LEVEL>10)_DBG_<<"Adding CDC stereo hit: ring="<<mytrkhit.hit->wire->ring<<" straw="<<mytrkhit.hit->wire->straw<<endl;
	
	if(DEBUG_LEVEL>15){
	  _DBG_<<" --- wire->udir X, Y, Z = "<<new_wire->udir.X()
	       <<", "<<new_wire->udir.Y()<<", "<<new_wire->udir.Z()<<endl;
	  _DBG_<<" -- ring="<<new_wire->ring<<" trkhit->z_stereo="
	       <<mytrkhit.z_stereo<<" trkhit->y_stereo="
	       <<mytrkhit.y_stereo<<" trkhit->x_stereo="
	       <<mytrkhit.x_stereo<<endl;
	  _DBG_<<" -- phi_stereo="<<mytrkhit.phi_stereo
	       <<"  R.Phi()="<<R.Phi()<<"  (X,Y)=("<<R.X()<<", "
	       <<R.Y()<<")"<<endl;
	  _DBG__;
	}	
      }
      else{
	// If the circle does not intersect the straw but the next ring is 
	// adjacent to the last ring included in the seed, also add this straw
	// to the seed...  This is to compensate for poor circle fits due to 
	// insufficient axial data.
	
	const DCDCWire *wire=seed.stereo_hits[0].hit->wire;
	if (abs(wire->ring-trkhit->hit->wire->ring)==1){	  
	  DCDCTrkHit mytrkhit = *trkhit; // we need to make a copy since x_stereo,... are unique based on the candidate
	  seed.stereo_hits.push_back(mytrkhit);
	  
	  trkhit->flags |= USED;
	  got_hits=true;
	  
	  if(DEBUG_LEVEL>10)
	    _DBG_<<"Adding CDC stereo hit: ring="<<mytrkhit.hit->wire->ring<<" straw="<<mytrkhit.hit->wire->straw<<endl;
	}
      }
    }
  }
}

// Code to determine tan(lambda)
void DTrackCandidate_factory_CDC::DoLineFit(DCDCSeed &seed){
  // First try the linear regression method
  if (FindThetaZRegression(seed)!=NOERROR){
    // If the linear regression doesn't work try the histogramming method
    // Fit stereo hits to get theta and vertex z position
    FindThetaZ(seed);
    if(!seed.valid){
      //_DBG_ << endl;
      // reset the valid flag
      seed.valid=true;
      
      // Assume that the track came from one end or the other 
      // of the target and use a point in one of the stereo 
      // layers to estimate tanl
      if (seed.z_vertex>TARGET_Z_MAX)
	seed.z_vertex=TARGET_Z_MAX;
      else
	seed.z_vertex=TARGET_Z_MIN;
      double x=seed.stereo_hits[0].x_stereo;
      double y=seed.stereo_hits[0].y_stereo;
      double tworc=2.0*seed.fit.r0;
      double ratio=sqrt(x*x+y*y)/tworc;
      if (ratio<1.){
	double tanl=(seed.stereo_hits[0].z_stereo-seed.z_vertex)/
	  (tworc*asin(ratio));
	seed.theta=M_PI_2-atan(tanl);
      }
      else seed.valid=false;
    }
    
  }
}

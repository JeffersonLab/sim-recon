// $Id$
//
//    File: DTrackCandidate_factory_CDC.cc
// Created: Thu Sep  6 14:47:48 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#include "DTrackCandidate_factory_CDC.h"
#include <cmath>
#include <JANA/JCalibration.h>

#define BeamRMS 0.5
#define EPS 1e-3

#define TWO(c)     (0x1u << (c))
#define MASK(c) ((unsigned int)(-1)) / (TWO(TWO(c)) + 1u)
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (TWO(c))) & MASK(c))

#ifndef M_TWO_PI
#define M_TWO_PI 6.28318530717958647692
#endif

using namespace std;

inline int bitcount(unsigned int n)
{
	n = COUNT(n, 0);
	n = COUNT(n, 1);
	n = COUNT(n, 2);
	n = COUNT(n, 3);
	n = COUNT(n, 4);
	//n = COUNT(n, 5);    for 64-bit integers
	return n;
}

inline bool CDCSortByRdecreasing(const DTrackCandidate_factory_CDC::DCDCTrkHit* hit1, const DTrackCandidate_factory_CDC::DCDCTrkHit* hit2)
{
	// use the ring number to sort by R(decreasing) and then straw(increasing)
	if(hit1->hit->wire->ring == hit2->hit->wire->ring)
		return hit1->hit->wire->straw < hit2->hit->wire->straw;
	return hit1->hit->wire->ring > hit2->hit->wire->ring;
}

inline bool CDCSortByChiSqPerNDFDecreasing(const DTrackCandidate_factory_CDC::DCDCTrackCircle* locTrackCircle1, const DTrackCandidate_factory_CDC::DCDCTrackCircle* locTrackCircle2)
{
	// largest weighted fit chisq/ndf is first
	return (locTrackCircle1->dWeightedChiSqPerDF > locTrackCircle2->dWeightedChiSqPerDF);
}

inline bool CDCSortByStereoChiSqPerNDFIncreasing(const DTrackCandidate_factory_CDC::DCDCTrackCircle* locTrackCircle1, const DTrackCandidate_factory_CDC::DCDCTrackCircle* locTrackCircle2)
{
	// smallest weighted theta/z chisq/ndf is first
	return (locTrackCircle1->dWeightedChiSqPerDF_Stereo < locTrackCircle2->dWeightedChiSqPerDF_Stereo);
}

inline bool CDCSort_Intersections(const DTrackCandidate_factory_CDC::intersection_t& locIntersection1, const DTrackCandidate_factory_CDC::intersection_t& locIntersection2)
{
	return (locIntersection1.perp2 < locIntersection2.perp2);
}

inline bool CDCSort_DeltaPhis(const pair<DTrackCandidate_factory_CDC::DCDCTrkHit*, double>& locDeltaPhiPair1, const pair<DTrackCandidate_factory_CDC::DCDCTrkHit*, double>& locDeltaPhiPair2)
{
	//smallest delta-phi is first
	return (locDeltaPhiPair1.second < locDeltaPhiPair2.second);
}

DTrackCandidate_factory_CDC::~DTrackCandidate_factory_CDC(){
  
}


//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_CDC::init(void)
{
	DEBUG_LEVEL = 0;
	MAX_DCDCTrkHitPoolSize = 200;
	MAX_DCDCSuperLayerSeedPoolSize = 50;
	MAX_HelicalFitPoolSize = 100;
	MAX_DCDCTrackCirclePoolSize = 100;

	MAX_ALLOWED_CDC_HITS = 10000;
	MAX_ALLOWED_TRACK_CIRCLES = 5000;
	MAX_HIT_DIST = 4.0; // cm //each straw is 5/8in (1.5875 cm) in diameter
	MAX_HIT_DIST2 = MAX_HIT_DIST*MAX_HIT_DIST;

	//when linking DCDCRingSeed's together to form DCDCSuperLayerSeed's, allow skipping of rings
		//for example, say there are hits in ring 1 and ring 3, but none in ring 2 because the particle didn't deposit enough energy.
			//setting MAX_NUM_RINGSEED_RINGS_SKIPABLE >= 1 will recover these cases
	MAX_NUM_RINGSEED_RINGS_SKIPABLE = 1; 

	MIN_SEED_HITS = 2;

	// to be labeled as a potential/definite spiral turn (in/out-wards), the following conditions need to be met (where appropriate):
	MIN_STRAWS_POTENTIAL_SPIRAL_TURN = 4; //minimum number of straws in a DCDCRingSeed necessary to be labeled as a potential spiral turn (in/out-wards)
	MIN_STRAWS_DEFINITE_SPIRAL_TURN = 6; //minimum number of straws in a DCDCRingSeed necessary to be labeled as a definite spiral turn (in/out-wards)
	MIN_STRAWS_ADJACENT_TO_SPIRAL_TURN = 3; // in some cases, this is the minimum number of straws in a DCDCRingSeed that is adjacent to the spiral turn (in/out-wards) 
	// if a spiral turn occurs between rings or outside the CDC, the below is the max # of straws that can be between two DCDCRingSeed's in the adjacent ring
		// the two DCDCRingSeed's would be in different super layer seeds
	MAX_STRAWS_BETWEEN_LINK_SPIRAL_TURN = 6;

	//within a super layer, if the number of seeds in a region of width SEED_DENSITY_BIN_STRAW_WIDTH is greater than MAX_SEEDS_IN_STRAW_BIN, reject all seeds passing through this region
		//note that since the # of straws increases with ring #, the region width is technically computed as a phi region: 
			//from phi of straw 'N' in the first ring of a given super layer, to the phi of straw 'N + DENSITY_BIN_STRAW_WIDTH - 1' 
	SEED_DENSITY_BIN_STRAW_WIDTH = 8;
	MAX_SEEDS_IN_STRAW_BIN = 15;

	//when true, will allow linking of super layer seeds to skip a super layer (in case there is a dead HV board there)
		//note that this feature is not fully tested!!
	ENABLE_DEAD_HV_BOARD_LINKING = false;

	//don't allow new track seeds to start after this super layer 
		//track seeds could start late if track is a decay product, or HV board is dead
		//this is to prevent tracks forming that are complete garbage (e.g. a spiral craziness, knockout electrons re-entering the CDC from the BCAL, etc.)
	MAX_SUPERLAYER_NEW_TRACK = 4;

	// the maximum # of hits allowed to be shared between the axial super layer seeds of DCDCTrackCircle's
		// if > than this amount, the track with the larger circle-fit weighted-chisq/ndf will be rejected
	MAX_COMMON_HIT_FRACTION = 0.49; //reject if exactly half

	MAX_DRIFT_TIME = 1000.0; // ns
	MAX_SEED_TIME_DIFF = 1000.0; // ns

	// used for identifying arms of a spiral: can be true if the centers of the fit circles are close together (relative to their difference from the origin)
	MIN_CIRCLE_ASYMMETRY = 0.10;

	// when calculating the final theta/z, include at least this many stereo hits for the calculation
		// don't want to include stereo hits whose projections onto the track circle are too far away
	MIN_PRUNED_STEREO_HITS = 4; 

	// when searching for unused axial hits to add to the track, require that the hit be within this #-degrees in phi to the circle fit
	MAX_UNUSED_HIT_LINK_ANGLE = 10.0; //degrees

	TARGET_Z = 65.0;
	VERTEX_Z_MIN = -100.0;
	VERTEX_Z_MAX = 200.0;

	cdchits_by_superlayer.resize(7);
	dSuperLayerSeeds.resize(7);
	for(unsigned int loc_i = 0; loc_i < 7; ++loc_i)
		superlayer_boundaries.push_back(4*(1 + loc_i));

	dNumStrawsPerRing.resize(28);

	dNumSeedDensityPhiBins = 360;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_CDC::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
	gPARMS->SetDefaultParameter("TRKFIND:DEBUG_LEVEL", DEBUG_LEVEL);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_ALLOWED_CDC_HITS", MAX_ALLOWED_CDC_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_ALLOWED_TRACK_CIRCLES", MAX_ALLOWED_TRACK_CIRCLES);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_HIT_DIST", MAX_HIT_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_NUM_RINGSEED_RINGS_SKIPABLE", MAX_NUM_RINGSEED_RINGS_SKIPABLE);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_SEED_HITS", MIN_SEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_STRAWS_POTENTIAL_SPIRAL_TURN", MIN_STRAWS_POTENTIAL_SPIRAL_TURN);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_STRAWS_DEFINITE_SPIRAL_TURN", MIN_STRAWS_DEFINITE_SPIRAL_TURN);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_STRAWS_ADJACENT_TO_SPIRAL_TURN", MIN_STRAWS_ADJACENT_TO_SPIRAL_TURN);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_STRAWS_BETWEEN_LINK_SPIRAL_TURN", MAX_STRAWS_BETWEEN_LINK_SPIRAL_TURN);
	gPARMS->SetDefaultParameter("TRKFIND:SEED_DENSITY_BIN_STRAW_WIDTH", SEED_DENSITY_BIN_STRAW_WIDTH);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEEDS_IN_STRAW_BIN", MAX_SEEDS_IN_STRAW_BIN);
	gPARMS->SetDefaultParameter("TRKFIND:ENABLE_DEAD_HV_BOARD_LINKING", ENABLE_DEAD_HV_BOARD_LINKING);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SUPERLAYER_NEW_TRACK", MAX_SUPERLAYER_NEW_TRACK);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_COMMON_HIT_FRACTION", MAX_COMMON_HIT_FRACTION);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_CIRCLE_ASYMMETRY", MIN_CIRCLE_ASYMMETRY);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_DRIFT_TIME", MAX_DRIFT_TIME);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_TIME_DIFF", MAX_SEED_TIME_DIFF);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_PRUNED_STEREO_HITS", MIN_PRUNED_STEREO_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_UNUSED_HIT_LINK_ANGLE", MAX_UNUSED_HIT_LINK_ANGLE);
	gPARMS->SetDefaultParameter("TRKFIND:VERTEX_Z_MIN", VERTEX_Z_MIN);
	gPARMS->SetDefaultParameter("TRKFIND:VERTEX_Z_MAX", VERTEX_Z_MAX);

	MAX_HIT_DIST2 = MAX_HIT_DIST*MAX_HIT_DIST;

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	dMagneticField = locApplication->GetBfield(runnumber);
	dFactorForSenseOfRotation=(dMagneticField->GetBz(0.,0.,65.)>0.)?-1.:1.;

	const DGeometry *locGeometry = locApplication->GetDGeometry(runnumber);
	JCalibration *jcalib = locApplication->GetJCalibration(runnumber);
	map<string, double> targetparms;
	if (jcalib->Get("TARGET/target_parms",targetparms)==false){
	  TARGET_Z = targetparms["TARGET_Z_POSITION"];
	}
	else{
	  locGeometry->GetTargetZ(TARGET_Z);
	}

	// Get the CDC wire table from the XML
	vector<vector<DCDCWire*> > locCDCWires;
	locGeometry->GetCDCWires(locCDCWires);
	for(size_t loc_i = 0; loc_i < locCDCWires.size(); ++loc_i)
		dNumStrawsPerRing[loc_i] = locCDCWires[loc_i].size();

	// Clean up after using wire map
	for (size_t i=0;i<locCDCWires.size();i++){
	  for (size_t j=0;j<locCDCWires[i].size();j++){
	    delete locCDCWires[i][j];
	  }
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_CDC::evnt(JEventLoop *locEventLoop, uint64_t eventnumber)
{
	// Reset
	dRejectedPhiRegions.clear();

	// Reset in case didn't clear before exiting event evaluation on last event. 
	Reset_Pools();

	// Get CDC hits
	if(Get_CDCHits(locEventLoop) != NOERROR)
	{
		Reset_Pools();
		return RESOURCE_UNAVAILABLE;
	}

	// Build Super Layer Seeds
	for(unsigned int loc_i = 0; loc_i < 7; ++loc_i)
	{
		if(DEBUG_LEVEL > 3)
			cout << "Find Seeds, Super Layer = " << loc_i + 1 << endl;
		Find_SuperLayerSeeds(cdchits_by_superlayer[loc_i], loc_i + 1);
		Reject_SuperLayerSeeds_HighSeedDensity(loc_i + 1);
	}
	// Search for Spiral Links in and between super layer seeds
	Set_SpiralLinkParams();
	if(DEBUG_LEVEL > 5)
	{
		cout << "init super layers" << endl;
		Print_SuperLayerSeeds();
	}

	// Build DCDCTrackCircle objects (each corresponds to (at most) one track):
	deque<DCDCTrackCircle*> locCDCTrackCircles;
	bool locStatusFlag = Build_TrackCircles(locCDCTrackCircles);
	if(!locStatusFlag)
	{
		//SHOULD SET JEVENT STATUS BIT HERE!!!
		Reset_Pools();
		return OBJECT_NOT_AVAILABLE;
	}
	if(locCDCTrackCircles.empty())
	{
		Reset_Pools();
		return NOERROR;
	}
	Handle_StereoAndFilter(locCDCTrackCircles, false); //false: not final pass: filter seeds, don't reject seeds with no stereo hits (will add unused below), etc.

	// If the last super layer of a track is not 7, search for lone, unused hits on the next super layer and add them to the track
		// Only add them if they are close enough to the track circle fit
	Add_UnusedHits(locCDCTrackCircles);

	// using axial on the track, redo the circle fit and re-calc theta-z
		// This is for when extra hits are picked up by Add_UnusedHits, and for including midpoints between SL2/SL3 & SL5/SL6
	//fit circles will reject fits if they aren't very good //false: fit all circles //true: add intersections between stereo layers
	Fit_Circles(locCDCTrackCircles, false, true); //will reject fits if they aren't very good
	if(DEBUG_LEVEL > 5)
	{
		cout << "final fit track circles" << endl;
		Print_TrackCircles(locCDCTrackCircles);
	}
	sort(locCDCTrackCircles.begin(), locCDCTrackCircles.end(), CDCSortByChiSqPerNDFDecreasing); //sort by circle-fit weighted chisq/ndf (largest first)

	Handle_StereoAndFilter(locCDCTrackCircles, true); //true: final pass: don't need to filter any more, get improved theta-z (although will reject if bad/no theta/z)
	if(DEBUG_LEVEL > 5)
	{
		cout << "final track circles" << endl;
		Print_TrackCircles(locCDCTrackCircles);
	}

	// Create track candidates (as long as p > 0!!)
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
		Create_TrackCandidiate(locCDCTrackCircles[loc_i]);

	// Reset memory before exiting event evaluation. 
	Reset_Pools();

	return NOERROR;
}

//------------------
// Reset_Pools
//------------------
void DTrackCandidate_factory_CDC::Reset_Pools(void)
{
	// delete pool contents if too large, preventing memory-leakage-like behavor.
	if(dCDCTrkHitPool_All.size() > MAX_DCDCTrkHitPoolSize)
	{
		for(size_t loc_i = MAX_DCDCTrkHitPoolSize; loc_i < dCDCTrkHitPool_All.size(); ++loc_i)
			delete dCDCTrkHitPool_All[loc_i];
		dCDCTrkHitPool_All.resize(MAX_DCDCTrkHitPoolSize);
	}
	dCDCTrkHitPool_Available = dCDCTrkHitPool_All;

	if(dCDCSuperLayerSeedPool_All.size() > MAX_DCDCSuperLayerSeedPoolSize)
	{
		for(size_t loc_i = MAX_DCDCSuperLayerSeedPoolSize; loc_i < dCDCSuperLayerSeedPool_All.size(); ++loc_i)
			delete dCDCSuperLayerSeedPool_All[loc_i];
		dCDCSuperLayerSeedPool_All.resize(MAX_DCDCSuperLayerSeedPoolSize);
	}
	dCDCSuperLayerSeedPool_Available = dCDCSuperLayerSeedPool_All;

	if(dHelicalFitPool_All.size() > MAX_HelicalFitPoolSize)
	{
		for(unsigned int loc_i = MAX_HelicalFitPoolSize; loc_i < dHelicalFitPool_All.size(); ++loc_i)
			delete dHelicalFitPool_All[loc_i];
		dHelicalFitPool_All.resize(MAX_HelicalFitPoolSize);
	}
	dHelicalFitPool_Available = dHelicalFitPool_All;

	if(dCDCTrackCirclePool_All.size() > MAX_DCDCTrackCirclePoolSize)
	{
		for(size_t loc_i = MAX_DCDCTrackCirclePoolSize; loc_i < dCDCTrackCirclePool_All.size(); ++loc_i)
			delete dCDCTrackCirclePool_All[loc_i];
		dCDCTrackCirclePool_All.resize(MAX_DCDCTrackCirclePoolSize);
	}
	dCDCTrackCirclePool_Available = dCDCTrackCirclePool_All;
}

DTrackCandidate_factory_CDC::DCDCTrkHit* DTrackCandidate_factory_CDC::Get_Resource_CDCTrkHit(void)
{
	DCDCTrkHit* locCDCTrkHit;
	if(dCDCTrkHitPool_Available.empty())
	{
		locCDCTrkHit = new DCDCTrkHit;
		dCDCTrkHitPool_All.push_back(locCDCTrkHit);
	}
	else
	{
		locCDCTrkHit = dCDCTrkHitPool_Available.back();
		dCDCTrkHitPool_Available.pop_back();
	}
	locCDCTrkHit->Reset();
	return locCDCTrkHit;
}

DTrackCandidate_factory_CDC::DCDCSuperLayerSeed* DTrackCandidate_factory_CDC::Get_Resource_CDCSuperLayerSeed(void)
{
	DCDCSuperLayerSeed* locCDCSuperLayerSeed;
	if(dCDCSuperLayerSeedPool_Available.empty())
	{
		locCDCSuperLayerSeed = new DCDCSuperLayerSeed;
		dCDCSuperLayerSeedPool_All.push_back(locCDCSuperLayerSeed);
	}
	else
	{
		locCDCSuperLayerSeed = dCDCSuperLayerSeedPool_Available.back();
		dCDCSuperLayerSeedPool_Available.pop_back();
	}
	locCDCSuperLayerSeed->Reset();
	return locCDCSuperLayerSeed;
}

DHelicalFit* DTrackCandidate_factory_CDC::Get_Resource_HelicalFit(void)
{
	DHelicalFit* locHelicalFit;
	if(dHelicalFitPool_Available.empty())
	{
		locHelicalFit = new DHelicalFit;
		dHelicalFitPool_All.push_back(locHelicalFit);
	}
	else
	{
		locHelicalFit = dHelicalFitPool_Available.back();
		dHelicalFitPool_Available.pop_back();
	}
	locHelicalFit->Reset();
	return locHelicalFit;
}

DTrackCandidate_factory_CDC::DCDCTrackCircle* DTrackCandidate_factory_CDC::Get_Resource_CDCTrackCircle(void)
{
	DCDCTrackCircle* locCDCTrackCircle;
	if(dCDCTrackCirclePool_Available.empty())
	{
		locCDCTrackCircle = new DCDCTrackCircle;
		dCDCTrackCirclePool_All.push_back(locCDCTrackCircle);
	}
	else
	{
		locCDCTrackCircle = dCDCTrackCirclePool_Available.back();
		dCDCTrackCirclePool_Available.pop_back();
	}
	locCDCTrackCircle->Reset();
	return locCDCTrackCircle;
}

//------------------
// Get_CDCHits
//------------------
jerror_t DTrackCandidate_factory_CDC::Get_CDCHits(JEventLoop* loop)
{
	// Get the "raw" hits. These already have the wire associated with them.
	vector<const DCDCTrackHit*> cdctrackhits;
	loop->Get(cdctrackhits);
	dNumCDCHits = cdctrackhits.size();

	// If there are no hits, then bail now
	if(cdctrackhits.empty())
		return RESOURCE_UNAVAILABLE;
	
	// If there are too many hits, bail with a warning message
	if(cdctrackhits.size() > MAX_ALLOWED_CDC_HITS)
	{
		cout << "Too many hits in CDC (" <<cdctrackhits.size() << ", max = " << MAX_ALLOWED_CDC_HITS << ")! Track finding in CDC bypassed for event " << loop->GetJEvent().GetEventNumber() << endl;
		cdctrackhits.clear();
		return UNRECOVERABLE_ERROR;
	}

	// clear old hits
	cdctrkhits.clear();
	for(unsigned int i = 0; i < cdchits_by_superlayer.size(); i++)
		cdchits_by_superlayer[i].clear();

	// Create DCDCTrkHit objects out of these.	
	int oldwire = -1;
	for(size_t i = 0; i< cdctrackhits.size(); ++i)
	{
		// Add to "master" list
		// ONLY FIRST HIT OF A WIRE
		int newwire = cdctrackhits[i]->wire->ring*1000 + cdctrackhits[i]->wire->straw;
		if(newwire == oldwire)
			continue;
		oldwire = newwire;

		if(DEBUG_LEVEL > 40)
			cout << "adding ring, straw = " << cdctrackhits[i]->wire->ring << ", " << cdctrackhits[i]->wire->straw << endl;

		// Add to "master" list
		DCDCTrkHit* cdctrkhit = Get_Resource_CDCTrkHit();
		cdctrkhit->index = i;
		cdctrkhit->hit = cdctrackhits[i];
		cdctrkhit->flags = NONE;
		cdctrkhit->flags |= NOISE; // (see below)
		cdctrkhits.push_back(cdctrkhit);
    
		// Sort into list of hits by superlayer
		for(size_t j = 0; j < superlayer_boundaries.size(); ++j)
		{
			if(cdctrkhit->hit->wire->ring <= int(superlayer_boundaries[j]))
			{
				cdchits_by_superlayer[j].push_back(cdctrkhit);
				break;
			}
		}
	}

	// Sort the individual superlayer lists by decreasing values of R
	for(size_t i = 0; i < cdchits_by_superlayer.size(); ++i)
		sort(cdchits_by_superlayer[i].begin(), cdchits_by_superlayer[i].end(), CDCSortByRdecreasing);

	// Filter out noise hits. All hits are initially flagged as "noise".
		// Hits with a neighbor within MAX_HIT_DIST have their noise flags cleared.
	// Also flag hits as out-of-time if their drift time is too large
	for(size_t i = 0; i < cdctrkhits.size(); ++i)
	{
		DCDCTrkHit *trkhit1 = cdctrkhits[i];
		if(trkhit1->hit->tdrift > MAX_DRIFT_TIME)
			trkhit1->flags |= OUT_OF_TIME;
		if(!(trkhit1->flags & NOISE))
			continue; // this hit already not marked for noise
		for(size_t j = 0; j < cdctrkhits.size(); ++j)
		{
			if(j == i)
				continue;
			double d2 = trkhit1->Dist2(cdctrkhits[j]);
			if(d2 > 9.0*MAX_HIT_DIST2)
				continue;
			trkhit1->flags &= ~NOISE;
			cdctrkhits[j]->flags &= ~NOISE;
			break;
		}
	}

	return NOERROR;
}

/*********************************************************************************************************************************************************************/
/********************************************************************** BUILD SUPER LAYER SEEDS **********************************************************************/
/*********************************************************************************************************************************************************************/

//---------------------
// Find_SuperLayerSeeds
//---------------------
void DTrackCandidate_factory_CDC::Find_SuperLayerSeeds(deque<DCDCTrkHit*>& locSuperLayerHits, unsigned int locSuperLayer)
{
	// Sort through hits ring by ring to find DCDCRingSeed's from neighboring wires in the same ring. 
	// What we want is a list of DCDCRingSeed's for each ring, which will then be combined to form DCDCSuperLayerSeeds
	// Each DCDCRingSeed is a list of adjacent hits ordered by straw number.
		// If a DCDCRingSeed crosses the straw = 1 barrier:
			// Then the first hit is the smallest straw with phi > pi, and the last hit is the largest straw with phi < pi

	// Clear DCDCSuperLayerSeed's
	deque<DCDCSuperLayerSeed*>& locSuperLayerSeeds = dSuperLayerSeeds[locSuperLayer - 1];
	locSuperLayerSeeds.clear();

	DCDCRingSeed locCDCRingSeed;
	deque<DCDCRingSeed> locCDCRingSeeds; //list of DCDCRingSeed's in a given ring
	deque<deque<DCDCRingSeed> > rings; //1st dimension is ring, 2nd dimension is DCDCRingSeed's in that ring
	int last_ring = -1;

	for(size_t i = 0; i < locSuperLayerHits.size(); ++i)
	{
		DCDCTrkHit *trkhit = locSuperLayerHits[i];
		if(DEBUG_LEVEL > 20)
			cout << "track hit ring, straw = " << trkhit->hit->wire->ring << ", " << trkhit->hit->wire->straw << endl;

		// Check if ring number has changed. 
		if(trkhit->hit->wire->ring != last_ring)
		{
			if(DEBUG_LEVEL > 20)
				cout << "new ring, last ring = " << trkhit->hit->wire->ring << ", " << last_ring << endl;
			//ring # has changed: save the current DCDCRingSeed (from the previous ring) (if not empty)
			if(!locCDCRingSeed.hits.empty())
				locCDCRingSeeds.push_back(locCDCRingSeed);
			//if > 1 DCDCRingSeed on the previous ring: compare first and last DCDCRingSeeds
				//if they are adjacent (extending through the straw = 1 boundary): merge DCDCRingSeeds
			if(locCDCRingSeeds.size() > 1)
			{
				unsigned int locMinStraw = locCDCRingSeeds[0].hits[0]->hit->wire->straw;
				unsigned int locMaxStraw = locCDCRingSeeds[locCDCRingSeeds.size() - 1].hits[locCDCRingSeeds[locCDCRingSeeds.size() - 1].hits.size() - 1]->hit->wire->straw;
				unsigned int locNumStrawsInRing = dNumStrawsPerRing[locCDCRingSeeds[0].hits[0]->hit->wire->ring - 1];
				if((locMinStraw + locNumStrawsInRing - locMaxStraw) <= 1) //merge ringseeds
				{
					if(DEBUG_LEVEL > 20)
						cout << "straw boundary: merge ringseeds" << endl;
					locCDCRingSeeds[0].hits.insert(locCDCRingSeeds[0].hits.begin(), locCDCRingSeeds[locCDCRingSeeds.size() - 1].hits.begin(), locCDCRingSeeds[locCDCRingSeeds.size() - 1].hits.end());
					//insert at beginning: straws now arranged as: ..., N - 2, N - 1, N, 1, 2, 3, ...
					locCDCRingSeeds.pop_back();
				}
			}
			//if there was at least one DCDCRingSeed found on the previous ring, save it in the 2d deque
			if(!locCDCRingSeeds.empty())
				rings.push_back(locCDCRingSeeds);
			if(DEBUG_LEVEL > 3)
				cout << "  ringseed hits:" << locCDCRingSeed.hits.size() << "  locCDCRingSeeds:" << locCDCRingSeeds.size() << endl;
			//reset for finding the next group of hits
			locCDCRingSeeds.clear();
			locCDCRingSeed.hits.clear();
			//save this hit and continue
			locCDCRingSeed.hits.push_back(trkhit);
			locCDCRingSeed.ring = trkhit->hit->wire->ring;
			locCDCRingSeed.linked = false;
			last_ring = trkhit->hit->wire->ring;
			continue;
		}
		
		// Check if this hit is a neighbor of the last hit added to the ringseed
		if((unsigned int)abs(locCDCRingSeed.hits[locCDCRingSeed.hits.size() - 1]->hit->wire->straw - trkhit->hit->wire->straw) > 1)
		{
			//not a neighbor: save old and create new ringseed
			if(DEBUG_LEVEL > 20)
				cout << "straw diff" << endl;
			if(!locCDCRingSeed.hits.empty())
				locCDCRingSeeds.push_back(locCDCRingSeed);
			if(DEBUG_LEVEL > 3)
				cout << "ringseed hits: " << locCDCRingSeed.hits.size() << endl;
			locCDCRingSeed.hits.clear();
			locCDCRingSeed.linked = false;
		}

		locCDCRingSeed.hits.push_back(trkhit);
		if(DEBUG_LEVEL > 20)
			cout << "push back hit straw = " << trkhit->hit->wire->straw << endl;
	}
	//save final seeds, check if need to merge
	if(!locCDCRingSeed.hits.empty())
		locCDCRingSeeds.push_back(locCDCRingSeed);
	if(locCDCRingSeeds.size() > 1)
	{
		//compare first and last ringseeds: if ringseed extends through straw boundary, merge ringseeds
		unsigned int locMinStraw = locCDCRingSeeds[0].hits[0]->hit->wire->straw;
		unsigned int locMaxStraw = locCDCRingSeeds[locCDCRingSeeds.size() - 1].hits[locCDCRingSeeds[locCDCRingSeeds.size() - 1].hits.size() - 1]->hit->wire->straw;
		unsigned int locNumStrawsInRing = dNumStrawsPerRing[locCDCRingSeeds[0].hits[0]->hit->wire->ring - 1];
		if((locMinStraw + locNumStrawsInRing - locMaxStraw) <= 1) //merge ringseeds
		{
			if(DEBUG_LEVEL > 20)
				cout << "straw boundary: merge ringseeds" << endl;
			locCDCRingSeeds[0].hits.insert(locCDCRingSeeds[0].hits.begin(), locCDCRingSeeds[locCDCRingSeeds.size() - 1].hits.begin(), locCDCRingSeeds[locCDCRingSeeds.size() - 1].hits.end());
			//insert at beginning: straws now arranged as: ..., N - 2, N - 1, N, 1, 2, 3, ...
			locCDCRingSeeds.pop_back();
		}
	}
	if(!locCDCRingSeeds.empty())
		rings.push_back(locCDCRingSeeds);
	if(DEBUG_LEVEL > 3)
		cout << "  ringseed hits:" << locCDCRingSeed.hits.size() << "  ringseeds:" << locCDCRingSeeds.size() << endl;
	if(DEBUG_LEVEL > 3)
		cout << "rings: " << rings.size() << endl;

	// Print all DCDCRingSeeds to screen
	if(DEBUG_LEVEL > 45)
	{
		for(size_t i = 0; i < rings.size(); ++i)
		{
			for(size_t k = 0; k < rings[i].size(); ++k)
			{
				cout << "hits for ringseed ring, seed indices " << i << ", " << k << ":" << endl;
				for(size_t j = 0; j < rings[i][k].hits.size(); ++j)
					cout << "wire ring, straw = " << rings[i][k].hits[j]->hit->wire->ring << ", " << rings[i][k].hits[j]->hit->wire->straw << endl;
			}
		}
	}

	// If we have no rings, then there must be no super layer seeds. Bail now. 
	if(rings.empty())
		return;

	// Loop over rings, creating DCDCSuperLayerSeed's from adjacent rings
	for(ringiter ring = rings.begin(); ring != rings.end(); ++ring)
	{
		deque<DCDCRingSeed>& locCDCRingSeeds = *ring;
		ringiter next_ring = ring;
		++next_ring;
		
		// Loop over ringseeds of this ring
		for(size_t j = 0; j < locCDCRingSeeds.size(); ++j)
		{
			if(locCDCRingSeeds[j].linked)
				continue;

			// This ringseed hasn't been used in a DCDCSuperLayerSeed yet. Start a new seed with it.
			deque<DCDCRingSeed*> parent;
			parent.push_back(&locCDCRingSeeds[j]);
			Link_RingSeeds(parent, next_ring, rings.end(), locSuperLayer, 0);
		}
	}

	// Set wire orientations 
	for(size_t loc_i = 0; loc_i < locSuperLayerSeeds.size(); ++loc_i)
	{
		locSuperLayerSeeds[loc_i]->dSuperLayer = locSuperLayer;
		locSuperLayerSeeds[loc_i]->dSeedIndex = loc_i;
		if((locSuperLayer == 7) || (locSuperLayer == 4) || (locSuperLayer == 1))
			locSuperLayerSeeds[loc_i]->dWireOrientation = WIRE_DIRECTION_AXIAL;
		else if((locSuperLayer == 2) || (locSuperLayer == 6))
			locSuperLayerSeeds[loc_i]->dWireOrientation = WIRE_DIRECTION_STEREOLEFT;
		else
			locSuperLayerSeeds[loc_i]->dWireOrientation = WIRE_DIRECTION_STEREORIGHT;
	}
}

//---------------
// Link_RingSeeds
//---------------
void DTrackCandidate_factory_CDC::Link_RingSeeds(deque<DCDCRingSeed*>& parent, ringiter ring, ringiter ringend, unsigned int locSuperLayer, unsigned int locNumPreviousRingsWithoutHit)
{
	/// Combine DCDCRingSeed's from rings into DCDCSuperLayerSeed's
	///
	/// This a a re-entrant routine (i.e. it calls itself recursively). Upon
	/// entry, <i>parent</i> contains a list of pointers to all of the ringseeds
	/// from the rings outside of <i>ring</i> that are to be combined into
	/// a seed. This will search through all ringseeds of <i>ring</i> and if
	/// any are found that can extend the parent, a copy of parent is made,
	/// the current ringseed of this ring is added to it, and then it is
	/// passed on to another call to this routine. If no matches are found, 
	/// then it will try skipping a ring to find matches (if enabled). 
	/// If still no matches are found (which will be the case for the outer-most ring), then
	/// the ringseeds in <i>parent</i> will be combined into a single DCDCSuperLayerSeed.

	// Make sure parent has at least one ringseed
	if(parent.empty())
	{
		cout << "parent has no ringseeds!!" << endl;
		return;
	}

	// Set flag to keep track of whether this is the end of the seed or not
	bool seed_extended = false;
	if(ring != ringend)
	{
		// Last ringseed in parent list is the one we need to compare to
		DCDCRingSeed *parent_ringseed = parent[parent.size() - 1];
		double r_parent = parent_ringseed->hits[0]->hit->wire->origin.Perp();
	
		// Loop over ringseeds in this ring
		deque<DCDCRingSeed> &locCDCRingSeeds = (*ring);
		++ring; // increment ring iterator to point to next level down in case we recall ouself below
		
		for(size_t i = 0; i < locCDCRingSeeds.size(); ++i)
		{
			// Calculate the transverse (to the beamline) distance between the two DCDCRingSeed's
			double dr = r_parent - locCDCRingSeeds[i].hits[0]->hit->wire->origin.Perp();
			double locTransverseDist2 = fabs(MinDist2(locCDCRingSeeds[i], *parent_ringseed) - dr*dr);
			// Check if this ringseed is close enough to the parent's to link them together
			if(DEBUG_LEVEL > 20)
				cout << "ring1, ring2, locTransverseDist2, MAX_HIT_DIST2, dr*dr = " << locCDCRingSeeds[i].hits[0]->hit->wire->ring << ", " << parent_ringseed->hits[0]->hit->wire->ring << ", " << locTransverseDist2 << ", " << MAX_HIT_DIST2 << ", " << dr*dr << endl;
			if(locTransverseDist2 < MAX_HIT_DIST2)
			{
				// link them together
				deque<DCDCRingSeed*> myparent = parent;
				myparent.push_back(&locCDCRingSeeds[i]);
				locCDCRingSeeds[i].linked = true;
				// recursive call: try to link this grouping of DCDCRingSeed's to a DCDCRingSeed in the next ring
				Link_RingSeeds(myparent, ring, ringend, locSuperLayer, 0);
				seed_extended = true;
			}
		}
		if((!seed_extended) && (locNumPreviousRingsWithoutHit < MAX_NUM_RINGSEED_RINGS_SKIPABLE))
		{
			// no link was found, but try to skip this ring and find a match in the next ring
			Link_RingSeeds(parent, ring, ringend, locSuperLayer, locNumPreviousRingsWithoutHit + 1);
			seed_extended = true;
		}
	}

	// Check if this is the end of the line.
	if(!seed_extended)
	{
		// This is the end of this seed. 
		// Check if this seed contains a spiral. 
		int locSpiralRingIndex = -1;
		for(size_t i = 0; i < parent.size(); ++i)
		{
			if(parent[i]->hits.size() < MIN_STRAWS_DEFINITE_SPIRAL_TURN)
				continue;
			locSpiralRingIndex = i;
			break;
		}

		// Check whether this seed is an obvious intersection of two tracks (with one of them a spiral)
		bool locSeparateSeedsFlag = false;
		if((locSpiralRingIndex != -1) && (parent.size() > 1))
		{
			double locAvgNumHitsInNonSpiralRing = 0.0;
			for(size_t i = 0; i < parent.size(); ++i)
			{
				if(int(i) == locSpiralRingIndex)
					continue;
				locAvgNumHitsInNonSpiralRing += parent[i]->hits.size();
			}
			locAvgNumHitsInNonSpiralRing /= double(parent.size() - 1);
			if(locAvgNumHitsInNonSpiralRing < 2)
				locSeparateSeedsFlag = true;
		}

		//if obvious intersection: save spiral ring-seed separately
		if(locSeparateSeedsFlag && (parent[locSpiralRingIndex]->hits.size() >= MIN_SEED_HITS))
		{
			DCDCSuperLayerSeed* locSuperLayerSeed = Get_Resource_CDCSuperLayerSeed();
			for(size_t loc_i = 0; loc_i < parent[locSpiralRingIndex]->hits.size(); ++loc_i)
				parent[locSpiralRingIndex]->hits[loc_i]->flags |= USED;
			dSuperLayerSeeds[locSuperLayer - 1].push_back(locSuperLayerSeed);
			locSuperLayerSeed->dCDCRingSeeds.push_back(*(parent[locSpiralRingIndex]));
		}

		//save normal seed if enough hits
		unsigned int locTotalNumHits = 0;
		for(size_t i = 0; i < parent.size(); ++i)
		{
			if((int(i) == locSpiralRingIndex) && locSeparateSeedsFlag)
				continue;
			locTotalNumHits += parent[i]->hits.size();
		}
		if(locTotalNumHits < MIN_SEED_HITS)
		{
			if(DEBUG_LEVEL > 10)
				cout << "rejecting seed due to too few hits (have " << locTotalNumHits << " need " << MIN_SEED_HITS << ")" << endl;
			return;
		}

		DCDCSuperLayerSeed* locSuperLayerSeed = Get_Resource_CDCSuperLayerSeed();
		for(size_t i = 0; i < parent.size(); ++i)
		{
			if((int(i) == locSpiralRingIndex) && locSeparateSeedsFlag)
				continue;
			locSuperLayerSeed->dCDCRingSeeds.push_front(*(parent[i])); //input rings were in reverse order!!
			for(size_t loc_j = 0; loc_j < parent[i]->hits.size(); ++loc_j)
				parent[i]->hits[loc_j]->flags |= USED;
		}
		dSuperLayerSeeds[locSuperLayer - 1].push_back(locSuperLayerSeed);
	}
}

//---------
// MinDist2
//---------
double DTrackCandidate_factory_CDC::MinDist2(const DCDCRingSeed& locInnerRingSeed, const DCDCRingSeed& locOuterRingSeed)
{
	const deque<DCDCTrkHit*>& locInnerSeedHits = locInnerRingSeed.hits;
	const deque<DCDCTrkHit*>& locOuterSeedHits = locOuterRingSeed.hits;
	return MinDist2(locInnerSeedHits, locOuterSeedHits);
}

//---------
// MinDist2
//---------
double DTrackCandidate_factory_CDC::MinDist2(const deque<DCDCTrkHit*>& locInnerSeedHits, const deque<DCDCTrkHit*>& locOuterSeedHits)
{
	/// Returns the minimum distance squared between the two seeds. Assumes all of the hits in a given set are on the same ring. 
	/// First it checks if the two seeds overlap in phi: if so, then the radial difference between the rings (squared) is returned. 
	/// Otherwise, only the first and last hits of the adjacent rings between each seed's hit list are used. 
	/// to calculate a maximum of 4 distances (minimum of 1), of which the smallest is returned.
	if(locInnerSeedHits.empty() || locOuterSeedHits.empty())
	{
		cout << "Number of seed hits 0! (Ninner = " << locInnerSeedHits.size() << " ,Nouter = " << locOuterSeedHits.size() << ")" << endl;
		return 1.0E10;
	}

	DCDCTrkHit* locInnermostRingFirstStrawHit = locInnerSeedHits.front();
	DCDCTrkHit* locInnermostRingLastStrawHit = locInnerSeedHits.back();
	DCDCTrkHit* locOutermostRingFirstStrawHit = locOuterSeedHits.front();
	DCDCTrkHit* locOutermostRingLastStrawHit = locOuterSeedHits.back();

	//see if seeds overlap in phi
	float locInnermostRingFirstStrawPhi = locInnermostRingFirstStrawHit->hit->wire->phi;
	float locInnermostRingLastStrawPhi = locInnermostRingLastStrawHit->hit->wire->phi;
	float locOutermostRingFirstStrawPhi = locOutermostRingFirstStrawHit->hit->wire->phi;
	float locOutermostRingLastStrawPhi = locOutermostRingLastStrawHit->hit->wire->phi;
	if(DEBUG_LEVEL > 100)
		cout << "inner ring: ring, first/last straws & phis = " << locInnermostRingFirstStrawHit->hit->wire->ring << ", " << locInnermostRingFirstStrawHit->hit->wire->straw << ", " << locInnermostRingLastStrawHit->hit->wire->straw << ", " << locInnermostRingFirstStrawPhi << ", " << locInnermostRingLastStrawPhi << endl;
	if(DEBUG_LEVEL > 100)
		cout << "outer ring: ring, first/last straws & phis = " << locOutermostRingFirstStrawHit->hit->wire->ring << ", " << locOutermostRingFirstStrawHit->hit->wire->straw << ", " << locOutermostRingLastStrawHit->hit->wire->straw << ", " << locOutermostRingFirstStrawPhi << ", " << locOutermostRingLastStrawPhi << endl;

	//account for phi = 0/2pi boundary
	bool locInnerRingCrossesBoundaryFlag = (locInnermostRingLastStrawPhi < locInnermostRingFirstStrawPhi);
	bool locOuterRingCrossesBoundaryFlag = (locOutermostRingLastStrawPhi < locOutermostRingFirstStrawPhi);
	if(DEBUG_LEVEL > 100)
		cout << "in/out boundary flags = " << locInnerRingCrossesBoundaryFlag << ", " << locOuterRingCrossesBoundaryFlag << endl;
	if(locOuterRingCrossesBoundaryFlag)
		locOutermostRingLastStrawPhi += M_TWO_PI;
	if(locInnerRingCrossesBoundaryFlag)
		locInnermostRingLastStrawPhi += M_TWO_PI;
	if(locOuterRingCrossesBoundaryFlag & (!locInnerRingCrossesBoundaryFlag) && ((locOutermostRingLastStrawPhi - locInnermostRingLastStrawPhi) > M_PI))
	{
		locInnermostRingFirstStrawPhi += M_TWO_PI;
		locInnermostRingLastStrawPhi += M_TWO_PI;
	}
	if(locInnerRingCrossesBoundaryFlag & (!locOuterRingCrossesBoundaryFlag) && ((locInnermostRingLastStrawPhi - locOutermostRingLastStrawPhi) > M_PI))
	{
		locOutermostRingFirstStrawPhi += M_TWO_PI;
		locOutermostRingLastStrawPhi += M_TWO_PI;
	}

	if(DEBUG_LEVEL > 100)
		cout << "final inner ring: ring, first/last straws & phis = " << locInnermostRingFirstStrawHit->hit->wire->ring << ", " << locInnermostRingFirstStrawHit->hit->wire->straw << ", " << locInnermostRingLastStrawHit->hit->wire->straw << ", " << locInnermostRingFirstStrawPhi << ", " << locInnermostRingLastStrawPhi << endl;
	if(DEBUG_LEVEL > 100)
		cout << "final outer ring: ring, first/last straws & phis = " << locOutermostRingFirstStrawHit->hit->wire->ring << ", " << locOutermostRingFirstStrawHit->hit->wire->straw << ", " << locOutermostRingLastStrawHit->hit->wire->straw << ", " << locOutermostRingFirstStrawPhi << ", " << locOutermostRingLastStrawPhi << endl;

	//finally check for overlaps
	double dr = locOutermostRingLastStrawHit->hit->wire->origin.Perp() - locInnermostRingFirstStrawHit->hit->wire->origin.Perp();
	if((locOutermostRingFirstStrawPhi >= locInnermostRingFirstStrawPhi) && (locOutermostRingFirstStrawPhi <= locInnermostRingLastStrawPhi))
		return dr*dr;
	if((locOutermostRingLastStrawPhi >= locInnermostRingFirstStrawPhi) && (locOutermostRingLastStrawPhi <= locInnermostRingLastStrawPhi))
		return dr*dr;
	if((locInnermostRingFirstStrawPhi >= locOutermostRingFirstStrawPhi) && (locInnermostRingFirstStrawPhi <= locOutermostRingLastStrawPhi))
		return dr*dr; //4th case not needed.  this case only needed if innermost ring is one wire across

	//no overlap, make all 4 comparisons between hits
	double d2, d2min;
	d2min = locInnermostRingFirstStrawHit->Dist2(locOutermostRingFirstStrawHit);
	if(locOutermostRingFirstStrawHit != locOutermostRingLastStrawHit)
	{
		d2 = locInnermostRingFirstStrawHit->Dist2(locOutermostRingLastStrawHit);
		if(d2 < d2min)
			d2min = d2;
	}
	if(locInnermostRingFirstStrawHit == locInnermostRingLastStrawHit)
		return d2min;

	d2 = locInnermostRingLastStrawHit->Dist2(locOutermostRingFirstStrawHit);
	if(d2 < d2min)
		d2min = d2;
	if(locOutermostRingFirstStrawHit != locOutermostRingLastStrawHit)
	{
		d2 = locInnermostRingLastStrawHit->Dist2(locOutermostRingLastStrawHit);
		if(d2 < d2min)
			d2min = d2;
	}

	return d2min;
}

//---------------------------------------
// Reject_SuperLayerSeeds_HighSeedDensity
//---------------------------------------
void DTrackCandidate_factory_CDC::Reject_SuperLayerSeeds_HighSeedDensity(unsigned int locSuperLayer)
{
	//The track finding algorithm will grind to a halt if there are too many super layer seeds in a given area
		//e.g. A slow spiral pion that loses energy very slowly will loop around many, many times before stopping, 
		//sometimes resulting in 40+ super layer seeds in a given super layer.
	//Therefore, if the density of super layer seeds in a given region is too high, reject all super layer seeds passing through that region of space, except for the ones on the edges
		//It will be impossible to determine track parameters in these regions anyway, so may as well reject them. 

	if(SEED_DENSITY_BIN_STRAW_WIDTH == 0)
		return; //rejection disabled

	deque<DCDCSuperLayerSeed*>& locSuperLayerSeeds = dSuperLayerSeeds[locSuperLayer - 1];
	if(locSuperLayerSeeds.size() <= MAX_SEEDS_IN_STRAW_BIN)
		return; //not enough for there to even possibly be an issue

	// histogram phis of all seeds within this super layer
	// We use a simple array to store our histogram here. We don't want to use ROOT histograms because they are not thread safe.
	// Setup histogram
	unsigned int hist[dNumSeedDensityPhiBins];
	for(unsigned int i = 0; i < dNumSeedDensityPhiBins; ++i)
		hist[i] = 0; // clear histogram
	double bin_width = M_TWO_PI/(double)dNumSeedDensityPhiBins;
	double hist_low_limit = 0.0; // lower edge of histogram limits

	//Find phi ranges of super layer seeds: store in map & histogram them
	map<unsigned int, pair<double, double> > locMapPhiRanges; //key is seed index, pair is phi range (first/last) //first > last if passes through phi = 0 boundary
	for(size_t loc_i = 0; loc_i < locSuperLayerSeeds.size(); ++loc_i)
	{
		double locSeedFirstPhi, locSeedLastPhi;
		Calc_SuperLayerPhiRange(locSuperLayerSeeds[loc_i], locSeedFirstPhi, locSeedLastPhi);
		locMapPhiRanges[locSuperLayerSeeds[loc_i]->dSeedIndex] = pair<double, double>(locSeedFirstPhi, locSeedLastPhi);
		if(DEBUG_LEVEL > 20)
			cout << "super layer, seed index, first phi, last phi = " << locSuperLayer << ", " << locSuperLayerSeeds[loc_i]->dSeedIndex << ", " << locSeedFirstPhi << ", " << locSeedLastPhi << endl;

		unsigned int locFirstPhiBin = (unsigned int)((locSeedFirstPhi - hist_low_limit)/bin_width);
		unsigned int locLastPhiBin = (unsigned int)((locSeedLastPhi - hist_low_limit)/bin_width);
		for(unsigned int locPhiBin = locFirstPhiBin; locPhiBin <= locLastPhiBin; ++locPhiBin)
			++hist[locPhiBin];
	}

	//determine search window size:
	unsigned int locOuterRing = superlayer_boundaries[locSuperLayer - 1];
	unsigned int locAverageNumStrawsInRing = (dNumStrawsPerRing[locOuterRing - 4] + dNumStrawsPerRing[locOuterRing - 1])/2; //I know it floored; it's close enough
	double locSearchBinPhiSize = double(SEED_DENSITY_BIN_STRAW_WIDTH)*M_TWO_PI/double(locAverageNumStrawsInRing);

	//find the search start point: try to start somewhere where there are no seeds for a range that is at least as large as the window; else start at 0
	int locStartPhiBin = -1;
	for(unsigned int locPhiBin = 0; locPhiBin < dNumSeedDensityPhiBins; ++locPhiBin)
	{
		if(hist[locPhiBin] > 0)
		{
			locStartPhiBin = -1;
			continue;
		}
		if(locStartPhiBin == -1)
			locStartPhiBin = locPhiBin;
		else if((locPhiBin - locStartPhiBin) > locSearchBinPhiSize)
			break;
	}
	if(locStartPhiBin == -1)
		locStartPhiBin = 0;

	// loop over phi range, scanning with a window of size "locSearchBinPhiSize," finding where the density of seeds is too high and marking those seeds for rejection
	set<unsigned int> locSeedsToReject;
	for(unsigned int locPhiBin = 0; locPhiBin < dNumSeedDensityPhiBins; ++locPhiBin)
	{
		unsigned int locReadPhiBin = locStartPhiBin + locPhiBin;
		if(locReadPhiBin >= dNumSeedDensityPhiBins)
			locReadPhiBin -= dNumSeedDensityPhiBins;

		double locWindowPhiRangeMin = hist_low_limit + bin_width*(double(locReadPhiBin));
		if(locWindowPhiRangeMin >= M_TWO_PI)
			locWindowPhiRangeMin -= M_TWO_PI;
		double locWindowPhiRangeMax = locWindowPhiRangeMin + locSearchBinPhiSize;

		map<unsigned int, pair<double, double> >::const_iterator locSeedIterator = locMapPhiRanges.begin();
		vector<unsigned int> locSeedsInThisRange;
		for(; locSeedIterator != locMapPhiRanges.end(); ++locSeedIterator)
		{
			double locSeedPhiRangeMin = locSeedIterator->second.first;
			double locSeedPhiRangeMax = locSeedIterator->second.second;
			if(Check_IfPhiRangesOverlap(locSeedPhiRangeMin, locSeedPhiRangeMax, locWindowPhiRangeMin, locWindowPhiRangeMax))
				locSeedsInThisRange.push_back(locSeedIterator->first);
		}
		if(locSeedsInThisRange.size() <= MAX_SEEDS_IN_STRAW_BIN)
			continue;

		//mark these seeds for deletion
		for(size_t loc_i = 0; loc_i < locSeedsInThisRange.size(); ++loc_i)
			locSeedsToReject.insert(locSeedsInThisRange[loc_i]);
	}

	//loop over seeds marked for rejection: keep the first and last seeds of each group
		//ignore phi = 0 barrier for now
		//also, mark phi regions as bad
	set<unsigned int>::const_iterator locRejectIterator = locSeedsToReject.begin();
	int locBeginIndex = -1, locPreviousIndex = 0, locFirstIndex = 0, locEndIndex = 0;
	set<unsigned int> locSeedsToNotReject;
	for(; locRejectIterator != locSeedsToReject.end(); ++locRejectIterator)
	{
		unsigned int locCurrentIndex = *locRejectIterator;
		if(DEBUG_LEVEL > 15)
			cout << "Marked for rejection: super layer, seed index = " << locSuperLayer << ", " << locCurrentIndex << endl;
		if(locBeginIndex == -1)
		{
			locFirstIndex = locCurrentIndex;
			locBeginIndex = locCurrentIndex;
			locSeedsToNotReject.insert(locBeginIndex);
			if(DEBUG_LEVEL > 15)
				cout << "don't reject: " << locBeginIndex << endl;
		}
		if(DEBUG_LEVEL > 30)
			cout << "window size, current phis, previous phis = " << locSearchBinPhiSize << ", " << locMapPhiRanges[locCurrentIndex].first << ", " << locMapPhiRanges[locCurrentIndex].second << ", " << locMapPhiRanges[locPreviousIndex].first << ", " << locMapPhiRanges[locPreviousIndex].second << endl;
		if(Check_IfPhiRangesOverlap(locMapPhiRanges[locCurrentIndex].first, locMapPhiRanges[locCurrentIndex].second, locMapPhiRanges[locPreviousIndex].first, locMapPhiRanges[locPreviousIndex].second + locSearchBinPhiSize))
		{
			locPreviousIndex = locCurrentIndex;
			locEndIndex = locCurrentIndex;
			continue;
		}
		if(Check_IfPhiRangesOverlap(locMapPhiRanges[locCurrentIndex].first, locMapPhiRanges[locCurrentIndex].second, locMapPhiRanges[locPreviousIndex].first - locSearchBinPhiSize, locMapPhiRanges[locPreviousIndex].second))
		{
			locPreviousIndex = locCurrentIndex;
			locEndIndex = locCurrentIndex;
			continue;
		}

		//seed isn't within locSearchBinPhiSize of previous seed: edge
		if(DEBUG_LEVEL > 15)
			cout << "Unmark for rejection: super layer, seed indexes = " << locSuperLayer << ", " << locBeginIndex << ", " << locPreviousIndex << endl;

		//too many seeds in this region: in future super layers, don't assume unmatched seeds in this region are new tracks
		dRejectedPhiRegions[locSuperLayer - 1].push_back(pair<double, double>(locMapPhiRanges[locBeginIndex].first, locMapPhiRanges[locPreviousIndex].second));

		locSeedsToNotReject.insert(locPreviousIndex);
		locSeedsToNotReject.insert(locCurrentIndex); //beginning of new section
		locBeginIndex = locCurrentIndex;
		locPreviousIndex = locCurrentIndex;
		locEndIndex = locCurrentIndex;
	}
	locSeedsToNotReject.insert(locEndIndex); //don't reject the last one //will check overlap below

	if(DEBUG_LEVEL > 30)
		cout << "first phis, last phis = " << locMapPhiRanges[locFirstIndex].first << ", " << locMapPhiRanges[locFirstIndex].second << ", " << locMapPhiRanges[locEndIndex].first << ", " << locMapPhiRanges[locEndIndex].second << endl;

	//check to see if last overlaps with first //overlap across phi = 0 barrier
	if(Check_IfPhiRangesOverlap(locMapPhiRanges[locFirstIndex].first, locMapPhiRanges[locFirstIndex].second, locMapPhiRanges[locEndIndex].first, locMapPhiRanges[locEndIndex].second))
	{
		if(DEBUG_LEVEL > 15)
			cout << "re-reject first/last = " << locFirstIndex << ", " << locEndIndex << endl;
		//do reject them
		locSeedsToNotReject.erase(locFirstIndex);
		locSeedsToNotReject.erase(locEndIndex);
	}


	// reject (and recycle) the seeds in areas where the seeds were too dense
	for(deque<DCDCSuperLayerSeed*>::iterator locDequeIterator = locSuperLayerSeeds.begin(); locDequeIterator != locSuperLayerSeeds.end();)
	{
		unsigned int locSeedIndex = (*locDequeIterator)->dSeedIndex;
		if((locSeedsToReject.find(locSeedIndex) != locSeedsToReject.end()) && (locSeedsToNotReject.find(locSeedIndex) == locSeedsToNotReject.end()))
		{
			if(DEBUG_LEVEL > 10)
				cout << "seed density too high, reject seed: " << locSeedIndex << endl;
			dCDCSuperLayerSeedPool_Available.push_back(*locDequeIterator); //recycle memory
		   locDequeIterator = locSuperLayerSeeds.erase(locDequeIterator);
		}
		else
			++locDequeIterator;
	}
}

//------------------------
// Calc_SuperLayerPhiRange
//------------------------
void DTrackCandidate_factory_CDC::Calc_SuperLayerPhiRange(DCDCSuperLayerSeed* locSuperLayerSeed, double& locSeedFirstPhi, double& locSeedLastPhi)
{
	// Calculate the phi-range across which the DCDCSuperLayerSeed extends. 
	locSeedFirstPhi = 9.9E9;
	locSeedLastPhi = -9.9E9;
	for(size_t loc_j = 0; loc_j < locSuperLayerSeed->dCDCRingSeeds.size(); ++loc_j)
	{
	  if (locSuperLayerSeed->dCDCRingSeeds[loc_j].hits.empty()) continue;
		DCDCTrkHit *locFirstStrawHit = locSuperLayerSeed->dCDCRingSeeds[loc_j].hits.front();
		DCDCTrkHit *locLastStrawHit = locSuperLayerSeed->dCDCRingSeeds[loc_j].hits.back();

		double locRingFirstPhi = locFirstStrawHit->hit->wire->phi;
		if(locRingFirstPhi < 0.0)
			locRingFirstPhi += M_TWO_PI;
		double locRingLastPhi = locLastStrawHit->hit->wire->phi;
		if(locRingLastPhi < 0.0)
			locRingLastPhi += M_TWO_PI;
		if(loc_j == 0)
		{
			locSeedFirstPhi = locRingFirstPhi;
			locSeedLastPhi = locRingLastPhi;
			continue;
		}

		if(locSeedFirstPhi > locSeedLastPhi) //seed goes across phi = 0 boundary
		{
			if(locRingFirstPhi > locRingLastPhi) //ring goes across phi = 0 boundary
			{
				if(locRingFirstPhi < locSeedFirstPhi)
					locSeedFirstPhi = locRingFirstPhi;
				if(locRingLastPhi > locSeedLastPhi)
					locSeedLastPhi = locRingLastPhi;
			}
			else
			{
				if((locRingFirstPhi > M_PI) && (locRingFirstPhi < locSeedFirstPhi))
					locSeedFirstPhi = locRingFirstPhi;
				else if((locRingLastPhi < M_PI) && (locRingLastPhi > locSeedLastPhi))
					locSeedLastPhi = locRingLastPhi;
			}
		}
		else //seed does not (so far) go across phi = 0 boundary
		{
			if(locRingFirstPhi > locRingLastPhi) //ring goes across phi = 0 boundary
			{
				locSeedFirstPhi = locRingFirstPhi;
				if(locRingLastPhi > locSeedLastPhi)
					locSeedLastPhi = locRingLastPhi;
			}
			else
			{
				if(locRingFirstPhi < locSeedFirstPhi)
					locSeedFirstPhi = locRingFirstPhi;
				if(locRingLastPhi > locSeedLastPhi)
					locSeedLastPhi = locRingLastPhi;
			}
		}
	}
}

//-------------------------
// Check_IfPhiRangesOverlap
//-------------------------
bool DTrackCandidate_factory_CDC::Check_IfPhiRangesOverlap(double locFirstSeedPhi, double locLastSeedPhi, double locTargetFirstPhi, double locTargetLastPhi)
{
	//if first phi > last phi 2nd, then it extends through the phi = 0 boundary
	if(locFirstSeedPhi > locLastSeedPhi) //seed extends through phi = 0 boundary
	{
		if(locTargetFirstPhi > locTargetLastPhi) //hist range extends through phi = 0 boundary
			return true; //clearly both intersect: both cross phi = 0
		else
		{
//S:      F---|---L
//H:    FXXXL |  FXXXL
			if(locTargetLastPhi > locFirstSeedPhi)
				return true;
			else if(locTargetFirstPhi < locLastSeedPhi)
				return true;
		}
	}
	else //seed does not extend through phi = 0 boundary
	{
		if(locTargetFirstPhi > locTargetLastPhi) //hist range extends through phi = 0 boundary
		{
//H:      F---|---L
//S:    FXXXL |  FXXXL
			if(locLastSeedPhi > locTargetFirstPhi)
				return true;
			else if(locFirstSeedPhi < locTargetLastPhi)
				return true;
		}
		else
		{
			if((locTargetFirstPhi > locFirstSeedPhi) && (locTargetFirstPhi < locLastSeedPhi))
				return true;
			else if((locTargetLastPhi > locFirstSeedPhi) && (locTargetLastPhi < locLastSeedPhi))
				return true;
			else if((locFirstSeedPhi > locTargetFirstPhi) && (locFirstSeedPhi < locTargetLastPhi))
				return true;
		}
	}
	return false;
}

/*********************************************************************************************************************************************************************/
/********************************************************************** SEARCH FOR SPIRAL LINKS **********************************************************************/
/*********************************************************************************************************************************************************************/

//---------------------
// Set_SpiralLinkParams
//---------------------
void DTrackCandidate_factory_CDC::Set_SpiralLinkParams(void)
{
	// Search for Spiral Links within and between super layer seeds
		//don't worry about getting every case: it's better to have extra particles reconstructed than to be overzealous and lose some!!
	for(size_t loc_i = 0; loc_i < dSuperLayerSeeds.size(); ++loc_i)
	{
		deque<DCDCSuperLayerSeed*>& locSuperLayerSeeds = dSuperLayerSeeds[loc_i];
		for(size_t loc_j = 0; loc_j < locSuperLayerSeeds.size(); ++loc_j)
		{
			for(size_t loc_k = loc_j + 1; loc_k < locSuperLayerSeeds.size(); ++loc_k)
			{
				if(SearchFor_SpiralTurn_TwoSeedsSharingManyHits(locSuperLayerSeeds[loc_j], locSuperLayerSeeds[loc_k]))
					continue;
				if(SearchFor_SpiralTurn_TwoSeedsSharingFewHits(locSuperLayerSeeds[loc_j], locSuperLayerSeeds[loc_k]))
					continue;
				if(SearchFor_SpiralTurn_MissingOrBetweenRings(locSuperLayerSeeds[loc_j], locSuperLayerSeeds[loc_k]))
					continue;
			}
			// if this super layer seed has not yet been linked to any other seed, check to see if itself is the spiral turn
				// e.g. only (and many) hits in one ring of super-layer 7.
			if(locSuperLayerSeeds[loc_j]->dSpiralLinkParams.empty()) //empty if not identified as a spiral yet
				SearchFor_SpiralTurn_SingleSeed(locSuperLayerSeeds[loc_j]);
		}
	}
}

//---------------------------------------------
// SearchFor_SpiralTurn_TwoSeedsSharingManyHits
//---------------------------------------------
bool DTrackCandidate_factory_CDC::SearchFor_SpiralTurn_TwoSeedsSharingManyHits(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2)
{
	//check if two seeds share at least MIN_STRAWS_POTENTIAL_SPIRAL_TURN on the same ring

	int locMaxSpiralNumHits = 0; //can have more than one potential spiral turn on a seed: could have two tracks (or two spiral arms) turning near each other
	//search all rings: spiral turn may not be on outermost ring if two tracks are crossing
	for(size_t loc_i = 0; loc_i < locSuperLayerSeed1->dCDCRingSeeds.size(); ++loc_i)
	{
		int locRing = locSuperLayerSeed1->dCDCRingSeeds[loc_i].ring;
		if(!locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locRing))
			continue; //the hits on this ring aren't the same in both DCDCSuperLayerSeed's
		int locLocalRingNumber = (locRing - 1)%4 + 1; //ranges from 1 -> 4 //ring of super layer

		// find how many straws are in this ring
		int locNumHits = locSuperLayerSeed1->dCDCRingSeeds[loc_i].hits.size();
		if(locNumHits < locMaxSpiralNumHits)
			continue; //already found a potential spiral link in this seed with a > # of straws

		if(locNumHits >= int(MIN_STRAWS_POTENTIAL_SPIRAL_TURN))
		{
			// check to make sure the hits on the first and last rings (in this super layer) of both DCDCSuperLayerSeed's aren't identical
				//if so, then if there truly is a spiral, this is the wrong combination of DCDCSuperLayerSeed's for it
			int locFirstRing1 = locSuperLayerSeed1->dCDCRingSeeds.front().ring;
			int locFirstRing2 = locSuperLayerSeed2->dCDCRingSeeds.front().ring;
			int locLastRing1 = locSuperLayerSeed1->dCDCRingSeeds.back().ring;
			int locLastRing2 = locSuperLayerSeed2->dCDCRingSeeds.back().ring;
			if((locFirstRing1 == locFirstRing2) && (locLastRing1 == locLastRing2))
			{
				if(locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locFirstRing1) && locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locLastRing1))
					continue;
			}

			int locTempSpiralNumHits = 0;
			if((locLocalRingNumber != 1) && (locLocalRingNumber != 4))
			{
				//spiral turn is not on an edge ring, it is on one of the two middle rings instead: there must also be MIN_STRAWS_ADJACENT_TO_SPIRAL_TURN hits in an adjacent ring
					//if spiral turn was on an edge ring, there may be many hits in the adjacent ring that is in another super layer, so can't do this check
				if(!SearchFor_SpiralTurn_ManyHitsAdjacentRing(locSuperLayerSeed1, locSuperLayerSeed2, locRing + 1, MIN_STRAWS_ADJACENT_TO_SPIRAL_TURN, locTempSpiralNumHits))
				{
					if(!SearchFor_SpiralTurn_ManyHitsAdjacentRing(locSuperLayerSeed1, locSuperLayerSeed2, locRing - 1, MIN_STRAWS_ADJACENT_TO_SPIRAL_TURN, locTempSpiralNumHits))
						continue; //not a spiral turn: there was only a few straws on the adjacent rings //perhaps this was two nearby tracks instead
				}
			}

			// is potential spiral turn, save results
			locMaxSpiralNumHits = locNumHits;
			DSpiralParams_t locSpiralTurnParams;
			//try to determine whether the track is turning inwards or outwards in this ring
			int locSpiralTurnRingFlag = 0;
			if(locSuperLayerSeed1->dCDCRingSeeds.size() == 1)
				locSpiralTurnRingFlag = ((locRing - 1)%4 > 1) ? 1 : -1;
			else if(loc_i == 0)
				locSpiralTurnRingFlag = -1; //turn on innermost ring
			else if(loc_i == (locSuperLayerSeed1->dCDCRingSeeds.size() - 1))
				locSpiralTurnRingFlag = 1; //turn on outermost ring
			else if(locSuperLayerSeed1->dCDCRingSeeds[loc_i + 1].hits.size() > locSuperLayerSeed1->dCDCRingSeeds[loc_i - 1].hits.size())
				locSpiralTurnRingFlag = -1; //more hits on outer ring than inner ring: turn on innermost ring
			else
				locSpiralTurnRingFlag = 1; //turn on outermost ring (if hit vector sizes are equal, assume outermost)
			locSpiralTurnParams.dSpiralTurnRingFlag = locSpiralTurnRingFlag;
			locSpiralTurnParams.dSpiralTurnRing = locRing;
			locSpiralTurnParams.dDefiniteSpiralTurnRingFlag = (locNumHits >= int(MIN_STRAWS_DEFINITE_SPIRAL_TURN)) ? locSpiralTurnRingFlag : 0;
			locSuperLayerSeed1->dSpiralLinkParams[locSuperLayerSeed2->dSeedIndex] = locSpiralTurnParams;
			locSuperLayerSeed2->dSpiralLinkParams[locSuperLayerSeed1->dSeedIndex] = locSpiralTurnParams;
			if(DEBUG_LEVEL > 10)
				cout << "SL" << locSuperLayerSeed1->dSuperLayer << " Seed" << locSuperLayerSeed1->dSeedIndex << " Spiral-linked to SL" << locSuperLayerSeed2->dSuperLayer << " Seed" << locSuperLayerSeed2->dSeedIndex << ": Share Many Hits, Ring = " << locRing << endl;
		}
	}
	return (locMaxSpiralNumHits > 0);
}

//--------------------------------------------
// SearchFor_SpiralTurn_TwoSeedsSharingFewHits
//--------------------------------------------
bool DTrackCandidate_factory_CDC::SearchFor_SpiralTurn_TwoSeedsSharingFewHits(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2)
{
	// check if: two seeds share a few hits on a given ring, and at least one has very many (unshared) hits on an adjacent ring: spiral turn

	int locMaxSpiralNumHits = 0; //can have more than one potential spiral turn on a seed: could have two tracks (or two spiral arms) turning near each other
	for(size_t loc_i = 0; loc_i < locSuperLayerSeed1->dCDCRingSeeds.size(); ++loc_i)
	{
		int locRing = locSuperLayerSeed1->dCDCRingSeeds[loc_i].ring;
		if(!locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locRing))
			continue;

		//check locRing - 1 for many hits: check if track is turning back inwards
		int locTempSpiralNumHits = -2;
		if(SearchFor_SpiralTurn_ManyHitsAdjacentRing(locSuperLayerSeed1, locSuperLayerSeed2, locRing - 1, MIN_STRAWS_POTENTIAL_SPIRAL_TURN, locTempSpiralNumHits))
		{
			//is potential spiral turn
			if(locTempSpiralNumHits > locMaxSpiralNumHits)
			{
				locMaxSpiralNumHits = locTempSpiralNumHits;
				DSpiralParams_t locSpiralTurnParams;
				int locSpiralTurnRingFlag = 1; //turning inwards
				locSpiralTurnParams.dSpiralTurnRingFlag = locSpiralTurnRingFlag;
				locSpiralTurnParams.dSpiralTurnRing = locRing;
				bool locIsDefiniteSpiralTurn = (locTempSpiralNumHits > int(MIN_STRAWS_DEFINITE_SPIRAL_TURN));
				locSpiralTurnParams.dDefiniteSpiralTurnRingFlag = locIsDefiniteSpiralTurn ? locSpiralTurnRingFlag : 0;
				locSuperLayerSeed1->dSpiralLinkParams[locSuperLayerSeed2->dSeedIndex] = locSpiralTurnParams;
				locSuperLayerSeed2->dSpiralLinkParams[locSuperLayerSeed1->dSeedIndex] = locSpiralTurnParams;
				if(DEBUG_LEVEL > 10)
					cout << "SL" << locSuperLayerSeed1->dSuperLayer << " Seed" << locSuperLayerSeed1->dSeedIndex << " Spiral-linked to SL" << locSuperLayerSeed2->dSuperLayer << " Seed" << locSuperLayerSeed2->dSeedIndex << ": Share Few Hits, Ring = " << locRing << endl;
			}
		}

		//check locRing + 1 for many hits: check if track is turning back outwards
		locTempSpiralNumHits = -2;
		if(SearchFor_SpiralTurn_ManyHitsAdjacentRing(locSuperLayerSeed1, locSuperLayerSeed2, locRing + 1, MIN_STRAWS_POTENTIAL_SPIRAL_TURN, locTempSpiralNumHits))
		{
			//is potential spiral turn
			if(locTempSpiralNumHits > locMaxSpiralNumHits)
			{
				locMaxSpiralNumHits = locTempSpiralNumHits;
				DSpiralParams_t locSpiralTurnParams;
				int locSpiralTurnRingFlag = -1; //turning outwards
				locSpiralTurnParams.dSpiralTurnRingFlag = locSpiralTurnRingFlag;
				locSpiralTurnParams.dSpiralTurnRing = locRing;
				bool locIsDefiniteSpiralTurn = (locTempSpiralNumHits > int(MIN_STRAWS_DEFINITE_SPIRAL_TURN));
				locSpiralTurnParams.dDefiniteSpiralTurnRingFlag = locIsDefiniteSpiralTurn ? locSpiralTurnRingFlag : 0;
				locSuperLayerSeed1->dSpiralLinkParams[locSuperLayerSeed2->dSeedIndex] = locSpiralTurnParams;
				locSuperLayerSeed2->dSpiralLinkParams[locSuperLayerSeed1->dSeedIndex] = locSpiralTurnParams;
				if(DEBUG_LEVEL > 10)
					cout << "SL" << locSuperLayerSeed1->dSuperLayer << " Seed" << locSuperLayerSeed1->dSeedIndex << " Spiral-linked to SL" << locSuperLayerSeed2->dSuperLayer << " Seed" << locSuperLayerSeed2->dSeedIndex << ": Share Few Hits, Ring = " << locRing << endl;
			}
		}
	}

	return (locMaxSpiralNumHits > 0);
}

//------------------------------------------
// SearchFor_SpiralTurn_ManyHitsAdjacentRing
//------------------------------------------
bool DTrackCandidate_factory_CDC::SearchFor_SpiralTurn_ManyHitsAdjacentRing(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2, int locRingToCheck, int locMinStrawsAdjacentRing, int& locMaxSpiralNumHits)
{
	//utility function, used to confirm if there are many hits on both seeds in the given ring
	if(locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locRingToCheck))
		return false; //these hits are the same

	deque<DCDCTrkHit*> locHits1;
	for(size_t loc_i = 0; loc_i < locSuperLayerSeed1->dCDCRingSeeds.size(); ++loc_i)
	{
		if(locSuperLayerSeed1->dCDCRingSeeds[loc_i].ring != locRingToCheck)
			continue;
		locHits1 = locSuperLayerSeed1->dCDCRingSeeds[loc_i].hits;
		break;
	}
	if(locHits1.empty())
		return false;

	deque<DCDCTrkHit*> locHits2;
	for(size_t loc_j = 0; loc_j < locSuperLayerSeed2->dCDCRingSeeds.size(); ++loc_j)
	{
		if(locSuperLayerSeed2->dCDCRingSeeds[loc_j].ring != locRingToCheck)
			continue;
		locHits2 = locSuperLayerSeed2->dCDCRingSeeds[loc_j].hits;
		break;
	}
	if(locHits2.empty())
		return false;

	int locNumHits1 = locHits1.size();
	int locNumHits2 = locHits2.size();

	if((locNumHits1 < int(locMinStrawsAdjacentRing)) || (locNumHits2 < int(locMinStrawsAdjacentRing)))
		return false; //neither seed has enough hits on the adjacent ring: return false

	// check to make sure the hits on the inner & outer rings of both seeds aren't identical
		//if there truly is a spiral, then this is the wrong combo of seeds for it
	int locFirstRing1 = locSuperLayerSeed1->dCDCRingSeeds.front().ring;
	int locFirstRing2 = locSuperLayerSeed2->dCDCRingSeeds.front().ring;
	int locLastRing1 = locSuperLayerSeed1->dCDCRingSeeds.back().ring;
	int locLastRing2 = locSuperLayerSeed2->dCDCRingSeeds.back().ring;
	if((locFirstRing1 == locFirstRing2) && (locLastRing1 == locLastRing2))
	{
		if(locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locFirstRing1) && locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locLastRing1))
			return false;
	}

	locMaxSpiralNumHits = locNumHits1;
	if(locNumHits2 > locMaxSpiralNumHits)
		locMaxSpiralNumHits = locNumHits2;

	return true;
}

//-------------------------------------------
// SearchFor_SpiralTurn_MissingOrBetweenRings
//-------------------------------------------
bool DTrackCandidate_factory_CDC::SearchFor_SpiralTurn_MissingOrBetweenRings(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2)
{
	// sometimes the spiral turn occurs between rings (or in a region with a dead high-voltage board). this function attempts to find these
		// the signature is that both seeds have many hits in a given ring.  these hits must be different but nearby, and the majority of the hits in the other rings must be further away from the other seed
	int locMaxSpiralNumHits = 0; //can have more than one potential spiral turn on a seed: could have two tracks (or two spiral arms) turning near each other
	for(size_t loc_i = 0; loc_i < locSuperLayerSeed1->dCDCRingSeeds.size(); ++loc_i)
	{
		int locRing = locSuperLayerSeed1->dCDCRingSeeds[loc_i].ring;

		//make sure there are enough straws on this ring in seed1
		deque<DCDCTrkHit*>& locHits1 = locSuperLayerSeed1->dCDCRingSeeds[loc_i].hits;
		int locNumHits1 = locHits1.size();
		if(locNumHits1 < int(MIN_STRAWS_POTENTIAL_SPIRAL_TURN))
			continue; //not enough straws

		for(size_t loc_j = 0; loc_j < locSuperLayerSeed2->dCDCRingSeeds.size(); ++loc_j)
		{
			if(locSuperLayerSeed2->dCDCRingSeeds[loc_j].ring != locRing)
				continue; //select the same ring as in seed1

			deque<DCDCTrkHit*>& locHits2 = locSuperLayerSeed2->dCDCRingSeeds[loc_j].hits;
			if(locHits2 == locHits1)
				break; //hits are shared: if truly was a spiral, will have been picked up earlier

			//make sure there are enough straws on this ring in seed2
			int locNumHits2 = locHits2.size();
			if(locNumHits2 < int(MIN_STRAWS_POTENTIAL_SPIRAL_TURN))
				break; //not enough straws

			bool locAtLeastOneSeedDefiniteTurnFlag = (locNumHits1 >= int(MIN_STRAWS_DEFINITE_SPIRAL_TURN));
			if(locNumHits2 >= int(MIN_STRAWS_DEFINITE_SPIRAL_TURN))
				locAtLeastOneSeedDefiniteTurnFlag = true;

			//check to make sure the seeds are nearby
			int locStrawNumDiff = abs(locHits2.front()->hit->wire->straw - locHits1.back()->hit->wire->straw);
			int locNumStrawsInRing = dNumStrawsPerRing[locRing - 1];
			if(locStrawNumDiff > locNumStrawsInRing/2)
				locStrawNumDiff = locNumStrawsInRing - locStrawNumDiff; //crosses straw count edge
			int locNumHits = locStrawNumDiff - 1;
			if(locNumHits > int(MAX_STRAWS_BETWEEN_LINK_SPIRAL_TURN))
				break; //straw groups are too far apart

			// check to make sure the hits on the inner & outer rings of both seeds aren't identical
				//if there truly is a spiral, then this is the wrong combo of seeds for it
			int locFirstRing1 = locSuperLayerSeed1->dCDCRingSeeds.front().ring;
			int locFirstRing2 = locSuperLayerSeed2->dCDCRingSeeds.front().ring;
			int locLastRing1 = locSuperLayerSeed1->dCDCRingSeeds.back().ring;
			int locLastRing2 = locSuperLayerSeed2->dCDCRingSeeds.back().ring;
			if((locFirstRing1 == locFirstRing2) && (locLastRing1 == locLastRing2))
			{
				if(locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locFirstRing1) && locSuperLayerSeed1->Are_AllHitsOnRingShared(locSuperLayerSeed2, locLastRing1))
					continue;
			}

			if((locNumHits1 < locMaxSpiralNumHits) && (locNumHits2 < locMaxSpiralNumHits))
				continue; //a seed with more hits was found earlier
			if(locNumHits1 > locMaxSpiralNumHits)
				locMaxSpiralNumHits = locNumHits1;
			if(locNumHits2 > locMaxSpiralNumHits)
				locMaxSpiralNumHits = locNumHits2;

			//will use results from circle fits to make sure the seeds are turning in the correct direction.
			DSpiralParams_t locSpiralTurnParams;
			// try to determine if the spiral is turning inwards or outwards
			int locSpiralTurnRingFlag = 0;
			if(locSuperLayerSeed1->dCDCRingSeeds.size() == 1)
				locSpiralTurnRingFlag = ((locRing - 1)%4 > 1) ? 1 : -1;
			else if(loc_i == 0)
				locSpiralTurnRingFlag = -1; //turn on innermost ring
			else if(loc_i == (locSuperLayerSeed1->dCDCRingSeeds.size() - 1))
				locSpiralTurnRingFlag = 1; //turn on outermost ring
			else if(locSuperLayerSeed1->dCDCRingSeeds[loc_i + 1].hits.size() > locSuperLayerSeed1->dCDCRingSeeds[loc_i - 1].hits.size())
				locSpiralTurnRingFlag = -1; //more hits on outer ring than inner ring: turn on innermost ring
			else
				locSpiralTurnRingFlag = 1; //turn on outermost ring (if hit vector sizes are equal, assume outermost)
			locSpiralTurnParams.dSpiralTurnRingFlag = locSpiralTurnRingFlag;
			locSpiralTurnParams.dSpiralTurnRing = locRing;
			locSpiralTurnParams.dDefiniteSpiralTurnRingFlag = locAtLeastOneSeedDefiniteTurnFlag ? locSpiralTurnRingFlag : 0;
			locSuperLayerSeed1->dSpiralLinkParams[locSuperLayerSeed2->dSeedIndex] = locSpiralTurnParams;
			locSuperLayerSeed2->dSpiralLinkParams[locSuperLayerSeed1->dSeedIndex] = locSpiralTurnParams;

			if(DEBUG_LEVEL > 10)
				cout << "SL" << locSuperLayerSeed1->dSuperLayer << " Seed" << locSuperLayerSeed1->dSeedIndex << " Spiral-linked to SL" << locSuperLayerSeed2->dSuperLayer << " Seed" << locSuperLayerSeed2->dSeedIndex << ": Share No Hits, Ring = " << locRing << endl;
		}
	}

	return (locMaxSpiralNumHits > 0);
}

//--------------------------------
// SearchFor_SpiralTurn_SingleSeed
//--------------------------------
bool DTrackCandidate_factory_CDC::SearchFor_SpiralTurn_SingleSeed(DCDCSuperLayerSeed* locSuperLayerSeed)
{
	// search for a spiral turn contained within this DCDCSuperLayerSeed
		// signature is a ring with many hits
	int locMaxSpiralNumHits = 0; //can have more than one potential spiral turn on a seed: could have two tracks (or two spiral arms) turning near each other
	for(size_t loc_i = 0; loc_i < locSuperLayerSeed->dCDCRingSeeds.size(); ++loc_i)
	{
		int locRing = locSuperLayerSeed->dCDCRingSeeds[loc_i].ring;
		int locNumHits = locSuperLayerSeed->dCDCRingSeeds[loc_i].hits.size();
		if(locNumHits < locMaxSpiralNumHits)
			continue; //already found a potential spiral link in this seed with a > # of straws

		if(locNumHits < int(MIN_STRAWS_POTENTIAL_SPIRAL_TURN))
			continue; //not enough straws in seed
		locMaxSpiralNumHits = locNumHits;

		DSpiralParams_t locSpiralTurnParams;
		int locSpiralTurnRingFlag = ((locRing - 1)%4 > 1) ? 1 : -1;
		locSpiralTurnParams.dSpiralTurnRingFlag = locSpiralTurnRingFlag;
		locSpiralTurnParams.dSpiralTurnRing = locRing;
		locSpiralTurnParams.dDefiniteSpiralTurnRingFlag = (locNumHits >= int(MIN_STRAWS_DEFINITE_SPIRAL_TURN)) ? locSpiralTurnRingFlag : 0;
		locSuperLayerSeed->dSpiralLinkParams[locSuperLayerSeed->dSeedIndex] = locSpiralTurnParams;
		if(DEBUG_LEVEL > 10)
			cout << "SL" << locSuperLayerSeed->dSuperLayer << " Seed" << locSuperLayerSeed->dSeedIndex << " Self-spiral-linked: Ring = " << locRing << endl;
	}
	return (locMaxSpiralNumHits > 0);
}

//----------------------
// Print_SuperLayerSeeds
//----------------------
void DTrackCandidate_factory_CDC::Print_SuperLayerSeeds(void)
{
	cout << "HITS BY SUPER LAYER & SEED INDEX:" << endl;
	for(unsigned int loc_i = 0; loc_i < 7; ++loc_i)
	{
		cout << "   SUPER LAYER " << loc_i + 1 << ":" << endl;
		for(unsigned int loc_j = 0; loc_j < dSuperLayerSeeds[loc_i].size(); ++loc_j)
		{
			cout << "      Seed Index " << loc_j << ":" << endl;
			DCDCSuperLayerSeed* locSuperLayerSeed = dSuperLayerSeeds[loc_i][loc_j];
			deque<DCDCTrkHit*> locHits;
			locSuperLayerSeed->Get_Hits(locHits);
			for(unsigned int loc_k = 0; loc_k < locHits.size(); ++loc_k)
				cout << "ring, straw, E (keV) = " << locHits[loc_k]->hit->wire->ring << ", " << locHits[loc_k]->hit->wire->straw << ", " << 1.0E6*locHits[loc_k]->hit->dE << endl;
			for(map<int, DSpiralParams_t>::iterator locIterator = locSuperLayerSeed->dSpiralLinkParams.begin(); locIterator != locSuperLayerSeed->dSpiralLinkParams.end(); ++locIterator)
				cout << "dSpiralLinkSeedIndex, dSpiralTurnRingFlag, dSpiralTurnRing, dDefiniteSpiralTurnRingFlag = " << locIterator->first << ", " << locIterator->second.dSpiralTurnRingFlag << ", " << locIterator->second.dSpiralTurnRing << ", " << locIterator->second.dDefiniteSpiralTurnRingFlag << endl;
		}
	}
}

/*********************************************************************************************************************************************************************/
/************************************************************************ Build Track Circles ************************************************************************/
/*********************************************************************************************************************************************************************/

//-------------------
// Build_TrackCircles
//-------------------
bool DTrackCandidate_factory_CDC::Build_TrackCircles(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	//return false if had to bail due to failure (i.e. too many track circles)
	//return true even if no track circles: it worked, it's just there aren't any

	//link DCDCSuperLayers together to form track circles
	locCDCTrackCircles.clear();
	for(unsigned int locOuterSuperLayer = 2; locOuterSuperLayer <= 7; ++locOuterSuperLayer)
	{
		unsigned int locInnerSuperLayer = locOuterSuperLayer - 1;
		if(locCDCTrackCircles.empty())
		{
			// no track circles yet: create new track circle objects using the super layer seeds in locInnerSuperLayer
			if(locInnerSuperLayer > MAX_SUPERLAYER_NEW_TRACK)
				return true; // don't create track circles beyond super layer MAX_SUPERLAYER_NEW_TRACK (e.g. knockout electrons from BCAL), return: no track circles!!!
			if(dSuperLayerSeeds[locInnerSuperLayer - 1].empty())
				continue; // no inner seeds, try next super layer
			//build new DCDCTrackCircle's: one for each super layer seed in this (inner) super layer
			for(size_t loc_j = 0; loc_j < dSuperLayerSeeds[locInnerSuperLayer - 1].size(); ++loc_j)
			{
				DCDCSuperLayerSeed* locSuperLayerSeed = dSuperLayerSeeds[locInnerSuperLayer - 1][loc_j];
				DCDCTrackCircle* locCDCTrackCircle = Get_Resource_CDCTrackCircle();
				if((locInnerSuperLayer == 1) || (locInnerSuperLayer == 4) || (locInnerSuperLayer == 7))
					locCDCTrackCircle->dSuperLayerSeeds_Axial.push_back(locSuperLayerSeed);
				else if((locInnerSuperLayer == 2) || (locInnerSuperLayer == 3))
					locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.push_back(deque<DCDCSuperLayerSeed*>(1, locSuperLayerSeed));
				else
					locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.push_back(deque<DCDCSuperLayerSeed*>(1, locSuperLayerSeed));
				locCDCTrackCircles.push_back(locCDCTrackCircle);
			}
		}
		//link existing track circles to DCDCSuperLayerSeed's in the outer super layer.  Any unused DCDCSuperLayerSeed's will be used to make new track circles
		if(!Link_SuperLayers(locCDCTrackCircles, locOuterSuperLayer))
			return false; //too many track circles
	}
	//if still no tracks, make DCDCSuperLayerSeed's in super layer 7 into new tracks (if allowed (probably shouldn't be))
	if(locCDCTrackCircles.empty() && (MAX_SUPERLAYER_NEW_TRACK >= 7))
	{
		for(size_t loc_j = 0; loc_j < dSuperLayerSeeds[6].size(); ++loc_j)
		{
			DCDCSuperLayerSeed* locSuperLayerSeed = dSuperLayerSeeds[6][loc_j];
			DCDCTrackCircle* locCDCTrackCircle = Get_Resource_CDCTrackCircle();
			locCDCTrackCircle->dSuperLayerSeeds_Axial.push_back(locSuperLayerSeed);
			locCDCTrackCircles.push_back(locCDCTrackCircle);
		}
	}
	if(locCDCTrackCircles.empty())
		return true;

	//shouldn't be possible, but check to see if too many (should've failed sooner)
	if(locCDCTrackCircles.size() >= MAX_ALLOWED_TRACK_CIRCLES)
	{
		if(DEBUG_LEVEL > 10)
			cout << "Too many track circles; bailing" << endl;
		locCDCTrackCircles.clear();
		return true;
	}

	// Reject track circles if they are definitely a spiral arm that turns back outwards in its inner super layer
	Reject_DefiniteSpiralArms(locCDCTrackCircles);
	if(locCDCTrackCircles.empty())
		return true;

	// Reject track circles that don't contain at least one axial and one stereo super layer, unless that axial is super layer 1
	Drop_IncompleteGroups(locCDCTrackCircles);
	if(locCDCTrackCircles.empty())
		return true;

	if(DEBUG_LEVEL > 5)
	{
		cout << "LINKED TRACK CIRCLES" << endl;
		Print_TrackCircles(locCDCTrackCircles);
	}

	//fit circles to the DCDCTrackCircle's.  This will reject fits if they aren't very good
		//1st false: fit all circles //2nd false: don't add intersections between stereo layers (wait until specific stereo super layers have been selected)
	Fit_Circles(locCDCTrackCircles, false, false); 
	if(locCDCTrackCircles.empty())
		return true;
	sort(locCDCTrackCircles.begin(), locCDCTrackCircles.end(), CDCSortByChiSqPerNDFDecreasing); //sort by circle-fit weighted chisq/ndf (largest first)

	if(DEBUG_LEVEL > 5)
	{
		cout << "post circle fit" << endl;
		Print_TrackCircles(locCDCTrackCircles);
	}

	return true;
}

//-----------------
// Link_SuperLayers
//-----------------
bool DTrackCandidate_factory_CDC::Link_SuperLayers(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer)
{
	/// Loop over the DCDCTrackCircles and the next super layer seeds and compare the positions of
		/// their first and last hits to see if we should link them together. 
	/// If at <= MAX_SUPERLAYER_NEW_TRACK: Any seeds from the outer super layer that are not linked will be added to the <i>DCDCTrackCircle</i> list.	
	/// Will try to link by skipping a super layer if a link was not previously performed and ENABLE_DEAD_HV_BOARD_LINKING is set to true (default false)

	if(DEBUG_LEVEL > 3)
		cout << "Linking seeds, locOuterSuperLayer = " << locOuterSuperLayer << endl;

	//Link locCDCTrackCircles's to the next super layer seeds
	unsigned int locInnerSuperLayer = locOuterSuperLayer - 1;
	if((locInnerSuperLayer == 1) || (locInnerSuperLayer == 4)) //else linking from stereo
		Link_SuperLayers_FromAxial(locCDCTrackCircles, locOuterSuperLayer, locInnerSuperLayer);
	else if((locOuterSuperLayer == 4) || (locOuterSuperLayer == 7))
	{
		if(!Link_SuperLayers_FromStereo_ToAxial(locCDCTrackCircles, locOuterSuperLayer, locInnerSuperLayer))
			return false; //linking failed: too many track circles
	}
	else
		Link_SuperLayers_FromStereo_ToStereo(locCDCTrackCircles, locOuterSuperLayer, locInnerSuperLayer);
	if(DEBUG_LEVEL > 25)
	{
		cout << "Link-to SL" << locOuterSuperLayer << ", v1" << endl;
		Print_TrackCircles(locCDCTrackCircles);
	}

	//For DCDCTrackCircle's that failed to link, try checking to see if can link to the previous super layer (e.g. a HV board was dead)
	if(ENABLE_DEAD_HV_BOARD_LINKING && (locInnerSuperLayer != 1))
	{
		--locInnerSuperLayer;
		if(DEBUG_LEVEL > 3)
			cout << "Try to skip super layer (e.g. HV board dead); new inner super layer = " << locInnerSuperLayer << endl;
		if((locInnerSuperLayer == 1) || (locInnerSuperLayer == 4))
			Link_SuperLayers_FromAxial(locCDCTrackCircles, locOuterSuperLayer, locInnerSuperLayer);
		else if((locOuterSuperLayer == 4) || (locOuterSuperLayer == 7))
		{
			if(!Link_SuperLayers_FromStereo_ToAxial(locCDCTrackCircles, locOuterSuperLayer, locInnerSuperLayer))
				return false; //linking failed: too many track circles
		}
		else
			Link_SuperLayers_FromStereo_ToStereo(locCDCTrackCircles, locOuterSuperLayer, locInnerSuperLayer);
		if(DEBUG_LEVEL > 25)
		{
			cout << "Link-to SL" << locOuterSuperLayer << ", v2" << endl;
			Print_TrackCircles(locCDCTrackCircles);
		}
	}

	// For outer seeds that failed to link, save as new track circles (even if they are stereo layers), UNLESS they are after MAX_SUPERLAYER_NEW_TRACK
		// UNLESS: a region at about the same phi in a previous super layer had a very high density of hits (such that all seeds in it were rejected)
			// Don't want to create new tracks in this case!
	if(locOuterSuperLayer <= MAX_SUPERLAYER_NEW_TRACK)
	{
		if(DEBUG_LEVEL > 3)
			cout << "Save any unlinked as new track circles" << endl;
		for(size_t loc_i = 0; loc_i < dSuperLayerSeeds[locOuterSuperLayer - 1].size(); ++loc_i)
		{
			DCDCSuperLayerSeed* locSuperLayerSeed = dSuperLayerSeeds[locOuterSuperLayer - 1][loc_i];
			if(locSuperLayerSeed->linked)
				continue; //previously linked, don't create new object
			//make sure this seed doesn't correspond to a region in phi where all of the seeds were deleted earlier
			//calc phi range of seed
			double locSeedFirstPhi, locSeedLastPhi;
			Calc_SuperLayerPhiRange(locSuperLayerSeed, locSeedFirstPhi, locSeedLastPhi);
			bool locHitDensityPreviouslyTooHighFlag = false;
			for(size_t loc_j = 1; loc_j < locOuterSuperLayer; ++loc_j)
			{
				for(size_t loc_k = 0; loc_k < dRejectedPhiRegions[loc_j].size(); ++loc_k)
				{
					if(!Check_IfPhiRangesOverlap(locSeedFirstPhi, locSeedLastPhi, dRejectedPhiRegions[loc_j][loc_k].first, dRejectedPhiRegions[loc_j][loc_k].second))
						continue;
					locHitDensityPreviouslyTooHighFlag = true;
					break;
				}
				if(locHitDensityPreviouslyTooHighFlag)
					break;
			}
			if(locHitDensityPreviouslyTooHighFlag)
				continue;
			//create new track circle
			if(DEBUG_LEVEL > 10)
				cout << "Unlinked seed; create new track circle: super layer & seed index = " << locSuperLayerSeed->dSuperLayer << ", " << locSuperLayerSeed->dSeedIndex << endl;
			if(locCDCTrackCircles.size() >= MAX_ALLOWED_TRACK_CIRCLES)
			{
				if(DEBUG_LEVEL > 10)
					cout << "Too many track circles; bailing" << endl;
				locCDCTrackCircles.clear();
				return false;
			}
			DCDCTrackCircle* locCDCTrackCircle = Get_Resource_CDCTrackCircle();
			if((locOuterSuperLayer == 4) || (locOuterSuperLayer == 7))
				locCDCTrackCircle->dSuperLayerSeeds_Axial.push_back(locSuperLayerSeed);
			else if((locOuterSuperLayer == 2) || (locOuterSuperLayer == 3))
				locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.push_back(deque<DCDCSuperLayerSeed*>(1, locSuperLayerSeed));
			else
				locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.push_back(deque<DCDCSuperLayerSeed*>(1, locSuperLayerSeed));
			locCDCTrackCircles.push_back(locCDCTrackCircle);
		}
		if(DEBUG_LEVEL > 25)
		{
			cout << "Link-to SL" << locOuterSuperLayer << ", v3" << endl;
			Print_TrackCircles(locCDCTrackCircles);
		}
	}

	return true;
}

//---------------------------
// Link_SuperLayers_FromAxial
//---------------------------
void DTrackCandidate_factory_CDC::Link_SuperLayers_FromAxial(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer, unsigned int locInnerSuperLayer)
{
	//Link locCDCTrackCircles from an axial super layer to a stereo super layer (must be stereo: cannot skip two (both stereo) super layers)
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		DCDCTrackCircle* locCDCTrackCircle = locCDCTrackCircles[loc_i];
		if(locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
			continue; // no axial super layer seeds to link from
		DCDCSuperLayerSeed* locSuperLayerSeed1 = locCDCTrackCircle->dSuperLayerSeeds_Axial.back();
		if(locSuperLayerSeed1->dSuperLayer != locInnerSuperLayer)
			continue; // e.g. the track ended earlier
		if(DEBUG_LEVEL > 10)
			cout << "Seed 1 Super Layer, Seed Index = " << locSuperLayerSeed1->dSuperLayer << ", " << locSuperLayerSeed1->dSeedIndex << endl;
		if(!Check_IfShouldAttemptLink(locSuperLayerSeed1, true))
			continue; //don't attempt seed link if the inner seed was a definite spiral turn that was turning back inwards, etc. 

		// if trying to skip a layer, first check to make sure this seed wasn't already linked (to layer "locInnerSuperLayer + 1")
		if((locInnerSuperLayer == 1) && (!locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.empty()))
			continue; //already linked to some inner stereo seeds
		else if((locInnerSuperLayer == 4) && (!locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.empty()))
			continue; //already linked to some outer stereo seeds

		for(size_t loc_j = 0; loc_j < dSuperLayerSeeds[locOuterSuperLayer - 1].size(); ++loc_j)
		{
			DCDCSuperLayerSeed* locSuperLayerSeed2 = dSuperLayerSeeds[locOuterSuperLayer - 1][loc_j];
			if(DEBUG_LEVEL > 10)
				cout << "Seed 2 Super Layer, Seed Index = " << locSuperLayerSeed2->dSuperLayer << ", " << locSuperLayerSeed2->dSeedIndex << endl;
			if(!Check_IfShouldAttemptLink(locSuperLayerSeed2, false))
				continue; //don't attempt seed link if the outer seed was a definite spiral turn that was turning back outwards, etc. 
			if(DEBUG_LEVEL > 10)
				cout << "Attempting Seed Link" << endl;
			if(!Attempt_SeedLink(locSuperLayerSeed1, locSuperLayerSeed2)) //see if seeds are nearby enough to link
				continue;
			//LINK SUCCESSFUL!!
			if(DEBUG_LEVEL > 10)
				cout << "LINK SUCCESSFUL" << endl;
			locSuperLayerSeed1->linked = true;
			locSuperLayerSeed2->linked = true;

			//create new group of stereo seeds (assumes linking to stereo: cannot axial (would skip too many))
			if(locOuterSuperLayer < 4)
				locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.push_back(deque<DCDCSuperLayerSeed*>(1, locSuperLayerSeed2));
			else
				locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.push_back(deque<DCDCSuperLayerSeed*>(1, locSuperLayerSeed2));
		}
	}
}

//------------------------------------
// Link_SuperLayers_FromStereo_ToAxial
//------------------------------------
bool DTrackCandidate_factory_CDC::Link_SuperLayers_FromStereo_ToAxial(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer, unsigned int locInnerSuperLayer)
{
	// Link from a stereo super layer to an axial super layer.  Unfortunately, this is extremely messy.  
		//If you're editing this, I stronly suggest reading the simpler Link_SuperLayers_FromStereo_ToStereo function to make sure you understand it first.  Good luck! :)

	//since linking to axial, will create new track circles (one for each unique combo of axial seeds)
	deque<DCDCTrackCircle*> locNewCDCTrackCircles; //will eventually return this deque (by setting locCDCTrackCircles to it at the end)

	// first save all track circles for output (into locNewCDCTrackCircles) that definitely will NOT be linked from: 
		//those that already have an axial layer greater than the stereo we're linking from
			//this happens when trying to skip a super layer (e.g. a dead hv board), but this circle link was linked successfully
		//do this first, so that if an axial combo is created later that somehow is identical (e.g. due to skipping a super layer due to a dead HV board)
			//the stereo seeds will just be merged with the already existing one
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		DCDCTrackCircle* locCDCTrackCircle = locCDCTrackCircles[loc_i];
		if(!locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
		{
			if(locCDCTrackCircle->dSuperLayerSeeds_Axial.back()->dSuperLayer > locInnerSuperLayer)
			{
				locNewCDCTrackCircles.push_back(locCDCTrackCircle);
				continue;
			}
		}
	}

	// loop over track circles
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		DCDCTrackCircle* locCDCTrackCircle = locCDCTrackCircles[loc_i];

		// if trying to skip a layer, first check to make sure this seed wasn't already linked (to layer "locInnerSuperLayer + 1")
			// axial checked here, stereo checked for each stereo series below
		if(!locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
		{
			// circle already saved in the output (above), so just skip linking here
			if(locCDCTrackCircle->dSuperLayerSeeds_Axial.back()->dSuperLayer > locInnerSuperLayer)
				continue;
		}

		// determine if linking from the stereo seeds stored in the dSuperLayerSeeds_InnerStereo or dSuperLayerSeeds_OuterStereo deques
		bool locFromInnerStereoFlag = (locInnerSuperLayer < 4);
		if(!locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
		{
			if(locCDCTrackCircle->dSuperLayerSeeds_Axial.back()->dSuperLayer != 4)
				locFromInnerStereoFlag = true; //middle axial layer is missing, store stereo layer as inner
		}

		//get a deque containing all of the series's of seeds to loop over
		deque<deque<DCDCSuperLayerSeed*> > locPreLinkedCDCSuperLayerSeeds;
		if(locFromInnerStereoFlag) //if on outer stereo seeds but axial seed 4 is missing, use inner stereo seeds
			locPreLinkedCDCSuperLayerSeeds = locCDCTrackCircle->dSuperLayerSeeds_InnerStereo;
		else //use outer stereo seeds			
			locPreLinkedCDCSuperLayerSeeds = locCDCTrackCircle->dSuperLayerSeeds_OuterStereo; //for looping over

		bool locTrackCircleSavedFlag = false; //make sure that this track circle is saved in the output even if no future link is found
		//loop over each series of stereo seeds
		for(size_t loc_j = 0; loc_j < locPreLinkedCDCSuperLayerSeeds.size(); ++loc_j)
		{
			//get the last DCDCSuperLayerSeed in this stereo seed series
			deque<DCDCSuperLayerSeed*> locStereoSeedSeries = locPreLinkedCDCSuperLayerSeeds[loc_j];
			DCDCSuperLayerSeed* locSuperLayerSeed1 = locPreLinkedCDCSuperLayerSeeds[loc_j].back();
			if(DEBUG_LEVEL > 10)
				cout << "Seed 1 Super Layer, Seed Index = " << locSuperLayerSeed1->dSuperLayer << ", " << locSuperLayerSeed1->dSeedIndex << endl;

			bool locInnerSeedLinkSuccessfulFlag = false;
			//don't attempt seed link if wrong super layer, or if the inner seed was a definite spiral turn that was turning back inwards, etc. 
			if((locSuperLayerSeed1->dSuperLayer == locInnerSuperLayer) && Check_IfShouldAttemptLink(locSuperLayerSeed1, true))
			{
				//attempt to link to the outer superlayer
				for(size_t loc_k = 0; loc_k < dSuperLayerSeeds[locOuterSuperLayer - 1].size(); ++loc_k)
				{
					DCDCSuperLayerSeed* locSuperLayerSeed2 = dSuperLayerSeeds[locOuterSuperLayer - 1][loc_k];
					if(DEBUG_LEVEL > 10)
						cout << "Seed 2 Super Layer, Seed Index = " << locSuperLayerSeed2->dSuperLayer << ", " << locSuperLayerSeed2->dSeedIndex << endl;
					if(!Check_IfShouldAttemptLink(locSuperLayerSeed2, false))
						continue; //don't attempt seed link if the outer seed was a definite spiral turn that was turning back outwards, etc. 
					if(DEBUG_LEVEL > 10)
						cout << "Attempting Seed Link" << endl;
					if(!Attempt_SeedLink(locSuperLayerSeed1, locSuperLayerSeed2)) //see if seeds are nearby enough to link
						continue;
					//LINK SUCCESSFUL!!
					if(DEBUG_LEVEL > 10)
						cout << "LINK SUCCESSFUL" << endl;
					locSuperLayerSeed1->linked = true;
					locSuperLayerSeed2->linked = true;
					locInnerSeedLinkSuccessfulFlag = true;

					//make new deque of axial seeds
					deque<DCDCSuperLayerSeed*> locNewCDCSuperLayerSeeds_Axial = locCDCTrackCircle->dSuperLayerSeeds_Axial;
					locNewCDCSuperLayerSeeds_Axial.push_back(locSuperLayerSeed2);

					//see if should create new DCDCTrackCircle object (axial combo is unique) or use a (probably recently created) existing one (it is not unique)
					DCDCTrackCircle* locNewCDCTrackCircle = NULL;
					for(size_t loc_l = 0; loc_l < locNewCDCTrackCircles.size(); ++loc_l)
					{
						if(locNewCDCSuperLayerSeeds_Axial != locNewCDCTrackCircles[loc_l]->dSuperLayerSeeds_Axial)
							continue;
						//all axial layers are identical: use previous object
						locNewCDCTrackCircle = locNewCDCTrackCircles[loc_l];
					}
					if(locNewCDCTrackCircle == NULL)
					{
						//new combination of axial layers, create new object
						if(locNewCDCTrackCircles.size() >= MAX_ALLOWED_TRACK_CIRCLES)
						{
							if(DEBUG_LEVEL > 10)
								cout << "Too many track circles; bailing" << endl;
							locCDCTrackCircles.clear();
							return false;
						}
						locNewCDCTrackCircle = Get_Resource_CDCTrackCircle();
						locNewCDCTrackCircle->dSuperLayerSeeds_Axial = locNewCDCSuperLayerSeeds_Axial;
						if(!locFromInnerStereoFlag) //if on outer stereo, keep inner stereo results (but not outer stereo!!: will save below)
							locNewCDCTrackCircle->dSuperLayerSeeds_InnerStereo = locCDCTrackCircle->dSuperLayerSeeds_InnerStereo;
						locNewCDCTrackCircles.push_back(locNewCDCTrackCircle); //store new circle for return
						locTrackCircleSavedFlag = true;
					}
					//save this current series of stereo seeds
					if(locFromInnerStereoFlag)
						locNewCDCTrackCircle->dSuperLayerSeeds_InnerStereo.push_back(locStereoSeedSeries);
					else
						locNewCDCTrackCircle->dSuperLayerSeeds_OuterStereo.push_back(locStereoSeedSeries);
				}
			}
			if(!locInnerSeedLinkSuccessfulFlag)
			{
				//this inner seed series was not linked to an axial seed: need to make sure this group is still saved in the "new" deque
				//see if should create new DCDCTrackCircle object (axial combo is unique) or use a (probably recently created) existing one (it is not unique)
				DCDCTrackCircle* locNewCDCTrackCircle = NULL;
				for(size_t loc_l = 0; loc_l < locNewCDCTrackCircles.size(); ++loc_l)
				{
					if(locCDCTrackCircle->dSuperLayerSeeds_Axial != locNewCDCTrackCircles[loc_l]->dSuperLayerSeeds_Axial)
						continue;
					//all axial layers are identical: use previous object
					locNewCDCTrackCircle = locNewCDCTrackCircles[loc_l];
				}
				if(locNewCDCTrackCircle == NULL)
				{
					//new combination of axial layers, create new object
					if(locNewCDCTrackCircles.size() >= MAX_ALLOWED_TRACK_CIRCLES)
					{
						if(DEBUG_LEVEL > 10)
							cout << "Too many track circles; bailing" << endl;
						locCDCTrackCircles.clear();
						return false;
					}
					locNewCDCTrackCircle = Get_Resource_CDCTrackCircle();
					locNewCDCTrackCircle->dSuperLayerSeeds_Axial = locCDCTrackCircle->dSuperLayerSeeds_Axial;
					if(!locFromInnerStereoFlag) //if on outer stereo, keep inner stereo results (but not outer stereo!!: will save below)
						locNewCDCTrackCircle->dSuperLayerSeeds_InnerStereo = locCDCTrackCircle->dSuperLayerSeeds_InnerStereo;
					locNewCDCTrackCircles.push_back(locNewCDCTrackCircle); //store new circle for return
					locTrackCircleSavedFlag = true;
				}
				//save this current series of stereo seeds
				if(locFromInnerStereoFlag)
					locNewCDCTrackCircle->dSuperLayerSeeds_InnerStereo.push_back(locStereoSeedSeries);
				else
					locNewCDCTrackCircle->dSuperLayerSeeds_OuterStereo.push_back(locStereoSeedSeries);
			}
		}
		//if true: no combination of stereo seeds from this track circle was linked as a new object, keep as a unique combination
		if(!locTrackCircleSavedFlag)
			locNewCDCTrackCircles.push_back(locCDCTrackCircle); 
	}

	locCDCTrackCircles = locNewCDCTrackCircles; //"return" "new" track circle objects
	return true;
}

//-------------------------------------
// Link_SuperLayers_FromStereo_ToStereo
//-------------------------------------
void DTrackCandidate_factory_CDC::Link_SuperLayers_FromStereo_ToStereo(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer, unsigned int locInnerSuperLayer)
{
	// Link from a stereo super layer to a stereo super layer. 

	// loop over track circles
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		DCDCTrackCircle* locCDCTrackCircle = locCDCTrackCircles[loc_i];

		// if trying to skip a layer, first check to make sure this seed wasn't already linked (to layer "locInnerSuperLayer + 1")
			// axial checked here, stereo checked for each stereo series below
		if(!locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
		{
			if(locCDCTrackCircle->dSuperLayerSeeds_Axial.back()->dSuperLayer > locInnerSuperLayer)
				continue;
		}

		// determine if linking from the stereo seeds stored in the dSuperLayerSeeds_InnerStereo or dSuperLayerSeeds_OuterStereo deques
		bool locFromInnerStereoFlag = (locInnerSuperLayer < 4);
		if(!locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
		{
			if(locCDCTrackCircle->dSuperLayerSeeds_Axial.back()->dSuperLayer != 4)
				locFromInnerStereoFlag = true; //middle axial layer is missing, store stereo layer as inner
		}

		//get series of seeds to loop over
		deque<deque<DCDCSuperLayerSeed*> > locPreLinkedCDCSuperLayerSeeds;
		if(locFromInnerStereoFlag)
		{
			//if on outer stereo seeds but axial seed 4 is missing, use inner stereo seeds
			locPreLinkedCDCSuperLayerSeeds = locCDCTrackCircle->dSuperLayerSeeds_InnerStereo;
			locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.clear(); //reset the existing combos: will re-fill below with the new combos
		}
		else
		{
			//use outer stereo seeds
			locPreLinkedCDCSuperLayerSeeds = locCDCTrackCircle->dSuperLayerSeeds_OuterStereo; //for looping over
			locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.clear(); //reset the existing combos: will re-fill below with the new combos
		}

		//loop over each series of stereo seeds
		for(size_t loc_j = 0; loc_j < locPreLinkedCDCSuperLayerSeeds.size(); ++loc_j)
		{
			//get the last DCDCSuperLayerSeed in this stereo seed series
			deque<DCDCSuperLayerSeed*> locStereoSeedSeries = locPreLinkedCDCSuperLayerSeeds[loc_j];
			DCDCSuperLayerSeed* locSuperLayerSeed1 = locPreLinkedCDCSuperLayerSeeds[loc_j].back();
			if(DEBUG_LEVEL > 10)
				cout << "Seed 1 Super Layer, Seed Index = " << locSuperLayerSeed1->dSuperLayer << ", " << locSuperLayerSeed1->dSeedIndex << endl;

			bool locInnerSeedLinkSuccessfulFlag = false;
			//don't attempt seed link if wrong super layer, or if the inner seed was a definite spiral turn that was turning back inwards, etc. 
			if((locSuperLayerSeed1->dSuperLayer == locInnerSuperLayer) && Check_IfShouldAttemptLink(locSuperLayerSeed1, true))
			{
				//attempt to link to the outer superlayer
				for(size_t loc_k = 0; loc_k < dSuperLayerSeeds[locOuterSuperLayer - 1].size(); ++loc_k)
				{
					DCDCSuperLayerSeed* locSuperLayerSeed2 = dSuperLayerSeeds[locOuterSuperLayer - 1][loc_k];
					if(DEBUG_LEVEL > 10)
						cout << "Seed 2 Super Layer, Seed Index = " << locSuperLayerSeed2->dSuperLayer << ", " << locSuperLayerSeed2->dSeedIndex << endl;
					if(!Check_IfShouldAttemptLink(locSuperLayerSeed2, false))
						continue; //don't attempt seed link if the outer seed was a definite spiral turn that was turning back outwards, etc. 
					if(DEBUG_LEVEL > 10)
						cout << "Attempting Seed Link" << endl;
					if(!Attempt_SeedLink(locSuperLayerSeed1, locSuperLayerSeed2)) //see if seeds are nearby enough to link
						continue;
					//LINK SUCCESSFUL!!
					if(DEBUG_LEVEL > 10)
						cout << "LINK SUCCESSFUL" << endl;
					locSuperLayerSeed1->linked = true;
					locSuperLayerSeed2->linked = true;
					locInnerSeedLinkSuccessfulFlag = true;

					//create deque of new stereo seed series
					deque<DCDCSuperLayerSeed*> locNewStereoSeedSeries = locStereoSeedSeries;
					locNewStereoSeedSeries.push_back(locSuperLayerSeed2);

					//save new combination of stereo seeds
					if(locFromInnerStereoFlag)
						locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.push_back(locNewStereoSeedSeries);
					else //outer and middle axial layer isn't missing
						locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.push_back(locNewStereoSeedSeries);
				}
			}
			if(!locInnerSeedLinkSuccessfulFlag)
			{
				//failed to match, but keep seed series
				if(locFromInnerStereoFlag)
					locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.push_back(locStereoSeedSeries);
				else //outer and middle axial layer isn't missing
					locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.push_back(locStereoSeedSeries);
			}
		}
	}
}

//--------------------------
// Check_IfShouldAttemptLink
//--------------------------
bool DTrackCandidate_factory_CDC::Check_IfShouldAttemptLink(const DCDCSuperLayerSeed* locSuperLayerSeed, bool locInnerSeedFlag)
{
	// don't attempt linking if a definite spiral turn on one of the seeds is turning in the wrong direction
	if(locSuperLayerSeed->dSpiralLinkParams.empty())
		return true; //should definitely attempt it if not a spiral turn

	//get spiral link information
	const map<int, DSpiralParams_t>& locSpiralLinkParams = locSuperLayerSeed->dSpiralLinkParams;
	bool locSelfLinkFlag = (locSpiralLinkParams.find(locSuperLayerSeed->dSeedIndex) != locSpiralLinkParams.end()); //true if spiral-link is within itself

	// determine if turning inwards, outwards, or indeterminate
	bool locIsTurningOutwardsFlag = false;
	bool locIsTurningInwardsFlag = false;
	map<int, DSpiralParams_t>::const_iterator locIterator;
	for(locIterator = locSuperLayerSeed->dSpiralLinkParams.begin(); locIterator != locSuperLayerSeed->dSpiralLinkParams.end(); ++locIterator)
	{
		if(locIterator->second.dDefiniteSpiralTurnRingFlag == -1)
			locIsTurningOutwardsFlag = true;
		else if(locIterator->second.dDefiniteSpiralTurnRingFlag == 1)
			locIsTurningInwardsFlag = true;
	}
	if(locIsTurningOutwardsFlag == locIsTurningInwardsFlag)
		return true; //either not a definite spiral turn (both are false), or indeterminate turning direction (both are true)

	//check spiral turn direction
	if(locInnerSeedFlag)
	{
		if(!locSelfLinkFlag)
		{
			if(locIsTurningInwardsFlag)
			{
				//spiral turn not linked to itself, inner seed is turning inwards : don't link to next super layer
				if(DEBUG_LEVEL > 10)
					cout << "Seed1 part of spiral turn in different direction, don't link it here." << endl;
				return false;
			}
		}
		else if(locIsTurningOutwardsFlag)
		{
			//spiral turn linked to itself (single ring), on ring in inner-half of the inner super layer used for matching : don't link to next super layer
			if(DEBUG_LEVEL > 10)
				cout << "Seed1 part of spiral turn in different direction, don't link it here." << endl;
			return false;
		}
	}
	else
	{
		if(!locSelfLinkFlag)
		{
			if(locIsTurningOutwardsFlag)
			{
				//spiral turn not linked to itself, outer seed is turning outwards : don't link to previous super layer
				if(DEBUG_LEVEL > 10)
					cout << "Seed2 part of spiral turn in different direction, don't link it here." << endl;
				return false;
			}
		}
		else if(locIsTurningInwardsFlag)
		{
			//spiral turn linked to itself (single ring), on ring in outer-half of the outer super layer used for matching : don't link to previous super layer
			if(DEBUG_LEVEL > 10)
				cout << "Seed2 part of spiral turn in different direction, don't link it here." << endl;
			return false;
		}
	}
	return true;
}

//-----------------
// Attempt_SeedLink
//-----------------
bool DTrackCandidate_factory_CDC::Attempt_SeedLink(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2)
{
	//locSuperLayerSeed1 should be the inner of the two
	wire_direction_t locWireDirection1 = locSuperLayerSeed1->dWireOrientation;
	DCDCRingSeed& locRingSeed1 = locSuperLayerSeed1->dCDCRingSeeds.back(); //largest ring of inner seed selected

	//locSuperLayerSeed2 should be the outer of the two
	wire_direction_t locWireDirection2 = locSuperLayerSeed2->dWireOrientation;
	DCDCRingSeed& locRingSeed2 = locSuperLayerSeed2->dCDCRingSeeds.front(); //smallest ring of outer seed selected

	return Attempt_SeedLink(locRingSeed1, locRingSeed2, locWireDirection1, locWireDirection2);
}

//-----------------
// Attempt_SeedLink
//-----------------
bool DTrackCandidate_factory_CDC::Attempt_SeedLink(DCDCRingSeed& locRingSeed1, DCDCRingSeed& locRingSeed2, wire_direction_t locWireDirection1, wire_direction_t locWireDirection2)
{
	//attempt seed link between ring-seeds
		//should only be called when linking super layers together, and when attempting to link to unused hits
	deque<DCDCTrkHit*> &hits1 = locRingSeed1.hits; 
	if(hits1.empty())
		return false;

	deque<DCDCTrkHit*> &hits2 = locRingSeed2.hits; 
	if(hits2.empty())
		return false;

	const DCDCWire* wire1 = hits1[0]->hit->wire;
	const DCDCWire* wire2 = hits2[0]->hit->wire;

	// Determine the minimum distance between the two sets of hits
	double locMinDist2 = MinDist2(locRingSeed1, locRingSeed2);

	// Determine the maximum-allowed transverse distance for linking: is dependent on the orientation of the wires
		//if axial, distance is MAX_HIT_DIST
		//if stereo, distance is MAX_HIT_DIST + 1/2 of the length of the projection of the straw onto the X-Y plane
			//why 1/2? because the wire position (origin) is reported at the midpoint of the straw
	double locMaxDist2;
	if((locWireDirection1 == WIRE_DIRECTION_AXIAL) && (locWireDirection2 == WIRE_DIRECTION_AXIAL))
		locMaxDist2 = MAX_HIT_DIST2;
	else if((locWireDirection1 == WIRE_DIRECTION_AXIAL) && (locWireDirection2 != WIRE_DIRECTION_AXIAL))
	{
		locMaxDist2 = MAX_HIT_DIST + 0.5*wire2->L*fabs(sin(wire2->stereo));
		locMaxDist2 *= locMaxDist2;
	}
	else if((locWireDirection1 != WIRE_DIRECTION_AXIAL) && (locWireDirection2 == WIRE_DIRECTION_AXIAL))
	{
		locMaxDist2 = MAX_HIT_DIST + 0.5*wire1->L*fabs(sin(wire1->stereo));
		locMaxDist2 *= locMaxDist2;
	}
	else if((locWireDirection1 == WIRE_DIRECTION_STEREORIGHT) && (locWireDirection2 == WIRE_DIRECTION_STEREOLEFT))
	{
		locMaxDist2 = MAX_HIT_DIST + 0.5*wire1->L*fabs(sin(wire1->stereo)) + 0.5*wire2->L*fabs(sin(wire2->stereo));
		locMaxDist2 *= locMaxDist2;
	}
	else if((locWireDirection1 == WIRE_DIRECTION_STEREOLEFT) && (locWireDirection2 == WIRE_DIRECTION_STEREORIGHT))
	{
		locMaxDist2 = MAX_HIT_DIST + 0.5*wire1->L*fabs(sin(wire1->stereo)) + 0.5*wire2->L*fabs(sin(wire2->stereo));
		locMaxDist2 *= locMaxDist2;
	}
	else //both stereo left or right
		locMaxDist2 = MAX_HIT_DIST2;

	double dr = wire2->origin.Perp() - wire1->origin.Perp();
	if(DEBUG_LEVEL > 20)
		cout << "locMinDist2, locMaxDist2, dr = " << locMinDist2 << ", " << locMaxDist2 << ", " << dr << endl;

	//calculate transverse distance, return whether or not it is too large
	double locTransverseDist2 = locMinDist2 - dr*dr;
	return (locTransverseDist2 < locMaxDist2);
}

//-------------------
// Print_TrackCircles
//-------------------
void DTrackCandidate_factory_CDC::Print_TrackCircles(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	cout << "Track Circle Super Layer Seeds (Axial, Inner/Outer Stereo):" << endl;
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		cout << "Track Circle Index = " << loc_i << endl;
		Print_TrackCircle(locCDCTrackCircles[loc_i]);
	}
}

//------------------
// Print_TrackCircle
//------------------
void DTrackCandidate_factory_CDC::Print_TrackCircle(DCDCTrackCircle* locCDCTrackCircle)
{
	for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_Axial.size(); ++loc_j)
		cout << "Axial Super Layer, Seed Index = " << locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_j]->dSuperLayer << ", " << locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_j]->dSeedIndex << endl;
	for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.size(); ++loc_j)
	{
		cout << "Inner Stereo Series " << loc_j << ":" << endl;
		for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_j].size(); ++loc_k)
			cout << "Super Layer, Seed Index = " << locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_j][loc_k]->dSuperLayer << ", " << locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_j][loc_k]->dSeedIndex << endl;
	}
	for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.size(); ++loc_j)
	{
		cout << "Outer Stereo Series " << loc_j << ":" << endl;
		for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_j].size(); ++loc_k)
			cout << "Super Layer, Seed Index = " << locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_j][loc_k]->dSuperLayer << ", " << locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_j][loc_k]->dSeedIndex << endl;
	}
	const DHelicalFit* locFit = locCDCTrackCircle->fit;
	if(locFit != NULL)
	{
		cout << "Fit h, x0, y0, r0, phi = " << locFit->h << ", " << locFit->x0 << ", " << locFit->y0 << ", " << locFit->r0 << ", " << locFit->phi << endl;
		cout << "Fit Weighted chisq/ndf = " << locCDCTrackCircle->dWeightedChiSqPerDF << endl;
		cout << "Stereo Weighted chisq/ndf = " << locCDCTrackCircle->dWeightedChiSqPerDF_Stereo << endl;
	}
}

//--------------------------
// Reject_DefiniteSpiralArms
//--------------------------
void DTrackCandidate_factory_CDC::Reject_DefiniteSpiralArms(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	//disregard all DCDCTrackCircle objects where the innermost super layer is definitely a spiral turn
		//e.g. the track has already turned backed inwards towards the target, but then turns back outward again towards the BCAL
	deque<DCDCTrackCircle*>::iterator locIterator;
	for(locIterator = locCDCTrackCircles.begin(); locIterator != locCDCTrackCircles.end();)
	{
		DCDCTrackCircle* locCDCTrackCircle = *locIterator;

		// get the innermost super layer seed
		DCDCSuperLayerSeed* locInnermostSuperLayerSeed = NULL;
		if(!locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
			locInnermostSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_Axial.front();
		if(!locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.empty())
		{
			DCDCSuperLayerSeed* locTempSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[0].front();
			if(locInnermostSuperLayerSeed == NULL)
				locInnermostSuperLayerSeed = locTempSuperLayerSeed;
			else if(locTempSuperLayerSeed->dSuperLayer < locInnermostSuperLayerSeed->dSuperLayer)
				locInnermostSuperLayerSeed = locTempSuperLayerSeed;
		}
		else if(!locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.empty())
		{
			DCDCSuperLayerSeed* locTempSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[0].front();
			if(locInnermostSuperLayerSeed == NULL)
				locInnermostSuperLayerSeed = locTempSuperLayerSeed;
			else if(locTempSuperLayerSeed->dSuperLayer < locInnermostSuperLayerSeed->dSuperLayer)
				locInnermostSuperLayerSeed = locTempSuperLayerSeed;
		}
		if (locInnermostSuperLayerSeed == NULL)
			continue; //is impossible, but clears the warning from the static analyzer

		//loop over spiral links, see if one of them is a definite spiral link
		bool locIsDefinitelyTurningFlag = false;
		map<int, DSpiralParams_t>::iterator locSpiralIterator;
		for(locSpiralIterator = locInnermostSuperLayerSeed->dSpiralLinkParams.begin(); locSpiralIterator != locInnermostSuperLayerSeed->dSpiralLinkParams.end(); ++locSpiralIterator)
		{
			if(locSpiralIterator->second.dDefiniteSpiralTurnRingFlag == 0)
				continue;
			locIsDefinitelyTurningFlag = true;
			break;
		}

		//if definitely a spiral turn in its innermost super layer: then delete the track circle: does not originate from the target
		if(locIsDefinitelyTurningFlag)
		{
			Recycle_DCDCTrackCircle(locCDCTrackCircle); //recycle
		   locIterator = locCDCTrackCircles.erase(locIterator);
		}
		else
			++locIterator;
	}
}

//----------------------
// Drop_IncompleteGroups
//----------------------
void DTrackCandidate_factory_CDC::Drop_IncompleteGroups(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	// Only interested in groups that contain either: 
		// At least super layer 1 (forward going particles)
		// Or at least one stereo layer and one axial layer 
	deque<DCDCTrackCircle*>::iterator locIterator;
	for(locIterator = locCDCTrackCircles.begin(); locIterator != locCDCTrackCircles.end();)
	{
		DCDCTrackCircle* locCDCTrackCircle = *locIterator;
		if(locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
		{
			Recycle_DCDCTrackCircle(locCDCTrackCircle); //recycle
		   locIterator = locCDCTrackCircles.erase(locIterator); //no axial super layers
			continue;
		}

		if(locCDCTrackCircle->dSuperLayerSeeds_Axial.front()->dSuperLayer == 1)
		{
			++locIterator; //has super layer 1: group is OK
			continue;
		}

		if(locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.empty() && locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.empty())
		{
			Recycle_DCDCTrackCircle(locCDCTrackCircle); //recycle
		   locIterator = locCDCTrackCircles.erase(locIterator); //no stereo super layers
		}
		else
			++locIterator; //has at least one axial & one stereo super layer: group is OK
	}
}

//------------
// Fit_Circles
//------------
void DTrackCandidate_factory_CDC::Fit_Circles(deque<DCDCTrackCircle*>& locCDCTrackCircles, bool locFitOnlyIfNullFitFlag, bool locAddStereoLayerIntersectionsFlag, bool locFitDuringLinkingFlag)
{
	/// Do a quick fit of the 2-D projection of the axial hits in the seed to a circle
		/// Include intersection between stereo layer seeds if locAddStereoLayerIntersectionsFlag = true (assumes only one stereo seed series for each inner/outer)
	/// Determine the sign of the charge (and correspondingly the initial phi angle)
	/// assuming that the majority of hits come from the outgoing part of the track.

	/// If the resulting circle passes within MAX_HIT_DIST the majority of the hits, 
	/// then the fit was a success. Otherwise it is a failure and the DCDCTrackCircle is discarded. 

	/// locFitOnlyIfNullFitFlag should be false unless just truncated the input set of circles
		/// if true, this will skip fitting circles that have DCDCTrackCircle::fit != NULL (i.e. were not truncated)

	double locAxialStrawVariance = 0.214401; //[d/sqrt(12)]^2, d = 1.604 = straw diameter
	deque<DCDCTrackCircle*>::iterator locIterator;
	for(locIterator = locCDCTrackCircles.begin(); locIterator != locCDCTrackCircles.end();)
	{
		DCDCTrackCircle* locCDCTrackCircle = *locIterator;
		size_t locNumAxialSuperLayers = locCDCTrackCircle->dSuperLayerSeeds_Axial.size();

		if(locFitDuringLinkingFlag && (locNumAxialSuperLayers < 2))
		{
			//just trying to reduce number of circles during super layer linking; don't try to fit & don't reject yet
			++locIterator;
			continue;
		}

		if(locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
		{
			//no axial hits: cannot fit circle: reject
			Recycle_DCDCTrackCircle(locCDCTrackCircle); //recycle
			locIterator = locCDCTrackCircles.erase(locIterator);
			continue;
		}

		if(locFitOnlyIfNullFitFlag && (locCDCTrackCircle->fit != NULL))
		{
			++locIterator;
			continue; //fit is not null: & locFitOnlyIfNullFitFlag is true: circle not truncated, don't refit
		}

		// Setup the fit (add hits)
		DHelicalFit* locFit = Get_Resource_HelicalFit();
		unsigned int locNumHitsInFit = 0;
		for(size_t loc_i = 0; loc_i < locNumAxialSuperLayers; ++loc_i)
		{
			deque<DCDCTrkHit*> hits;
			locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_i]->Get_Hits(hits);
			for(size_t k = 0; k < hits.size(); ++k)
			{
				if(hits[k]->flags & OUT_OF_TIME)
					continue;
				DVector3 pos = hits[k]->hit->wire->origin;
				locFit->AddHitXYZ(pos.x(), pos.y(), pos.z(), locAxialStrawVariance, locAxialStrawVariance, 0.0);
				++locNumHitsInFit;
			}
		}

		//add intersections between stereo super layers if desired
			//only uses the first combination of stereo super layers (series)
		if(locAddStereoLayerIntersectionsFlag)
		{
			DCDCSuperLayerSeed* locInnerSuperLayerSeed = locCDCTrackCircle->Get_SuperLayerSeed(2);
			DCDCSuperLayerSeed* locOuterSuperLayerSeed = locCDCTrackCircle->Get_SuperLayerSeed(3);
			if((locInnerSuperLayerSeed != NULL) && (locOuterSuperLayerSeed != NULL))
			{
				// Add intersection between super layers 2 & 3
				DVector3 locIntersection = Find_IntersectionBetweenSuperLayers(locInnerSuperLayerSeed, locOuterSuperLayerSeed);
				//are these the correct uncertainties?
				locFit->AddHitXYZ(locIntersection.x(), locIntersection.y(), locIntersection.z(), locAxialStrawVariance, locAxialStrawVariance, 0.0);
			}
			locInnerSuperLayerSeed = locCDCTrackCircle->Get_SuperLayerSeed(5);
			locOuterSuperLayerSeed = locCDCTrackCircle->Get_SuperLayerSeed(6);
			if((locInnerSuperLayerSeed != NULL) && (locOuterSuperLayerSeed != NULL))
			{
				// Add intersection between super layers 5 & 6
				DVector3 locIntersection = Find_IntersectionBetweenSuperLayers(locInnerSuperLayerSeed, locOuterSuperLayerSeed);
				//are these the correct uncertainties?
				locFit->AddHitXYZ(locIntersection.x(), locIntersection.y(), locIntersection.z(), locAxialStrawVariance, locAxialStrawVariance, 0.0);
			}
		}

		//place a tighter constraint on the beam center if fewer hits: tracks with detached vertices may not go through the center
		double locBeamRMS = BeamRMS*locNumAxialSuperLayers; //1sigma = 0.5, 1.0, 1.5 //3sigma = 1.5, 3.0, 4.5

		// Perform the fit
		if(locFit->FitCircleRiemann(TARGET_Z, locBeamRMS) != NOERROR)
		{
			if(DEBUG_LEVEL > 3)
				cout << "Riemann fit failed. Attempting regression fit..." << endl;
			if(locFit->FitCircle() != NOERROR)
			{
				if(DEBUG_LEVEL > 3)
					cout << "Regression circle fit failed. Trying straight line." << endl;
				if(locFit->FitCircleStraightTrack() != NOERROR)
				{
					if(DEBUG_LEVEL > 3)
						cout << "Trying straight line fit failed. Giving up." << endl;
					Recycle_DCDCTrackCircle(locCDCTrackCircle); //recycle
				   locIterator = locCDCTrackCircles.erase(locIterator);
					continue;
				}
			}
			else
				locFit->FindSenseOfRotation();
		}
		else
			locFit->GuessChargeFromCircleFit(); // for Riemann fit

		// Check if majority of hits are close to circle. Also calculate the avg drift time for the hits close to the circle.
		double x0 = locFit->x0;
		double y0 = locFit->y0;
		double r0 = locFit->r0;
		size_t locNumHitsCloseToCircle = 0;
		double locAverageDriftTime = 0.0;
		unsigned int locTotalNumHits = 0;
		for(size_t loc_i = 0; loc_i < locNumAxialSuperLayers; ++loc_i)
		{
			deque<DCDCTrkHit*> hits;
			locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_i]->Get_Hits(hits);
			locTotalNumHits += hits.size();
			for(size_t k = 0; k < hits.size(); ++k)
			{
				if(hits[k]->flags & OUT_OF_TIME)
					continue;
				double dx = hits[k]->hit->wire->origin.X() - x0;
				double dy = hits[k]->hit->wire->origin.Y() - y0;
				double d = sqrt(dx*dx + dy*dy);
				if(DEBUG_LEVEL > 15)
					cout << "dist = " << d - r0 << endl;
				if(fabs(d - r0) > MAX_HIT_DIST)
					continue;
				++locNumHitsCloseToCircle;
				locAverageDriftTime += hits[k]->hit->tdrift;
			}
		}
		locAverageDriftTime /= ((double)locNumHitsCloseToCircle);

		if(DEBUG_LEVEL > 3)
			cout << "Circle fit: Nhits=" << locFit->GetNhits() << " h=" << locFit->h << " N=" << locNumHitsCloseToCircle << " phi=" << locFit->phi << endl;
		if(locNumHitsCloseToCircle < MIN_SEED_HITS)
		{
			if(DEBUG_LEVEL > 3)
				cout << "Rejected circle fit due to too few hits on track (N=" << locNumHitsCloseToCircle << " MIN_SEED_HITS=" << MIN_SEED_HITS << ")" << endl;
			Recycle_DCDCTrackCircle(locCDCTrackCircle); //recycle
		   locIterator = locCDCTrackCircles.erase(locIterator);
			continue;
		}
	
		if(locNumHitsCloseToCircle < locTotalNumHits/2)
		{
			if(DEBUG_LEVEL > 3)
				cout << "Rejected circle fit due to minority of hits on track (N=" << locNumHitsCloseToCircle << " locTotalNumHits/2=" << locTotalNumHits/2 << ")" << endl;
			Recycle_DCDCTrackCircle(locCDCTrackCircle); //recycle
		   locIterator = locCDCTrackCircles.erase(locIterator);
			continue;
		}

		// Fit is good, save fit results
		locCDCTrackCircle->fit = locFit;
		double locWeightedChiSqPerDF = ((fabs(locFit->chisq) > 0.0) && (locFit->ndof > 0)) ? locFit->chisq/(float(locFit->ndof*locNumAxialSuperLayers*locNumAxialSuperLayers)) : 9.9E50;
		if(DEBUG_LEVEL > 10)
			cout << "chisq, ndof, numaxial, weightedchisq = " << locFit->chisq << ", " << locFit->ndof << ", " << locNumAxialSuperLayers << ", " << locWeightedChiSqPerDF << endl;
		locCDCTrackCircle->dWeightedChiSqPerDF = locWeightedChiSqPerDF;
		locCDCTrackCircle->dAverageDriftTime = locAverageDriftTime;

		++locIterator;
	}
}

//------------------------------------
// Find_IntersectionBetweenSuperLayers
//------------------------------------
DVector3 DTrackCandidate_factory_CDC::Find_IntersectionBetweenSuperLayers(const DCDCSuperLayerSeed* locInnerSuperLayerSeed, const DCDCSuperLayerSeed* locOuterSuperLayerSeed)
{
	const DCDCRingSeed& locInnerSuperLayerRing = locInnerSuperLayerSeed->dCDCRingSeeds.back(); //last ring of inner super layer
	const DCDCRingSeed& locOuterSuperLayerRing = locOuterSuperLayerSeed->dCDCRingSeeds.front(); //first ring of outer super layer

	const DCDCWire* first_wire = locInnerSuperLayerRing.hits.front()->hit->wire;
	const DCDCWire* second_wire = locOuterSuperLayerRing.hits.front()->hit->wire;

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

	if(DEBUG_LEVEL > 10)
		cout << "XYZ intersection between SL" << locInnerSuperLayerSeed->dSuperLayer << " and SL" << locOuterSuperLayerSeed->dSuperLayer << ": " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << endl;
	return pos;
}

/*********************************************************************************************************************************************************************/
/*************************************************************** Filter Track Circles and Stereo Wires ***************************************************************/
/*********************************************************************************************************************************************************************/

//-----------------------
// Handle_StereoAndFilter
//-----------------------
void DTrackCandidate_factory_CDC::Handle_StereoAndFilter(deque<DCDCTrackCircle*>& locCDCTrackCircles, bool locFinalPassFlag)
{
	// If not on final (refinement) pass: First (potentially) truncate and then filter track circles based on track circle fit results
	if(!locFinalPassFlag)
	{
		Truncate_TrackCircles(locCDCTrackCircles);
		Set_HitBitPattern_Axial(locCDCTrackCircles);
		Filter_TrackCircles_Axial(locCDCTrackCircles);
		if(DEBUG_LEVEL > 5)
		{
			cout << "post filter clone seeds" << endl;
			Print_TrackCircles(locCDCTrackCircles);
		}
	}

	// Create new stereo super layer seeds and select the best ones. 
		//This is done one track at a time to prevent memory spikes (memory is recycled as "new" stereo hits are rejected)
	deque<DCDCTrackCircle*>::iterator locIterator;
	size_t locCircleCounter = 0;
	for(locIterator = locCDCTrackCircles.begin(); locIterator != locCDCTrackCircles.end();)
	{
		if(DEBUG_LEVEL > 5)
			cout << "Create new Super Layer Seeds for Track Circle " << locCircleCounter << endl;
		// Create new super layer seeds: find the intersection between the stereo hits and the circle fit and create new seeds with that information
		Create_NewCDCSuperLayerSeeds(*locIterator);
		if(DEBUG_LEVEL > 5)
			cout << "Select Super Layer Seeds for Track Circle " << locCircleCounter << endl;
		// Select the combination of super layer seeds that give the best determination of theta/z for this track. 
			// Rejects the track if on the final pass and no valid theta/z can be calculated. 
		if(Select_CDCSuperLayerSeeds(*locIterator, locFinalPassFlag))
			++locIterator;
		else
		{
			Recycle_DCDCTrackCircle(*locIterator); //recycle
			locIterator = locCDCTrackCircles.erase(locIterator);
		}
		++locCircleCounter;
	}
	Set_HitBitPattern_All(locCDCTrackCircles);

	// If not on final (refinement) pass: Filter track circles based on stereo results
	if(!locFinalPassFlag)
	{
		if(DEBUG_LEVEL > 5)
		{
			cout << "stereo-selected track circles" << endl;
			Print_TrackCircles(locCDCTrackCircles);
		}
		// Filter out duplicates seeds and definite spiral arms. 
		Filter_TrackCircles_Stereo(locCDCTrackCircles);
		if(DEBUG_LEVEL > 5)
		{
			cout << "stereo-filtered track circles" << endl;
			Print_TrackCircles(locCDCTrackCircles);
		}
	}
}

//------------------------
// Set_HitBitPattern_Axial
//------------------------
void DTrackCandidate_factory_CDC::Set_HitBitPattern_Axial(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	unsigned int locNumBits = 8*sizeof(unsigned int);
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		DCDCTrackCircle* locCDCTrackCircle = locCDCTrackCircles[loc_i];
		locCDCTrackCircle->HitBitPattern.clear();
		locCDCTrackCircle->HitBitPattern.resize(dNumCDCHits/(8*sizeof(unsigned int)) + 1);
		for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_Axial.size(); ++loc_j)
		{
			deque<DCDCTrkHit*> locHits;
			locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_j]->Get_Hits(locHits);
			for(size_t loc_k = 0; loc_k < locHits.size(); ++loc_k)
				locCDCTrackCircle->HitBitPattern[locHits[loc_k]->index/locNumBits] |= 1 << locHits[loc_k]->index % locNumBits;
		}
	}
}

//----------------------
// Truncate_TrackCircles
//----------------------
void DTrackCandidate_factory_CDC::Truncate_TrackCircles(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	if(locCDCTrackCircles.empty())
		return;

	// This routine is designed to save tracks that exit the drift chamber early (e.g. before SL7) when there are other tracks nearby (in phi)
		// The problem is when these tracks get mixed together, the track that goes most-radially-outward in the CDC generally wins
			// e.g. A track not reaching SL7 can match to the other track's SL7, so its fit is worse than it should be
			// And, if the two tracks share many hits (such as SL7), then the track not hitting SL7 tends to get rejected

	// It should almost be impossible for two tracks to have identical DCDCSuperLayerSeed's AND for all hits to belong to both tracks. 
	// Two crossing tracks should almost always have different DCDCSuperLayerSeed's that share hits between them.  
	// If the DCDCSuperLayerSeed truly ought to belong to both tracks, stripping hits won't kill the track anyway: it should still be reconstructable

	// To save these tracks, look for pairs of tracks that have identical axial DCDCSuperLayerSeed's: first in super layer 7
		// If the seeds in SL7 are identical, strip SL6 & SL7 from the track with the worse circle-fit weighted chisq/ndf.  
		// If the seeds in SL4 are also identical in this pair, strip everything but SL1 & SL2 from the track with the worse circle-fit weighted chisq/ndf and refit it.  
	// If the axial super layer seeds in the truncated circle are not unique, merge it with the other existing DCDCTrackCircle
	// If this truncated DCDCTrackCircle is merely a subset of a different DCDCTrackCircle, delete it. 
	// Refit the newly-truncated track circles. 

	// Next look for pairs of tracks that have identical DCDCSuperLayerSeed's in super layer 4
		// If the track with the worse circle-fit weighted chisq/ndf has an SL7, DO NOT truncate it (it is a unique SL7 (else would have been rejected earlier))
		// Otherwise, if the seeds in SL4 are identical, strip SL3+ from the track with the worse circle-fit weighted chisq/ndf and refit it.  
	// If the axial super layer seeds in the truncated circle are not unique, merge it with the other existing DCDCTrackCircle
	// If this truncated DCDCTrackCircle is merely a subset of a different DCDCTrackCircle, delete it. 
	// Refit the newly-truncated track circles. 

	//first initialize "DCDCTrackCircle::dHasNonTruncatedSeedsFlag" variables
		//these are useful for determining whether or not any of the stereo seeds in the track circle are unique, or whether they all came from a truncation
		//this is necessary because after truncating a circle, it may have the same axial super layers as another track circle, with which it will be merged. 
		//You need to know if any of the stereo combinations in a track circle are unique when determining whether the truncation result is merely a subset of another track
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		if(!locCDCTrackCircles[loc_i]->dSuperLayerSeeds_InnerStereo.empty())
			locCDCTrackCircles[loc_i]->dHasNonTruncatedSeedsFlag_InnerStereo = true;
		if(!locCDCTrackCircles[loc_i]->dSuperLayerSeeds_OuterStereo.empty())
			locCDCTrackCircles[loc_i]->dHasNonTruncatedSeedsFlag_OuterStereo = true;
		else if(locCDCTrackCircles[loc_i]->Get_LastSuperLayerSeed()->dSuperLayer > 4)
			locCDCTrackCircles[loc_i]->dHasNonTruncatedSeedsFlag_OuterStereo = true; //e.g. SL4 is missing, so outers are grouped with inners
	}

	//assumes input circles are sorted, with largest chisq/ndf first and smallest last
		//want to compare the worst fit to the best fit
	deque<DCDCTrackCircle*>::iterator locIterator_Validating, locIterator_ToCompareTo;

	// FIRST CHECK FOR SAME SUPER LAYER SEED IN SL7
	bool locTruncationPerformedFlag = false;
	for(locIterator_Validating = locCDCTrackCircles.begin(); locIterator_Validating != locCDCTrackCircles.end(); ++locIterator_Validating)
	{
		DCDCTrackCircle* locCDCTrackCircle_Validating = *locIterator_Validating; //this has the worst circle-fit weighted chisq/ndf (and decreasing)
		if(DEBUG_LEVEL > 10)
		{
			cout << "validating: search for SL7 truncation for circle with axial seeds:" << endl;
			for(size_t loc_j = 0; loc_j < locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial.size(); ++loc_j)
				cout << "Axial Super Layer, Seed Index = " << locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial[loc_j]->dSuperLayer << ", " << locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial[loc_j]->dSeedIndex << endl;
		}

		DCDCSuperLayerSeed* locSuperLayerSeed_Validating = locCDCTrackCircle_Validating->Get_SuperLayerSeed(7);
		if(locSuperLayerSeed_Validating == NULL)
			continue; //this track circle has no SL7

		for(locIterator_ToCompareTo = --(locCDCTrackCircles.end()); locIterator_ToCompareTo != locIterator_Validating; --locIterator_ToCompareTo)
		{
			DCDCTrackCircle* locCDCTrackCircle_ToCompareTo = *locIterator_ToCompareTo; //this has the best circle-fit weighted chisq/ndf (and increasing)
			if(DEBUG_LEVEL > 10)
			{
				cout << "to-compare-to: search for SL7 truncation for circle with axial seeds:" << endl;
				for(size_t loc_j = 0; loc_j < locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial.size(); ++loc_j)
					cout << "Axial Super Layer, Seed Index = " << locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial[loc_j]->dSuperLayer << ", " << locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial[loc_j]->dSeedIndex << endl;
			}

			DCDCSuperLayerSeed* locSuperLayerSeed_ToCompareTo = locCDCTrackCircle_ToCompareTo->Get_SuperLayerSeed(7);
			if(locSuperLayerSeed_Validating != locSuperLayerSeed_ToCompareTo)
				continue; //different seed for SL7

			//track circles have identical seeds in SL7: should almost never be possible for both to be good tracks AND for all of the hits to belong to both
				// if two tracks are crossing they should almost always have different super layer seeds (that share hits)
				// if they are somehow both good tracks, stripping hits won't be the end of the world: track should still be reconstructable with a subset of them

			//check to see if SL4 is identical also
			locSuperLayerSeed_Validating = locCDCTrackCircle_Validating->Get_SuperLayerSeed(4);
			locSuperLayerSeed_ToCompareTo = locCDCTrackCircle_ToCompareTo->Get_SuperLayerSeed(4);
			if(locSuperLayerSeed_Validating == locSuperLayerSeed_ToCompareTo)
			{
				// SL4 is identical (or missing from both): strip everything from (and including) SL3 and outwards from the track with lower chisq/ndf
					//don't trust SL3 since it led to SL4: hit selector can (in theory) pick up the hits later if needed
				//impossible for SL1 to also be identical: would be same object because all axials would be the same
				if(DEBUG_LEVEL > 5)
					cout << "SL7's and SL4's identical: truncate to SL2" << endl;
				locCDCTrackCircle_Validating->Truncate_Circle(2); //2: new last super layer
				if(DEBUG_LEVEL > 20)
				{
					cout << "post truncate" << endl;
					Print_TrackCircle(locCDCTrackCircle_Validating);
				}
			}
			else
			{
				// SL4 is different: strip SL6 & SL7 from the track with lower chisq/ndf
					//don't trust SL6 since it led to SL7: hit selector can (in theory) pick up the hits later if needed
				if(DEBUG_LEVEL > 5)
					cout << "SL7's identical (but not SL4's): truncate to SL5" << endl;
				locCDCTrackCircle_Validating->Truncate_Circle(5); //5: new last super layer
				if(DEBUG_LEVEL > 20)
				{
					cout << "post truncate" << endl;
					Print_TrackCircle(locCDCTrackCircle_Validating);
				}
			}

			//reset fit
			dHelicalFitPool_Available.push_back(locCDCTrackCircle_Validating->fit); //will redo the fit: recycle the memory
			locCDCTrackCircle_Validating->fit = NULL; //indicate need to reperform fit

			//update truncation sources
			locCDCTrackCircle_Validating->dTruncationSourceCircles.push_back(locCDCTrackCircle_ToCompareTo);

			locTruncationPerformedFlag = true;
			break; //truncation successful
		}
	}

	if(locTruncationPerformedFlag)
	{
		//now merge any track circles that have identical axial super layers
			//e.g. two tracks have SL1 & SL4 but different SL7, and their SL7's are rejected by other tracks who have the same SL7s
			//loop order doesn't really matter here, keeping consistent anyway
		for(locIterator_Validating = locCDCTrackCircles.begin(); locIterator_Validating != locCDCTrackCircles.end();)
		{
			DCDCTrackCircle* locCDCTrackCircle_Validating = *locIterator_Validating;
			bool locMergedTrackCircleFlag = false;
			for(locIterator_ToCompareTo = --(locCDCTrackCircles.end()); locIterator_ToCompareTo != locIterator_Validating; --locIterator_ToCompareTo)
			{
				DCDCTrackCircle* locCDCTrackCircle_ToCompareTo = *locIterator_ToCompareTo;
				if(locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial != locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial)
					continue; //not identical
				//identical axial super layers: merge track circles (absorb the "Validating" one into the "ToCompareTo" one //which is which doesn't matter)
				if(DEBUG_LEVEL > 20)
				{
					cout << "sl7 merging circles: validating = " << endl;
					Print_TrackCircle(locCDCTrackCircle_Validating);
					cout << "sl7 merging circles: to-compare-to = " << endl;
					Print_TrackCircle(locCDCTrackCircle_ToCompareTo);
				}
				locCDCTrackCircle_ToCompareTo->Absorb_TrackCircle(locCDCTrackCircle_Validating);
				if(DEBUG_LEVEL > 20)
				{
					cout << "sl7 merge circles: output = " << endl;
					Print_TrackCircle(locCDCTrackCircle_ToCompareTo);
				}
				locMergedTrackCircleFlag = true;
				break;
			}
			if(locMergedTrackCircleFlag)
			{
				Recycle_DCDCTrackCircle(locCDCTrackCircle_Validating); //recycle
				locIterator_Validating = locCDCTrackCircles.erase(locIterator_Validating);
			}
			else
				++locIterator_Validating;
		}

		//now check to see if any of the newly-truncated track circles is merely a subset of a different circle (regardless of which chisq/ndf is lower)
			//is a subset if the axials are a subset AND the outermost stereo layers dHasNonTruncatedSeedsFlag is false
		//if one is a subset of another, delete it
		for(locIterator_Validating = locCDCTrackCircles.begin(); locIterator_Validating != locCDCTrackCircles.end();)
		{
			DCDCTrackCircle* locCDCTrackCircle_Validating = *locIterator_Validating;
			if(locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial.size() == 3)
			{
				++locIterator_Validating;
				continue; //not going to be a subset
			}
			if(DEBUG_LEVEL > 20)
			{
				cout << "sl7 checking-for-subsets: validating = " << endl;
				Print_TrackCircle(locCDCTrackCircle_Validating);
			}
			bool locRejectTrackFlag = false;
			for(locIterator_ToCompareTo = locCDCTrackCircles.begin(); locIterator_ToCompareTo != locCDCTrackCircles.end(); ++locIterator_ToCompareTo)
			{
				if(locIterator_ToCompareTo == locIterator_Validating)
					continue;
				DCDCTrackCircle* locCDCTrackCircle_ToCompareTo = *locIterator_ToCompareTo;
				if(DEBUG_LEVEL > 20)
				{
					cout << "sl7 checking-for-subsets: to-compare-to = " << endl;
					Print_TrackCircle(locCDCTrackCircle_ToCompareTo);
				}

				if(!locCDCTrackCircle_ToCompareTo->Check_IfInputIsSubset(locCDCTrackCircle_Validating))
					continue; //not a subset

				if(DEBUG_LEVEL > 10)
					cout << "rejecting subset circle" << endl;
				locRejectTrackFlag = true; //is a subset
				break;
			}
			if(locRejectTrackFlag)
			{
				Recycle_DCDCTrackCircle(locCDCTrackCircle_Validating); //recycle
				locIterator_Validating = locCDCTrackCircles.erase(locIterator_Validating);
			}
			else
				++locIterator_Validating;
		}

		//now fit the newly-truncated track circles
		Fit_Circles(locCDCTrackCircles, true, false); //true: fit only truncated circles //false: don't add stereo intersections
		sort(locCDCTrackCircles.begin(), locCDCTrackCircles.end(), CDCSortByChiSqPerNDFDecreasing); //sort by fit chisq/ndf
		if(DEBUG_LEVEL > 5)
		{
			cout << "Post-SL7-turncation track circles" << endl;
			Print_TrackCircles(locCDCTrackCircles);
		}
	}

	// NOW CHECK FOR SAME SUPER LAYER SEED IN SL4
	locTruncationPerformedFlag = false;
	for(locIterator_Validating = locCDCTrackCircles.begin(); locIterator_Validating != locCDCTrackCircles.end();)
	{
		DCDCTrackCircle* locCDCTrackCircle_Validating = *locIterator_Validating;
		if(DEBUG_LEVEL > 10)
		{
			cout << "validating: search for SL4 truncation for circle with axial seeds:" << endl;
			for(size_t loc_j = 0; loc_j < locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial.size(); ++loc_j)
				cout << "Axial Super Layer, Seed Index = " << locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial[loc_j]->dSuperLayer << ", " << locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial[loc_j]->dSeedIndex << endl;
		}

		DCDCSuperLayerSeed* locSuperLayerSeed_Validating = locCDCTrackCircle_Validating->Get_SuperLayerSeed(4);
		if(locSuperLayerSeed_Validating == NULL)
		{
			++locIterator_Validating;
			continue; //this track circle has no SL4
		}

		if(locCDCTrackCircle_Validating->Get_SuperLayerSeed(7) != NULL)
		{
			//this track circle has a unique SL7 (if not unique would have been truncated earlier): do not truncate it
				//rejecting SL4 from the track with the SL7 would leave an unphysical combination
				//this track may have a low chisq/ndf because there is a kink in the track, or maybe the eloss is large and the circle doesn't quite catch SL7 very well
				//if the SL7 is truly bad (e.g. slightly different from a different, true SL7), then the 50% common-hit requirement should kill it (if SL4 is indeed identical)
			++locIterator_Validating;
			continue;
		}

		DCDCSuperLayerSeed* locSuperLayerSeed1_Validating = locCDCTrackCircle_Validating->Get_SuperLayerSeed(1);
		bool locRejectTrackFlag = false;
		for(locIterator_ToCompareTo = --(locCDCTrackCircles.end()); locIterator_ToCompareTo != locIterator_Validating; --locIterator_ToCompareTo)
		{
			DCDCTrackCircle* locCDCTrackCircle_ToCompareTo = *locIterator_ToCompareTo;
			if(DEBUG_LEVEL > 10)
			{
				cout << "to-compare-to: search for SL4 truncation for circle with axial seeds:" << endl;
				for(size_t loc_j = 0; loc_j < locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial.size(); ++loc_j)
					cout << "Axial Super Layer, Seed Index = " << locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial[loc_j]->dSuperLayer << ", " << locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial[loc_j]->dSeedIndex << endl;
			}

			DCDCSuperLayerSeed* locSuperLayerSeed_ToCompareTo = locCDCTrackCircle_ToCompareTo->Get_SuperLayerSeed(4);
			if(locSuperLayerSeed_Validating != locSuperLayerSeed_ToCompareTo)
				continue; //different seed for SL4

			//SL4 is identical, check status of SL1
			DCDCSuperLayerSeed* locSuperLayerSeed1_ToCompareTo = locCDCTrackCircle_ToCompareTo->Get_SuperLayerSeed(1);
			if(locSuperLayerSeed1_Validating == locSuperLayerSeed1_ToCompareTo)
			{
				//SL1 & SL4 are identical: only difference is that locCDCTrackCircle_ToCompareTo has SL7
					//locCDCTrackCircle_Validating is merely a subset of locCDCTrackCircle_ToCompareTo (which has a better chisq/ndf): reject it
				if(DEBUG_LEVEL > 5)
					cout << "SL1 and SL4 are identical, while SL7 of validating is missing: reject validating track" << endl;
				locRejectTrackFlag = true;
				break;
			}

			//track circles have identical seeds in SL4, and at least one track does not have an SL7 (although both may not):
				//should almost never be possible for both to be good tracks AND for all of the hits to belong to both
					// if two tracks are crossing they should almost always have different super layer seeds (that share hits)
					// if they are somehow both good tracks, stripping hits won't be the end of the world: track should still be reconstructable with a subset of them

			//truncate the circle //don't trust SL3 since it led to SL4: hit selector can (in theory) pick up the hits later if needed
			if(DEBUG_LEVEL > 5)
				cout << "SL4 is identical, SL1 is not, while SL7 of validating is missing: truncate to SL2" << endl;
			locCDCTrackCircle_Validating->Truncate_Circle(2); //2: new last super layer

			//reset fit if not already done
			if(locCDCTrackCircle_Validating->fit != NULL)
			{
				dHelicalFitPool_Available.push_back(locCDCTrackCircle_Validating->fit); //will redo the fit: recycle the memory
				locCDCTrackCircle_Validating->fit = NULL; //indicate need to reperform fit
			}

			//update truncation sources
			bool locIsAlreadyTruncationSourceFlag = false;
			for(size_t loc_i = 0; loc_i < locCDCTrackCircle_Validating->dTruncationSourceCircles.size(); ++loc_i)
			{
				if(locCDCTrackCircle_Validating->dTruncationSourceCircles[loc_i] != locCDCTrackCircle_ToCompareTo)
					continue;
				locIsAlreadyTruncationSourceFlag = true;
				break;
			}
			if(!locIsAlreadyTruncationSourceFlag)
				locCDCTrackCircle_Validating->dTruncationSourceCircles.push_back(locCDCTrackCircle_ToCompareTo);

			locTruncationPerformedFlag = true;
			break; //truncation successful
		}

		if(locRejectTrackFlag)
		{
			Recycle_DCDCTrackCircle(locCDCTrackCircle_Validating); //recycle
			locIterator_Validating = locCDCTrackCircles.erase(locIterator_Validating);
		}
		else
			++locIterator_Validating;
	}

	if(locTruncationPerformedFlag)
	{
		//now merge any track circles that have identical axial super layers
		for(locIterator_Validating = locCDCTrackCircles.begin(); locIterator_Validating != locCDCTrackCircles.end();)
		{
			DCDCTrackCircle* locCDCTrackCircle_Validating = *locIterator_Validating;
			bool locMergedTrackCircleFlag = false;
			for(locIterator_ToCompareTo = --(locCDCTrackCircles.end()); locIterator_ToCompareTo != locIterator_Validating; --locIterator_ToCompareTo)
			{
				DCDCTrackCircle* locCDCTrackCircle_ToCompareTo = *locIterator_ToCompareTo;
				if(locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial != locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial)
					continue; //not identical
				//identical axial super layers: merge track circles (absorb the "Validating" one into the "ToCompareTo" one //which is which doesn't matter)
				locCDCTrackCircle_ToCompareTo->Absorb_TrackCircle(locCDCTrackCircle_Validating);
				locMergedTrackCircleFlag = true;
				break;
			}
			if(locMergedTrackCircleFlag)
			{
				Recycle_DCDCTrackCircle(locCDCTrackCircle_Validating); //recycle
				locIterator_Validating = locCDCTrackCircles.erase(locIterator_Validating);
			}
			else
				++locIterator_Validating;
		}

		//now check to see if any of the newly-truncated track circles is merely a subset of a different circle (regardless of which chisq/ndf is lower)
			//is a subset if the axials are a subset AND the outermost stereo layers dHasNonTruncatedSeedsFlag is false
		//if one is a subset of another, delete it
		for(locIterator_Validating = locCDCTrackCircles.begin(); locIterator_Validating != locCDCTrackCircles.end();)
		{
			DCDCTrackCircle* locCDCTrackCircle_Validating = *locIterator_Validating;
			bool locRejectTrackFlag = false;
			if(DEBUG_LEVEL > 20)
			{
				cout << "sl4 checking-for-subsets: validating = " << endl;
				Print_TrackCircle(locCDCTrackCircle_Validating);
			}
			for(locIterator_ToCompareTo = locCDCTrackCircles.begin(); locIterator_ToCompareTo != locCDCTrackCircles.end(); ++locIterator_ToCompareTo)
			{
				if(locIterator_ToCompareTo == locIterator_Validating)
					continue;
				DCDCTrackCircle* locCDCTrackCircle_ToCompareTo = *locIterator_ToCompareTo;
				if(DEBUG_LEVEL > 20)
				{
					cout << "sl4 checking-for-subsets: to-compare-to = " << endl;
					Print_TrackCircle(locCDCTrackCircle_ToCompareTo);
				}

				if(!locCDCTrackCircle_ToCompareTo->Check_IfInputIsSubset(locCDCTrackCircle_Validating))
					continue; //not a subset

				locRejectTrackFlag = true; //is a subset
				if(DEBUG_LEVEL > 10)
					cout << "rejecting subset circle" << endl;
				break;
			}
			if(locRejectTrackFlag)
			{
				Recycle_DCDCTrackCircle(locCDCTrackCircle_Validating); //recycle
				locIterator_Validating = locCDCTrackCircles.erase(locIterator_Validating);
			}
			else
				++locIterator_Validating;
		}

		//now fit the newly-truncated track circles
		Fit_Circles(locCDCTrackCircles, true, false); //true: fit only truncated circles //false: don't add stereo intersections
		sort(locCDCTrackCircles.begin(), locCDCTrackCircles.end(), CDCSortByChiSqPerNDFDecreasing); //sort by fit chisq/ndf
		if(DEBUG_LEVEL > 5)
		{
			cout << "Post-SL4-turncation track circles" << endl;
			Print_TrackCircles(locCDCTrackCircles);
		}
	}
}

//--------------------------
// Filter_TrackCircles_Axial
//--------------------------
void DTrackCandidate_factory_CDC::Filter_TrackCircles_Axial(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	if(locCDCTrackCircles.empty())
		return;

	//assumes input circles are sorted, with largest chisq/ndf first and smallest last
		//want to compare the worst fit to the best fit
	deque<DCDCTrackCircle*>::iterator locIterator_Validating, locIterator_ToCompareTo;

	//FIRST: If circles share > MAX_COMMON_HIT_FRACTION of axial hits, reject circle with larger chisq/ndf
	for(locIterator_Validating = locCDCTrackCircles.begin(); locIterator_Validating != locCDCTrackCircles.end();)
	{
		DCDCTrackCircle* locCDCTrackCircle_Validating = *locIterator_Validating;
		if(locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial.empty())
			continue; //no axial hits to share (somehow)

		size_t locNumHits_Validating = 0;
		deque<DCDCTrkHit*> hits;
		for(size_t loc_i = 0; loc_i < locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial.size(); ++loc_i)
		{
			locCDCTrackCircle_Validating->dSuperLayerSeeds_Axial[loc_i]->Get_Hits(hits);
			locNumHits_Validating += hits.size();
		}

		bool locRejectTrackCircleFlag = false;
		for(locIterator_ToCompareTo = --(locCDCTrackCircles.end()); locIterator_ToCompareTo != locIterator_Validating; --locIterator_ToCompareTo)
		{
			DCDCTrackCircle* locCDCTrackCircle_ToCompareTo = *locIterator_ToCompareTo;
			if(locCDCTrackCircle_ToCompareTo->dSuperLayerSeeds_Axial.empty())
				continue; //no axial hits to share (somehow)

			//If seeds share > MAX_COMMON_HIT_FRACTION of hits, reject seed with larger chisq/ndf
			size_t locNumCommonHits = 0;
			size_t locNumWords = locCDCTrackCircle_Validating->HitBitPattern.size();
			for(size_t loc_i = 0; loc_i < locNumWords; ++loc_i)
				locNumCommonHits += bitcount(locCDCTrackCircle_Validating->HitBitPattern[loc_i] & locCDCTrackCircle_ToCompareTo->HitBitPattern[loc_i]);
			double locHitFraction = double(locNumCommonHits)/double(locNumHits_Validating);

			if(locHitFraction > MAX_COMMON_HIT_FRACTION)
			{
				locRejectTrackCircleFlag = true;
				break;
			}
      }
		if(locRejectTrackCircleFlag)
		{
			//free up some memory by clearing the seed deques via reset
			Recycle_DCDCTrackCircle(locCDCTrackCircle_Validating); //recycle
			locIterator_Validating = locCDCTrackCircles.erase(locIterator_Validating);
		}
		else
			++locIterator_Validating; //track circle is valid (for now)
	}
}

//-----------------------------
// Create_NewCDCSuperLayerSeeds
//-----------------------------
void DTrackCandidate_factory_CDC::Create_NewCDCSuperLayerSeeds(DCDCTrackCircle* locCDCTrackCircle)
{
	// Create new stereo DCDCSuperLayerSeed objects, finding the intersections of each stereo wire with the fit circle

	map<DCDCSuperLayerSeed*, DCDCSuperLayerSeed*> locConvertedSuperLayerSeeds;
	map<DCDCSuperLayerSeed*, DCDCSuperLayerSeed*>::iterator locMapIterator; //map from orig super layer to new super layer
	map<DCDCTrkHit*, DCDCTrkHit*> locProjectedStereoHitMap; //map from orig (super layer, non-circle-projected) stereo hit to circle-projected (hit-z-group) hit

	//inner stereo
	for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i].size(); ++loc_j)
		{
			DCDCSuperLayerSeed* locCDCSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i][loc_j];
			locMapIterator = locConvertedSuperLayerSeeds.find(locCDCSuperLayerSeed);
			if(locMapIterator != locConvertedSuperLayerSeeds.end())
			{
				locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i][loc_j] = locMapIterator->second;
				continue; //already converted
			}
			DCDCSuperLayerSeed* locNewCDCSuperLayerSeed = Create_NewStereoSuperLayerSeed(locCDCSuperLayerSeed, locCDCTrackCircle, locProjectedStereoHitMap);
			locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i][loc_j] = locNewCDCSuperLayerSeed;
			locConvertedSuperLayerSeeds[locCDCSuperLayerSeed] = locNewCDCSuperLayerSeed;
		}
	}

	//outer stereo
	for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_i].size(); ++loc_j)
		{
			DCDCSuperLayerSeed* locCDCSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_i][loc_j];
			locMapIterator = locConvertedSuperLayerSeeds.find(locCDCSuperLayerSeed);
			if(locMapIterator != locConvertedSuperLayerSeeds.end())
			{
				locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_i][loc_j] = locMapIterator->second;
				continue; //already converted
			}
			DCDCSuperLayerSeed* locNewCDCSuperLayerSeed = Create_NewStereoSuperLayerSeed(locCDCSuperLayerSeed, locCDCTrackCircle, locProjectedStereoHitMap);
			locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_i][loc_j] = locNewCDCSuperLayerSeed;
			locConvertedSuperLayerSeeds[locCDCSuperLayerSeed] = locNewCDCSuperLayerSeed;
		}
	}
}

//-------------------------------
// Create_NewStereoSuperLayerSeed
//-------------------------------
DTrackCandidate_factory_CDC::DCDCSuperLayerSeed* DTrackCandidate_factory_CDC::Create_NewStereoSuperLayerSeed(DCDCSuperLayerSeed* locCDCSuperLayerSeed, const DCDCTrackCircle* locCDCTrackCircle, map<DCDCTrkHit*, DCDCTrkHit*>& locProjectedStereoHitMap)
{
	//locProjectedStereoHitMap is map from orig (super layer, non-circle-projected) stereo hit to circle-projected hit
	DCDCSuperLayerSeed* locNewCDCSuperLayerSeed = Get_Resource_CDCSuperLayerSeed();
	*locNewCDCSuperLayerSeed = *locCDCSuperLayerSeed;

	// Project all stereo hits in this super layer onto the circle fit

	for(size_t loc_i = 0; loc_i < locCDCSuperLayerSeed->dCDCRingSeeds.size(); ++loc_i)
	{
		deque<DCDCTrkHit*>& hits = locCDCSuperLayerSeed->dCDCRingSeeds[loc_i].hits;
		deque<DCDCTrkHit*> locProjectedHits;
		for(size_t loc_l = 0; loc_l < hits.size(); ++loc_l)
		{
			if(fabs(hits[loc_l]->hit->tdrift - locCDCTrackCircle->dAverageDriftTime) > MAX_SEED_TIME_DIFF)
				continue; // Ignore hits that are out of time with the group

			// If haven't done so already, Calculate intersection points between circle and stereo wire
			DCDCTrkHit* locProjectedCDCTrkHit = NULL;
			map<DCDCTrkHit*, DCDCTrkHit*>::iterator locHitIterator = locProjectedStereoHitMap.find(hits[loc_l]);
			if(locHitIterator == locProjectedStereoHitMap.end())
			{
				// Clone the hit and set it's stereo-hit-position
				DVector3 locStereoHitPos;
				double locPhiStereo = 0.0, var_z = 9.9E9;
				locProjectedCDCTrkHit = Get_Resource_CDCTrkHit();
				*locProjectedCDCTrkHit = *hits[loc_l];
				//if below is false: wire doesn't intersect the circle: in this case, don't reject hit outright: ignore for theta/z, but let wire-based tracking try to use it
					//e.g., track has a spiral turn in SL6, and since there is no SL7, the circle-fit is extrapolated into SL6 and doesn't quite catch the spiral turn
				if(Calc_StereoPosition(locProjectedCDCTrkHit->hit->wire, locCDCTrackCircle->fit, locStereoHitPos, var_z, locPhiStereo))
					locProjectedCDCTrkHit->dValidStereoHitPosFlag = true; 
				locProjectedCDCTrkHit->dStereoHitPos = locStereoHitPos;
				locProjectedCDCTrkHit->var_z = var_z;
				locProjectedCDCTrkHit->dPhiStereo = locPhiStereo;
				dStereoHitNumUsedMap[locProjectedCDCTrkHit] = 1;
			}
			else
			{
				//hit was already projected onto this circle (e.g. in a different super layer seed), just reuse the results/memory
				locProjectedCDCTrkHit = locHitIterator->second;
				++dStereoHitNumUsedMap[locProjectedCDCTrkHit];
			}

			// Save the hit
			locProjectedHits.push_back(locProjectedCDCTrkHit);
		}
		locNewCDCSuperLayerSeed->dCDCRingSeeds[loc_i].hits = locProjectedHits;
	}

	return locNewCDCSuperLayerSeed;
}

//--------------------
// Calc_StereoPosition
//--------------------
bool DTrackCandidate_factory_CDC::Calc_StereoPosition(const DCDCWire *wire, const DHelicalFit* fit, DVector3 &pos, double &var_z, double& locPhiStereo, double d)
{
	// Calculate intersection point between circle and stereo wire
	DVector3 origin = wire->origin;
	DVector3 dir = (1./wire->udir.z())*wire->udir;
	double dx = origin.x() - fit->x0;
	double dy = origin.y() - fit->y0;
	double ux = dir.x();
	double uy = dir.y();
	double temp1 = ux*ux + uy*uy;
	double temp2 = ux*dy - uy*dx;
	double b = -ux*dx - uy*dy;
	double dr = fit->r0 - d;
	double r0_sq = dr*dr;
	double A = r0_sq*temp1 - temp2*temp2;

	// Check that this wire intersects this circle
	if(A < 0.0)
		return false; // line along wire does not intersect circle, ever.

	// Guess for variance for z: assume straw cell size??
	double temp = 1.6/sin(wire->stereo);
	var_z = temp*temp/12.;
	
	// Calculate intersection points for the two roots 
	double B = sqrt(A);
	double dz1 = (b - B)/temp1;
	double dz2 = (b + B)/temp1;

	if(DEBUG_LEVEL > 15)
		cout<<"dz1="<<dz1<<" dz2="<<dz2<<endl;
	
	// At this point we must decide which value of alpha to use. 
	// For now, we just use the value closest to zero (i.e. closest to
	// the center of the wire).
	double dz = dz1;
	if(fabs(dz2) < fabs(dz1))
		dz = dz2;
		
	// Compute the position for this hit
	pos = origin + dz*dir;
 
	// distance along wire relative to origin
	double s = dz/cos(wire->stereo);
	
	if(DEBUG_LEVEL > 15)
		cout<<"s="<<s<<" ring="<<wire->ring<<" straw="<<wire->straw<<" stereo="<<wire->stereo<<endl;

	// Compute phi for the stereo wire
	DVector2 R(fit->x0, fit->y0);
	locPhiStereo = atan2(pos.Y() - R.Y(), pos.X() - R.X());
	R *= -1.0; // make R point from center of circle to beamline instead of other way around
	locPhiStereo -= R.Phi(); // make angle relative to beamline

	// We want this to go either from 0 to +2pi for positive charge, or 0 to -2pi for negative.
	double phi_hi = fit->h > 0.0 ? +M_TWO_PI : 0.0;
	double phi_lo = fit->h > 0.0 ? 0.0 : -M_TWO_PI;
	while(locPhiStereo < phi_lo)
		locPhiStereo += M_TWO_PI;
	while(locPhiStereo > phi_hi)
		locPhiStereo -= M_TWO_PI;

	return true;
}

//--------------------------
// Select_CDCSuperLayerSeeds
//--------------------------
bool DTrackCandidate_factory_CDC::Select_CDCSuperLayerSeeds(DCDCTrackCircle* locCDCTrackCircle, bool locFinalPassFlag)
{
	// If on initial pass (locFinalPassFlag = false):
		// Calculates the theta/z of the track for each possible combination of stereo super layer seeds. 
		// The combination with the smallest chisq/ndf is selected, and the rest of the super layer seeds are deleted. 
		// If there are no stereo hits, DO NOT reject the track, as Add_UnusedHits hasn't been called yet, and may find some.

	// If on final pass (locFinalPassFlag = true): 
		// Calculate the theta/z of the track using a subset of the remaining stereo hits
		// A subset is used to give the best-possible calculation of theta/z by ignoring hits where the circle-fit may be inaccurate
		// This is because the projection of the stereo hits onto the circle may be bad in this region, giving a bad theta/z result
		// If there are no stereo hits, reject the track. 

	// Select the best combinations of stereo hits and determine theta/z
	double locBestTheta = 0.0, locBestZ = TARGET_Z, locBestChiSqPerNDF = 9.9E99;
	deque<DCDCSuperLayerSeed*> locBestSuperLayerSeeds_Inner;
	deque<DCDCSuperLayerSeed*> locBestSuperLayerSeeds_Outer;
	bool locGoodStereoComboFoundFlag = false;

	if((!locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.empty()) && (!locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.empty()))
	{
		//hits in both the inner & outer stereo super layers
		for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.size(); ++loc_i)
		{
			for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.size(); ++loc_j)
			{
				//get the hits from this inner combination
				deque<DCDCTrkHit*> locComboHits;
				Select_ThetaZStereoHits(locCDCTrackCircle, loc_i, loc_j, locFinalPassFlag, locComboHits);

				//Evaluate theta & z for this combination, returning the chisq/ndf from the fit
				double locTheta = 0.0, locZ = TARGET_Z, locChiSqPerNDF = 9.9E50;
				if(DEBUG_LEVEL > 5)
					cout << "in/out try theta/z, num stereo hits = " << locComboHits.size() << endl;
				if(!Find_ThetaZ(locCDCTrackCircle->fit, locComboHits, locTheta, locZ, locChiSqPerNDF))
					continue; //combo didn't work for some reason, try a different one
				if(!((locChiSqPerNDF < 1.0) || (locChiSqPerNDF > -1.0)))
					continue; // NaN
				locGoodStereoComboFoundFlag = true;
				if(DEBUG_LEVEL > 5)
					cout << "in/out good theta/z: theta, z, chisq, best-chisq = " << locTheta << ", " << locZ << ", " << locChiSqPerNDF << ", " << locBestChiSqPerNDF << endl;
				if(locChiSqPerNDF >= locBestChiSqPerNDF)
					continue;
				// This is the best combination of stereo seeds so far, save the results
				locBestSuperLayerSeeds_Inner.clear();
				for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i].size(); ++loc_k)
					locBestSuperLayerSeeds_Inner.push_back(locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i][loc_k]);
				locBestSuperLayerSeeds_Outer.clear();
				for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_j].size(); ++loc_k)
					locBestSuperLayerSeeds_Outer.push_back(locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_j][loc_k]);
				locBestTheta = locTheta;
				locBestZ = locZ;
				locBestChiSqPerNDF = locChiSqPerNDF;
			}
		}
	}
	else if(!locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.empty())
	{
		//no hits in the outer super layers (or super layer 4 was missing)
		for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.size(); ++loc_i)
		{
			//get the hits from this inner combination
			deque<DCDCTrkHit*> locComboHits;
			Select_ThetaZStereoHits(locCDCTrackCircle, loc_i, -1, locFinalPassFlag, locComboHits);

			//Evaluate theta & z for this combination, returning the chisq/ndf from the fit
			double locTheta = 0.0, locZ = TARGET_Z, locChiSqPerNDF = 9.9E50;
			if(DEBUG_LEVEL > 5)
				cout << "in-only try theta/z, num stereo hits = " << locComboHits.size() << endl;
			if(!Find_ThetaZ(locCDCTrackCircle->fit, locComboHits, locTheta, locZ, locChiSqPerNDF))
				continue; //combo didn't work for some reason, try a different one
			if(!((locChiSqPerNDF < 1.0) || (locChiSqPerNDF > -1.0)))
				continue; // NaN
			locGoodStereoComboFoundFlag = true;
			if(DEBUG_LEVEL > 5)
				cout << "in-only good theta/z: theta, z, chisq, best-chisq = " << locTheta << ", " << locZ << ", " << locChiSqPerNDF << ", " << locBestChiSqPerNDF << endl;
			if(locChiSqPerNDF >= locBestChiSqPerNDF)
				continue;
			// This is the best combination of stereo seeds so far, save the results
			locBestSuperLayerSeeds_Inner.clear();
			for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i].size(); ++loc_k)
				locBestSuperLayerSeeds_Inner.push_back(locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i][loc_k]);
			locBestSuperLayerSeeds_Outer.clear();
			locBestTheta = locTheta;
			locBestZ = locZ;
			locBestChiSqPerNDF = locChiSqPerNDF;
		}
	}
	else if(!locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.empty())
	{
		// no hits in the inner super layers: e.g. decay product (e.g. p) of a long-lived decaying neutral particle (e.g. lambda)
		for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.size(); ++loc_j)
		{
			//get the hits from this inner combination
			deque<DCDCTrkHit*> locComboHits;
			Select_ThetaZStereoHits(locCDCTrackCircle, -1, loc_j, locFinalPassFlag, locComboHits);

			//Evaluate theta & z for this combination, returning the chisq/ndf from the fit
			double locTheta = 0.0, locZ = TARGET_Z, locChiSqPerNDF = 9.9E50;
			if(DEBUG_LEVEL > 5)
				cout << "out-only try theta/z, num stereo hits = " << locComboHits.size() << endl;
			if(!Find_ThetaZ(locCDCTrackCircle->fit, locComboHits, locTheta, locZ, locChiSqPerNDF))
				continue; //combo didn't work for some reason, try a different one
			if(!((locChiSqPerNDF < 1.0) || (locChiSqPerNDF > -1.0)))
				continue; // NaN
			locGoodStereoComboFoundFlag = true;
			if(DEBUG_LEVEL > 5)
					cout << "out-only good theta/z: theta, z, chisq, best-chisq = " << locTheta << ", " << locZ << ", " << locChiSqPerNDF << ", " << locBestChiSqPerNDF << endl;
			if(locChiSqPerNDF >= locBestChiSqPerNDF)
				continue;
			// This is the best combination of stereo seeds so far, save the results
			locBestSuperLayerSeeds_Inner.clear();
			locBestSuperLayerSeeds_Outer.clear();
			for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_j].size(); ++loc_k)
				locBestSuperLayerSeeds_Outer.push_back(locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_j][loc_k]);
			locBestTheta = locTheta;
			locBestZ = locZ;
			locBestChiSqPerNDF = locChiSqPerNDF;
		}
	}
	else //no stereo hits
		return (!locFinalPassFlag); //return true if don't wan't to filter, false if you do

	//save the results to the track circle
	locCDCTrackCircle->dTheta = locBestTheta;
	locCDCTrackCircle->dVertexZ = locBestZ;
	double locNumStereoSuperLayers = double(locBestSuperLayerSeeds_Inner.size() + locBestSuperLayerSeeds_Outer.size());
	locCDCTrackCircle->dWeightedChiSqPerDF_Stereo = locBestChiSqPerNDF/(locNumStereoSuperLayers*locNumStereoSuperLayers);

	set<DCDCSuperLayerSeed*> locAlreadyRecycledSuperLayerSeeds; //a seed can appear in more than one combo: don't recycle the same memory more than once!!

	// recycle the unused super layer seeds and store only the best ones: inner
	if(!locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.empty())
	{
		for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.size(); ++loc_i)
		{
			for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i].size(); ++loc_j)
			{
				DCDCSuperLayerSeed* locCDCSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i][loc_j];
				if(locAlreadyRecycledSuperLayerSeeds.find(locCDCSuperLayerSeed) != locAlreadyRecycledSuperLayerSeeds.end())
					continue; //already recycled!
				bool locKeepSuperLayerSeed = false;
				for(size_t loc_k = 0; loc_k < locBestSuperLayerSeeds_Inner.size(); ++loc_k)
				{
					if(locBestSuperLayerSeeds_Inner[loc_k] != locCDCSuperLayerSeed)
						continue;
					locKeepSuperLayerSeed = true; //one of the best, don't recycle!!
					break;
				}
				if(locKeepSuperLayerSeed)
					continue;
				Recycle_DCDCSuperLayerSeed(locCDCSuperLayerSeed); //no longer in use
				locAlreadyRecycledSuperLayerSeeds.insert(locCDCSuperLayerSeed);
			}
		}
		locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.clear();
		if(locGoodStereoComboFoundFlag)
			locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.push_back(locBestSuperLayerSeeds_Inner);
	}
	locAlreadyRecycledSuperLayerSeeds.clear();

	// recycle the unused super layer seeds and store only the best ones: outer
	if(!locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.empty())
	{
		for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.size(); ++loc_i)
		{
			for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_i].size(); ++loc_j)
			{
				DCDCSuperLayerSeed* locCDCSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[loc_i][loc_j];
				if(locAlreadyRecycledSuperLayerSeeds.find(locCDCSuperLayerSeed) != locAlreadyRecycledSuperLayerSeeds.end())
					continue; //already recycled!
				bool locKeepSuperLayerSeed = false;
				for(size_t loc_k = 0; loc_k < locBestSuperLayerSeeds_Outer.size(); ++loc_k)
				{
					if(locBestSuperLayerSeeds_Outer[loc_k] != locCDCSuperLayerSeed)
						continue;
					locKeepSuperLayerSeed = true; //one of the best, don't recycle!!
					break;
				}
				if(locKeepSuperLayerSeed)
					continue;
				Recycle_DCDCSuperLayerSeed(locCDCSuperLayerSeed); //no longer in use
				locAlreadyRecycledSuperLayerSeeds.insert(locCDCSuperLayerSeed);
			}
		}
		locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.clear();
		if(locGoodStereoComboFoundFlag)
			locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.push_back(locBestSuperLayerSeeds_Outer);
	}

	return locGoodStereoComboFoundFlag;
}

//------------------------
// Select_ThetaZStereoHits
//------------------------
void DTrackCandidate_factory_CDC::Select_ThetaZStereoHits(const DCDCTrackCircle* locCDCTrackCircle, int locInnerSeedSeriesIndex, int locOuterSeedSeriesIndex, bool locFinalPassFlag, deque<DCDCTrkHit*>& locComboHits)
{
	//tracks at the edges of the CDC's phase space OFTEN fail because the circle-fit is not perfect
		//this causes the projected-hit-position on the circle of the stereo hits to be incorrect
		//which then causes the calculated theta/z to be bad, which causes the momentum magnitude to be bad
		//this is especially a concern for tracks without hits in Super Layer (SL) 7: e.g. spiraling tracks in SL6 or tracks leaving the CDC
			//this is because the circle fit is extrapolated out to SL5 & SL6 and the errors are much larger
	//therefore, you must be VERY careful when selecting which stereo hits to use for the final calculation of theta & z
		//for the initial calculation: use all stereo hits: the initial calculation is intended to select which stereo super layer seeds are best, nothing more

	locComboHits.clear();

	// Get super layer seeds for this combination
	deque<DCDCSuperLayerSeed*> locSuperLayerSeeds;
	if(locInnerSeedSeriesIndex >= 0)
	{
		for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[locInnerSeedSeriesIndex].size(); ++loc_k)
			locSuperLayerSeeds.push_back(locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[locInnerSeedSeriesIndex][loc_k]);
	}
	if(locOuterSeedSeriesIndex >= 0)
	{
		for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[locOuterSeedSeriesIndex].size(); ++loc_k)
			locSuperLayerSeeds.push_back(locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[locOuterSeedSeriesIndex][loc_k]);
	}

	//select initial hits: prune hits that don't intersect the circle
	deque<DCDCTrkHit*> locHits;
	map<unsigned int, deque<DCDCTrkHit*> > locHitsBySuperLayer; //key is super layer
	unsigned int locTotalNumStereoHits = 0;
	for(size_t loc_k = 0; loc_k < locSuperLayerSeeds.size(); ++loc_k)
	{
		locSuperLayerSeeds[loc_k]->Get_Hits(locHits);
		for(deque<DCDCTrkHit*>::iterator locIterator = locHits.begin(); locIterator != locHits.end();)
		{
			if((*locIterator)->dValidStereoHitPosFlag)
				++locIterator;
			else
				locIterator = locHits.erase(locIterator);
		}
		locHitsBySuperLayer[locSuperLayerSeeds[loc_k]->dSuperLayer] = locHits;
		locTotalNumStereoHits += locHits.size();
	}

	// see if no more pruning is necessary
		// don't prune anymore if on initial pass, or if very few stereo hits
	if((!locFinalPassFlag) || (locTotalNumStereoHits <= MIN_PRUNED_STEREO_HITS))
	{
		// no more pruning, add hits to locComboHits and return
		map<unsigned int, deque<DCDCTrkHit*> >::iterator locMapIterator;
		for(locMapIterator = locHitsBySuperLayer.begin(); locMapIterator != locHitsBySuperLayer.end(); ++locMapIterator)
			locComboHits.insert(locComboHits.end(), locMapIterator->second.begin(), locMapIterator->second.end());
		if(DEBUG_LEVEL > 10)
			cout << "no more pruning, total num stereo hits = " << locTotalNumStereoHits << endl;
		return; //either on first pass (eval'ing which stereo seed is best), or not enough hits to prune: return
	}

	// Now, select the stereo hits whose projection onto the cirlce-fit are closest to the axial hits
		// the closer the stereo hits are to the axial hits, the better the estimation of theta/z will be

	//calc delta-phi for all remaining hits and sort them
		//delta-phi: between the intersection-of-the-stereo-hit-with-the-circle and the nearest axial hits
	deque<pair<DCDCTrkHit*, double> > locDeltaPhis;
	for(size_t loc_k = 0; loc_k < locSuperLayerSeeds.size(); ++loc_k)
	{
		unsigned int locSuperLayer = locSuperLayerSeeds[loc_k]->dSuperLayer;
		Calc_StereoHitDeltaPhis(locSuperLayer, locHitsBySuperLayer[locSuperLayer], locCDCTrackCircle, locDeltaPhis);
	}
	sort(locDeltaPhis.begin(), locDeltaPhis.end(), CDCSort_DeltaPhis);

	if(locDeltaPhis.size() <= MIN_PRUNED_STEREO_HITS)
	{
		for(size_t loc_k = 0; loc_k < locDeltaPhis.size(); ++loc_k)
			locComboHits.push_back(locDeltaPhis[loc_k].first);
		return;
	}

	// take at least the "MIN_PRUNED_STEREO_HITS" hits with the smallest delta_phi
	double locMaxHitDeltaPhi = locDeltaPhis[MIN_PRUNED_STEREO_HITS - 1].second;

	// now, potentially expand the max-delta-phi range in certain cases:
	double locDeltaPhiRangeExtension = 0.0;
	//if first statement below is true, want all the stereo hits:
		//track was matched to another track circle as a spiral turn, and it turns in its last axial layer (4 or 7)
		//because it is turning sharply, as much info/hits as possible is needed to give an accurate theta
		//however, if it was turning in a stereo layer, then the sharpest part of the turn was not included in the circle fit
			//in this case the hit-projections onto the circle are probably way off, so only expand to pi if turning on axial layer
	//if second statement below is true, then the track was very unlikely to spiral, because the radius of the circle indicates it (likely) passed through the bcal
		//in this case, all stereo hits within a reasonable distance (10 degrees) from the track circle are probably OK. 
		//note: expanding beyond 10-15 degrees kills candidates at theta >= 120:
			//these tracks leave the CDC at SL6 and sooner, and extrapolating the circle fit out to SL6 gives spurious results, so truncate the theta-search
	//if neither statement below is true, then the track likely spiraled in a stereo super layer:
		//don't trust the stereo information as much: only take stereo hits very close to the axial hits

	unsigned int locLastSuperLayer = locCDCTrackCircle->Get_LastSuperLayerSeed()->dSuperLayer;
	if(((locLastSuperLayer == 4) || (locLastSuperLayer == 7)) && (locCDCTrackCircle->dSpiralTurnRing != -1))
		locDeltaPhiRangeExtension = M_PI; //spiral turn on axial layer: all stereo hits good
	else if(locCDCTrackCircle->fit->r0 > 65.0/2.0) //BCAL is at r = ~65
		locDeltaPhiRangeExtension = 10.0; //unlikely to spiral, all stereo hits reasonably close are good
	else
		locDeltaPhiRangeExtension = 5.0; //likely spiraled, trust stereo information less
	locMaxHitDeltaPhi += locDeltaPhiRangeExtension;

	if(DEBUG_LEVEL > 10)
		cout << "prune with max delta-phi, fit circle r0 = " << locMaxHitDeltaPhi << ", " << locCDCTrackCircle->fit->r0 << endl;
	size_t locDeltaPhiIndex = 0;
	// add stereo hits to locComboHits whose projection onto the circle fit are within locMaxHitDeltaPhi of the nearest axial hit phi
	for(locDeltaPhiIndex = 0; locDeltaPhiIndex < locDeltaPhis.size(); ++locDeltaPhiIndex)
	{
		if(locDeltaPhis[locDeltaPhiIndex].second > locMaxHitDeltaPhi)
			break;
		locComboHits.push_back(locDeltaPhis[locDeltaPhiIndex].first);
	}

	//make sure to have at least one wire from each super layer
		// this is especially useful for spiral turns in SL6, where the majority of the theta-constraining information comes from SL5 and SL6
		// you don't want very many of these hits because they can throw you off, but having at least one really helps
		// this is experimental: it may not be necessary if the above locMaxHitDeltaPhi section is better fine-tuned

	// to keep track of which super layers need hits, delete keys from locHitsBySuperLayer map if hits from them are already selected
	map<unsigned int, deque<DCDCTrkHit*> >::iterator locMapIterator;
	for(size_t loc_i = 0; loc_i < locComboHits.size(); ++loc_i)
	{
		unsigned int locSuperLayer = (locComboHits[loc_i]->hit->wire->ring - 1)/4 + 1;
		locMapIterator = locHitsBySuperLayer.find(locSuperLayer);
		if(locMapIterator == locHitsBySuperLayer.end())
			continue; //already have a hit of that type
		locHitsBySuperLayer.erase(locMapIterator); //have a hit of this type
		if(locHitsBySuperLayer.empty())
			break;
	}

	//locHitsBySuperLayer now only contains the keys of super layers that don't have hits yet: loop over them
	for(locMapIterator = locHitsBySuperLayer.begin(); locMapIterator != locHitsBySuperLayer.end(); ++locMapIterator)
	{
		// search for the hit with the lowest delta-phi from this (locMapIterator->first) super layer
		for(size_t loc_i = locDeltaPhiIndex; loc_i < locDeltaPhis.size(); ++loc_i) //everything before locDeltaPhiIndex was already included
		{
			unsigned int locSuperLayer = (locDeltaPhis[loc_i].first->hit->wire->ring - 1)/4 + 1;
			if(locSuperLayer != locMapIterator->first)
				continue;
			locComboHits.push_back(locDeltaPhis[loc_i].first); //use best hit on this super layer
			if(DEBUG_LEVEL > 10)
				cout << "Add stereo hit on SL" << locSuperLayer << ": ring, straw = " << locDeltaPhis[loc_i].first->hit->wire->ring << ", " << locDeltaPhis[loc_i].first->hit->wire->straw << endl;
			break;
		}
	}

	if(DEBUG_LEVEL > 10)
		cout << "final #hits = " << locComboHits.size() << endl;
}

//------------------------
// Calc_StereoHitDeltaPhis
//------------------------
void DTrackCandidate_factory_CDC::Calc_StereoHitDeltaPhis(unsigned int locSuperLayer, deque<DCDCTrkHit*>& locHits, const DCDCTrackCircle* locCDCTrackCircle, deque<pair<DCDCTrkHit*, double> >& locDeltaPhis)
{
	// calc delta-phi for hits: distance from projected-stereo-hit-position to the nearest axial hits
		// because the circle fit is not perfect, the stereo hit position gets progessively worse the sharper the track is turning
		// this screws up the theta/z calculation, so just ignore the hits that are the farthest away
		// don't do this unless on the final pass though: even a moderately bad theta/z is (in theory) still good enough for vetoing totally wrong hit combinations

	if(DEBUG_LEVEL > 10)
		cout << "Calc_StereoHitDeltaPhis, SL = " << locSuperLayer << ", #hits = " << locHits.size() << endl;

	//get the axial super layer seeds nearest this seed:
	DCDCSuperLayerSeed* locPriorAxialSuperLayerSeed = NULL;
	DCDCSuperLayerSeed* locNextAxialSuperLayerSeed = NULL;
	for(size_t loc_k = 0; loc_k < locCDCTrackCircle->dSuperLayerSeeds_Axial.size(); ++loc_k)
	{
		if(locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_k]->dSuperLayer < locSuperLayer)
			locPriorAxialSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_k];
		else if((locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_k]->dSuperLayer > locSuperLayer) && (locNextAxialSuperLayerSeed == NULL))
			locNextAxialSuperLayerSeed = locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_k];
	}

	//get the rings in the axial super layer seeds nearest this stereo seed:
	DCDCRingSeed* locPriorAxialRingSeed = NULL;
	if(locPriorAxialSuperLayerSeed != NULL)
	{
		locPriorAxialRingSeed = &(locPriorAxialSuperLayerSeed->dCDCRingSeeds.back());
		if(DEBUG_LEVEL > 10)
			cout << "Prior axial ring = " << locPriorAxialRingSeed->ring << endl;
	}
	DCDCRingSeed* locNextAxialRingSeed = NULL;
	if(locNextAxialSuperLayerSeed != NULL)
	{
		locNextAxialRingSeed = &(locNextAxialSuperLayerSeed->dCDCRingSeeds.front());
		if(DEBUG_LEVEL > 10)
			cout << "Next axial ring = " << locNextAxialRingSeed->ring << endl;
	}
	if((locPriorAxialRingSeed == NULL) && (locNextAxialRingSeed == NULL))
		return; //no axial hits: shoudldn't be possible ...

	// calculate min-delta-phi for each hit (to the nearest axial ring seed)
	for(size_t loc_i = 0; loc_i < locHits.size(); ++loc_i)
	{
		deque<DCDCTrkHit*> locTempDeque(1, locHits[loc_i]);
		double locMinDeltaPhi = M_PI;
		//compare to prior axial ring
		if(locPriorAxialRingSeed != NULL)
		{
			double locDeltaPhi = MinDeltaPhi(locPriorAxialRingSeed->hits, locTempDeque);
			if(locDeltaPhi < locMinDeltaPhi)
				locMinDeltaPhi = locDeltaPhi;
		}
		//compare to next axial ring
		if(locNextAxialRingSeed != NULL)
		{
			double locDeltaPhi = MinDeltaPhi(locTempDeque, locNextAxialRingSeed->hits);
			if(locDeltaPhi < locMinDeltaPhi)
				locMinDeltaPhi = locDeltaPhi;
		}
		locMinDeltaPhi *= 180.0/M_PI;
		if(DEBUG_LEVEL > 10)
			cout << "Ring, Straw, min delta phi = " << locHits[loc_i]->hit->wire->ring << ", " << locHits[loc_i]->hit->wire->straw << ", " << locMinDeltaPhi << endl;
		//store the minimum
		locDeltaPhis.push_back(pair<DCDCTrkHit*, double>(locHits[loc_i], locMinDeltaPhi));
	}
}

//------------
// MinDeltaPhi
//------------
double DTrackCandidate_factory_CDC::MinDeltaPhi(const deque<DCDCTrkHit*>& locInnerSeedHits, const deque<DCDCTrkHit*>& locOuterSeedHits)
{
	/// Returns the minimum delta-phi between the two groups of hits. Assumes all of the hits in a given set are on the same ring. 
	/// First it checks if the two seeds overlap in phi: if so, then return 0
	/// Otherwise, only the first and last hits of the adjacent rings between each seed's hit list are used. 
	/// to calculate a maximum of 4 delta-phis (minimum of 1) of which, the smallest is returned.
	if(locInnerSeedHits.empty() || locOuterSeedHits.empty())
	{
		cout << "Number of seed hits 0! (Ninner = " << locInnerSeedHits.size() << " ,Nouter = " << locOuterSeedHits.size() << ")" << endl;
		return M_PI;
	}

	DCDCTrkHit* locInnermostRingFirstStrawHit = locInnerSeedHits.front();
	DCDCTrkHit* locInnermostRingLastStrawHit = locInnerSeedHits.back();
	DCDCTrkHit* locOutermostRingFirstStrawHit = locOuterSeedHits.front();
	DCDCTrkHit* locOutermostRingLastStrawHit = locOuterSeedHits.back();

	//see if seeds overlap in phi
	float locInnermostRingFirstStrawPhi = locInnermostRingFirstStrawHit->hit->wire->phi;
	float locInnermostRingLastStrawPhi = locInnermostRingLastStrawHit->hit->wire->phi;
	float locOutermostRingFirstStrawPhi = locOutermostRingFirstStrawHit->hit->wire->phi;
	float locOutermostRingLastStrawPhi = locOutermostRingLastStrawHit->hit->wire->phi;
	if(DEBUG_LEVEL > 100)
		cout << "inner ring: ring, first/last straws & phis = " << locInnermostRingFirstStrawHit->hit->wire->ring << ", " << locInnermostRingFirstStrawHit->hit->wire->straw << ", " << locInnermostRingLastStrawHit->hit->wire->straw << ", " << locInnermostRingFirstStrawPhi << ", " << locInnermostRingLastStrawPhi << endl;
	if(DEBUG_LEVEL > 100)
		cout << "outer ring: ring, first/last straws & phis = " << locOutermostRingFirstStrawHit->hit->wire->ring << ", " << locOutermostRingFirstStrawHit->hit->wire->straw << ", " << locOutermostRingLastStrawHit->hit->wire->straw << ", " << locOutermostRingFirstStrawPhi << ", " << locOutermostRingLastStrawPhi << endl;

	//account for phi = 0/2pi boundary
	bool locInnerRingCrossesBoundaryFlag = (locInnermostRingLastStrawPhi < locInnermostRingFirstStrawPhi);
	bool locOuterRingCrossesBoundaryFlag = (locOutermostRingLastStrawPhi < locOutermostRingFirstStrawPhi);
	if(DEBUG_LEVEL > 100)
		cout << "in/out boundary flags = " << locInnerRingCrossesBoundaryFlag << ", " << locOuterRingCrossesBoundaryFlag << endl;
	if(locOuterRingCrossesBoundaryFlag)
		locOutermostRingLastStrawPhi += M_TWO_PI;
	if(locInnerRingCrossesBoundaryFlag)
		locInnermostRingLastStrawPhi += M_TWO_PI;
	if(locOuterRingCrossesBoundaryFlag & (!locInnerRingCrossesBoundaryFlag) && ((locOutermostRingLastStrawPhi - locInnermostRingLastStrawPhi) > M_PI))
	{
		locInnermostRingFirstStrawPhi += M_TWO_PI;
		locInnermostRingLastStrawPhi += M_TWO_PI;
	}
	if(locInnerRingCrossesBoundaryFlag & (!locOuterRingCrossesBoundaryFlag) && ((locInnermostRingLastStrawPhi - locOutermostRingLastStrawPhi) > M_PI))
	{
		locOutermostRingFirstStrawPhi += M_TWO_PI;
		locOutermostRingLastStrawPhi += M_TWO_PI;
	}

	if(DEBUG_LEVEL > 100)
		cout << "final inner ring: ring, first/last straws & phis = " << locInnermostRingFirstStrawHit->hit->wire->ring << ", " << locInnermostRingFirstStrawHit->hit->wire->straw << ", " << locInnermostRingLastStrawHit->hit->wire->straw << ", " << locInnermostRingFirstStrawPhi << ", " << locInnermostRingLastStrawPhi << endl;
	if(DEBUG_LEVEL > 100)
		cout << "final outer ring: ring, first/last straws & phis = " << locOutermostRingFirstStrawHit->hit->wire->ring << ", " << locOutermostRingFirstStrawHit->hit->wire->straw << ", " << locOutermostRingLastStrawHit->hit->wire->straw << ", " << locOutermostRingFirstStrawPhi << ", " << locOutermostRingLastStrawPhi << endl;

	if((locOutermostRingFirstStrawPhi >= locInnermostRingFirstStrawPhi) && (locOutermostRingFirstStrawPhi <= locInnermostRingLastStrawPhi))
		return 0.0;
	if((locOutermostRingLastStrawPhi >= locInnermostRingFirstStrawPhi) && (locOutermostRingLastStrawPhi <= locInnermostRingLastStrawPhi))
		return 0.0;
	if((locInnermostRingFirstStrawPhi >= locOutermostRingFirstStrawPhi) && (locInnermostRingFirstStrawPhi <= locOutermostRingLastStrawPhi))
		return 0.0; //4th case not needed.  this case only needed if innermost ring is one wire across

	//make all 4 comparisons between hits
	double locDeltaPhi, locMinDeltaPhi;
	locMinDeltaPhi = fabs(locInnermostRingFirstStrawPhi - locOutermostRingFirstStrawPhi);
	if(locMinDeltaPhi > M_PI)
		locMinDeltaPhi = fabs(locMinDeltaPhi - M_TWO_PI);
	if(locOutermostRingFirstStrawHit != locOutermostRingLastStrawHit)
	{
		locDeltaPhi = fabs(locInnermostRingFirstStrawPhi - locOutermostRingLastStrawPhi);
		if(locDeltaPhi > M_PI)
			locDeltaPhi = fabs(locDeltaPhi - M_TWO_PI);
		if(locDeltaPhi < locMinDeltaPhi)
			locMinDeltaPhi = locDeltaPhi;
	}
	if(locInnermostRingFirstStrawHit == locInnermostRingLastStrawHit)
		return locMinDeltaPhi;

	locDeltaPhi = fabs(locInnermostRingLastStrawPhi - locOutermostRingFirstStrawPhi);
	if(locDeltaPhi > M_PI)
		locDeltaPhi = fabs(locDeltaPhi - M_TWO_PI);
	if(locDeltaPhi < locMinDeltaPhi)
		locMinDeltaPhi = locDeltaPhi;
	if(locOutermostRingFirstStrawHit != locOutermostRingLastStrawHit)
	{
		locDeltaPhi = fabs(locInnermostRingLastStrawPhi - locOutermostRingLastStrawPhi);
		if(locDeltaPhi > M_PI)
			locDeltaPhi = fabs(locDeltaPhi - M_TWO_PI);
		if(locDeltaPhi < locMinDeltaPhi)
			locMinDeltaPhi = locDeltaPhi;
	}

	return locMinDeltaPhi;
}

//------------
// Find_ThetaZ
//------------
bool DTrackCandidate_factory_CDC::Find_ThetaZ(const DHelicalFit* locFit, const deque<DCDCTrkHit*>& locStereoHits, double& locTheta, double& locZ, double& locChiSqPerNDF)
{
	// Calculate theta/z for the input stereo hits
	if(locStereoHits.empty())
		return false;

	if(Find_ThetaZ_Regression(locFit, locStereoHits, locTheta, locZ, locChiSqPerNDF))
		return true;
	double locThetaMin, locThetaMax;

	// Regression fit failed, try using histogram methods
	bool locThetaOKFlag = Find_Theta(locFit, locStereoHits, locTheta, locThetaMin, locThetaMax, locChiSqPerNDF);
	if(locThetaOKFlag)
	{
		if(Find_Z(locFit, locStereoHits, locThetaMin, locThetaMax, locZ))
			return true;
	}

	// Histogram methods failed. 

	// Assume that the track came from the center of the target
	locChiSqPerNDF = 9.9E8;
	locZ = TARGET_Z;
	if(locThetaOKFlag)
		return true;

	// Use a point in one of the stereo layers to estimate tanl
	double x = locStereoHits[0]->dStereoHitPos.X();
	double y = locStereoHits[0]->dStereoHitPos.Y();
	double tworc = 2.0*locFit->r0;
	double ratio = sqrt(x*x + y*y)/tworc;
	if(ratio >= 1.0)
		return false;

	double tanl = (locStereoHits[0]->dStereoHitPos.Z() - locZ)/(tworc*asin(ratio));
	locTheta = M_PI_2 - atan(tanl);
	return true;
}

//-----------------------
// Find_ThetaZ_Regression
//-----------------------
// Linear regression to find tan(lambda) and z_vertex.
// This method assumes that there are errors in both the z positions and 
// the arc lengths.
// Algorithm from Numerical Recipes in C (2nd. ed.), pp. 668-669.
bool DTrackCandidate_factory_CDC::Find_ThetaZ_Regression(const DHelicalFit* locFit, const deque<DCDCTrkHit*>& locStereoHits, double& locTheta, double& locZ, double& locChiSqPerNDF)
{
	if(DEBUG_LEVEL > 3)
		cout<<"Finding theta and z via linear regression method."<<endl;

	if(locStereoHits.empty() || (!(locFit->normal.Mag() > 0.0)))
		return false;

	// Vector of intersections between the circle and the stereo wires
	vector<intersection_t> intersections;
	for(size_t m = 0; m < locStereoHits.size(); ++m)
	{
		DCDCTrkHit* trkhit = locStereoHits[m];

		//DVector3_with_perp intersection;
		intersection_t intersection;
		intersection.x = trkhit->dStereoHitPos.X();
		intersection.y = trkhit->dStereoHitPos.Y();
		intersection.perp2 = intersection.x*intersection.x + intersection.y*intersection.y;
		intersection.z = trkhit->dStereoHitPos.Z();
		intersection.var_z = trkhit->var_z;
		intersections.push_back(intersection);
	}

	// Now, sort the entries
	sort(intersections.begin(), intersections.end(), CDCSort_Intersections);
 
	// Compute the arc lengths between the origin in x and y and (xi,yi)
	vector<double> arclengths(intersections.size()); 
	vector<double> ratios(intersections.size());
	double xc = locFit->x0;
	double yc = locFit->y0;
	double rc = locFit->r0;
	double two_rc = 2.*rc;

	// Find POCA to beam line
	double myphi = atan2(yc, xc);
	double y0 = yc - rc*sin(myphi);
	double x0 = xc - rc*cos(myphi);

	// Arc length to first measurement
	double diffx = intersections[0].x - x0;
	double diffy = intersections[0].y - y0;
	double chord = sqrt(diffx*diffx + diffy*diffy);
	double ratio = chord/two_rc;
	double s = (ratio < 1.) ? two_rc*asin(ratio) : M_PI_2*two_rc;
	arclengths[0] = s;
	ratios[0] = ratio;

	// Find arc lengths for the rest of the stereo hits
	for(size_t m = 1; m < arclengths.size(); ++m)
	{
		diffx = intersections[m].x - intersections[m - 1].x;
		diffy = intersections[m].y - intersections[m - 1].y;
		chord = sqrt(diffx*diffx + diffy*diffy);
		ratio = chord/two_rc;
		if(ratio > 0.999)
			return false;
		double ds = two_rc*asin(ratio);
		s += ds;
		arclengths[m] = s;
		ratios[m] = ratio;
	}

	//Linear regression to find z0, tanl
	double tanl = 0.,z0 = 0.;
	if(arclengths.size() > 1) // Do fit only if have more than one measurement
	{
		DCDCLineFit fit;
		size_t n = fit.n = intersections.size();
		fit.s.resize(n);
		fit.var_s.resize(n);
		fit.z.resize(n);
		fit.var_z.resize(n);
		fit.w.resize(n);

		// Find average variances for z and s
		double avg_var_s = 0., avg_var_z = 0.;
		double var_r = 1.6*1.6/12.; // assume cell size
		for (size_t m = 0; m < n; ++m)
		{
			fit.s[m] = arclengths[m];
			fit.var_s[m] = var_r/(1. - ratios[m]*ratios[m]);

			avg_var_s += fit.var_s[m];
			avg_var_z += intersections[m].var_z;

			if(DEBUG_LEVEL>5)
				cout<<"Using CDC hit "<<m<<" z="<<intersections[m].z << " s=" << arclengths[m] <<endl;
		}

		// Scale z errors according to the ratio of the average variances
		double scale2 = avg_var_s/avg_var_z;
		double scale = sqrt(scale2);
		vector<double> weight(n);
		for (size_t m = 0; m < n; ++m)
		{
			fit.z[m] = scale*intersections[m].z;
			fit.var_z[m] = scale2*intersections[m].var_z;
			weight[m] = fit.var_s[m] + fit.var_z[m];
		}

		// Perform preliminary fit to find the (scaled) slope tanl
		double sumv=0., sumx=0.;
		double sumy=0., sumxx=0., sumxy=0.;
		for(size_t m = 0; m < n; ++m)
		{
			//double temp = 1./var_z[m];
			double temp = 1./weight[m];
			sumv += temp;
			sumx += arclengths[m]*temp;
			sumy += fit.z[m]*temp;
			sumxx += arclengths[m]*arclengths[m]*temp;
			sumxy += arclengths[m]*fit.z[m]*temp;
		}
		double Delta = sumv*sumxx - sumx*sumx;
		if(!(fabs(Delta) > 0.0))
			return false;

		tanl = (sumv*sumxy - sumx*sumy)/Delta;
		fit.z0 = (sumxx*sumy - sumx*sumxy)/Delta;

		// Convert tanl to an angle and create two other reference angles
		double angle[3];
		angle[0] = 0.;
		angle[1] = atan(tanl);
		angle[2] = 1.571;
		// Compute chi^2 values for line fits with these three angles
		double ch[3];
		for (unsigned int m = 0; m < 3; ++m)
			ch[m] = fit.ChiXY(angle[m]);

		// Bracket the minimum chi^2
		fit.BracketMinimumChisq(angle[0], angle[1], angle[2], ch[0], ch[1], ch[2]);
		// Find the minimum chi^2 using Brent's method and compute the best value for lambda
		double lambda = 0.;
		locChiSqPerNDF = fit.FindMinimumChisq(angle[0], angle[1], angle[2], lambda)/2.0; //2 degrees of freedom
		// Undo the scaling
		z0 = fit.z0/scale;
		tanl = tan(lambda)/scale;
	}
	else
	{
		z0 = TARGET_Z;
		tanl = (intersections[0].z - z0)/arclengths[0];
		locChiSqPerNDF = 9.9E9; //only two hits: technically is zero, but if possible want to pick a group of stereo hits with more hits
	}

	locTheta = M_PI_2 - atan(tanl);
	locZ = z0;

	return true;
}

//-----------
// Find_Theta
//-----------
bool DTrackCandidate_factory_CDC::Find_Theta(const DHelicalFit* locFit, const deque<DCDCTrkHit*>& locStereoHits, double& locTheta, double& locThetaMin, double& locThetaMax, double& locChiSqPerNDF)
{
	if(locStereoHits.empty())
		return false;
	/// Find the theta value using the input stereo hits.
	/// The values for dPhiStereo and dStereoHitPos.Z() are assumed to be valid. 
	/// The value of locFit.r0 is also used to calculate theta.
	///
	/// This uses a histogramming technique that looks at the overlaps of the
	/// angle ranges subtended by each hit between the given target limits.
	/// The overlaps usually lead to a range of values for theta. The limits
	/// of these are stored in locThetaMin and locThetaMax. 
	/// The centroid of the range is stored in the theta field.
		
	// We use a simple array to store our histogram here. We don't want to use
	// ROOT histograms because they are not thread safe.
	unsigned int Nbins = 1000;
	unsigned int hist[Nbins];
	for(unsigned int i=0; i<Nbins; ++i)
		hist[i] = 0; // clear histogram
	double bin_width = M_TWO_PI/(double)Nbins;
	double hist_low_limit = -M_PI; // lower edge of histogram limits
	
	// Loop over CDC hits, filling the histogram
	double r0 = locFit->r0;
	for(unsigned int i=0; i < locStereoHits.size(); ++i)
	{
		DCDCTrkHit *trkhit = locStereoHits[i];

		// Calculate upper and lower limits in theta
		double alpha = r0*trkhit->dPhiStereo;
		if(locFit->h < 0.0)
			alpha = -alpha;
		double tmin = atan2(alpha, trkhit->dStereoHitPos.Z() - VERTEX_Z_MIN);
		double tmax = atan2(alpha, trkhit->dStereoHitPos.Z() - VERTEX_Z_MAX);
		if(tmin>tmax)
		{
			double tmp = tmin;
			tmin=tmax;
			tmax=tmp;
		}
		if(DEBUG_LEVEL>3)
			cout<<" -- phi_stereo="<<trkhit->dPhiStereo<<" z_stereo="<<trkhit->dStereoHitPos.Z()<<"  alpha="<<alpha<<endl;
		if(DEBUG_LEVEL>3)
			cout<<" -- tmin="<<tmin<<"  tmax="<<tmax<<endl;
		
		// Find index of bins corresponding to tmin and tmax
		unsigned int imin = (unsigned int)floor((tmin-hist_low_limit)/bin_width);
		unsigned int imax = (unsigned int)floor((tmax-hist_low_limit)/bin_width);
		
		// If entire range of this hit is outside of the histogram limit
		// then discard this hit.
		if(imin>=Nbins)
			continue;
		
		// Clip limits of bin range to our histogram limits
		if(imin>=Nbins)
			imin=Nbins-1;
		if(imax>=Nbins)
			imax=Nbins-1;
		
		// Increment all bins between imin and imax
		for(unsigned int j=imin; j<=imax; ++j)
			++hist[j];
	}
		
	// Look for the indexes of the plateau
	unsigned int istart=0;
	unsigned int iend=0;
	for(unsigned int i=1; i<Nbins; ++i)
	{
		if(hist[i]>hist[istart])
		{
			istart = i;
			if(DEBUG_LEVEL>3)
				cout<<" -- istart="<<istart<<" (theta="<<hist_low_limit + bin_width*(0.5+(double)istart)<<" , N="<<hist[i]<<")"<<endl;
		}
		if(hist[i] == hist[istart])
			iend = i;
	}
	
	// If there are no entries in the histogram, then return false
	if(hist[istart] == 0)
		return false;
	
	// Calculate theta limits
	locThetaMin = hist_low_limit + bin_width*(0.5+(double)istart);
	locThetaMax = hist_low_limit + bin_width*(0.5+(double)iend);
	locTheta = (locThetaMax + locThetaMin)/2.0;
	if(DEBUG_LEVEL>3)
		cout<<"istart="<<istart<<" iend="<<iend<<" theta_min="<<locThetaMin<<" theta_max="<<locThetaMax<<endl;
	locChiSqPerNDF = 9.9E5; //NEED TO CALCULATE THIS!!!

	return true;
}

//-------
// Find_Z
//-------
bool DTrackCandidate_factory_CDC::Find_Z(const DHelicalFit* locFit, const deque<DCDCTrkHit*>& locStereoHits, double locThetaMin, double locThetaMax, double& locZ)
{
	if(locStereoHits.empty())
		return false;

	/// Find the z value of the vertex using the stereo hits. 
	/// The values for dPhiStereo and dStereoHitPos.Z() are assumed to be valid. 
	///
	/// This uses a histogramming technique that looks at the overlaps of the
	/// z ranges subtended by each hit between the given theta limits.
	/// The overlaps usually lead to a range of values for z_vertex. 
	/// The centroid of the range is returned as locZ.
	
	// We use a simple array to store our histogram here. We don't want to use
	// ROOT histograms because they are not thread safe.
	unsigned int Nbins = 300;
	unsigned int hist[Nbins];
	for(unsigned int i=0; i<Nbins; ++i)
		hist[i] = 0; // clear histogram
	double bin_width = 0.5; // bins are 5mm
	double hist_low_limit = 0.0; // lower edge of histogram limits
	
	// Loop over CDC hits, filling the histogram
	double r0 = locFit->r0;
	double tan_alpha_min = tan(locThetaMin)/r0;
	double tan_alpha_max = tan(locThetaMax)/r0;
	for(unsigned int i=0; i< locStereoHits.size(); ++i)
	{
		DCDCTrkHit* trkhit = locStereoHits[i];
		
		// Calculate upper and lower limits in z
		double q_sign = locFit->h > 0.0 ? +1.0:-1.0;
		double zmin = trkhit->dStereoHitPos.Z() - q_sign*trkhit->dPhiStereo/tan_alpha_min;
		double zmax = trkhit->dStereoHitPos.Z() - q_sign*trkhit->dPhiStereo/tan_alpha_max;
		if(zmin>zmax)
		{
			double tmp = zmin;
			zmin=zmax;
			zmax=tmp;
		}
		if(DEBUG_LEVEL>3)
			cout<<" -- phi_stereo="<<trkhit->dPhiStereo<<" z_stereo="<<trkhit->dStereoHitPos.Z()<<endl;
		if(DEBUG_LEVEL>3)
			cout<<" -- zmin="<<zmin<<"  zmax="<<zmax<<endl;
		
		// Find index of bins corresponding to tmin and tmax
		unsigned int imin = (unsigned int)floor((zmin-hist_low_limit)/bin_width);
		unsigned int imax = (unsigned int)floor((zmax-hist_low_limit)/bin_width);
		
		// If entire range of this hit is outside of the histogram limit
		// then discard this hit.
		if(imax<=0 || imin>=Nbins)
			continue;
		
		// Clip limits of bin range to our histogram limits
		if(imin>=Nbins)
			imin=Nbins-1;
		if(imax>=Nbins)
			imax=Nbins-1;
		
		// Increment all bins between imin and imax
		for(unsigned int j=imin; j<=imax; ++j)
			++hist[j];
	}
	
	// Look for the indexes of the plateau
	unsigned int istart=0;
	unsigned int iend=0;
	for(unsigned int i=1; i<Nbins; ++i)
	{
		if(hist[i]>hist[istart])
		{
			istart = i;
			if(DEBUG_LEVEL>3)
				cout<<" -- istart="<<istart<<" (z="<<hist_low_limit + bin_width*(0.5+(double)istart)<<" , N="<<hist[i]<<")"<<endl;
		}
		if(hist[i] == hist[istart])
			iend = i;
	}

	// If there are no entries in the histogram, then return false
	if(hist[istart] == 0)
		return false;
	
	// Calculate z limits
	double locZMin = hist_low_limit + bin_width*(0.5+(double)istart);
	double locZMax = hist_low_limit + bin_width*(0.5+(double)iend);
	locZ = (locZMax + locZMin)/2.0;
	if(DEBUG_LEVEL>3)
		cout<<"istart="<<istart<<" iend="<<iend<<" z_min="<<locZMin<<" z_max="<<locZMax<<" hits[istart]="<<hist[istart]<<endl;
	return true;
}

//---------------------------
// Recycle_DCDCSuperLayerSeed
//---------------------------
void DTrackCandidate_factory_CDC::Recycle_DCDCSuperLayerSeed(DCDCSuperLayerSeed* locCDCSuperLayerSeed)
{
	// this function should ONLY be called for stereo super layers AFTER the hits have been projected onto the circle (new hit objects were made)
		// this is useful for recycling the memory used by the projected hits and seeds
	// first loop over the hits: see if can recycle those as well (each hit can be used in multiple seeds)

	deque<DCDCTrkHit*> locHits;
	locCDCSuperLayerSeed->Get_Hits(locHits);
	for(size_t loc_i = 0; loc_i < locHits.size(); ++loc_i)
	{
		DCDCTrkHit* locCDCTrkHit = locHits[loc_i];
		if(dStereoHitNumUsedMap[locCDCTrkHit] == 1) //this is the last remaining DCDCSuperLayerSeed this hit is used in: recycle it
			dCDCTrkHitPool_Available.push_back(locCDCTrkHit);
		else
			--dStereoHitNumUsedMap[locCDCTrkHit];
	}
	dCDCSuperLayerSeedPool_Available.push_back(locCDCSuperLayerSeed); //recycle
}

//---------------------------
// Recycle_DCDCTrackCircle
//---------------------------
void DTrackCandidate_factory_CDC::Recycle_DCDCTrackCircle(DCDCTrackCircle* locCDCTrackCircle)
{
	if(locCDCTrackCircle->fit != NULL)
	{
		dHelicalFitPool_Available.push_back(locCDCTrackCircle->fit);
		locCDCTrackCircle->fit = NULL;
	}
	locCDCTrackCircle->Reset();
	dCDCTrackCirclePool_Available.push_back(locCDCTrackCircle); //recycle
}

//----------------------
// Set_HitBitPattern_All
//----------------------
void DTrackCandidate_factory_CDC::Set_HitBitPattern_All(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	unsigned int locNumBits = 8*sizeof(unsigned int);
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		DCDCTrackCircle* locCDCTrackCircle = locCDCTrackCircles[loc_i];
		locCDCTrackCircle->HitBitPattern.clear();
		locCDCTrackCircle->HitBitPattern.resize(dNumCDCHits/(8*sizeof(unsigned int)) + 1);
		deque<DCDCTrkHit*> locHits;
		for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_Axial.size(); ++loc_j)
		{
			locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_j]->Get_Hits(locHits);
			for(size_t loc_k = 0; loc_k < locHits.size(); ++loc_k)
				locCDCTrackCircle->HitBitPattern[locHits[loc_k]->index/locNumBits] |= 1 << locHits[loc_k]->index % locNumBits;
		}

		if(!locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.empty())
		{
			for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[0].size(); ++loc_j)
			{
				locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[0][loc_j]->Get_Hits(locHits);
				for(size_t loc_k = 0; loc_k < locHits.size(); ++loc_k)
					locCDCTrackCircle->HitBitPattern[locHits[loc_k]->index/locNumBits] |= 1 << locHits[loc_k]->index % locNumBits;
			}
		}
		if(!locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.empty())
		{
			for(size_t loc_j = 0; loc_j < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[0].size(); ++loc_j)
			{
				locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[0][loc_j]->Get_Hits(locHits);
				for(size_t loc_k = 0; loc_k < locHits.size(); ++loc_k)
					locCDCTrackCircle->HitBitPattern[locHits[loc_k]->index/locNumBits] |= 1 << locHits[loc_k]->index % locNumBits;
			}
		}
	}
}

//---------------------------
// Filter_TrackCircles_Stereo
//---------------------------
void DTrackCandidate_factory_CDC::Filter_TrackCircles_Stereo(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	if(locCDCTrackCircles.empty())
		return;

	//sort track circles so that the ones with the best stereo chisq/ndf are first
	sort(locCDCTrackCircles.begin(), locCDCTrackCircles.end(), CDCSortByStereoChiSqPerNDFIncreasing);

	if(DEBUG_LEVEL > 5)
	{
		cout << "filter stereo, circles sorted by stereo chisq/ndf" << endl;
		Print_TrackCircles(locCDCTrackCircles);
	}

	//FIRST: Delete super layer seeds if they are shared between circles: delete it from the one with the worst dWeightedChiSqPerDF_Stereo
		// This is typically a problem when a track exits the CDC before reaching SL5 or SL6, and picks up the SL5 & SL6 hits from a nearby track instead
		// However, DO NOT delete the super layer seed if it is the last one on the track
	for(size_t loc_i = 0; loc_i < (locCDCTrackCircles.size() - 1); ++loc_i)
	{
		//assume this track circle has no wrong seeds (since it has the best chisq/ndf of all circles containing any of these seeds)
		DCDCTrackCircle* locTrackCircle_ToCompareTo = locCDCTrackCircles[loc_i];
		deque<DCDCSuperLayerSeed*> locStereoSuperLayerSeeds;
		locTrackCircle_ToCompareTo->Get_AllStereoSuperLayerSeeds(locStereoSuperLayerSeeds);
		//if any of the following track circles has any of the stereo super layer seeds in this track circle, then delete those seeds from those circles
			//should recycle...
		bool locSeedsStrippedFlag_AnyCircle = false;
		for(size_t loc_j = loc_i + 1; loc_j < locCDCTrackCircles.size(); ++loc_j)
		{
			DCDCTrackCircle* locTrackCircle_Validating = locCDCTrackCircles[loc_j];
			if(locTrackCircle_Validating->Get_NumStereoSuperLayerSeeds() <= 1)
				continue; //don't strip the last stereo hits!! (stripping isn't perfrect, maybe they are correct, or maybe hit selector will pick the correct ones)
			bool locSeedsStrippedFlag = false;
			for(size_t loc_k = 0; loc_k < locStereoSuperLayerSeeds.size(); ++loc_k)
			{
				if(locTrackCircle_Validating->Get_NumStereoSuperLayerSeeds() <= 1)
					break;
				//remember, can't compare pointers directly (made new objects for projection to track circle), must compare dSuperLayer and dSeedIndex
				unsigned int locSuperLayer = locStereoSuperLayerSeeds[loc_k]->dSuperLayer;
				unsigned int locSeedIndex = locStereoSuperLayerSeeds[loc_k]->dSeedIndex;
				DCDCSuperLayerSeed* locSuperLayerSeed = locTrackCircle_Validating->Get_SuperLayerSeed(locSuperLayer);
				if(locSuperLayerSeed == NULL)
					continue;
				if((locSuperLayerSeed->dSuperLayer != locSuperLayer) || (locSuperLayerSeed->dSeedIndex != locSeedIndex))
					continue;
				//stereo super layer is shared by both track circles, reject it from the one with the worse stereo chisq/ndf (weighted)
				if(DEBUG_LEVEL > 10)
					cout << "strip SL" << locSuperLayer << " from track circle " << loc_j << endl;
				locTrackCircle_Validating->Strip_StereoSuperLayerSeed(locSuperLayer);
				locSeedsStrippedFlag = true; //may want to recalc theta/z
				locSeedsStrippedFlag_AnyCircle = true;
			}
			if(locSeedsStrippedFlag) //recalc theta/z to get new chisq/ndf
				Select_CDCSuperLayerSeeds(locTrackCircle_Validating, false);
		}
		
		if(locSeedsStrippedFlag_AnyCircle) //re-sort if fits have been re-performed
			sort(locCDCTrackCircles.begin() + loc_i + 1, locCDCTrackCircles.end(), CDCSortByStereoChiSqPerNDFIncreasing);
	}

	//restore sort by circle-fit chisq/ndf (weighted)
	sort(locCDCTrackCircles.begin(), locCDCTrackCircles.end(), CDCSortByChiSqPerNDFDecreasing);

	if(DEBUG_LEVEL > 5)
	{
		cout << "filter stereo, circles sorted by axial chisq/ndf" << endl;
		Print_TrackCircles(locCDCTrackCircles);
	}
}

//---------------
// Add_UnusedHits
//---------------
void DTrackCandidate_factory_CDC::Add_UnusedHits(deque<DCDCTrackCircle*>& locCDCTrackCircles)
{
	if(DEBUG_LEVEL > 5)
		cout << "Add unused hits" << endl;

	// If the last super layer of a track is not 7, search for lone, unused hits on the next super layer and add them to the track
		// Only add them if they can link to the previous super layer, and wouldn't skip too many rings
		// Add at most one hit: prefer the hit on the innermost ring
			//If more than one hit on the innermost axial ring: choose the one closest to the circle fit
			//If more than one hit on the innermost stereo ring: ambiguous, don't add any hits
	for(size_t loc_i = 0; loc_i < locCDCTrackCircles.size(); ++loc_i)
	{
		DCDCSuperLayerSeed* locLastSuperLayerSeed = locCDCTrackCircles[loc_i]->Get_LastSuperLayerSeed();
		unsigned int locLastSuperLayer = locLastSuperLayerSeed->dSuperLayer;
		if(DEBUG_LEVEL > 5)
			cout << "i, last super layer = " << loc_i << ", " << locLastSuperLayer << endl;
		if(locLastSuperLayer == 7)
			continue;
		unsigned int locSearchSuperLayer = locLastSuperLayer + 1;
		const DHelicalFit* locFit = locCDCTrackCircles[loc_i]->fit;

		//make sure the next super layer doesn't correspond to a region where all of the seeds were deleted earlier (density too high)
		double locSeedFirstPhi, locSeedLastPhi;
		Calc_SuperLayerPhiRange(locLastSuperLayerSeed, locSeedFirstPhi, locSeedLastPhi);
		bool locHitDensityTooHighFlag = false;
		for(size_t loc_k = 0; loc_k < dRejectedPhiRegions[locSearchSuperLayer].size(); ++loc_k)
		{
			if(!Check_IfPhiRangesOverlap(locSeedFirstPhi, locSeedLastPhi, dRejectedPhiRegions[locSearchSuperLayer][loc_k].first, dRejectedPhiRegions[locSearchSuperLayer][loc_k].second))
				continue;
			locHitDensityTooHighFlag = true;
			break;
		}
		if(locHitDensityTooHighFlag)
		{
			if(DEBUG_LEVEL > 5)
				cout << "hit density too high, don't add unused hits" << endl;
			continue; //hit density in this region is too high, ignore all hits here
		}

		//make sure we don't skip too many rings
		int locLastHitRing = locLastSuperLayerSeed->dCDCRingSeeds.back().ring;
		if((4*locLastSuperLayer - locLastHitRing) > MAX_NUM_RINGSEED_RINGS_SKIPABLE)
		{
			if(DEBUG_LEVEL > 5)
				cout << "too many rings missing at the end of the last super layer: actual, max = " << 4*locLastSuperLayer - locLastHitRing << ", " << MAX_NUM_RINGSEED_RINGS_SKIPABLE << endl;
			continue; //too many rings missing at the end of the last super layer
		}

		wire_direction_t locWireDirection;
		if((locSearchSuperLayer == 1) || (locSearchSuperLayer == 4) || (locSearchSuperLayer == 7))
			locWireDirection = WIRE_DIRECTION_AXIAL;
		else if((locSearchSuperLayer == 2) || (locSearchSuperLayer == 6))
			locWireDirection = WIRE_DIRECTION_STEREOLEFT;
		else
			locWireDirection = WIRE_DIRECTION_STEREORIGHT;

		//ok, look for unused hits within a certain angle/distance from the circle fit
		DCDCTrkHit* locBestTrkHit = NULL;
		int locAmbiguousHitRing = -1; //if != -1, ignore all hits in or above this ring //set if > 1 stereo hit that matches on this ring
		double locBestDeltaPhi = 9.9E9; //for axial if > 1 hit on innermost ring
		for(size_t loc_j = 0; loc_j < cdchits_by_superlayer[locSearchSuperLayer - 1].size(); ++loc_j)
		{
			DCDCTrkHit* locTrkHit = cdchits_by_superlayer[locSearchSuperLayer - 1][loc_j];
			if(locTrkHit->flags & USED)
				continue; //not a lone hit: used in a super layer seed
			if(locTrkHit->hit->wire->ring == locAmbiguousHitRing)
				continue; // already > 1 stereo hit on this ring: ambiguous as to which hit is best

			//make sure we don't skip too many rings
			int locNumRingsSkipped = locTrkHit->hit->wire->ring - locLastHitRing - 1;
			if(locNumRingsSkipped > int(MAX_NUM_RINGSEED_RINGS_SKIPABLE))
			{
				if(DEBUG_LEVEL > 5)
					cout << "would skip too many rings: actual, max = " << locNumRingsSkipped << ", " << MAX_NUM_RINGSEED_RINGS_SKIPABLE << endl;
				continue; //would skip too many rings
			}

			// see if the hit has small-enough transverse distance to the previous super layer
			DCDCRingSeed locRingSeed;
			locRingSeed.hits.push_back(locTrkHit);
			locRingSeed.ring = locTrkHit->hit->wire->ring;
			locRingSeed.linked = false;
			if(!Attempt_SeedLink(locLastSuperLayerSeed->dCDCRingSeeds.back(), locRingSeed, locLastSuperLayerSeed->dWireOrientation, locWireDirection))
			{
				if(DEBUG_LEVEL > 5)
					cout << "new hit isn't close to wires in previous seed" << endl;
				continue; //new hit isn't close to wires in previous seed
			}
			if(DEBUG_LEVEL > 5)
				cout << "hit is close to wires in previous seed, ring = " << locTrkHit->hit->wire->ring << endl;

			// If axial, require that the hit be near the circle fit. 
			double locDeltaPhi = 9.9E9;
			if((locSearchSuperLayer == 4) || (locSearchSuperLayer == 7))
			{
				// Find the position on the circle that is closest to locTrkHit
				const DVector3 locOrigin = locTrkHit->hit->wire->origin;
				double dx = locOrigin.x() - locFit->x0;
				double dy = locOrigin.y() - locFit->y0;
				double one_over_denom = 1.0/sqrt(dx*dx + dy*dy);
				double x = locFit->x0 + locFit->r0*dx*one_over_denom;
				double y = locFit->y0 + locFit->r0*dy*one_over_denom;
				DVector2 locCirclePosition(x, y);

				// Compare phi values to see if the seeds are close enough to link
				locDeltaPhi = fabs(locCirclePosition.Phi() - locOrigin.Phi());
				while(locDeltaPhi > M_PI)
					locDeltaPhi -= M_TWO_PI;
				locDeltaPhi *= 180.0/M_PI;
				if(DEBUG_LEVEL > 5)
					cout << "hit is axial, check if delta phi is close enough: " << fabs(locDeltaPhi) << ", " << MAX_UNUSED_HIT_LINK_ANGLE << endl;
				if(fabs(locDeltaPhi) > MAX_UNUSED_HIT_LINK_ANGLE)
					continue; //hit is too far away from the current circle fit
			}

			//Have a matching unused hit. Now make sure it is the best one so far. 

			if(locBestTrkHit == NULL)
			{
				//No best hit before, save this result
				locBestTrkHit = locTrkHit;
				locBestDeltaPhi = locDeltaPhi;
				if(DEBUG_LEVEL > 5)
					cout << "brand new track hit, delta phi = " << locBestDeltaPhi << endl;
				continue;
			}

			if(locTrkHit->hit->wire->ring > locBestTrkHit->hit->wire->ring)
			{
				//ring is larger: new hit not as good as the current best one
				if(DEBUG_LEVEL > 5)
					cout << "new hit not as good as the current best one" << endl;
				continue;
			}

			if(locTrkHit->hit->wire->ring < locBestTrkHit->hit->wire->ring)
			{
				//ring is smaller: new hit better than the one we had before
				locAmbiguousHitRing = -1;
				locBestTrkHit = locTrkHit;
				locBestDeltaPhi = locDeltaPhi;
				if(DEBUG_LEVEL > 5)
					cout << "new best track hit (better ring), delta phi = " << locBestDeltaPhi << endl;
				continue;
			}

			//The new hit is on the same ring as the previous best hit. 

			if((locSearchSuperLayer != 4) && (locSearchSuperLayer != 7))
			{
				//stereo: can't tell which hit is best, label ring as ambiguous
				locAmbiguousHitRing = locBestTrkHit->hit->wire->ring;
				locBestTrkHit = NULL;
				if(DEBUG_LEVEL > 5)
					cout << "stereo, can't tell which hit is best, label ring " << locAmbiguousHitRing << " as ambiguous" << endl;
				continue;
			}

			//Axial, see if closest to the track circle (smallest delta phi)
			if(locDeltaPhi >= locBestDeltaPhi)
			{
				//delta-phi is larger: new hit not as good as the current best one
				if(DEBUG_LEVEL > 5)
					cout << "axial, not the best hit, phis = " << locDeltaPhi << ", " << locBestDeltaPhi << endl;
				continue;
			}

			//delta-phi is smaller: new hit better than the current best one: save it
			locBestTrkHit = locTrkHit;
			locBestDeltaPhi = locDeltaPhi;
			if(DEBUG_LEVEL > 5)
				cout << "new best track hit (same ring), delta phi = " << locBestDeltaPhi << endl;
		}
		if(locBestTrkHit == NULL)
			continue; // no hit found for this track (or hits were ambiguous)

		// add best hit to track: create new DCDCSuperLayerSeed for it, add to circle fit
		locBestTrkHit->flags |= USED;
		DCDCRingSeed locRingSeed;
		locRingSeed.hits.push_back(locBestTrkHit);
		locRingSeed.ring = locBestTrkHit->hit->wire->ring;
		locRingSeed.linked = true;

		DCDCSuperLayerSeed* locNewSuperLayerSeed = Get_Resource_CDCSuperLayerSeed();
		locNewSuperLayerSeed->dCDCRingSeeds.push_back(locRingSeed);
		locNewSuperLayerSeed->dSuperLayer = locSearchSuperLayer;
		locNewSuperLayerSeed->dSeedIndex = dSuperLayerSeeds[locSearchSuperLayer - 1].size();
		locNewSuperLayerSeed->linked = true;
		locNewSuperLayerSeed->dWireOrientation = locWireDirection;
		dSuperLayerSeeds[locSearchSuperLayer - 1].push_back(locNewSuperLayerSeed);
		locCDCTrackCircles[loc_i]->Add_LastSuperLayerSeed(locNewSuperLayerSeed);
	}
}

//-----------------------
// Create_TrackCandidiate
//-----------------------
void DTrackCandidate_factory_CDC::Create_TrackCandidiate(DCDCTrackCircle* locCDCTrackCircle)
{
	DVector3 pos, mom;
	if(!Calc_PositionAndMomentum(locCDCTrackCircle, pos, mom))
	{
		if(DEBUG_LEVEL > 5)
			cout << "Track momentum not greater than zero (or NaN), DTrackCandidate object not created." << endl;
		return; //don't create object!!
	}

	DTrackCandidate *locTrackCandidate = new DTrackCandidate;
	locTrackCandidate->setCharge(locCDCTrackCircle->fit->h*dFactorForSenseOfRotation);

	locTrackCandidate->chisq = locCDCTrackCircle->fit->chisq;
	locTrackCandidate->Ndof = locCDCTrackCircle->fit->ndof;
	locTrackCandidate->setPosition(pos);
	locTrackCandidate->setMomentum(mom);

	// Add axial hits (if any)
	deque<DCDCTrkHit*> locHits;
	for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_Axial.size(); ++loc_i)
	{
		locCDCTrackCircle->dSuperLayerSeeds_Axial[loc_i]->Get_Hits(locHits);
		for(size_t loc_j = 0; loc_j < locHits.size(); ++loc_j)
		{
			locTrackCandidate->AddAssociatedObject(locHits[loc_j]->hit);
			locTrackCandidate->used_cdc_indexes.push_back(locHits[loc_j]->index);
		}
	}

	// Add inner stereo hits (if any)
	if(!locCDCTrackCircle->dSuperLayerSeeds_InnerStereo.empty())
	{
		for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[0].size(); ++loc_i) //only one seed series: the "best" one
		{
			locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[0][loc_i]->Get_Hits(locHits);
			for(size_t loc_j = 0; loc_j < locHits.size(); ++loc_j)
			{
				locTrackCandidate->AddAssociatedObject(locHits[loc_j]->hit);
				locTrackCandidate->used_cdc_indexes.push_back(locHits[loc_j]->index);
			}
		}
	}

	// Add outer stereo hits (if any)
	if(!locCDCTrackCircle->dSuperLayerSeeds_OuterStereo.empty())
	{
		for(size_t loc_i = 0; loc_i < locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[0].size(); ++loc_i) //only one seed series: the "best" one
		{
			locCDCTrackCircle->dSuperLayerSeeds_OuterStereo[0][loc_i]->Get_Hits(locHits);
			for(size_t loc_j = 0; loc_j < locHits.size(); ++loc_j)
			{
				locTrackCandidate->AddAssociatedObject(locHits[loc_j]->hit);
				locTrackCandidate->used_cdc_indexes.push_back(locHits[loc_j]->index);
			}
		}
	}

	if(DEBUG_LEVEL>3)
		cout<<"Final Candidate parameters: p="<<mom.Mag()<<" theta="<<mom.Theta()<<"  phi="<<mom.Phi()<<" z="<<pos.Z()<<endl;
	
	_data.push_back(locTrackCandidate);
}

//-------------------------
// Calc_PositionAndMomentum
//-------------------------
bool DTrackCandidate_factory_CDC::Calc_PositionAndMomentum(DCDCTrackCircle* locCDCTrackCircle, DVector3 &pos, DVector3 &mom)
{
	// Get the position and momentum at a fixed radius from the beam line

	// Direction tangent
	double tanl = tan(M_PI_2 - locCDCTrackCircle->dTheta);

	// Squared radius of cylinder outside start counter but inside CDC inner 
	// radius
	double r2 = 90.0;
 
	// Circle parameters
	double xc = locCDCTrackCircle->fit->x0;
	double yc = locCDCTrackCircle->fit->y0;
	double rc = locCDCTrackCircle->fit->r0;
	double rc2 = rc*rc;
	double xc2 = xc*xc;
	double yc2 = yc*yc;
	double xc2_plus_yc2 = xc2 + yc2;

	// variables needed for intersection of circles
	double a = (r2 - xc2_plus_yc2 - rc2)/(2.*rc);
	double temp1 = yc*sqrt(xc2_plus_yc2 - a*a);
	double temp2 = xc*a;
	double cosphi_plus = (temp2 + temp1)/xc2_plus_yc2;
	double cosphi_minus = (temp2 - temp1)/xc2_plus_yc2;

	// Check for intersections
	if(!isfinite(temp1) || !isfinite(temp2))
	{
		// We did not find an intersection between the two circles, so return 
		// sensible defaults for pos and mom
		double my_seed_phi = locCDCTrackCircle->fit->phi;
		double my_center_phi = atan2(yc,xc);
		double xv = xc - rc*cos(my_center_phi);
		double yv = yc - rc*sin(my_center_phi);
		pos.SetXYZ(xv, yv, locCDCTrackCircle->dVertexZ);

		double pt = 0.003*fabs(dMagneticField->GetBz(pos.x(), pos.y(), pos.z()))*rc;
		mom.SetXYZ(pt*cos(my_seed_phi), pt*sin(my_seed_phi), pt*tanl);
		return (mom.Mag() > 0.0);
	}

	// if we have intersections, find both solutions
	double phi_plus = -acos(cosphi_plus);
	double phi_minus = -acos(cosphi_minus);
	double x_plus = xc + rc*cosphi_plus;
	double x_minus = xc + rc*cosphi_minus;
	double y_plus = yc + rc*sin(phi_plus);
	double y_minus = yc + rc*sin(phi_minus);

	// if the resulting radial position on the circle from the fit does not agree
	// with the radius to which we are matching, we have the wrong sign for phi+ 
	// or phi-
	double r2_plus = x_plus*x_plus + y_plus*y_plus;
	double r2_minus = x_minus*x_minus + y_minus*y_minus;
	if(fabs(r2 - r2_plus) > EPS)
	{
		phi_plus *= -1.;
		y_plus = yc + rc*sin(phi_plus);
	}
	if(fabs(r2 - r2_minus) > EPS)
	{
		phi_minus *= -1.;
		y_minus = yc + rc*sin(phi_minus);
	}

	// Choose phi- or phi+ depending on proximity to one of the cdc hits
	DCDCTrkHit* locFirstHit = NULL;
	if(locCDCTrackCircle->dSuperLayerSeeds_Axial.empty())
		locFirstHit = locCDCTrackCircle->dSuperLayerSeeds_InnerStereo[0][0]->dCDCRingSeeds[0].hits[0];
	else
		locFirstHit = locCDCTrackCircle->dSuperLayerSeeds_Axial[0]->dCDCRingSeeds[0].hits[0];

	double xwire = locFirstHit->hit->wire->origin.x();
	double ywire = locFirstHit->hit->wire->origin.y();
	double dx = x_minus - xwire;
	double dy = y_minus - ywire;
	double d2_minus = dx*dx + dy*dy;
	dx = x_plus - xwire;
	dy = y_plus - ywire;
	double d2_plus = dx*dx + dy*dy;
	if(d2_plus > d2_minus)
	{
		phi_minus *= -1.;
		if(locCDCTrackCircle->fit->h < 0)
			phi_minus += M_PI;
		double myphi = atan2(yc, xc);
		double xv = xc - rc*cos(myphi);
		double yv = yc - rc*sin(myphi);
		double dx = x_minus - xv;
		double dy = y_minus - yv;
		double chord = sqrt(dx*dx + dy*dy);
		double two_rc = 2.*rc;
		double ratio = chord/two_rc;
		double ds = (ratio < 1.) ? (two_rc*asin(ratio)) : (two_rc*M_PI_2);
		pos.SetXYZ(x_minus, y_minus, locCDCTrackCircle->dVertexZ + ds*tanl);

		double pt = 0.003*fabs(dMagneticField->GetBz(pos.x(), pos.y(), pos.z()))*rc;
		mom.SetXYZ(pt*sin(phi_minus), pt*cos(phi_minus), pt*tanl);
	}
	else
	{
		phi_plus *= -1.;
		if(locCDCTrackCircle->fit->h < 0)
			phi_plus += M_PI;
		double myphi = atan2(yc, xc);
		double xv = xc - rc*cos(myphi);
		double yv = yc - rc*sin(myphi);
		double dx = x_plus - xv;
		double dy = y_plus - yv;
		double chord = sqrt(dx*dx + dy*dy);
		double two_rc = 2.*rc;
		double ratio = chord/two_rc;
		double ds = (ratio < 1.) ? (two_rc*asin(ratio)) : (two_rc*M_PI_2);
		pos.SetXYZ(x_plus, y_plus, locCDCTrackCircle->dVertexZ + ds*tanl); 

		double pt =0.003*fabs(dMagneticField->GetBz(pos.x(), pos.y(), pos.z()))*rc;
		mom.SetXYZ(pt*sin(phi_plus), pt*cos(phi_plus), pt*tanl);
	}
	return (mom.Mag() > 0.0);
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
	// Delete memory in resource pools
	for(size_t loc_i = 0; loc_i < dCDCTrkHitPool_All.size(); ++loc_i)
		delete dCDCTrkHitPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dCDCSuperLayerSeedPool_All.size(); ++loc_i)
		delete dCDCSuperLayerSeedPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dHelicalFitPool_All.size(); ++loc_i)
		delete dHelicalFitPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dCDCTrackCirclePool_All.size(); ++loc_i)
		delete dCDCTrackCirclePool_All[loc_i];

	return NOERROR;
}

void DTrackCandidate_factory_CDC::DCDCTrkHit::Reset(void)
{
	hit = NULL;
	index = 0;
	flags = NONE;
	var_z = 0.0;
	dStereoHitPos.SetXYZ(0.0, 0.0, 0.0);
	dPhiStereo = 0.0;
	dValidStereoHitPosFlag = false;
}

void DTrackCandidate_factory_CDC::DCDCSuperLayerSeed::Reset(void)
{
	dSuperLayer = 0;
	dSeedIndex = 0;
	linked = false;
	dSpiralLinkParams.clear();
	dCDCRingSeeds.clear();
}

bool DTrackCandidate_factory_CDC::DCDCSuperLayerSeed::Are_AllHitsOnRingShared(const DCDCSuperLayerSeed* locCDCSuperLayerSeed, int locRing) const
{
	deque<DCDCTrkHit*> locRingHits;
	for(size_t loc_i = 0; loc_i < dCDCRingSeeds.size(); ++loc_i)
	{
		if(dCDCRingSeeds[loc_i].ring != locRing)
			continue;
		locRingHits = dCDCRingSeeds[loc_i].hits;
		break;
	}

	deque<DCDCTrkHit*> locRingHits_CompareTo;
	for(size_t loc_i = 0; loc_i < locCDCSuperLayerSeed->dCDCRingSeeds.size(); ++loc_i)
	{
		if(locCDCSuperLayerSeed->dCDCRingSeeds[loc_i].ring != locRing)
			continue;
		locRingHits_CompareTo = locCDCSuperLayerSeed->dCDCRingSeeds[loc_i].hits;
		break;
	}

	return (locRingHits == locRingHits_CompareTo);
}

void DTrackCandidate_factory_CDC::DCDCTrackCircle::Reset(void)
{
	dSuperLayerSeeds_Axial.clear();
	dSuperLayerSeeds_InnerStereo.clear();
	dSuperLayerSeeds_OuterStereo.clear();
	fit = NULL;
	dWeightedChiSqPerDF = 0.0;
	dAverageDriftTime = 0.0;
	HitBitPattern.clear();
	dTheta = 0.0;
	dVertexZ = 0.0;
	dSpiralTurnRing = -1;
	dTruncationSourceCircles.clear();
	dHasNonTruncatedSeedsFlag_InnerStereo = false;
	dHasNonTruncatedSeedsFlag_OuterStereo = false;
	dWeightedChiSqPerDF_Stereo = 0.0;
}

DTrackCandidate_factory_CDC::DCDCSuperLayerSeed* DTrackCandidate_factory_CDC::DCDCTrackCircle::Get_LastSuperLayerSeed(void) const
{
	//only checks the first combination of stereo seeds
	DCDCSuperLayerSeed* locLastAxialSuperLayerSeed = NULL;
	if(!dSuperLayerSeeds_Axial.empty())
	{
		locLastAxialSuperLayerSeed = dSuperLayerSeeds_Axial.back();
		if(locLastAxialSuperLayerSeed->dSuperLayer == 7)
			return locLastAxialSuperLayerSeed;
	}
	if(!dSuperLayerSeeds_OuterStereo.empty())
	{
		DCDCSuperLayerSeed* locLastOuterStereoSuperLayerSeed = dSuperLayerSeeds_OuterStereo[0].back();
		if(locLastOuterStereoSuperLayerSeed != NULL){
			if(locLastAxialSuperLayerSeed == NULL) return locLastOuterStereoSuperLayerSeed;
			if(locLastOuterStereoSuperLayerSeed->dSuperLayer > locLastAxialSuperLayerSeed->dSuperLayer){
				return locLastOuterStereoSuperLayerSeed;
			}else{
				return locLastAxialSuperLayerSeed;
			}
		}
	}
	if(!dSuperLayerSeeds_InnerStereo.empty())
	{
		DCDCSuperLayerSeed* locLastInnerStereoSuperLayerSeed = dSuperLayerSeeds_InnerStereo[0].back();
		if(locLastInnerStereoSuperLayerSeed != NULL){
			if(locLastAxialSuperLayerSeed == NULL) return locLastInnerStereoSuperLayerSeed;
			if(locLastInnerStereoSuperLayerSeed->dSuperLayer > locLastAxialSuperLayerSeed->dSuperLayer){
				return locLastInnerStereoSuperLayerSeed;
			}else{
				return locLastAxialSuperLayerSeed;
			}
		}
	}
	return locLastAxialSuperLayerSeed;
}

DTrackCandidate_factory_CDC::DCDCSuperLayerSeed* DTrackCandidate_factory_CDC::DCDCTrackCircle::Get_SuperLayerSeed(unsigned int locSuperLayer) const
{
	//only checks the first combination of stereo seeds
	if((locSuperLayer == 1) || (locSuperLayer == 4) || (locSuperLayer == 7))
	{
		for(size_t loc_i = 0; loc_i < dSuperLayerSeeds_Axial.size(); ++loc_i)
		{
			if(dSuperLayerSeeds_Axial[loc_i]->dSuperLayer == locSuperLayer)
				return dSuperLayerSeeds_Axial[loc_i];
		}
		return NULL;
	}
	else if((locSuperLayer == 2) || (locSuperLayer == 3))
	{
		if(dSuperLayerSeeds_InnerStereo.empty())
			return NULL;
		for(size_t loc_i = 0; loc_i < dSuperLayerSeeds_InnerStereo[0].size(); ++loc_i)
		{
			if(dSuperLayerSeeds_InnerStereo[0][loc_i]->dSuperLayer == locSuperLayer)
				return dSuperLayerSeeds_InnerStereo[0][loc_i];
		}
	}
	else //outer SL seeds could be listed as inner (e.g. no SL4)
	{
		if(!dSuperLayerSeeds_OuterStereo.empty())
		{
			for(size_t loc_i = 0; loc_i < dSuperLayerSeeds_OuterStereo[0].size(); ++loc_i)
			{
				if(dSuperLayerSeeds_OuterStereo[0][loc_i]->dSuperLayer == locSuperLayer)
					return dSuperLayerSeeds_OuterStereo[0][loc_i];
			}
		}
		if(!dSuperLayerSeeds_InnerStereo.empty())
		{
			for(size_t loc_i = 0; loc_i < dSuperLayerSeeds_InnerStereo[0].size(); ++loc_i)
			{
				if(dSuperLayerSeeds_InnerStereo[0][loc_i]->dSuperLayer == locSuperLayer)
					return dSuperLayerSeeds_InnerStereo[0][loc_i];
			}
		}
	}
	return NULL;
}

void DTrackCandidate_factory_CDC::DCDCTrackCircle::Strip_StereoSuperLayerSeed(unsigned int locSuperLayer)
{
	//assumes at most one stereo series of each type!!!
	deque<DCDCSuperLayerSeed*>::iterator locIterator;
	if(!dSuperLayerSeeds_OuterStereo.empty())
	{
		for(locIterator = dSuperLayerSeeds_OuterStereo[0].begin(); locIterator != dSuperLayerSeeds_OuterStereo[0].end(); ++locIterator)
		{
			if((*locIterator)->dSuperLayer != locSuperLayer)
				continue;
			dSuperLayerSeeds_OuterStereo[0].erase(locIterator);
			if(dSuperLayerSeeds_OuterStereo[0].empty())
				dSuperLayerSeeds_OuterStereo.clear();
			return;
		}
	}
	if(!dSuperLayerSeeds_InnerStereo.empty())
	{
		for(locIterator = dSuperLayerSeeds_InnerStereo[0].begin(); locIterator != dSuperLayerSeeds_InnerStereo[0].end(); ++locIterator)
		{
			if((*locIterator)->dSuperLayer != locSuperLayer)
				continue;
			dSuperLayerSeeds_InnerStereo[0].erase(locIterator);
			if(dSuperLayerSeeds_InnerStereo[0].empty())
				dSuperLayerSeeds_InnerStereo.clear();
			return;
		}
	}
}

void DTrackCandidate_factory_CDC::DCDCTrackCircle::Add_LastSuperLayerSeed(DCDCSuperLayerSeed* locSuperLayerSeed)
{
	//assumes at most one stereo series of each type!!!
	unsigned int locSuperLayer = locSuperLayerSeed->dSuperLayer;
	if((locSuperLayer == 1) || (locSuperLayer == 4) || (locSuperLayer == 7))
		dSuperLayerSeeds_Axial.push_back(locSuperLayerSeed);
	else if((locSuperLayer == 2) || (locSuperLayer == 3))
	{
		if(dSuperLayerSeeds_InnerStereo.empty())
			dSuperLayerSeeds_InnerStereo.resize(1);
		dSuperLayerSeeds_InnerStereo[0].push_back(locSuperLayerSeed);
	}
	else //is super layer 5 or 6
	{
		unsigned int locLastAxialLayer = dSuperLayerSeeds_Axial.empty() ? 0 : dSuperLayerSeeds_Axial.back()->dSuperLayer;
		if(dSuperLayerSeeds_InnerStereo.empty() || (locLastAxialLayer == 4))
		{
			if(dSuperLayerSeeds_OuterStereo.empty())
				dSuperLayerSeeds_OuterStereo.resize(1);
			dSuperLayerSeeds_OuterStereo[0].push_back(locSuperLayerSeed);
		}
		else //put in inner stereo to keep combinations correct
			dSuperLayerSeeds_InnerStereo[0].push_back(locSuperLayerSeed);
	}
}

void DTrackCandidate_factory_CDC::DCDCTrackCircle::Truncate_Circle(unsigned int locNewLastSuperLayer)
{
	//truncate this circle by rejecting all super layers down to locNewLastSuperLayer
		//DOES NOT recycle any memory (these seeds may be used by other circles!!)
	for(size_t loc_i = 0; loc_i < dSuperLayerSeeds_Axial.size(); ++loc_i)
	{
		if(dSuperLayerSeeds_Axial[loc_i]->dSuperLayer <= locNewLastSuperLayer)
			continue;
		((loc_i == 0) ? dSuperLayerSeeds_Axial.clear() : dSuperLayerSeeds_Axial.resize(loc_i));
		break;
	}

	//outer stereo
	deque<deque<DCDCSuperLayerSeed*> >::iterator locIterator, locIterator2;
	if(locNewLastSuperLayer < 5)
		dSuperLayerSeeds_OuterStereo.clear();
	else if(locNewLastSuperLayer < 7)
	{
		bool locClippedEveryStereoSeedFlag = true;
		for(locIterator = dSuperLayerSeeds_OuterStereo.begin(); locIterator != dSuperLayerSeeds_OuterStereo.end();)
		{
			deque<DCDCSuperLayerSeed*>& locSeedSeries = *locIterator;
			bool locClippedstereoSeriesFlag = false;
			for(size_t loc_j = 0; loc_j < locSeedSeries.size(); ++loc_j)
			{
				if(locSeedSeries[loc_j]->dSuperLayer <= locNewLastSuperLayer)
					continue;
				locClippedstereoSeriesFlag = true;
				if(loc_j == 0)
				{
					locSeedSeries.clear();
					break;
				}
				locSeedSeries.resize(loc_j);
				//now, check if another seed series is exactly like this one: if so, clear it
				for(locIterator2 = dSuperLayerSeeds_OuterStereo.begin(); locIterator2 != dSuperLayerSeeds_OuterStereo.end(); ++locIterator2)
				{
					if(locIterator2 == locIterator)
						continue;
					deque<DCDCSuperLayerSeed*>& locSeedSeries2 = *locIterator2;
					if(locSeedSeries == locSeedSeries2)
					{
						locSeedSeries.clear();
						break;
					}
				}
				break;
			}
			if(!locClippedstereoSeriesFlag)
				locClippedEveryStereoSeedFlag = false;
			(locSeedSeries.empty() ? (locIterator = dSuperLayerSeeds_OuterStereo.erase(locIterator)) : ++locIterator);
		}
		if(locClippedEveryStereoSeedFlag)
			dHasNonTruncatedSeedsFlag_OuterStereo = false;
	}

	//inner stereo
	if(locNewLastSuperLayer < 2)
	{
		dSuperLayerSeeds_InnerStereo.clear();
		return;
	}
	else if(locNewLastSuperLayer < 7) //7 instead of 4: could be outer super layers in inner (if SL4 is missing)
	{
		bool locClippedEveryStereoSeedFlag = true;
		for(locIterator = dSuperLayerSeeds_InnerStereo.begin(); locIterator != dSuperLayerSeeds_InnerStereo.end();)
		{
			deque<DCDCSuperLayerSeed*>& locSeedSeries = *locIterator;
			bool locClippedstereoSeriesFlag = false;
			for(size_t loc_j = 0; loc_j < locSeedSeries.size(); ++loc_j)
			{
				if(locSeedSeries[loc_j]->dSuperLayer <= locNewLastSuperLayer)
					continue;
				locClippedstereoSeriesFlag = true;
				if(loc_j == 0)
				{
					locSeedSeries.clear();
					break;
				}
				locSeedSeries.resize(loc_j);
				//now, check if another seed series is exactly like this one: if so, clear it
				for(locIterator2 = dSuperLayerSeeds_InnerStereo.begin(); locIterator2 != dSuperLayerSeeds_InnerStereo.end(); ++locIterator2)
				{
					if(locIterator2 == locIterator)
						continue;
					deque<DCDCSuperLayerSeed*>& locSeedSeries2 = *locIterator2;
					if(locSeedSeries == locSeedSeries2)
					{
						locSeedSeries.clear();
						break;
					}
				}
				break;
			}
			if(!locClippedstereoSeriesFlag)
				locClippedEveryStereoSeedFlag = false;
			(locSeedSeries.empty() ? (locIterator = dSuperLayerSeeds_InnerStereo.erase(locIterator)) : ++locIterator);
		}
		if(locClippedEveryStereoSeedFlag)
			dHasNonTruncatedSeedsFlag_InnerStereo = false;
	}
}

void DTrackCandidate_factory_CDC::DCDCTrackCircle::Absorb_TrackCircle(const DCDCTrackCircle* locTrackCircle)
{
	//used when merging track circles: this track circle will absorb the stereo combinations from the other one
	for(size_t loc_i = 0; loc_i < locTrackCircle->dSuperLayerSeeds_InnerStereo.size(); ++loc_i)
	{
		const deque<DCDCSuperLayerSeed*>& locSeedSeries = locTrackCircle->dSuperLayerSeeds_InnerStereo[loc_i];
		bool locComboAlreadyPresentFlag = false;
		for(size_t loc_j = 0; loc_j < dSuperLayerSeeds_InnerStereo.size(); ++loc_j)
		{
			if(locSeedSeries == dSuperLayerSeeds_InnerStereo[loc_j])
			{
				locComboAlreadyPresentFlag = true;
				break;
			}
			else if((locSeedSeries.size() == 1) && (locSeedSeries[0] == dSuperLayerSeeds_InnerStereo[loc_j][0]))
			{
				locComboAlreadyPresentFlag = true; //is a subset of an existing combo
				break;
			}
		}
		if(!locComboAlreadyPresentFlag)
			dSuperLayerSeeds_InnerStereo.push_back(locSeedSeries);
	}

	for(size_t loc_i = 0; loc_i < locTrackCircle->dSuperLayerSeeds_OuterStereo.size(); ++loc_i)
	{
		const deque<DCDCSuperLayerSeed*>& locSeedSeries = locTrackCircle->dSuperLayerSeeds_OuterStereo[loc_i];
		bool locComboAlreadyPresentFlag = false;
		for(size_t loc_j = 0; loc_j < dSuperLayerSeeds_OuterStereo.size(); ++loc_j)
		{
			if(locSeedSeries == dSuperLayerSeeds_OuterStereo[loc_j])
			{
				locComboAlreadyPresentFlag = true;
				break;
			}
			else if((locSeedSeries.size() == 1) && (locSeedSeries[0] == dSuperLayerSeeds_OuterStereo[loc_j][0]))
			{
				locComboAlreadyPresentFlag = true; //is a subset of an existing combo
				break;
			}
		}
		if(!locComboAlreadyPresentFlag)
			dSuperLayerSeeds_OuterStereo.push_back(locSeedSeries);
	}

	if(locTrackCircle->dHasNonTruncatedSeedsFlag_InnerStereo)
		dHasNonTruncatedSeedsFlag_InnerStereo = true;
	if(locTrackCircle->dHasNonTruncatedSeedsFlag_OuterStereo)
		dHasNonTruncatedSeedsFlag_OuterStereo = true;

	for(size_t loc_i = 0; loc_i < locTrackCircle->dTruncationSourceCircles.size(); ++loc_i)
	{
		bool locIsAlreadyTruncationSourceFlag = false;
		for(size_t loc_j = 0; loc_j < dTruncationSourceCircles.size(); ++loc_j)
		{
			if(locTrackCircle->dTruncationSourceCircles[loc_i] != dTruncationSourceCircles[loc_j])
				continue;
			locIsAlreadyTruncationSourceFlag = true;
			break;
		}
		if(!locIsAlreadyTruncationSourceFlag)
			dTruncationSourceCircles.push_back(locTrackCircle->dTruncationSourceCircles[loc_i]);
	}
}

bool DTrackCandidate_factory_CDC::DCDCTrackCircle::Check_IfInputIsSubset(const DCDCTrackCircle* locTrackCircle)
{
	//returns false if they are effectively identical
	if(locTrackCircle->dSuperLayerSeeds_Axial.empty())
		return false;
	if(locTrackCircle->Get_LastSuperLayerSeed()->dSuperLayer == 7)
		return false; //can't be subset if hasn't been truncated yet
	if(locTrackCircle->dSuperLayerSeeds_Axial.size() >= dSuperLayerSeeds_Axial.size())
		return false; //can't be subset if same # or greater of axial super layers (if identical should merge them instead)

	for(size_t loc_i = 0; loc_i < locTrackCircle->dSuperLayerSeeds_Axial.size(); ++loc_i)
	{
		DCDCSuperLayerSeed* locSuperLayerSeed = locTrackCircle->dSuperLayerSeeds_Axial[loc_i];
		if(Get_SuperLayerSeed(locSuperLayerSeed->dSuperLayer) != locSuperLayerSeed)
			return false;
	}
	//all axial seeds are a subset, now check the stereo dHasNonTruncatedSeedsFlag flags
	DCDCSuperLayerSeed* locLastSuperLayerSeed = locTrackCircle->Get_LastSuperLayerSeed();
	if((locLastSuperLayerSeed->dSuperLayer == 1) || (locLastSuperLayerSeed->dSuperLayer == 4))
		return true;

	DCDCSuperLayerSeed* locLastAxialSuperLayerSeed = locTrackCircle->dSuperLayerSeeds_Axial.back();
	if(locLastAxialSuperLayerSeed->dSuperLayer == 4)
		return locTrackCircle->dHasNonTruncatedSeedsFlag_OuterStereo;
	else
		return locTrackCircle->dHasNonTruncatedSeedsFlag_InnerStereo;
}

void DTrackCandidate_factory_CDC::DCDCTrackCircle::Get_AllStereoSuperLayerSeeds(deque<DCDCSuperLayerSeed*>& locStereoSuperLayerSeeds)
{
	//assumes there is only one seed series
	locStereoSuperLayerSeeds.clear();
	if(!dSuperLayerSeeds_InnerStereo.empty())
		locStereoSuperLayerSeeds = dSuperLayerSeeds_InnerStereo[0];
	if(!dSuperLayerSeeds_OuterStereo.empty())
		locStereoSuperLayerSeeds.insert(locStereoSuperLayerSeeds.begin(), dSuperLayerSeeds_OuterStereo[0].begin(), dSuperLayerSeeds_OuterStereo[0].end());
}

unsigned int DTrackCandidate_factory_CDC::DCDCTrackCircle::Get_NumStereoSuperLayerSeeds(void)
{
	//assumes there is only one seed series
	unsigned int locNumStereoSuperLayerSeeds = 0;
	if(!dSuperLayerSeeds_InnerStereo.empty())
		locNumStereoSuperLayerSeeds += dSuperLayerSeeds_InnerStereo[0].size();
	if(!dSuperLayerSeeds_OuterStereo.empty())
		locNumStereoSuperLayerSeeds += dSuperLayerSeeds_OuterStereo[0].size();
	return locNumStereoSuperLayerSeeds;
}

//-------------------------------------------------------------------------
// Routines for fitting a line to the stereo data
//-------------------------------------------------------------------------
// Compute the chi^2 for a line fit given errors in both s and z.  Also 
// computes current best guess for the scaled intercept z0. 
// Taken from Numerical Recipes in C (2nd ed.), p. 670.
double DTrackCandidate_factory_CDC::DCDCLineFit::ChiXY(double lambda)
{
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
void DTrackCandidate_factory_CDC::DCDCLineFit::BracketMinimumChisq(double &a,double &b,double &c,double &chi2a,double &chi2b,double &chi2c)
{
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
double DTrackCandidate_factory_CDC::DCDCLineFit::FindMinimumChisq(double ax,double bx, double cx, double &xmin)
{
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


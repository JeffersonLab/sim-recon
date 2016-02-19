// $Id$
//
//    File: DTrackCandidate_factory_CDC.h
// Created: Thu Sep  6 14:47:48 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#ifndef _DTrackCandidate_factory_CDC_
#define _DTrackCandidate_factory_CDC_

#include <map>
#include <deque>
using namespace std;

#include "TDirectory.h"

#include <JANA/JFactory.h>
using namespace jana;

#include "DTrackCandidate.h"
#include "DHelicalFit.h"
#include "CDC/DCDCTrackHit.h"
#include <DVector3.h>
#include <HDGEOMETRY/DGeometry.h>
#include <HDGEOMETRY/DMagneticFieldMap.h>

class DTrackCandidate_factory_CDC : public JFactory<DTrackCandidate>
{
	public:
		DTrackCandidate_factory_CDC(){};
		~DTrackCandidate_factory_CDC();
		const char* Tag(void){return "CDC";}

		enum trk_flags_t
		{
			NONE					= 0x000,
			NOISE					= 0x001,
			USED					= 0x002, //set if used in a super layer seed (whether that seed was rejected later or not)
			OUT_OF_TIME		   = 0x004
		};

		enum wire_direction_t
		{
			WIRE_DIRECTION_AXIAL = 0, 
			WIRE_DIRECTION_STEREOLEFT = 1, //rotated by 6 degrees counter-clockwise
			WIRE_DIRECTION_STEREORIGHT = 2 //rotated by 6 degrees clockwise
		};

		class DCDCTrkHit
		{
			public:
				void Reset(void);
				inline double Dist2(DCDCTrkHit* trkhit) const{return (trkhit->hit->wire->origin - this->hit->wire->origin).Mag2();}

				const DCDCTrackHit* hit;
				unsigned int index;
				unsigned int flags;
				float var_z;
				DVector3 dStereoHitPos;
				float dPhiStereo;
				bool dValidStereoHitPosFlag; //false if prior to calc, or if hit doesn't intersect circle
		};

		// Collection of adjacent DCDCTrkHit's on a ring. 
		class DCDCRingSeed
		{
			public:
				//hits are stored in order from the first straw to the last straw
					//generally, this in order from lowest straw # to highest straw #
					//unless it crosses straw = 1 barrier: then first straw is lowest-# straw with phi > pi, and last straw is highest-# straw with phi < pi
				deque<DCDCTrkHit*> hits;
				bool linked;
				int ring;
		};

		// Collection of information about a given potential spiral turn
		typedef struct
		{
			int dSpiralTurnRingFlag; //-1 if potential turn looks like the track is turning back outwards, 1 if turning back inwards
			int dDefiniteSpiralTurnRingFlag; //same as dSpiralTurnRingFlag if DEFINITE turn, 0 otherwise
			int dSpiralTurnRing; //is ring# of turn
		} DSpiralParams_t;

		// Collection of adjacent DCDCRingSeed's in a super layer
		class DCDCSuperLayerSeed
		{
			public:
				void Reset(void);
				bool Are_AllHitsOnRingShared(const DCDCSuperLayerSeed* locCDCSuperLayerSeed, int locRing) const;
				inline void Get_Hits(deque<DCDCTrkHit*>& locHits) const
				{
					locHits.clear();
					for(size_t loc_i = 0; loc_i < dCDCRingSeeds.size(); ++loc_i)
						locHits.insert(locHits.end(), dCDCRingSeeds[loc_i].hits.begin(), dCDCRingSeeds[loc_i].hits.end());
				}

				deque<DCDCRingSeed> dCDCRingSeeds; //stored in order from innermost (lowest) ring to outermost (highest) ring
				unsigned int dSuperLayer;
				unsigned int dSeedIndex;
				bool linked;
				wire_direction_t dWireOrientation;
				map<int, DSpiralParams_t> dSpiralLinkParams; //key is the dSeedIndex of the DCDCSuperLayerSeed in this super layer it is linked to (can point to itself (self-linked)!)
		};

		// Collection of adjacent DCDCSuperLayerSeed's in the CDC
			// Each unique combination of axial DCDCSuperLayerSeed's has it's own DCDCTrackCircle object
			// Every possible combination of stereo DCDCSuperLayerSeed's used to link these axial DCDCSuperLayerSeed's together is stored in the 2D-deques
				// This is until the best combination is found: then the unused combinations are cleared and only one combination will remain
		class DCDCTrackCircle
		{
			public:
				void Reset(void);
				DCDCSuperLayerSeed* Get_LastSuperLayerSeed(void) const;
				DCDCSuperLayerSeed* Get_SuperLayerSeed(unsigned int locSuperLayer) const;
				void Strip_StereoSuperLayerSeed(unsigned int locSuperLayer);
				void Add_LastSuperLayerSeed(DCDCSuperLayerSeed* locSuperLayerSeed);
				void Truncate_Circle(unsigned int locNewLastSuperLayer);
				void Absorb_TrackCircle(const DCDCTrackCircle* locTrackCircle);
				bool Check_IfInputIsSubset(const DCDCTrackCircle* locTrackCircle);
				void Get_AllStereoSuperLayerSeeds(deque<DCDCSuperLayerSeed*>& locStereoSuperLayerSeeds);
				unsigned int Get_NumStereoSuperLayerSeeds(void);

				deque<DCDCSuperLayerSeed*> dSuperLayerSeeds_Axial;
				//for dSuperLayerSeeds_InnerStereo and dSuperLayerSeeds_OuterStereo:
					//each deque<DCDCSuperLayerSeed*> is a series of adjacent stereo seeds that could potentially belong together
					//for example: 
						// dSuperLayerSeeds_InnerStereo[0] is likely (not necessarily) a deque with size = 2: 
							//the first in this deque is a DCDCSuperLayerSeed from super layer 2 (SL2), the second is from SL3 and is definitely adjacent to the one from SL2.
						// dSuperLayerSeeds_InnerStereo[1] is also likely a deque with size = 2, but a different combination of adjacent super layer seeds than the previous one. 
						// once the best seeds have been selected, dSuperLayerSeeds_InnerStereo will have size = 1. 
					//the reason that InnerStereo and OuterStereo are kept separately is to keep track of all possible combinations of adjacent super layers
						//that way they are evaluated together when determining which stereo super layer seeds are best
					//Note: if super layer 4 is missing (e.g. dead high voltage board), but inner (and/or outer) super layers are present, all super layers will be "inner"
						//this is so that the stereo super layer combinatorics are evaluated correctly in the Find_ThetaZ() function
					//Note: if super layers 1 through 3 are missing (e.g. a decay product), then dSuperLayerSeeds_InnerStereo will be empty
				deque<deque<DCDCSuperLayerSeed*> > dSuperLayerSeeds_InnerStereo;
				deque<deque<DCDCSuperLayerSeed*> > dSuperLayerSeeds_OuterStereo;

				DHelicalFit* fit; //the circle fit
				float dWeightedChiSqPerDF; //of circle fit //weighted: is (chisq/ndf)/(#axial_super_layers^2): prefer fits with more axial super layers (not necessarily more hits)
				float dWeightedChiSqPerDF_Stereo; //of theta/z determination //weighted: is (chisq/ndf)/(#axial_super_layers^2): prefer fits with more stereo super layers (not necessarily more hits)
				float dAverageDriftTime; //of axial hits close to the circle fit
				deque<unsigned int> HitBitPattern; //bit pattern of hits in the track circle (first is just axial, then is all)
				float dTheta;
				float dVertexZ;
				int dSpiralTurnRing; //is ring# of confirmed match spiral turn, -1 if no turn

				//dTruncationSourceCircles: if this object's (not this member's) circle was truncated, this member points to the track circles that indicated truncation was necessary
					//more than 1 possible if two circles were merged after truncation:
						//e.g. two tracks have SL1 & SL4 but different SL7, and their SL7's are rejected by other tracks who have the same SL7s
				deque<const DCDCTrackCircle*> dTruncationSourceCircles;
				bool dHasNonTruncatedSeedsFlag_InnerStereo; //true if any inner stereo seeds are present and are unique: not from a truncated source
				bool dHasNonTruncatedSeedsFlag_OuterStereo; //true if any outer stereo seeds are present and are unique: not from a truncated source
		};

		class DCDCLineFit
		{
			public:
				void BracketMinimumChisq(double &a,double &b, double &c, double &chi2a, double &chi2b, double &chi2c);
				double FindMinimumChisq(double a, double b, double c, double &lambda);
				double ChiXY(double lambda);

				unsigned int n;
				vector<float> s;
				vector<float> var_s;
				vector<float> z;
				vector<float> var_z;
				vector<float> w;
				float z0;
				float tanl;
		};

		typedef struct
		{
			float x, y, perp2;
			float z, var_z;
		} intersection_t;

		typedef deque<deque<DCDCRingSeed> >::iterator ringiter;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		// Utility Functions
		void Reset_Pools(void);
		DCDCTrkHit* Get_Resource_CDCTrkHit(void);
		DCDCSuperLayerSeed* Get_Resource_CDCSuperLayerSeed(void);
		DHelicalFit* Get_Resource_HelicalFit(void);
		DCDCTrackCircle* Get_Resource_CDCTrackCircle(void);

		// Make Super Layer Seeds
		jerror_t Get_CDCHits(JEventLoop* loop);
		void Find_SuperLayerSeeds(deque<DCDCTrkHit*>& locSuperLayerHits, unsigned int locSuperLayer);
		void Link_RingSeeds(deque<DCDCRingSeed*>& parent, ringiter ring, ringiter ringend, unsigned int locSuperLayer, unsigned int locNumPreviousRingsWithoutHit);
		double MinDist2(const DCDCRingSeed& locInnerRingSeed, const DCDCRingSeed& locOuterRingSeed);
		double MinDist2(const deque<DCDCTrkHit*>& locInnerSeedHits, const deque<DCDCTrkHit*>& locOuterSeedHits);
		void Reject_SuperLayerSeeds_HighSeedDensity(unsigned int locSuperLayer);
		void Calc_SuperLayerPhiRange(DCDCSuperLayerSeed* locSuperLayerSeed, double& locSeedFirstPhi, double& locSeedLastPhi);
		bool Check_IfPhiRangesOverlap(double locFirstSeedPhi, double locLastSeedPhi, double locTargetFirstPhi, double locTargetLastPhi);

		// Search for spirals
		void Set_SpiralLinkParams(void);
		bool SearchFor_SpiralTurn_TwoSeedsSharingManyHits(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2);
		bool SearchFor_SpiralTurn_TwoSeedsSharingFewHits(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2);
		bool SearchFor_SpiralTurn_ManyHitsAdjacentRing(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2, int locRingToCheck, int locMinStrawsAdjacentRing, int& locMaxSpiralNumHits);
		bool SearchFor_SpiralTurn_MissingOrBetweenRings(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2);
		bool SearchFor_SpiralTurn_SingleSeed(DCDCSuperLayerSeed* locSuperLayerSeed);
		void Print_SuperLayerSeeds(void);

		//Link Super Layers to Create DCDCTrackCircle Objects
		bool Build_TrackCircles(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		bool Link_SuperLayers(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer);
		void Link_SuperLayers_FromAxial(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer, unsigned int locInnerSuperLayer);
		void Link_SuperLayers_FromStereo(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer, unsigned int locInnerSuperLayer);
		bool Link_SuperLayers_FromStereo_ToAxial(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer, unsigned int locInnerSuperLayer);
		void Link_SuperLayers_FromStereo_ToStereo(deque<DCDCTrackCircle*>& locCDCTrackCircles, unsigned int locOuterSuperLayer, unsigned int locInnerSuperLayer);
		bool Check_IfShouldAttemptLink(const DCDCSuperLayerSeed* locSuperLayerSeed, bool locInnerSeedFlag);
		bool Attempt_SeedLink(DCDCSuperLayerSeed* locSuperLayerSeed1, DCDCSuperLayerSeed* locSuperLayerSeed2);
		bool Attempt_SeedLink(DCDCRingSeed& locRingSeed1, DCDCRingSeed& locRingSeed2, wire_direction_t locWireDirection1, wire_direction_t locWireDirection2);
		void Recycle_DCDCTrackCircle(DCDCTrackCircle* locCDCTrackCircle);

		// Continue DCDCTrackCircle Creation
		void Print_TrackCircles(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		void Print_TrackCircle(DCDCTrackCircle* locCDCTrackCircle);
		void Reject_DefiniteSpiralArms(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		void Drop_IncompleteGroups(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		void Fit_Circles(deque<DCDCTrackCircle*>& locCDCTrackCircles, bool locFitOnlyIfNullFitFlag, bool locAddStereoLayerIntersectionsFlag, bool locFitDuringLinkingFlag = false);
		DVector3 Find_IntersectionBetweenSuperLayers(const DCDCSuperLayerSeed* locInnerSuperLayerSeed, const DCDCSuperLayerSeed* locOuterSuperLayerSeed);

		// Filter Track Circles and Stereo Wires
		void Handle_StereoAndFilter(deque<DCDCTrackCircle*>& locCDCTrackCircles, bool locFinalPassFlag);
		void Truncate_TrackCircles(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		void Set_HitBitPattern_Axial(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		void Filter_TrackCircles_Axial(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		void Create_NewCDCSuperLayerSeeds(DCDCTrackCircle* locCDCTrackCircle);
		DCDCSuperLayerSeed* Create_NewStereoSuperLayerSeed(DCDCSuperLayerSeed* locCDCSuperLayerSeed, const DCDCTrackCircle* locCDCTrackCircle, map<DCDCTrkHit*, DCDCTrkHit*>& locProjectedStereoHitMap);
		bool Calc_StereoPosition(const DCDCWire *wire, const DHelicalFit* fit, DVector3 &pos, double &var_z, double& locPhiStereo, double d = 0.0);

		// Select Best Stereo Seeds and Calculate Theta & Z
		bool Select_CDCSuperLayerSeeds(DCDCTrackCircle* locCDCTrackCircle, bool locFinalPassFlag);
		void Select_ThetaZStereoHits(const DCDCTrackCircle* locCDCTrackCircle, int locInnerSeedSeriesIndex, int locOuterSeedSeriesIndex, bool locFinalPassFlag, deque<DCDCTrkHit*>& locComboHits);
		void Calc_StereoHitDeltaPhis(unsigned int locSuperLayer, deque<DCDCTrkHit*>& locHits, const DCDCTrackCircle* locCDCTrackCircle, deque<pair<DCDCTrkHit*, double> >& locDeltaPhis);
		double MinDeltaPhi(const deque<DCDCTrkHit*>& locInnerSeedHits, const deque<DCDCTrkHit*>& locOuterSeedHits);
		void Recycle_DCDCSuperLayerSeed(DCDCSuperLayerSeed* locCDCSuperLayerSeed);
		bool Find_ThetaZ(const DHelicalFit* locFit, const deque<DCDCTrkHit*>& locStereoHits, double& locTheta, double& locZ, double& locChiSqPerNDF);
		bool Find_ThetaZ_Regression(const DHelicalFit* locFit, const deque<DCDCTrkHit*>& locStereoHits, double& locTheta, double& locZ, double& locChiSqPerNDF);
		bool Find_Theta(const DHelicalFit* locFit, const deque<DCDCTrkHit*>& locStereoHits, double& locTheta, double& locThetaMin, double& locThetaMax, double& locChiSqPerNDF);
		bool Find_Z(const DHelicalFit* locFit, const deque<DCDCTrkHit*>& locStereoHits, double locThetaMin, double locThetaMax, double& locZ);
		void Set_HitBitPattern_All(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		void Filter_TrackCircles_Stereo(deque<DCDCTrackCircle*>& locCDCTrackCircles);

		// Finalize and Create DTrackCandidate Objects
		void Add_UnusedHits(deque<DCDCTrackCircle*>& locCDCTrackCircles);
		void Create_TrackCandidiate(DCDCTrackCircle* locCDCTrackCircle);
		bool Calc_PositionAndMomentum(DCDCTrackCircle* locCDCTrackCircle, DVector3 &pos, DVector3 &mom);

		int DEBUG_LEVEL;
		unsigned int MAX_ALLOWED_CDC_HITS;
		unsigned int MAX_NUM_RINGSEED_RINGS_SKIPABLE;
		double MAX_HIT_DIST; // cm
		double MAX_HIT_DIST2; // cm
		unsigned int MIN_SEED_HITS;
		unsigned int MAX_STRAWS_BETWEEN_LINK_SPIRAL_TURN;
		unsigned int MIN_STRAWS_POTENTIAL_SPIRAL_TURN;
		unsigned int MIN_STRAWS_DEFINITE_SPIRAL_TURN;
		unsigned int MIN_STRAWS_ADJACENT_TO_SPIRAL_TURN;
		unsigned int SEED_DENSITY_BIN_STRAW_WIDTH;
		unsigned int MAX_SEEDS_IN_STRAW_BIN;
		unsigned int MAX_ALLOWED_TRACK_CIRCLES; //bail if exceed this many
		bool ENABLE_DEAD_HV_BOARD_LINKING;
		unsigned int MAX_SUPERLAYER_NEW_TRACK; //don't allow new track seeds to start after this super layer //track seeds could start late if track is a decay product, or HV board is dead
		double MAX_COMMON_HIT_FRACTION;
		double MIN_CIRCLE_ASYMMETRY;
		double MAX_DRIFT_TIME;
		double MAX_SEED_TIME_DIFF;
		unsigned int MIN_PRUNED_STEREO_HITS;
		double MAX_UNUSED_HIT_LINK_ANGLE;

		size_t MAX_DCDCTrkHitPoolSize;
		size_t MAX_DCDCSuperLayerSeedPoolSize;
		size_t MAX_DCDCTrackCirclePoolSize;
		size_t MAX_HelicalFitPoolSize;

		double TARGET_Z;
		double VERTEX_Z_MIN;
		double VERTEX_Z_MAX;

		map<DCDCTrkHit*, unsigned int> dStereoHitNumUsedMap; //map of circle-projected (hit-z-group) hit to # DCDCSuperLayerSeeds it's in (when drops to 0, can recycle the hit memory)

		// Resource Pools for saving time and re-using memory
		deque<DCDCTrkHit*> dCDCTrkHitPool_Available;
		deque<DCDCTrkHit*> dCDCTrkHitPool_All;

		deque<DCDCSuperLayerSeed*> dCDCSuperLayerSeedPool_Available;
		deque<DCDCSuperLayerSeed*> dCDCSuperLayerSeedPool_All;

		deque<DCDCTrackCircle*> dCDCTrackCirclePool_Available;
		deque<DCDCTrackCircle*> dCDCTrackCirclePool_All;

		deque<DHelicalFit*> dHelicalFitPool_All;
		deque<DHelicalFit*> dHelicalFitPool_Available;

		deque<unsigned int> dNumStrawsPerRing; //index is ring index
		deque<unsigned int> superlayer_boundaries;

		unsigned int dNumCDCHits;
		deque<DCDCTrkHit*> cdctrkhits;
		deque<deque<DCDCTrkHit*> > cdchits_by_superlayer;
		deque<deque<DCDCSuperLayerSeed*> > dSuperLayerSeeds; //index 0 -> 6 is super layer 1 -> 7

		const DMagneticFieldMap* dMagneticField;
		double dFactorForSenseOfRotation;

		unsigned int dNumSeedDensityPhiBins;
		//dRejectedPhiRegions: due to hit density being too high in that region
		map<unsigned int, deque<pair<double, double> > > dRejectedPhiRegions; //key is super layer, pair is first/last phi (not necessarily min/max! (could pass thorugh phi = 0 barrier))
};

#endif // _DTrackCandidate_factory_CDC_


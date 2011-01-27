// $Id$
//
//    File: DTrackCandidate_factory_CDC.h
// Created: Thu Sep  6 14:47:48 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#ifndef _DTrackCandidate_factory_CDC_
#define _DTrackCandidate_factory_CDC_

#include <JANA/JFactory.h>
using namespace jana;

#include "DTrackCandidate.h"
#include "DHelicalFit.h"
#include "CDC/DCDCTrackHit.h"

class DFDCPseudo;

class DTrackCandidate_factory_CDC:public JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory_CDC(){};
		~DTrackCandidate_factory_CDC(){};
		const char* Tag(void){return "CDC";}
		
		enum trk_flags_t{
			NONE					= 0x000,
			NOISE					= 0x001,
			USED					= 0x002,
			IS_STEREO			= 0x004,
			CANT_BE_IN_SEED	= 0x008,
			IN_SEED				= 0x010,
			IN_LINE				= 0x020,
			IN_TRACK				= 0x040,
			VALID_STEREO		= 0x080,
			OUT_OF_TIME			= 0x100
		};
		
		enum ret_cond_t{
			FIT_OK			= 0,
			NO_SEED,
			BAD_SEED,
			FIND_FAILED
		};
		
		class DCDCTrkHit{
			public:
				const DCDCTrackHit *hit;
				unsigned int flags;
				double x_stereo;
				double y_stereo;
				double z_stereo;
				double phi_stereo;
				
				double Dist2(DCDCTrkHit *trkhit){
					DVector3 d = trkhit->hit->wire->origin - this->hit->wire->origin;
					return d.Mag2();
				}
		};
		
		class DCDCSeed{
			public:
				vector<DCDCTrkHit*> hits;
				vector<DCDCTrkHit> stereo_hits;
				vector<const DFDCPseudo*> fdchits;
				double phi_avg;
				double tdrift_avg;
				bool linked;
				bool valid;
				DHelicalFit fit;
				double theta;
				double z_vertex;
				double q;
				double theta_min, theta_max;
				double z_min, z_max;
				void Merge(DCDCSeed& seed);
				double MinDist2(DCDCSeed& seed);
				double FindAverageBz(JEventLoop *loop);
				
				DCDCSeed();
		};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		vector<DCDCTrkHit*> cdctrkhits;
		vector<vector<DCDCTrkHit*> > cdchits_by_superlayer;
		vector<DCDCTrkHit*> seedhits;
		
		typedef vector<vector<DCDCSeed > >::iterator ringiter;
		
		void GetCDCHits(JEventLoop *loop);
		void FindSeeds(vector<DCDCTrkHit*> &hits, vector<DCDCSeed> &seeds);
		void LinkSubSeeds(vector<DCDCSeed*> &parent, ringiter ring, ringiter ringend, vector<DCDCSeed> &seeds);
		void LinkSeeds(vector<DCDCSeed> &in_seeds1, vector<DCDCSeed> &in_seeds2, vector<DCDCSeed> &seeds, unsigned int max_linked_hits);
		bool FitCircle(DCDCSeed &seed);
		void PickupUnmatched(vector<DCDCSeed> &seeds);
		void DropIncompleteSeeds(vector<DCDCSeed> &seeds);
		void FilterCloneSeeds(vector<DCDCSeed> &seeds);
		void AddStereoHits(vector<DCDCTrkHit*> &stereo_hits, DCDCSeed &seed);
		void FindThetaZ(DCDCSeed &seed);
		jerror_t FindThetaZRegression(DCDCSeed &seed);
		void FindThetaZStraightTrack(DCDCSeed &seed);
		void FindTheta(DCDCSeed &seed, double target_z_min, double target_z_max);
		void FindZ(DCDCSeed &seed, double theta_min, double theta_max);
		int NumEligibleSeedHits(vector<DCDCTrkHit*> &hits);

		vector<int> superlayer_boundaries;

		unsigned int MAX_ALLOWED_CDC_HITS;
		unsigned int MAX_SUBSEED_STRAW_DIFF;
		unsigned int MIN_SEED_HITS;
		unsigned int MAX_SUBSEED_LINKED_HITS;
		unsigned int MAX_RING_SUBSEED_HITS;
		double MAX_HIT_DIST; // cm
		double MAX_HIT_DIST2; // cm
		double MAX_SEED_TIME_DIFF; // ns
		double MAX_CDC_MATCH_ANGLE; // degrees
		double MAX_FDC_MATCH_ANGLE; // degrees
		double MAX_SEED_LINK_ANGLE;
		double TARGET_Z_MIN;
		double TARGET_Z_MAX;
		int DEBUG_LEVEL;
};

#endif // _DTrackCandidate_factory_CDC_


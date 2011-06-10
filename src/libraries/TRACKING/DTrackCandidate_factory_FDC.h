// $Id$
//
//    File: DTrackCandidate_factory_FDC.h
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_factory_FDC_
#define _DTrackCandidate_factory_FDC_

#include <TH1.h>
#include <TH2F.h>
#include <TH3F.h>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
using namespace jana;

#include "DQuickFit.h"
#include "DHoughFind.h"
#include "DTrackCandidate.h"
#include "FDC/DFDCIntersection.h"
#include "FDC/DFDCWire.h"

class DMagneticFieldMap;

/// This is an alternate FDC track finder that is not used as part of the default
/// reconstruction. It is no longer maintained, but is kept around as an independent
/// check against the default FDC finder in DTrackCandiate_factory_FDCpseudo .
///
/// This was originally written to help study the possibility of a wires only
/// design of the FDC. As such, it uses the DFDCIntersection objects as inputs
/// which have a far worse position resolution than what is achievable with the
/// the cathode strips. Consequently, this tended to produce more ghost tracks. 

class DTrackCandidate_factory_FDC:public JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory_FDC();
		~DTrackCandidate_factory_FDC(){};
		virtual const char* Tag(void){return "FDC";}

		enum trk_flags_t{
			NONE					= 0x000,
			NOISE					= 0x001,
			USED					= 0x002,
			CANT_BE_IN_SEED	= 0x008,
			ON_CIRCLE			= 0x010,
			IN_THETA_RANGE		= 0x020,
			IN_Z_RANGE			= 0x040,
			VALID_HIT			= 0x080,
			OUT_OF_TIME			= 0x100
		};
		
		enum ret_cond_t{
			FIT_OK			= 0,
			NO_SEED,
			BAD_SEED,
			FIND_FAILED
		};
		
		class DFDCTrkHit{
			public:
				const DFDCIntersection *hit;
				double phi_hit;
				unsigned int flags;
				double theta_min;
				double theta_max;
				double zmin;
				double zmax;
				
				double Dist2(const DFDCTrkHit* trkhit){
					DVector3 delta = trkhit->hit->pos - this->hit->pos;
					//return delta.Mag2();
					double dx = trkhit->hit->pos.X() - this->hit->pos.X();
					double dy = trkhit->hit->pos.Y() - this->hit->pos.Y();
					return dx*dx + dy*dy;
				}
		};
		
		class DFDCSeed{
			public:
				vector<DFDCTrkHit*> hits;
				double p_trans;
				double tdrift_avg;
				bool valid;
				double r0, x0, y0;
				double phi;
				double theta;
				double z_vertex;
				double q;
				double theta_min, theta_max;
				double z_min, z_max;
				void Merge(DFDCSeed& seed);
				double MinDist2(DFDCSeed& seed);
		};


	protected:
		virtual jerror_t init(void);
		virtual jerror_t brun(JEventLoop *loop, int runnumber);
		virtual jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		virtual jerror_t fini(void);	///< Invoked via JEventProcessor virtual method

		DHoughFind hough;

		vector<DFDCTrkHit*> fdctrkhits;

		void GetTrkHits(JEventLoop *loop);
		void FindSeeds(vector<DFDCSeed> &seeds);
		void FillSeedHits(DFDCSeed &seed);
		unsigned int NumAvailableHits(void);
		void FindThetaZ(DFDCSeed &seed);
		void FindTheta(DFDCSeed &seed, double target_z_min, double target_z_max);
		void FindZ(DFDCSeed &seed, double theta_min, double theta_max);

		float TARGET_Z_MIN;
		float TARGET_Z_MAX;
		double MAX_HIT_DIST;
		double MAX_HIT_DIST2;

};

#endif // _DTrackCandidate_factory_FDC_


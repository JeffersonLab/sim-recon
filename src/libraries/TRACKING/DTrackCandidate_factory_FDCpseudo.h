// $Id: DTrackCandidate_factory_FDCpseudo.h 2955 2007-11-11 04:54:05Z davidl $
//
//    File: DTrackCandidate_factory_FDCpseudo.h
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_factory_FDCpseudo_
#define _DTrackCandidate_factory_FDCpseudo_

#include <TH1.h>
#include <TH2F.h>
#include <TH3F.h>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
using namespace jana;

#include "DQuickFit.h"
#include "DTrackCandidate.h"
#include "FDC/DFDCPseudo.h"
#include "FDC/DFDCWire.h"
#include "DHoughFind.h"

class DMagneticFieldMap;

class DTrackCandidate_factory_FDCpseudo:public JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory_FDCpseudo();
		~DTrackCandidate_factory_FDCpseudo(){};
		const string toString(void);
		virtual const char* Tag(void){return "FDCpseudo";}

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
		  const DFDCPseudo *hit;
		  double phi_hit;
		  unsigned int flags;
		  double theta_min;
		  double theta_max;
		  double zmin;
		  double zmax;
		  
		  double Dist2(const DFDCTrkHit* trkhit){
		    //double dx = trkhit->hit->x - this->hit->x;
		    //double dy = trkhit->hit->y - this->hit->y;
		    //return dx*dx + dy*dy;
		    return (trkhit->hit->xy-this->hit->xy).Mod2();
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

#endif // _DTrackCandidate_factory_FDCpseudo_


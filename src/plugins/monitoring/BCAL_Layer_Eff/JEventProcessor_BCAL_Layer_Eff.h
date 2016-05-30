// $Id$
//
// File: JEventProcessor_BCAL_Layer_Eff.h
//

#ifndef _JEventProcessor_BCAL_Layer_Eff_
#define _JEventProcessor_BCAL_Layer_Eff_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include "TH1I.h"
#include "TH2I.h"

#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALUnifiedHit.h"

#include "ANALYSIS/DCutActions.h"
#include "ANALYSIS/DTreeInterface.h"

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>

using namespace jana;
using namespace std;

class JEventProcessor_BCAL_Layer_Eff : public jana::JEventProcessor
{
	public:
		JEventProcessor_BCAL_Layer_Eff(){};
		~JEventProcessor_BCAL_Layer_Eff(){};
		const char* className(void){return "JEventProcessor_BCAL_Layer_Eff";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		pair<const DBCALPoint*, double> Find_NearestClusterPoint(double locProjectedSector, const set<const DBCALPoint*>& locClusterLayerBCALPoints);
		pair<const DBCALUnifiedHit*, double> Find_NearestClusterHit(double locProjectedSector, const set<const DBCALUnifiedHit*>& locBCALUnifiedHits);

		double Calc_AverageSector(const set<const DBCALPoint*>& locBCALPoints);
		double Calc_ProjectedSector(int locLayer, const map<int, set<const DBCALPoint*> >& locSortedPoints);

		//TRACK REQUIREMENTS
		double dMinTrackingFOM, dMinPIDFOM;
		int dMinNumTrackHits;
		int dMinHitRingsPerCDCSuperlayer, dMinHitPlanesPerFDCPackage;
		DCutAction_TrackHitPattern* locCutAction_TrackHitPattern;

		//HISTOGRAMS
		map<int, map<bool, TH1I*> > dHistMap_HitFound, dHistMap_HitTotal; //int = layer, bool = isUpstream
};

#endif // _JEventProcessor_BCAL_Layer_Eff_


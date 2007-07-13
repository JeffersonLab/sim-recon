// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceHDDM
//
/// Implements JEventSource for HDDM files

#ifndef _JEVENT_SOURCEHDDM_H_
#define _JEVENT_SOURCEHDDM_H_

#include <vector>
#include <string>
using namespace std;

#include <pthread.h>

#include <JANA/JEventSource.h>
#include <JANA/jerror.h>
#include "hddm_s.h"
#include "TRACKING/DMCTrackHit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DMCTrajectoryPoint.h"
#include "BCAL/DHDDMBCALHit.h"
#include "BCAL/DBCALTruthShower.h"
#include "CDC/DCDCHit.h"
#include "FDC/DFDCHit.h"
#include "FCAL/DFCALTruthShower.h"
#include "FCAL/DFCALHit.h"
#include "UPV/DUPVHit.h"
#include "UPV/DUPVTruthHit.h"
#include "TOF/DTOFTruth.h"
#include "TOF/DHDDMTOFHit.h"
#include "START_COUNTER/DSCTruthHit.h"
#include "START_COUNTER/DSCHit.h"

class DEventSourceHDDM:public JEventSource
{
	public:
		DEventSourceHDDM(const char* source_name);
		virtual ~DEventSourceHDDM();		
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventSourceHDDM";}

		jerror_t GetEvent(JEvent &event);
		void FreeEvent(JEvent &event);
		jerror_t GetObjects(JEvent &event, JFactory_base *factory);
		
		jerror_t Extract_DMCTrackHit(s_HDDM_t *hddm_s, JFactory<DMCTrackHit> *factory);
		jerror_t GetCDCTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetFDCTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetBCALTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetTOFTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetCherenkovTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetFCALTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetUPVTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);

		jerror_t Extract_DBCALTruthShower(s_HDDM_t *hddm_s, JFactory<DBCALTruthShower> *factory);
		jerror_t Extract_DHDDMBCALHit(s_HDDM_t *hddm_s, JFactory<DHDDMBCALHit> *factory);
		jerror_t Extract_DMCThrown(s_HDDM_t *hddm_s, JFactory<DMCThrown> *factory);
		jerror_t Extract_DCDCHit(s_HDDM_t *hddm_s, JFactory<DCDCHit> *factory);
		jerror_t Extract_DFDCHit(s_HDDM_t *hddm_s, JFactory<DFDCHit> *factory);
		jerror_t Extract_DFCALTruthShower(s_HDDM_t *hddm_s, JFactory<DFCALTruthShower> *factory);
		jerror_t Extract_DFCALHit(s_HDDM_t *hddm_s, JFactory<DFCALHit> *factory, JEventLoop* eventLoop );
		jerror_t Extract_DUPVHit(s_HDDM_t *hddm_s, JFactory<DUPVHit> *factory);
		jerror_t Extract_DUPVTruthHit(s_HDDM_t *hddm_s, JFactory<DUPVTruthHit> *factory);
		jerror_t Extract_DMCTrajectoryPoint(s_HDDM_t *hddm_s, JFactory<DMCTrajectoryPoint> *factory);
		jerror_t Extract_DTOFTruth(s_HDDM_t *hddm_s,  JFactory<DTOFTruth> *factory);
		jerror_t Extract_DHDDMTOFHit( s_HDDM_t *hddm_s,  JFactory<DHDDMTOFHit>* factory);
		jerror_t Extract_DSCHit(s_HDDM_t *hddm_s,  JFactory<DSCHit> *factory);
		jerror_t Extract_DSCTruthHit(s_HDDM_t *hddm_s,  JFactory<DSCTruthHit> *factory);
		s_iostream_t *fin;
		s_HDDM_t *hddm_s;
		bool flush_on_free;
		
	private:

};

#endif //_JEVENT_SOURCEHDDM_H_

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

#include "JANA/JEventSource.h"
#include "JANA/jerror.h"
#include "hddm_s.h"
#include "TRACKING/DMCTrackHit.h"
#include "TRACKING/DMCThrown.h"

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
		jerror_t GetCDCHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetFDCHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetBCALHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetTOFHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetCherenkovHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetFCALHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetUPVHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);

		jerror_t Extract_DMCThrown(s_HDDM_t *hddm_s, JFactory<DMCThrown> *factory);

		s_iostream_t *fin;
		s_HDDM_t *hddm_s;
		bool flush_on_free;
		
	private:

};

#endif //_JEVENT_SOURCEHDDM_H_

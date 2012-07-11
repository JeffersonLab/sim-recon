//
// Author: Richard Jones  June 29, 2012
//
// DEventSourceREST
//
/// Implements JEventSource for REST files

#ifndef _JEVENT_SOURCEREST_H_
#define _JEVENT_SOURCEREST_H_

#include <vector>
#include <string>

#include <pthread.h>

#include <JANA/JEventSource.h>
#include <JANA/jerror.h>
#include <JANA/JCalibration.h>

#include "hddm_r.hpp"

#include <PID/DMCReaction.h>
#include <PID/DBeamPhoton.h>
#include <TAGGER/DTagger.h>
#include "TRACKING/DMCThrown.h"
#include <FCAL/DFCALShower.h>
#include <BCAL/DBCALShower.h>
#include <PID/DNeutralShower.h>
#include <PID/DChargedTrackHypothesis.h>

class DEventSourceREST:public JEventSource
{
 public:
   DEventSourceREST(const char* source_name);
   virtual ~DEventSourceREST();		
   virtual const char* className(void) {
      return static_className();
   }
   static const char* static_className(void) {
      return "DEventSourceREST";
   }

   jerror_t GetEvent(JEvent &event);
   void FreeEvent(JEvent &event);
   jerror_t GetObjects(JEvent &event, JFactory_base *factory);
		
   jerror_t Extract_DMCReaction(hddm_r::HDDM *record,
                    JFactory<DMCReaction> *factory);
   jerror_t Extract_DBeamPhoton(hddm_r::HDDM *record,
                    JFactory<DBeamPhoton> *factory);
   jerror_t Extract_DMCThrown(hddm_r::HDDM *record,
                    JFactory<DMCThrown> *factory);
   jerror_t Extract_DTagger(hddm_r::HDDM *record,
                    JFactory<DTagger>* factory);
   jerror_t Extract_DFCALShower(hddm_r::HDDM *record,
                    JFactory<DFCALShower>* factory);
   jerror_t Extract_DBCALShower(hddm_r::HDDM *record,
                    JFactory<DBCALShower>* factory);
   jerror_t Extract_DNeutralShower(hddm_r::HDDM *record,
                    JFactory<DNeutralShower>* factory);
   jerror_t Extract_DChargedTrackHypothesis(hddm_r::HDDM *record,
                    JFactory<DChargedTrackHypothesis>* factory);

   Particle_t PDGtoPtype(int pdgtype);

 private:
   // Warning: Class JEventSource methods must be re-entrant, so do not
   // store any data here that might change from event to event.
   std::ifstream *ifs;
   hddm_r::istream *fin;
};

#endif //_JEVENT_SOURCEREST_H_

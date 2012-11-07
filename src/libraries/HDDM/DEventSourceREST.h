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
#include <TRACKING/DTrackTimeBased.h>
#include <FCAL/DFCALShower.h>
#include <BCAL/DBCALShower.h>
#include <START_COUNTER/DSCHit.h>
#include <TOF/DTOFPoint.h>
#include <TRIGGER/DMCTrigger.h>

#include <DMatrix.h>

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
   jerror_t Extract_DSCHit(hddm_r::HDDM *record,
                    JFactory<DSCHit>* factory);
   jerror_t Extract_DTOFPoint(hddm_r::HDDM *record,
                    JFactory<DTOFPoint>* factory);
   jerror_t Extract_DFCALShower(hddm_r::HDDM *record,
                    JFactory<DFCALShower>* factory);
   jerror_t Extract_DBCALShower(hddm_r::HDDM *record,
                    JFactory<DBCALShower>* factory);
   jerror_t Extract_DTrackTimeBased(hddm_r::HDDM *record,
                    JFactory<DTrackTimeBased>* factory);
   jerror_t Extract_DMCTrigger(hddm_r::HDDM *record,
                    JFactory<DMCTrigger>* factory);
#if 0
   jerror_t Extract_DRFTime(hddm_r::HDDM *record,
                    JFactory<DRFTime>* factory);
#endif

   Particle_t PDGtoPtype(int pdgtype);

   DMatrixDSym Get7x7ErrorMatrix(double mass, const double vec[5], const DMatrixDSym& C5x5);
 private:
   // Warning: Class JEventSource methods must be re-entrant, so do not
   // store any data here that might change from event to event.

   std::ifstream *ifs;		// input hddm file ifstream
   hddm_r::istream *fin;	// provides hddm layer on top of ifstream

   pthread_mutex_t rt_mutex;
   std::map<hddm_r::HDDM*, std::vector<DReferenceTrajectory*> > rt_by_event;
   std::list<DReferenceTrajectory*> rt_pool;

   DApplication *dapp;
   const DMagneticFieldMap *saved_bfield;
   const DGeometry *saved_geom;
   int saved_runnumber;
};

#endif //_JEVENT_SOURCEREST_H_

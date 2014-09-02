// $Id$
//
// Author: David Lawrence  June 24, 2004
//
// changes: Wed Jun 20 17:08:13 EDT 2007 B. Zihlmann
//          modify TOF section to add several new variables incuding the 
//          GEANT particle type to the Truth hits and the hit and track-hit
//          list.
//
// Oct 3, 2012 Yi Qiang: add functions for Cherenkov RICH detector
// Oct 11, 2012 Yi Qiang: complete functions for Cherenkov detector
// Oct 8, 2013 Yi Qiang: added dedicated object for RICH Truth Hit
// July 5, 2014 R.T.Jones: changed over from c to c++ API for hddm
//
//
// DEventSourceHDDM methods
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include <JANA/JFactory_base.h>
#include <JANA/JEventLoop.h>
#include <JANA/JEvent.h>

#include "BCAL/DBCALGeometry.h"

#include <DVector2.h>
#include <DEventSourceHDDM.h>
#include <FDC/DFDCGeometry.h>
#include <FCAL/DFCALGeometry.h>
#include <FCAL/DFCALHit.h>
#include <CCAL/DCCALGeometry.h>
#include <CCAL/DCCALHit.h>

//------------------------------------------------------------------
// Binary predicate used to sort hits
//------------------------------------------------------------------
class MCTrackHitSort{
   public:
      bool operator()(DMCTrackHit* const &thit1, 
                      DMCTrackHit* const &thit2) const {
         return thit1->z < thit2->z;
      }
};

bool MCTrackHitSort_C(DMCTrackHit* const &thit1, 
                      DMCTrackHit* const &thit2) {
   return thit1->z < thit2->z;
}


//----------------
// Constructor
//----------------
DEventSourceHDDM::DEventSourceHDDM(const char* source_name)
: JEventSource(source_name)
{
   /// Constructor for DEventSourceHDDM object
   ifs = new ifstream(source_name);
   if (ifs->is_open())
      fin = new hddm_s::istream(*ifs);
   else
      fin = NULL;
   initialized = false;
   dapp = NULL;
   bfield = NULL;
   geom = NULL;
   
   pthread_mutex_init(&rt_mutex, NULL);
}

//----------------
// Destructor
//----------------
DEventSourceHDDM::~DEventSourceHDDM()
{
   if (fin)
      delete fin;
   if (ifs)
      delete ifs;
}

//----------------
// GetEvent
//----------------
jerror_t DEventSourceHDDM::GetEvent(JEvent &event)
{
   /// Implementation of JEventSource virtual function

   if (!fin)
      return EVENT_SOURCE_NOT_OPEN;

   // Each open HDDM file takes up about 1M of memory so it's
   // worthwhile to close it as soon as we can.
   else if (!ifs->good()) {
      delete fin;
      fin = NULL;
      delete ifs;
      ifs = NULL;
      return NO_MORE_EVENTS_IN_SOURCE;
   }
   
   ++Nevents_read;
   
   hddm_s::HDDM *record = new hddm_s::HDDM();
   *fin >> *record;

   int event_number = -1;
   int run_number = -1;
   
   // Get event/run numbers from HDDM
   hddm_s::PhysicsEvent &pe = record->getPhysicsEvent(0);
   event_number = pe.getEventNo();     
   run_number = pe.getRunNo();

   // Copy the reference info into the JEvent object
   event.SetJEventSource(this);
   event.SetEventNumber(event_number);
   event.SetRunNumber(run_number);
   event.SetRef(record);
 
   return NOERROR;
}

//----------------
// FreeEvent
//----------------
void DEventSourceHDDM::FreeEvent(JEvent &event)
{
   hddm_s::HDDM *record = (hddm_s::HDDM*)event.GetRef();
   delete record;

   // Check for DReferenceTrajectory objects we need to delete
   pthread_mutex_lock(&rt_mutex);
   map<hddm_s::HDDM*, vector<DReferenceTrajectory*> >::iterator iter =
                                              rt_by_event.find(record);
   if(iter != rt_by_event.end()){
      vector<DReferenceTrajectory*> &rts = iter->second;
      for (unsigned int i=0; i<rts.size(); i++)
         rt_pool.push_back(rts[i]);
      rt_by_event.erase(iter);
   }
   pthread_mutex_unlock(&rt_mutex);
}

//----------------
// GetObjects
//----------------
jerror_t DEventSourceHDDM::GetObjects(JEvent &event, JFactory_base *factory)
{
   /// This gets called through the virtual method of the
   /// JEventSource base class. It creates the objects of the type
   /// on which factory is based. It uses the hddm_s::HDDM* object
   /// kept in the ref field of the JEvent object passed.

   // We must have a factory to hold the data
   if (!factory)
      throw RESOURCE_UNAVAILABLE;
   
   // HDDM doesn't exactly support tagged factories, but the tag
   // can be used to direct filling of the correct factory.
   string tag = (factory->Tag()==NULL)? "" : factory->Tag();
   
   // The ref field of the JEvent is just the HDDM object pointer.
   hddm_s::HDDM *record = (hddm_s::HDDM*)event.GetRef();
   if (!record)
      throw RESOURCE_UNAVAILABLE;

   // Get pointer to the B-field object and Geometry object
   JEventLoop *loop = event.GetJEventLoop();
   if (initialized == false && loop) {
      initialized = true;
      dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
      if (dapp) {
         jcalib = dapp->GetJCalibration(event.GetRunNumber());
         // Make sure jcalib is set
         if (!jcalib) {
            _DBG_ << "ERROR - no jcalib set!" <<endl;
            return RESOURCE_UNAVAILABLE;
         }
         // Get constants and do basic check on number of elements
         vector< map<string, float> > tvals;
         jcalib->Get("FDC/strip_calib", tvals);
 
         if (tvals.size() != 192) {
            _DBG_ << "ERROR - strip calibration vectors are not the right size!"
                  << endl;
            return VALUE_OUT_OF_RANGE;
         }
         // Notify user
         jout << "Read " << tvals.size()
              << " values from FDC/strip_calib in calibDB"
              << endl;
         jout << "   strip_calib columns (alphabetical): ";
         map<string,float>::iterator iter;
         for (iter=tvals[0].begin(); iter!=tvals[0].end(); iter++) {
            jout << iter->first << " ";
            jout << endl;

            // Copy values into tables. We preserve the order since
            // that is how it was originally done in hitFDC.c
            for (unsigned int i=0; i<tvals.size(); i++) {
               map<string, float> &row = tvals[i];
               uscale[i]=row["qru"];
               vscale[i]=row["qrv"];
            }
         }
      }
   }

   //Get target center
      //multiple reader threads can access this object: need lock
   bool locNewRunNumber = false;
   unsigned int locRunNumber = event.GetRunNumber();
   LockRead();
   {
      locNewRunNumber = (bTargetCenterZMap.find(locRunNumber) == bTargetCenterZMap.end());
   }
   UnlockRead();
   if(locNewRunNumber)
   {
      DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
      DGeometry* locGeometry = dapp->GetDGeometry(loop->GetJEvent().GetRunNumber());
      double locTargetCenterZ = 0.0;
      locGeometry->GetTargetZ(locTargetCenterZ);
      LockRead();
      {
         bTargetCenterZMap[locRunNumber] = locTargetCenterZ;
      }
      UnlockRead();
   }

   // Get name of data class we're trying to extract
   string dataClassName = factory->GetDataClassName();

   if (dataClassName == "DRFTime")
      return Extract_DRFTime(record, 
                     dynamic_cast<JFactory<DRFTime>*>(factory), loop);

   if (dataClassName == "DTAGMHit")
      return Extract_DTAGMHit(record, 
                     dynamic_cast<JFactory<DTAGMHit>*>(factory), tag);
 
   if (dataClassName == "DTAGHHit")
      return Extract_DTAGHHit(record, 
                     dynamic_cast<JFactory<DTAGHHit>*>(factory), tag);

   if (dataClassName == "DMCTrackHit")
      return Extract_DMCTrackHit(record,
                     dynamic_cast<JFactory<DMCTrackHit>*>(factory), tag);
 
   if (dataClassName == "DMCReaction")
      return Extract_DMCReaction(record,
                     dynamic_cast<JFactory<DMCReaction>*>(factory), tag, loop);

   if (dataClassName == "DBeamPhoton")
      return Extract_DBeamPhoton(record, 
                     dynamic_cast<JFactory<DBeamPhoton>*>(factory), tag, loop);
 
   if (dataClassName == "DMCThrown")
      return Extract_DMCThrown(record,
                     dynamic_cast<JFactory<DMCThrown>*>(factory), tag);
 
   if (dataClassName == "DBCALTruthShower")
      return Extract_DBCALTruthShower(record, 
                     dynamic_cast<JFactory<DBCALTruthShower>*>(factory), tag);
 
   if (dataClassName == "DBCALSiPMSpectrum")
      return Extract_DBCALSiPMSpectrum(record,
                     dynamic_cast<JFactory<DBCALSiPMSpectrum>*>(factory), tag);
 
   if (dataClassName == "DBCALTruthCell")
      return Extract_DBCALTruthCell(record,
                     dynamic_cast<JFactory<DBCALTruthCell>*>(factory), tag);
 
   if (dataClassName == "DBCALSiPMHit")
      return Extract_DBCALSiPMHit(record,
                     dynamic_cast<JFactory<DBCALSiPMHit>*>(factory), tag);
 
   if (dataClassName == "DBCALHit")
      return Extract_DBCALHit(record,
                     dynamic_cast<JFactory<DBCALHit>*>(factory), tag);

   if (dataClassName == "DBCALIncidentParticle")
      return Extract_DBCALIncidentParticle(record,
                     dynamic_cast<JFactory<DBCALIncidentParticle>*>(factory), tag);
 
   if (dataClassName == "DBCALTDCHit")
      return Extract_DBCALTDCHit(record,
                     dynamic_cast<JFactory<DBCALTDCHit>*>(factory), tag);
 
   if (dataClassName == "DCDCHit")
      return Extract_DCDCHit(record,
                     dynamic_cast<JFactory<DCDCHit>*>(factory) , tag);
 
   if (dataClassName == "DFDCHit")
      return Extract_DFDCHit(record, 
                     dynamic_cast<JFactory<DFDCHit>*>(factory), tag);
 
   if (dataClassName == "DFCALTruthShower")
      return Extract_DFCALTruthShower(record, 
                     dynamic_cast<JFactory<DFCALTruthShower>*>(factory), tag);
 
   if (dataClassName == "DFCALHit")
      return Extract_DFCALHit(record,
                     dynamic_cast<JFactory<DFCALHit>*>(factory), tag,
                     event.GetJEventLoop());
 
   if (dataClassName == "DCCALTruthShower")
      return Extract_DCCALTruthShower(record,
                     dynamic_cast<JFactory<DCCALTruthShower>*>(factory), tag);
 
   if (dataClassName == "DCCALHit")
      return Extract_DCCALHit(record,
                     dynamic_cast<JFactory<DCCALHit>*>(factory), tag,
                     event.GetJEventLoop());
 
   if (dataClassName == "DMCTrajectoryPoint" && tag == "")
      return Extract_DMCTrajectoryPoint(record,
                     dynamic_cast<JFactory<DMCTrajectoryPoint>*>(factory), tag);
 
   if (dataClassName == "DTOFTruth")
      return Extract_DTOFTruth(record, 
                     dynamic_cast<JFactory<DTOFTruth>*>(factory), tag);
 
   // TOF is a special case: TWO factories are needed at the same time
   // DTOFHit and DTOFHitMC
   if (dataClassName == "DTOFHit") {
      JFactory_base* factory2 = loop->GetFactory("DTOFHitMC", tag.c_str()); 
      return Extract_DTOFHit(record, 
                     dynamic_cast<JFactory<DTOFHit>*>(factory),
                     dynamic_cast<JFactory<DTOFHitMC>*>(factory2), tag);
   }
   if (dataClassName == "DTOFHitMC") {
      JFactory_base* factory2 = loop->GetFactory("DTOFHit", tag.c_str()); 
      return Extract_DTOFHit(record, 
                     dynamic_cast<JFactory<DTOFHit>*>(factory2),
                     dynamic_cast<JFactory<DTOFHitMC>*>(factory), tag);
   }

   if (dataClassName == "DSCHit")
      return Extract_DSCHit(record, 
                     dynamic_cast<JFactory<DSCHit>*>(factory), tag);

   if (dataClassName == "DSCTruthHit")
      return Extract_DSCTruthHit(record, 
                     dynamic_cast<JFactory<DSCTruthHit>*>(factory), tag);

   if (dataClassName == "DTrackTimeBased")
      return Extract_DTrackTimeBased(record,
                     dynamic_cast<JFactory<DTrackTimeBased>*>(factory), tag,
                     event.GetRunNumber());

   // extract CereTruth and CereRichHit hits, yqiang Oct 3, 2012
   // removed CereTruth (merged into MCThrown), added CereHit, yqiang Oct 10 2012
   if (dataClassName == "DCereHit")
      return Extract_DCereHit(record, 
                     dynamic_cast<JFactory<DCereHit>*>(factory), tag);
   if (dataClassName == "DRichHit")
      return Extract_DRichHit(record, 
                     dynamic_cast<JFactory<DRichHit>*>(factory), tag);

   // added DRichTruthHit
   if (dataClassName == "DRichTruthHit")
      return Extract_DRichTruthHit(record,
                     dynamic_cast<JFactory<DRichTruthHit>*>(factory), tag);

   return OBJECT_NOT_AVAILABLE;
}

//------------------
// Extract_DRFTime
//------------------
jerror_t DEventSourceHDDM::Extract_DRFTime(hddm_s::HDDM *record,
                                   JFactory<DRFTime> *factory, JEventLoop* locEventLoop)
{
   if (factory==NULL)
      return OBJECT_NOT_AVAILABLE;
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DRFTime*> locRFTimes;

   // loop over RF-time records
   const hddm_s::RFtimeList &rftimes = record->getRFtimes();
   hddm_s::RFtimeList::iterator iter;
   for (iter = rftimes.begin(); iter != rftimes.end(); ++iter)
   {
      if (iter->getJtag() != tag)
         continue;
      DRFTime *locRFTime = new DRFTime;
      locRFTime->dTime = iter->getTsync();
      locRFTime->dTimeVariance = 0.0015; //1.5ps
      locRFTimes.push_back(locRFTime);
   }

   if(locRFTimes.empty())
   {
      //See if MC data. If so, generate the DRFTime object here (not in input file)
      // https://halldweb1.jlab.org/wiki/index.php/How_HDGeant_defines_time-zero_for_physics_events
      vector<const DMCThrown*> locMCThrowns;
      locEventLoop->Get(locMCThrowns);
      if(!locMCThrowns.empty())
      {
         DRFTime *locRFTime = new DRFTime;
         locRFTime->dTime = 0.0;
         locRFTime->dTimeVariance = 0.0;
         locRFTimes.push_back(locRFTime);
      }
   }

   // Copy into factories
   factory->CopyTo(locRFTimes);

   return NOERROR;
}

//------------------
// Extract_DMCTrackHit
//------------------
jerror_t DEventSourceHDDM::Extract_DMCTrackHit(hddm_s::HDDM *record,
                                   JFactory<DMCTrackHit> *factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
   
   // The following routines will create DMCTrackHit objects and add them
   // to data.
   vector<DMCTrackHit*> data;
   GetCDCTruthHits(record, data);
   GetFDCTruthHits(record, data);
   GetBCALTruthHits(record, data);
   GetTOFTruthHits(record, data);
   GetCherenkovTruthHits(record, data);
   GetFCALTruthHits(record, data);
   GetSCTruthHits(record, data);
   // added RICH truth hits, yqiang, Oct 10 2012
   GetRichTruthHits(record, data);

   // It has happened that some CDC hits have "nan" for the drift time
   // in a peculiar event Alex Somov came across. This ultimately caused
   // a seg. fault in MCTrackHitSort_C. I hate doing this since it
   // is treating the symptom rather than the cause, but nonetheless,
   // it patches up the problem for now until there is time to revisit
   // it later.
   for (unsigned int i=0; i < data.size(); i++)
      if (!isfinite(data[i]->z))
         data[i]->z = -1000.0;
   
   // sort hits by z
   sort(data.begin(), data.end(), MCTrackHitSort_C);
   
   // Some systems will use negative phis. Force them all to
   // be in the 0 to 2pi range
   for (unsigned int i=0; i < data.size(); i++) {
      DMCTrackHit *mctrackhit = data[i];
      if (mctrackhit->phi < 0.0)
         mctrackhit->phi += 2.0*M_PI;
   }
   
   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//-------------------
// GetCDCTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetCDCTruthHits(hddm_s::HDDM *record, 
                                           vector<DMCTrackHit*>& data)
{
   const hddm_s::CdcTruthPointList &points = record->getCdcTruthPoints();
   hddm_s::CdcTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = iter->getR();
      mctrackhit->phi     = iter->getPhi();
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype();
      mctrackhit->system  = SYS_CDC;
      data.push_back(mctrackhit);
   }
      
   return NOERROR;
}

//-------------------
// GetFDCTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetFDCTruthHits(hddm_s::HDDM *record,
                                           vector<DMCTrackHit*>& data)
{
   const hddm_s::FdcTruthPointList &points = record->getFdcTruthPoints();
   hddm_s::FdcTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      float x = iter->getX();
      float y = iter->getY();
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = sqrt(x*x + y*y);
      mctrackhit->phi     = atan2(y,x);
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype();
      mctrackhit->system  = SYS_FDC;
      data.push_back(mctrackhit);
   }

   return NOERROR;
}

//-------------------
// GetBCALTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetBCALTruthHits(hddm_s::HDDM *record,
                                            vector<DMCTrackHit*>& data)
{
   const hddm_s::BcalTruthShowerList &showers = record->getBcalTruthShowers();
   hddm_s::BcalTruthShowerList::iterator iter;
   for (iter = showers.begin(); iter != showers.end(); ++iter) {
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = iter->getR();
      mctrackhit->phi     = iter->getPhi();
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype();
      mctrackhit->system  = SYS_BCAL;
      data.push_back(mctrackhit);
   }

   return NOERROR;
}

//-------------------
// GetTOFTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetTOFTruthHits(hddm_s::HDDM *record,
                                           vector<DMCTrackHit*>& data)
{
   const hddm_s::FtofTruthPointList &points = record->getFtofTruthPoints();
   hddm_s::FtofTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      float x = iter->getX();
      float y = iter->getY();
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = sqrt(x*x + y*y);
      mctrackhit->phi     = atan2(y,x);
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype(); // save GEANT particle type 
      mctrackhit->system  = SYS_TOF;
      data.push_back(mctrackhit);
   }
      
   return NOERROR;
}

//-------------------
// GetCherenkovTruthHits
// modified by yqiang, Oct 10 2012
//-------------------
jerror_t DEventSourceHDDM::GetCherenkovTruthHits(hddm_s::HDDM *record,
                                                 vector<DMCTrackHit*>& data)
{
   const hddm_s::CereTruthPointList &points = record->getCereTruthPoints();
   hddm_s::CereTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      float x = iter->getX();
      float y = iter->getY();
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = sqrt(x*x + y*y);
      mctrackhit->phi     = atan2(y,x);
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype();    // save GEANT particle typ()e
      mctrackhit->system  = SYS_CHERENKOV;
      data.push_back(mctrackhit);
   }

   return NOERROR;
}

//-------------------
// GetCherenkovRichTruthHits
// added by yqiang, Oct 11 2012
//-------------------
jerror_t DEventSourceHDDM::GetRichTruthHits(hddm_s::HDDM *record,
                                            vector<DMCTrackHit*>& data)
{
   const hddm_s::RichTruthPointList &points = record->getRichTruthPoints();
   hddm_s::RichTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      float x = iter->getX();
      float y = iter->getY();
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = sqrt(x*x + y*y);
      mctrackhit->phi     = atan2(y,x);
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype();    // save GEANT particle typ()e
      mctrackhit->system  = SYS_RICH;
      data.push_back(mctrackhit);
   }

   return NOERROR;
}

//-------------------
// GetFCALTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetFCALTruthHits(hddm_s::HDDM *record,
                                            vector<DMCTrackHit*>& data)
{
   const hddm_s::FcalTruthShowerList &points = record->getFcalTruthShowers();
   hddm_s::FcalTruthShowerList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      float x = iter->getX();
      float y = iter->getY();
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = sqrt(x*x + y*y);
      mctrackhit->phi     = atan2(y,x);
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype();
      mctrackhit->system  = SYS_FCAL;
      data.push_back(mctrackhit);
   }

   return NOERROR;
}

//-------------------
// GetCCALTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetCCALTruthHits(hddm_s::HDDM *record,
                                            vector<DMCTrackHit*>& data)
{
   const hddm_s::CcalTruthShowerList &points = record->getCcalTruthShowers();
   hddm_s::CcalTruthShowerList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      float x = iter->getX();
      float y = iter->getY();
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = sqrt(x*x + y*y);
      mctrackhit->phi     = atan2(y,x);
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype();
      mctrackhit->system  = SYS_CCAL;
      data.push_back(mctrackhit);
   }

   return NOERROR;
}


//-------------------
// GetSCTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetSCTruthHits(hddm_s::HDDM *record,
                                          vector<DMCTrackHit*>& data)
{
   const hddm_s::StcTruthPointList &points = record->getStcTruthPoints();
   hddm_s::StcTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      DMCTrackHit *mctrackhit = new DMCTrackHit;
      mctrackhit->r       = iter->getR();
      mctrackhit->phi     = iter->getPhi();
      mctrackhit->z       = iter->getZ();
      mctrackhit->track   = iter->getTrack();
      mctrackhit->primary = iter->getPrimary();
      mctrackhit->ptype   = iter->getPtype();    // save GEANT particle type
      mctrackhit->system  = SYS_START;
      data.push_back(mctrackhit);
   }
      
   return NOERROR;
}

//------------------
// Extract_DBCALSiPMHit
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALSiPMHit(hddm_s::HDDM *record,
                                   JFactory<DBCALSiPMHit> *factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DBCALSiPMHit*> data;

   const hddm_s::BcalSiPMUpHitList &uphits = record->getBcalSiPMUpHits();
   hddm_s::BcalSiPMUpHitList::iterator uiter;
   for (uiter = uphits.begin(); uiter != uphits.end(); ++uiter) {
      DBCALSiPMHit *response = new DBCALSiPMHit;
      response->module = uiter->getModule();
      response->layer  = uiter->getLayer();
      response->sector = uiter->getSector();
      response->E      = uiter->getE();
      response->t      = uiter->getT();
      response->end    = DBCALGeometry::kUpstream;
      response->cellId = DBCALGeometry::cellId(uiter->getModule(),
                                               uiter->getLayer(), 
                                               uiter->getSector());
      data.push_back(response);
   }
         
   const hddm_s::BcalSiPMDownHitList &downhits = record->getBcalSiPMDownHits();
   hddm_s::BcalSiPMDownHitList::iterator diter;
   for (diter = downhits.begin(); diter != downhits.end(); ++diter) {
      DBCALSiPMHit *response = new DBCALSiPMHit;
      response->module = diter->getModule();
      response->layer  = diter->getLayer();
      response->sector = diter->getSector();
      response->E      = diter->getE();
      response->t      = diter->getT();
      response->end    = DBCALGeometry::kDownstream;
      response->cellId = DBCALGeometry::cellId(diter->getModule(),
                                               diter->getLayer(), 
                                               diter->getSector());
      data.push_back(response);
   }
   
   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DBCALHit
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALHit(hddm_s::HDDM *record,
                                   JFactory<DBCALHit> *factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DBCALHit*> data;

   const hddm_s::BcalfADCHitList &hits = record->getBcalfADCHits();
   hddm_s::BcalfADCHitList::iterator iter;
   for (iter = hits.begin(); iter != hits.end(); ++iter) {
      DBCALHit *response = new DBCALHit;
      response->module = iter->getModule();
      response->layer  = iter->getLayer();
      response->sector = iter->getSector();
      response->E      = iter->getE();
      response->t      = iter->getT();
      response->end    = (iter->getEnd() == 0)? DBCALGeometry::kUpstream :
                                                DBCALGeometry::kDownstream;
      response->cellId = DBCALGeometry::cellId(iter->getModule(),
                                               iter->getLayer(),
                                               iter->getSector());
      data.push_back(response);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DBCALIncidentParticle
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALIncidentParticle(hddm_s::HDDM *record,
                                   JFactory<DBCALIncidentParticle> *factory,
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
  
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
  
   vector<DBCALIncidentParticle*> data;

   const hddm_s::BcalTruthIncidentParticleList &plist =
                          record->getBcalTruthIncidentParticles();
   hddm_s::BcalTruthIncidentParticleList::iterator iter;
   for (iter = plist.begin(); iter != plist.end(); ++iter) {
      DBCALIncidentParticle *part = new DBCALIncidentParticle;
      part->ptype = iter->getPtype();
      part->px    = iter->getPx();
      part->py    = iter->getPy();
      part->pz    = iter->getPz();
      part->x     = iter->getX();
      part->y     = iter->getY();
      part->z     = iter->getZ();
      data.push_back(part);
   }
  
   factory->CopyTo(data);
  
   return NOERROR;
}

//------------------
// Extract_DBCALSiPMSpectrum
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALSiPMSpectrum(hddm_s::HDDM *record,
                                   JFactory<DBCALSiPMSpectrum> *factory,
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DBCALSiPMSpectrum*> data;

   const hddm_s::BcalSiPMSpectrumList &specs = record->getBcalSiPMSpectrums();
   hddm_s::BcalSiPMSpectrumList::iterator iter;
   for (iter = specs.begin(); iter != specs.end(); ++iter) {
      DBCALSiPMSpectrum *dana_spectrum = new DBCALSiPMSpectrum;
      dana_spectrum->module      = iter->getModule();
      dana_spectrum->layer       = iter->getLayer();
      dana_spectrum->sector      = iter->getSector();
      dana_spectrum->end         = (iter->getEnd()==0)?
                                   DBCALGeometry::kUpstream :
                                   DBCALGeometry::kDownstream;
      if (tag == "") 
         dana_spectrum->incident_id = 0;
      else if (tag == "TRUTH") {
         const hddm_s::BcalSiPMTruthList &truths = iter->getBcalSiPMTruths();
         if (truths.size() > 0)
            dana_spectrum->incident_id = truths.begin()->getIncident_id();
         else
            dana_spectrum->incident_id = 0;
      }

      double t = iter->getTstart();
      double bin_width = iter->getBin_width();
      stringstream ss(iter->getVals());

      // Extract values and use them to fill histo
      double dE;
      while (ss >> dE) {
         dana_spectrum->spectrum.Fill(t, dE);
         t += bin_width;
      }

      data.push_back(dana_spectrum);
   }
         
   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DBCALTDCHit
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALTDCHit(hddm_s::HDDM *record,
                                   JFactory<DBCALTDCHit> *factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DBCALTDCHit*> data;

   const hddm_s::BcalTDCHitList &hits = record->getBcalTDCHits();
   hddm_s::BcalTDCHitList::iterator iter;
   for (iter = hits.begin(); iter != hits.end(); ++iter) {
      DBCALTDCHit *bcaltdchit = new DBCALTDCHit;
      bcaltdchit->module = iter->getModule();
      bcaltdchit->layer  = iter->getLayer();
      bcaltdchit->sector = iter->getSector();
      bcaltdchit->end    = (iter->getEnd() == 0)? DBCALGeometry::kUpstream :
                                                DBCALGeometry::kDownstream;
      bcaltdchit->cellId = DBCALGeometry::cellId(iter->getModule(),
                                                 iter->getLayer(),
                                                 iter->getSector());
      bcaltdchit->t      = iter->getT();
      data.push_back(bcaltdchit);
   }
         
   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DMCReaction
//------------------
jerror_t DEventSourceHDDM::Extract_DMCReaction(hddm_s::HDDM *record,
                                   JFactory<DMCReaction> *factory, string tag,
                                   JEventLoop *loop)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   
   double locTargetCenterZ = 0.0;
   int locRunNumber = loop->GetJEvent().GetRunNumber();
   LockRead();
   {
      locTargetCenterZ = bTargetCenterZMap[locRunNumber];
   }
   UnlockRead();
   DVector3 locPosition(0.0, 0.0, locTargetCenterZ);

   vector<DMCReaction*> dmcreactions;

   const hddm_s::ReactionList &reacts = record->getReactions();
   hddm_s::ReactionList::iterator iter;
   for (iter = reacts.begin(); iter != reacts.end(); ++iter) {
      DMCReaction *mcreaction = new DMCReaction;
      dmcreactions.push_back(mcreaction);
      mcreaction->type = iter->getType();
      mcreaction->weight = iter->getWeight();
      hddm_s::Origin &origin = iter->getVertex().getOrigin();
      double torig = origin.getT();
      double zorig = origin.getVz();
      
      const hddm_s::BeamList &beams = record->getBeams();
      if (beams.size() > 0) {
         hddm_s::Beam &beam = iter->getBeam();
         DVector3 mom(beam.getMomentum().getPx(),
                      beam.getMomentum().getPy(),
                      beam.getMomentum().getPz());
         mcreaction->beam.setPosition(locPosition);
         mcreaction->beam.setMomentum(mom);
         mcreaction->beam.setPID(Gamma);
         mcreaction->beam.setMass(beam.getProperties().getMass());
         mcreaction->beam.setCharge(beam.getProperties().getCharge());
         mcreaction->target.setPID(IDTrack(mcreaction->beam.charge(),
                                           mcreaction->beam.mass()));
         mcreaction->beam.clearErrorMatrix();
         mcreaction->beam.setTime(torig - (zorig - locTargetCenterZ)/30.0);
      }
      else {
         // fake values for DMCReaction
         DVector3 mom(0.0, 0.0, 0.0);
         mcreaction->beam.setPosition(locPosition);
         mcreaction->beam.setMomentum(mom);
         mcreaction->beam.setPID(Gamma);
         mcreaction->beam.setMass(0.0);
         mcreaction->beam.setCharge(0.0);
         mcreaction->beam.clearErrorMatrix();
         mcreaction->beam.setTime(0.0);
      }

      const hddm_s::TargetList &targets = record->getTargets();
      if (targets.size() > 0) {
         hddm_s::Target &target = iter->getTarget();
         DKinematicData target_kd;
         DVector3 mom(target.getMomentum().getPx(),
                      target.getMomentum().getPy(),
                      target.getMomentum().getPz());
         mcreaction->target.setPosition(locPosition);
         mcreaction->target.setMomentum(mom);
         mcreaction->target.setMass(target.getProperties().getMass());
         mcreaction->target.setCharge(target.getProperties().getCharge());
         mcreaction->target.setPID(IDTrack(mcreaction->target.charge(),
                                           mcreaction->target.mass()));
         mcreaction->target.clearErrorMatrix();
         mcreaction->target.setTime(0.0);
      }
      else {
         // fake values for DMCReaction
         DVector3 mom(0.0, 0.0, 0.0);
         mcreaction->target.setPosition(locPosition);
         mcreaction->target.setMomentum(mom);
         mcreaction->target.setPID(Unknown);
         mcreaction->target.setMass(0.0);
         mcreaction->target.setCharge(0.0);
         mcreaction->target.clearErrorMatrix();
         mcreaction->target.setTime(0.0);
      }
   }
   
   // Copy into factories
   //_DBG_<<"Creating "<<dmcreactions.size()<<" DMCReaction objects"<<endl;

   factory->CopyTo(dmcreactions);

   return NOERROR;
}


//------------------
// Extract_DBeamPhoton
//------------------
jerror_t DEventSourceHDDM::Extract_DBeamPhoton(hddm_s::HDDM *record,
                                   JFactory<DBeamPhoton> *factory, string tag,
                                   JEventLoop *loop)
{
   /// If tag="MCGEN" then defer to the Extract_DMCReaction method which
   /// extracts both the DMCReaction and DBeamPhoton objects at the same time.

   if (factory==NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "MCGEN")
      return OBJECT_NOT_AVAILABLE;

   vector<const DMCReaction*> dmcreactions;
   loop->Get(dmcreactions);

   vector<DBeamPhoton*> dbeam_photons;
   for(size_t loc_i = 0; loc_i < dmcreactions.size(); ++loc_i)
   {
      DBeamPhoton *beamphoton = new DBeamPhoton;
      *(DKinematicData*)beamphoton = dmcreactions[loc_i]->beam;
      dbeam_photons.push_back(beamphoton);
   }

   // Copy into factories
   factory->CopyTo(dbeam_photons);

   return NOERROR;
}

//------------------
// Extract_DMCThrown
//------------------
jerror_t DEventSourceHDDM::Extract_DMCThrown(hddm_s::HDDM *record,
                                   JFactory<DMCThrown> *factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DMCThrown*> data;

   const hddm_s::VertexList &verts = record->getVertices();
   hddm_s::VertexList::iterator iter;
   for (iter = verts.begin(); iter != verts.end(); ++iter) {
      const hddm_s::OriginList &origs = iter->getOrigins();
      const hddm_s::ProductList &prods = iter->getProducts();
      double vertex[4] = {0., 0., 0., 0.};
      if (origs.size() > 0) {
         vertex[0] = iter->getOrigin().getT();
         vertex[1] = iter->getOrigin().getVx();
         vertex[2] = iter->getOrigin().getVy();
         vertex[3] = iter->getOrigin().getVz();
      }
      hddm_s::ProductList::iterator piter;
      for (piter = prods.begin(); piter != prods.end(); ++piter) {
         double E  = piter->getMomentum().getE();
         double px = piter->getMomentum().getPx();
         double py = piter->getMomentum().getPy();
         double pz = piter->getMomentum().getPz();
         double mass = sqrt(E*E - (px*px + py*py + pz*pz));
         if (!isfinite(mass))
            mass = 0.0;
         DMCThrown *mcthrown = new DMCThrown;
         mcthrown->type     = piter->getType();
         mcthrown->myid     = piter->getId();
         mcthrown->parentid = piter->getParentid();
         mcthrown->mech     = piter->getMech();
         mcthrown->pdgtype  = piter->getPdgtype();
         mcthrown->setPID((Particle_t)mcthrown->type);
         mcthrown->setMass(mass);
         mcthrown->setMomentum(DVector3(px, py, pz));
         mcthrown->setPosition(DVector3(vertex[1], vertex[2], vertex[3]));
         mcthrown->setTime(vertex[0]);
         mcthrown->setCharge(ParticleCharge(piter->getType()));
         data.push_back(mcthrown);
      }
   }
   
   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DCDCHit
//------------------
jerror_t DEventSourceHDDM::Extract_DCDCHit(hddm_s::HDDM *record,
                                   JFactory<DCDCHit> *factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
       return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DCDCHit*> data;

   if (tag == "") {
      const hddm_s::CdcStrawHitList &hits = record->getCdcStrawHits();
      hddm_s::CdcStrawHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         DCDCHit *hit = new DCDCHit;
         hit->ring   = iter->getRing();
         hit->straw  = iter->getStraw();
         hit->q      = iter->getQ();
         hit->t      = iter->getT();
         hit->d      = 0.; // initialize to zero to avoid any NaN
         hit->itrack = 0;  // track information is in TRUTH tag
         hit->ptype  = 0;  // ditto
         data.push_back(hit);
      }
   }
   else if (tag == "TRUTH") {
      const hddm_s::CdcStrawTruthHitList &thits = record->getCdcStrawTruthHits();
      hddm_s::CdcStrawTruthHitList::iterator iter;
      for (iter = thits.begin(); iter != thits.end(); ++iter) {
         DCDCHit *hit = new DCDCHit;
         hit->ring   = iter->getRing();
         hit->straw  = iter->getStraw();
         hit->q      = iter->getQ();
         hit->t      = iter->getT();
         hit->d      = iter->getD();
         hit->itrack = iter->getItrack();
         hit->ptype  = iter->getPtype();
         data.push_back(hit);
      }
   }
   
   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}


//------------------
// Extract_DFDCHit
//------------------
jerror_t DEventSourceHDDM::Extract_DFDCHit(hddm_s::HDDM *record,
                                   JFactory<DFDCHit> *factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH" && tag != "CALIB")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DFDCHit*> data;

   if (tag == "") {
      const hddm_s::FdcAnodeHitList &ahits = record->getFdcAnodeHits();
      hddm_s::FdcAnodeHitList::iterator ahiter;
      for (ahiter = ahits.begin(); ahiter != ahits.end(); ++ahiter) {
         DFDCHit* newHit = new DFDCHit();
         newHit->layer   = ahiter->getLayer();
         newHit->module  = ahiter->getModule();
         newHit->element = ahiter->getWire();
         newHit->q       = ahiter->getDE();
         newHit->t       = ahiter->getT();
         newHit->d       = 0.; // initialize to zero to avoid any NaN
         newHit->itrack  = 0;  // track information is in TRUTH tag
         newHit->ptype   = 0;  // ditto
         newHit->plane   = 2;
         newHit->type    = 0;
         newHit->gPlane  = DFDCGeometry::gPlane(newHit);
         newHit->gLayer  = DFDCGeometry::gLayer(newHit);
         newHit->r       = DFDCGeometry::getWireR(newHit);
         data.push_back(newHit);
      }

      // Ditto for the cathodes.
      const hddm_s::FdcCathodeHitList &chits = record->getFdcCathodeHits();
      hddm_s::FdcCathodeHitList::iterator chiter;
      for (chiter = chits.begin(); chiter != chits.end(); ++chiter) {
         DFDCHit* newHit = new DFDCHit();
         newHit->layer   = chiter->getLayer();
         newHit->module  = chiter->getModule();
         newHit->element = chiter->getStrip();
         if (newHit->element > 1000)
            newHit->element -= 1000;
         newHit->plane   = chiter->getPlane();
         newHit->q       = chiter->getQ();
         newHit->t       = chiter->getT();
         newHit->d       = 0.; // initialize to zero to avoid any NaN
         newHit->itrack  = 0;  // track information is in TRUTH tag
         newHit->ptype   = 0;  // ditto
         newHit->type    = 1;
         newHit->gPlane  = DFDCGeometry::gPlane(newHit);    
         newHit->gLayer  = DFDCGeometry::gLayer(newHit);
         newHit->r       = DFDCGeometry::getStripR(newHit);
         data.push_back(newHit);
      }
   }

   else if (tag == "TRUTH"){
      const hddm_s::FdcAnodeTruthHitList &aths = record->getFdcAnodeTruthHits();
      hddm_s::FdcAnodeTruthHitList::iterator atiter;
      for (atiter = aths.begin(); atiter != aths.end(); ++atiter) {
         DFDCHit* newHit = new DFDCHit();
         newHit->layer   = atiter->getLayer();
         newHit->module  = atiter->getModule();
         newHit->element = atiter->getWire();
         newHit->q       = atiter->getDE();
         newHit->t       = atiter->getT();
         newHit->d       = atiter->getD();
         newHit->itrack  = atiter->getItrack();
         newHit->ptype   = atiter->getPtype();
         newHit->plane   = 2;
         newHit->type    = 0;
         newHit->gPlane  = DFDCGeometry::gPlane(newHit);
         newHit->gLayer  = DFDCGeometry::gLayer(newHit);
         newHit->r       = DFDCGeometry::getWireR(newHit);
         data.push_back(newHit);
      }
          
      // Ditto for the cathodes.
      const hddm_s::FdcCathodeTruthHitList &cths = 
                                           record->getFdcCathodeTruthHits();
      hddm_s::FdcCathodeTruthHitList::iterator ctiter;
      for (ctiter = cths.begin(); ctiter != cths.end(); ++ctiter) {
         DFDCHit* newHit = new DFDCHit();
         newHit->layer   = ctiter->getLayer();
         newHit->module  = ctiter->getModule();
         newHit->element = ctiter->getStrip();
         if (newHit->element > 1000)
            newHit->element -= 1000;
         newHit->plane   = ctiter->getPlane();
         newHit->q       = ctiter->getQ();
         newHit->t       = ctiter->getT();
         newHit->d       = 0.; // initialize to zero to avoid any NaN
         newHit->itrack  = ctiter->getItrack();
         newHit->ptype   = ctiter->getPtype();
         newHit->type    = 1;
         newHit->gPlane  = DFDCGeometry::gPlane(newHit);    
         newHit->gLayer  = DFDCGeometry::gLayer(newHit);
         newHit->r       = DFDCGeometry::getStripR(newHit);
         data.push_back(newHit);
      }
   }

   else if (tag == "CALIB") {
      // Deal with the wires
      const hddm_s::FdcAnodeHitList &ahits = record->getFdcAnodeHits();
      hddm_s::FdcAnodeHitList::iterator ahiter;
      for (ahiter = ahits.begin(); ahiter != ahits.end(); ++ahiter) {
         DFDCHit* newHit = new DFDCHit();
         newHit->layer   = ahiter->getLayer();
         newHit->module  = ahiter->getModule();
         newHit->element = ahiter->getWire();
         newHit->q       = ahiter->getDE();
         newHit->t       = ahiter->getT();
         newHit->d       = 0.; // initialize to zero to avoid any NaN
         newHit->itrack  = 0;  // track information is in TRUTH tag
         newHit->ptype   = 0;  // ditto
         newHit->plane   = 2;
         newHit->type    = 0;
         newHit->gPlane  = DFDCGeometry::gPlane(newHit);
         newHit->gLayer  = DFDCGeometry::gLayer(newHit);
         newHit->r       = DFDCGeometry::getWireR(newHit);
         data.push_back(newHit);
      }
          
      // Ditto for the cathodes.
      const hddm_s::FdcCathodeHitList &chits = record->getFdcCathodeHits();
      hddm_s::FdcCathodeHitList::iterator chiter;
      for (chiter = chits.begin(); chiter != chits.end(); ++chiter) {
         DFDCHit* newHit = new DFDCHit();
         newHit->layer   = chiter->getLayer();
         newHit->module  = chiter->getModule();
         newHit->element = chiter->getStrip();
         if (newHit->element > 1000)
            newHit->element-=1000;
         newHit->plane   = chiter->getPlane();
         if (newHit->plane == 1) // v 
            newHit->q = chiter->getQ()*vscale[newHit->element-1];
         else  // u
            newHit->q = chiter->getQ()*uscale[newHit->element-1];
         newHit->t       = chiter->getT();
         newHit->d       = 0.; // initialize to zero to avoid any NaN
         newHit->itrack  = 0;  // track information is in TRUTH tag
         newHit->ptype   = 0;  // ditto
         newHit->type    = 1;
         newHit->gPlane  = DFDCGeometry::gPlane(newHit);    
         newHit->gLayer  = DFDCGeometry::gLayer(newHit);
         newHit->r       = DFDCGeometry::getStripR(newHit);
         data.push_back(newHit);
      }
   }
   
   // Copy into factory
   factory->CopyTo(data);
   
   return NOERROR;
}

//------------------
// Extract_DBCALTruthShower
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALTruthShower(hddm_s::HDDM *record,
                                   JFactory<DBCALTruthShower> *factory,
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;

   vector<DBCALTruthShower*> data;

   const hddm_s::BcalTruthShowerList &shows = record->getBcalTruthShowers();
   hddm_s::BcalTruthShowerList::iterator iter;
   for (iter = shows.begin(); iter != shows.end(); ++iter) {
      DBCALTruthShower *bcaltruth = new DBCALTruthShower;
      bcaltruth->track   = iter->getTrack();
      bcaltruth->ptype   = iter->getPtype();
      bcaltruth->primary = (iter->getPrimary())? 1 : 0;
      bcaltruth->phi     = iter->getPhi();
      bcaltruth->r       = iter->getR();
      bcaltruth->z       = iter->getZ();
      bcaltruth->t       = iter->getT();
      bcaltruth->E       = iter->getE();
      bcaltruth->px      = iter->getPx();
      bcaltruth->py      = iter->getPy();
      bcaltruth->pz      = iter->getPz();
      data.push_back(bcaltruth);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DBCALTruthCell
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALTruthCell(hddm_s::HDDM *record,
                                   JFactory<DBCALTruthCell> *factory, 
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;

   vector<DBCALTruthCell*> data;

   const hddm_s::BcalTruthHitList &hits = record->getBcalTruthHits();
   hddm_s::BcalTruthHitList::iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      DBCALTruthCell *truthcell = new DBCALTruthCell();
      truthcell->module = hiter->getModule();
      truthcell->layer  = hiter->getLayer();
      truthcell->sector = hiter->getSector();
      truthcell->E      = hiter->getE();
      truthcell->t      = hiter->getT();
      truthcell->zLocal = hiter->getZLocal();
      data.push_back(truthcell);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DFCALTruthShower
//------------------
jerror_t DEventSourceHDDM::Extract_DFCALTruthShower(hddm_s::HDDM *record,
                                   JFactory<DFCALTruthShower> *factory,
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DFCALTruthShower*> data;
   JObject::oid_t id=1;

   const hddm_s::FcalTruthShowerList &shows = record->getFcalTruthShowers();
   hddm_s::FcalTruthShowerList::iterator iter;
   for (iter = shows.begin(); iter != shows.end(); ++iter) {
      DFCALTruthShower *dfcaltruthshower = new DFCALTruthShower(
            id++,
            iter->getX(),
            iter->getY(),
            iter->getZ(),
            iter->getPx(),
            iter->getPy(),
            iter->getPz(),
            iter->getE(),
            iter->getT(),
            iter->getPrimary(),
            iter->getTrack(),
            iter->getPtype()
            );
      data.push_back(dfcaltruthshower);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DFCALHit
//------------------
jerror_t DEventSourceHDDM::Extract_DFCALHit(hddm_s::HDDM *record,
                                   JFactory<DFCALHit> *factory, string tag,
                                   JEventLoop* eventLoop)
{
  /// Copies the data from the given hddm_s structure. This is called
  /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
  /// returs OBJECT_NOT_AVAILABLE immediately.
   
  if (factory == NULL)
     return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;

   // extract the FCAL Geometry (for isBlockActive() and positionOnFace())
   vector<const DFCALGeometry*> fcalGeomVect;
   eventLoop->Get( fcalGeomVect );
   if (fcalGeomVect.size() < 1)
      return OBJECT_NOT_AVAILABLE;
   const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);

   vector<DFCALHit*> data;
   int hitId = 0;

   if (tag == "") {
      const hddm_s::FcalHitList &hits = record->getFcalHits();
      hddm_s::FcalHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         int row = iter->getRow();
         int column = iter->getColumn();
 
         // Filter out non-physical blocks here
         if (!fcalGeom.isBlockActive(row, column))
            continue;

         // Get position of blocks on front face. (This should really come from
         // hdgeant directly so the positions can be shifted in mcsmear.)
         DVector2 pos = fcalGeom.positionOnFace(row, column);

         DFCALHit *mchit = new DFCALHit();
         mchit->row    = row;
         mchit->column = column;
         mchit->x      = pos.X();
         mchit->y      = pos.Y();
         mchit->E      = iter->getE();
         mchit->t      = iter->getT();
         mchit->id     = hitId++;
         data.push_back(mchit);
       }
    }
    else if (tag == "TRUTH") {
      const hddm_s::FcalTruthHitList &hits = record->getFcalTruthHits();
      hddm_s::FcalTruthHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         int row = iter->getRow();
         int column = iter->getColumn();
 
         // Filter out non-physical blocks here
         if (!fcalGeom.isBlockActive(row, column))
            continue;

         // Get position of blocks on front face. (This should really come from
         // hdgeant directly so the positions can be shifted in mcsmear.)
         DVector2 pos = fcalGeom.positionOnFace(row, column);

         DFCALHit *mchit = new DFCALHit();
         mchit->row    = row;
         mchit->column = column;
         mchit->x      = pos.X();
         mchit->y      = pos.Y();
         mchit->E      = iter->getE();
         mchit->t      = iter->getT();
         mchit->id     = hitId++;
         data.push_back(mchit);
      }
   }
  
  // Copy into factory
  factory->CopyTo(data);
  
  return NOERROR;
}

//------------------
// Extract_DCCALTruthShower
//------------------
jerror_t DEventSourceHDDM::Extract_DCCALTruthShower(hddm_s::HDDM *record,
                                   JFactory<DCCALTruthShower> *factory,
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DCCALTruthShower*> data;
   JObject::oid_t id=1;

   const hddm_s::CcalTruthShowerList &shows = record->getCcalTruthShowers();
   hddm_s::CcalTruthShowerList::iterator iter;
   for (iter = shows.begin(); iter != shows.end(); ++iter) {
      DCCALTruthShower *dccaltruthshower = new DCCALTruthShower(
            id++,
            iter->getX(),
            iter->getY(),
            iter->getZ(),
            iter->getPx(),
            iter->getPy(),
            iter->getPz(),
            iter->getE(),
            iter->getT(),
            iter->getPrimary(),
            iter->getTrack(),
            iter->getPtype()
            );
      data.push_back(dccaltruthshower);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DCCALHit
//------------------
jerror_t DEventSourceHDDM::Extract_DCCALHit(hddm_s::HDDM *record,
                                   JFactory<DCCALHit> *factory, string tag,
                                   JEventLoop* eventLoop)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returs OBJECT_NOT_AVAILABLE immediately.

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;

   // extract the CCAL Geometry (for isBlockActive() and positionOnFace())
   vector<const DCCALGeometry*> ccalGeomVect;
   eventLoop->Get( ccalGeomVect );
   if (ccalGeomVect.size() < 1)
      return OBJECT_NOT_AVAILABLE;
   const DCCALGeometry& ccalGeom = *(ccalGeomVect[0]);

   vector<DCCALHit*> data;
   int hitId = 0;

   if (tag == "") {
      const hddm_s::CcalHitList &hits = record->getCcalHits();
      hddm_s::CcalHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         int row = iter->getRow();
         int column = iter->getColumn();
 
         // Filter out non-physical blocks here
         if (!ccalGeom.isBlockActive(row, column))
            continue;

         // Get position of blocks on front face. (This should really come from
         // hdgeant directly so the poisitions can be shifted in mcsmear.)
         DVector2 pos = ccalGeom.positionOnFace(row, column);

         DCCALHit *mchit = new DCCALHit();
         mchit->row    = row;
         mchit->column = column;
         mchit->x      = pos.X();
         mchit->y      = pos.Y();
         mchit->E      = iter->getE();
         mchit->t      = iter->getT();
         mchit->id     = hitId++;
         data.push_back(mchit);
      }
   }

   else if (tag == "TRUTH") {
      const hddm_s::CcalTruthHitList &hits = record->getCcalTruthHits();
      hddm_s::CcalTruthHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         int row = iter->getRow();
         int column = iter->getColumn();
 
         // Filter out non-physical blocks here
         if (!ccalGeom.isBlockActive(row, column))
            continue;

         // Get position of blocks on front face. (This should really come from
         // hdgeant directly so the poisitions can be shifted in mcsmear.)
         DVector2 pos = ccalGeom.positionOnFace(row, column);

         DCCALHit *mchit = new DCCALHit();
         mchit->row    = row;
         mchit->column = column;
         mchit->x      = pos.X();
         mchit->y      = pos.Y();
         mchit->E      = iter->getE();
         mchit->t      = iter->getT();
         mchit->id     = hitId++;
         data.push_back(mchit);
      }
   }
  
  // Copy into factory
  factory->CopyTo(data);
  
  return NOERROR;
}

//------------------
// Extract_DMCTrajectoryPoint
//------------------
jerror_t DEventSourceHDDM::Extract_DMCTrajectoryPoint(hddm_s::HDDM *record,
                                   JFactory<DMCTrajectoryPoint> *factory, 
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
   
   vector<DMCTrajectoryPoint*> data;

   const hddm_s::McTrajectoryPointList &pts = record->getMcTrajectoryPoints();
   hddm_s::McTrajectoryPointList::iterator iter;
   for (iter = pts.begin(); iter != pts.end(); ++iter) {
      DMCTrajectoryPoint *p = new DMCTrajectoryPoint;
      p->x             = iter->getX();
      p->y             = iter->getY();
      p->z             = iter->getZ();
      p->t             = iter->getT();
      p->px            = iter->getPx();
      p->py            = iter->getPy();
      p->pz            = iter->getPz();
      p->E             = iter->getE();
      p->dE            = iter->getDE();
      p->primary_track = iter->getPrimary_track();
      p->track         = iter->getTrack();
      p->part          = iter->getPart();
      p->radlen        = iter->getRadlen();
      p->step          = iter->getStep();
      p->mech          = iter->getMech();
      data.push_back(p);
   }
   
   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DTOFTruth
//------------------
jerror_t DEventSourceHDDM::Extract_DTOFTruth(hddm_s::HDDM *record,
                                   JFactory<DTOFTruth>* factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;
  
   vector<DTOFTruth*> data;

   const hddm_s::FtofTruthPointList &points = record->getFtofTruthPoints();
   hddm_s::FtofTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      DTOFTruth *toftruth = new DTOFTruth;
      toftruth->primary = iter->getPrimary();
      toftruth->track   = iter->getTrack();
      toftruth->x       = iter->getX();
      toftruth->y       = iter->getY();
      toftruth->z       = iter->getZ();
      toftruth->t       = iter->getT();
      toftruth->px      = iter->getPx();
      toftruth->py      = iter->getPy();
      toftruth->pz      = iter->getPz();
      toftruth->E       = iter->getE();
      toftruth->ptype   = iter->getPtype();
      data.push_back(toftruth);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DTOFHit
//------------------
jerror_t DEventSourceHDDM::Extract_DTOFHit( hddm_s::HDDM *record,
                                   JFactory<DTOFHit>* factory,
                                   JFactory<DTOFHitMC> *factoryMC,
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;

   vector<DTOFHit*> data;
   vector<DTOFHitMC*> dataMC;

   const hddm_s::FtofCounterList &ctrs = record->getFtofCounters();
   hddm_s::FtofCounterList::iterator iter;
   for (iter = ctrs.begin(); iter != ctrs.end(); ++iter) {
      if (tag == "") {
         vector<DTOFHit*> north_hits;
         vector<DTOFHit*> south_hits;

         // Loop over north AND south hits
         const hddm_s::FtofHitList &hits = iter->getFtofHits();
         hddm_s::FtofHitList::iterator hiter;
         for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
            DTOFHit *tofhit = new DTOFHit;
            tofhit->bar   = hiter->getBar();
            tofhit->plane = hiter->getPlane();
            tofhit->end   = hiter->getEnd();
            tofhit->dE    = hiter->getDE();
            tofhit->t     = hiter->getT();
            data.push_back(tofhit);
            if (tofhit->end == 0)
               north_hits.push_back(tofhit);
            else 
               south_hits.push_back(tofhit);
         }

         // return truth hits in a different factory
         const hddm_s::FtofTruthHitList &truths = iter->getFtofTruthHits();
         hddm_s::FtofTruthHitList::iterator titer;
         unsigned int north_mchits = 0;
         unsigned int south_mchits = 0;
         for (titer = truths.begin(); titer != truths.end(); ++titer) {
            DTOFHitMC *tofmchit = new DTOFHitMC;
            tofmchit->bar    = titer->getBar();
            tofmchit->plane  = titer->getPlane();
            tofmchit->end    = titer->getEnd();
            tofmchit->itrack = titer->getFtofTruthExtra(0).getItrack();
            tofmchit->ptype  = titer->getFtofTruthExtra(0).getPtype();
            tofmchit->dist   = titer->getFtofTruthExtra(0).getDist();
            tofmchit->x      = titer->getFtofTruthExtra(0).getX();
            tofmchit->y      = titer->getFtofTruthExtra(0).getY();
            tofmchit->z      = titer->getFtofTruthExtra(0).getZ();
            tofmchit->px     = titer->getFtofTruthExtra(0).getPx();
            tofmchit->py     = titer->getFtofTruthExtra(0).getPy();
            tofmchit->pz     = titer->getFtofTruthExtra(0).getPz();
            tofmchit->E      = titer->getFtofTruthExtra(0).getE();
            dataMC.push_back(tofmchit);

            // best-guess at tofhit-tofMChit association, not exact
            if (tofmchit->end == 0) {
               if (north_mchits < north_hits.size()) {
                  north_hits[north_mchits]->AddAssociatedObject(tofmchit);
               }
               north_mchits++;
            }
            else {
               if (south_mchits < south_hits.size()) {
                  south_hits[south_mchits]->AddAssociatedObject(tofmchit);
               }
               south_mchits++;
            }
         }
     }

     else if (tag == "TRUTH") {
        const hddm_s::FtofTruthHitList &truths = iter->getFtofTruthHits();
        hddm_s::FtofTruthHitList::iterator titer;
        for (titer = truths.begin(); titer != truths.end(); ++titer) {
           DTOFHit *tofhit = new DTOFHit;
           tofhit->bar   = titer->getBar();
           tofhit->plane = titer->getPlane();
           tofhit->end   = titer->getEnd();
           tofhit->dE    = titer->getDE();
           tofhit->t     = titer->getT();
           data.push_back(tofhit);

           DTOFHitMC *tofmchit = new DTOFHitMC;
           tofmchit->bar    = tofhit->bar;
           tofmchit->plane  = tofhit->plane;
           tofmchit->end    = tofhit->end;
           tofmchit->itrack = titer->getFtofTruthExtra().getItrack();
           tofmchit->ptype  = titer->getFtofTruthExtra().getPtype();
           tofmchit->dist   = titer->getFtofTruthExtra().getDist();
           tofmchit->x      = titer->getFtofTruthExtra().getX();
           tofmchit->y      = titer->getFtofTruthExtra().getY();
           tofmchit->z      = titer->getFtofTruthExtra().getZ();
           tofmchit->px     = titer->getFtofTruthExtra().getPx();
           tofmchit->py     = titer->getFtofTruthExtra().getPy();
           tofmchit->pz     = titer->getFtofTruthExtra().getPz();
           tofmchit->E      = titer->getFtofTruthExtra().getE();
           dataMC.push_back(tofmchit);
           tofhit->AddAssociatedObject(tofmchit);
        }
     }
   }

   // Copy into factory
   factory->CopyTo(data);
   factoryMC->CopyTo(dataMC);

   return NOERROR;
}

//------------------
// Extract_DSCHit
//------------------
jerror_t DEventSourceHDDM::Extract_DSCHit(hddm_s::HDDM *record,
                                   JFactory<DSCHit>* factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;

   vector<DSCHit*> data;
  
   if (tag == "") {
      const hddm_s::StcHitList &hits = record->getStcHits();
      hddm_s::StcHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         DSCHit *hit = new DSCHit;
         hit->sector = iter->getSector();
         hit->dE = iter->getDE();
         hit->t = iter->getT();
         hit->sigma_t = 0.350;
         data.push_back(hit);
      }
   }
   else if (tag == "TRUTH") {
      const hddm_s::StcTruthHitList &hits = record->getStcTruthHits();
      hddm_s::StcTruthHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         DSCHit *hit = new DSCHit;
         hit->sector = iter->getSector();
         hit->dE = iter->getDE();
         hit->t = iter->getT();
         hit->sigma_t = 1e-6;
         data.push_back(hit);
      }
   }

  // Copy into factory
  factory->CopyTo(data);

  return NOERROR;
}

//------------------
// Extract_DSCTruthHit
//------------------
jerror_t DEventSourceHDDM::Extract_DSCTruthHit(hddm_s::HDDM *record,
                                   JFactory<DSCTruthHit>* factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;

   vector<DSCTruthHit*> data;

   const hddm_s::StcTruthPointList &points = record->getStcTruthPoints();
   hddm_s::StcTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      DSCTruthHit *hit = new DSCTruthHit;
      hit->dEdx    = iter->getDEdx();
      hit->phi     = iter->getPhi();
      hit->primary = iter->getPrimary();
      hit->ptype   = iter->getPtype();
      hit->r       = iter->getR();
      hit->t       = iter->getT();
      hit->z       = iter->getZ();
      hit->track   = iter->getTrack();
      hit->sector  = iter->getSector();
      data.push_back(hit);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DTrackTimeBased
//------------------
jerror_t DEventSourceHDDM::Extract_DTrackTimeBased(hddm_s::HDDM *record,
                                   JFactory<DTrackTimeBased> *factory, 
                                   string tag, int runnumber)
{
   // Note: Since this is a reconstructed factory, we want to generally return OBJECT_NOT_AVAILABLE
   // rather than NOERROR. The reason being that the caller interprets "NOERROR" to mean "yes I
   // usually can provide objects of that type, but this event has none." This will cause it to
   // skip any attempt at reconstruction. On the other hand, a value of "OBJECT_NOT_AVAILABLE" tells
   // it "I cannot provide those type of objects for this event."

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;

   vector<DTrackTimeBased*> data;
   vector<DReferenceTrajectory*> rts;

   const hddm_s::TracktimebasedList &ttbs = record->getTracktimebaseds();

   // Get enough DReferenceTrajectory objects for all of the DTrackTimeBased Objects
   // we're about to read in. This seems a little complicated, but that's because it
   // is expensive to allocate these things so we recycle as much as possible.
   list<DReferenceTrajectory*> my_rts;
   pthread_mutex_lock(&rt_mutex);
   while (my_rts.size() < (unsigned int)ttbs.size()) {
      if (rt_pool.size() > 0) {
         my_rts.push_back(rt_pool.back());
         rt_pool.pop_back();
      }
      else {
         if (dapp && !bfield)
            // delay getting the bfield object until we have to!
            bfield = dapp->GetBfield();
         if (dapp && !geom)
            geom = dapp->GetDGeometry(runnumber);
         my_rts.push_back(new DReferenceTrajectory(bfield));
      }
   }
   pthread_mutex_unlock(&rt_mutex);

   hddm_s::TracktimebasedList::iterator iter;
   for (iter = ttbs.begin(); iter != ttbs.end(); ++iter) {
      DVector3 mom(iter->getMomentum().getPx(),
                   iter->getMomentum().getPy(),
                   iter->getMomentum().getPz());
      DVector3 pos(iter->getOrigin().getVx(),
                   iter->getOrigin().getVy(),
                   iter->getOrigin().getVz());
      DTrackTimeBased *track = new DTrackTimeBased();
      track->setMomentum(mom);
      track->setPosition(pos);
      track->setCharge(iter->getProperties().getCharge());
      track->setMass(iter->getProperties().getMass());
      track->setPID(IDTrack(iter->getProperties().getCharge(),
                            iter->getProperties().getMass()));
      track->chisq = iter->getChisq();
      track->Ndof = iter->getNdof();
      track->FOM = iter->getFOM();
      track->candidateid = iter->getCandidateid();
      track->id = iter->getId();

      // Reconstitute errorMatrix
      string str_vals = iter->getErrorMatrix().getVals();
      DMatrixDSym errMatrix;
      StringToDMatrixDSym(str_vals, errMatrix,
                          iter->getErrorMatrix().getNrows(),
                          iter->getErrorMatrix().getNcols());
      track->setErrorMatrix(errMatrix);

      // Reconstitute TrackingErrorMatrix
      str_vals = iter->getTrackingErrorMatrix().getVals();
      DMatrixDSym TrackingErrorMatrix;
      StringToDMatrixDSym(str_vals, TrackingErrorMatrix,
                          iter->getTrackingErrorMatrix().getNrows(),
                          iter->getTrackingErrorMatrix().getNcols());
      track->setTrackingErrorMatrix(TrackingErrorMatrix);
 
      // Use DReferenceTrajectory objects (either recycled or new)
      DReferenceTrajectory *rt = my_rts.back();
      my_rts.pop_back();
      if (rt) {
         rt->SetMass(track->mass());
         rt->SetDGeometry(geom);
         rt->Swim(pos, mom, track->charge());
         rts.push_back(rt);
      }
      track->rt = rt;

      data.push_back(track);
   }

   // Copy into factory
   if (ttbs.size() > 0){
      factory->CopyTo(data);
      
      // Add DReferenceTrajectory objects to rt_by_event so they can be
      // deleted later. The rt_by_event maintains lists indexed by the
      // hddm record pointer since multiple threads may be calling us.
      // Note that we first look to see if a list already exists for
      // this event and append to it if it does. This is so we can the
      // same list for all objects that use DReferenceTrajectories.
      pthread_mutex_lock(&rt_mutex);
      map<hddm_s::HDDM*, vector<DReferenceTrajectory*> >::iterator iter =
                                                   rt_by_event.find(record);
      if (iter != rt_by_event.end()) {
         vector<DReferenceTrajectory*> &my_rts = iter->second;
         my_rts.insert(my_rts.end(), rts.begin(), rts.end());
      }
      else {
         rt_by_event[record] = rts;
      }
      pthread_mutex_unlock(&rt_mutex);

      // If the event had a s_Tracktimebased_t pointer, then report
      // back that we read them in from the file. Otherwise, report
      // OBJECT_NOT_AVAILABLE
      return NOERROR;
   }

   // If we get to here then there was not even a placeholder in the HDDM file.
   // Return OBJECT_NOT_AVAILABLE to indicate reconstruction should be tried.
   return OBJECT_NOT_AVAILABLE;
}


//-------------------------------
// StringToDMatrixDSym
//-------------------------------
string DEventSourceHDDM::StringToDMatrixDSym(string &str_vals, DMatrixDSym &mat,
                                             int Nrows, int Ncols)
{
   /// This is the inverse of the DMatrixDSymToString method in the
   /// danahddm plugin.

   // Convert the given string into a symmetric matrix
   mat.ResizeTo(Nrows, Ncols);
   stringstream ss(str_vals);
   for (int irow=0; irow<mat.GetNrows(); irow++) {
      for (int icol=irow; icol<mat.GetNcols(); icol++) {
         ss >> mat[irow][icol];
         mat[icol][irow] = mat[irow][icol];
      }
   }
   
   return ss.str();
}

//------------------
// Extract_DTAGMHit
//------------------
jerror_t DEventSourceHDDM::Extract_DTAGMHit(hddm_s::HDDM *record,
                                   JFactory<DTAGMHit>* factory,
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;

   vector<DTAGMHit*> data;

   // loop over microChannel/taggerHit records
   const hddm_s::MicroChannelList &tags = record->getMicroChannels();
   hddm_s::MicroChannelList::iterator iter;
   for (iter = tags.begin(); iter != tags.end(); ++iter) {
      if (tag == "") {
         const hddm_s::TaggerHitList &hits = iter->getTaggerHits();
         hddm_s::TaggerHitList::iterator hiter;
         for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
            DTAGMHit *taghit = new DTAGMHit();
            taghit->E = hiter->getE();
            taghit->t = hiter->getT();
            taghit->npix_fadc = hiter->getNpe();
            taghit->time_fadc = hiter->getTADC();
            taghit->column = hiter->getColumn();
            taghit->row = hiter->getRow();
            data.push_back(taghit);
         }
      }
      else if (tag == "TRUTH") {
         const hddm_s::TaggerTruthHitList &hits = iter->getTaggerTruthHits();
         hddm_s::TaggerTruthHitList::iterator hiter;
         for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
            DTAGMHit *taghit = new DTAGMHit();
            taghit->E = hiter->getE();
            taghit->t = hiter->getT();
            taghit->npix_fadc = hiter->getDE() * 1e5; // ~1e5 pixels/GeV
            taghit->time_fadc = hiter->getT();
            taghit->column = hiter->getColumn();
            taghit->row = hiter->getRow();
            data.push_back(taghit);
         }
      }
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DTAGHHit
//------------------
jerror_t DEventSourceHDDM::Extract_DTAGHHit( hddm_s::HDDM *record,
                                   JFactory<DTAGHHit>* factory,
                                   string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;

   vector<DTAGHHit*> data;
  
   // loop over hodoChannel/taggerHit records
   const hddm_s::HodoChannelList &tags = record->getHodoChannels();
   hddm_s::HodoChannelList::iterator iter;
   for (iter = tags.begin(); iter != tags.end(); ++iter) {
      if (tag == "") {
         const hddm_s::TaggerHitList &hits = iter->getTaggerHits();
         hddm_s::TaggerHitList::iterator hiter;
         for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
            DTAGHHit *taghit = new DTAGHHit();
            taghit->E = hiter->getE();
            taghit->t = hiter->getT();
            taghit->npe_fadc = hiter->getNpe();
            taghit->time_fadc = hiter->getTADC();
            taghit->counter_id = hiter->getCounterId();
            data.push_back(taghit);
         }
      }
      else if (tag == "TRUTH") {
         const hddm_s::TaggerTruthHitList &hits = iter->getTaggerTruthHits();
         hddm_s::TaggerTruthHitList::iterator hiter;
         for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
            DTAGHHit *taghit = new DTAGHHit();
            taghit->E = hiter->getE();
            taghit->t = hiter->getT();
            taghit->npe_fadc = hiter->getDE() * 5e5; // ~5e5 pe/GeV
            taghit->time_fadc = hiter->getT();
            taghit->counter_id = hiter->getCounterId();
            data.push_back(taghit);
         }
      }
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

Particle_t DEventSourceHDDM::IDTrack(float locCharge, float locMass) const
{
   float locMassTolerance = 0.010;
   if (locCharge > 0.1) // Positive particles
   {
      if (fabs(locMass - ParticleMass(Proton)) < locMassTolerance) return Proton;
      if (fabs(locMass - ParticleMass(PiPlus)) < locMassTolerance) return PiPlus;
      if (fabs(locMass - ParticleMass(KPlus)) < locMassTolerance) return KPlus;
      if (fabs(locMass - ParticleMass(Positron)) < locMassTolerance) return Positron;
      if (fabs(locMass - ParticleMass(MuonPlus)) < locMassTolerance) return MuonPlus;
   }
   else if(locCharge < -0.1) // Negative particles
   {
      if (fabs(locMass - ParticleMass(PiMinus)) < locMassTolerance) return PiMinus;
      if (fabs(locMass - ParticleMass(KMinus)) < locMassTolerance) return KMinus;
      if (fabs(locMass - ParticleMass(MuonMinus)) < locMassTolerance) return MuonMinus;
      if (fabs(locMass - ParticleMass(Electron)) < locMassTolerance) return Electron;
      if (fabs(locMass - ParticleMass(AntiProton)) < locMassTolerance) return AntiProton;
   }
   else //Neutral Track
   {
      if (fabs(locMass - ParticleMass(Gamma)) < locMassTolerance) return Gamma;
      if (fabs(locMass - ParticleMass(Neutron)) < locMassTolerance) return Neutron;
   }
   return Unknown;
}

//------------------
// Extract_DCereRichHit
// added by yqiang Oct 3, 2012
//------------------
jerror_t DEventSourceHDDM::Extract_DRichHit(hddm_s::HDDM *record,
                                   JFactory<DRichHit>* factory, string tag)
{
   /// Copies the data from the given hddm_s structure. This is called
   /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;

   vector<DRichHit*> data;

   const hddm_s::RichHitList &hits = record->getRichHits();
   hddm_s::RichHitList::iterator iter;
   for (iter = hits.begin(); iter != hits.end(); ++iter) {
      DRichHit *hit = new DRichHit;
      hit->x = iter->getX();
      hit->y = iter->getY();
      hit->z = iter->getZ();
      hit->t = iter->getT();
      data.push_back(hit);
   }

  // Copy into factory
  factory->CopyTo(data);

  return NOERROR;
}

//------------------
// Extract_DCereHit
// added by yqiang Oct 11, 2012
//------------------
jerror_t DEventSourceHDDM::Extract_DCereHit(hddm_s::HDDM *record,
                                   JFactory<DCereHit>* factory, string tag)
{
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "" && tag != "TRUTH")
      return OBJECT_NOT_AVAILABLE;

   vector<DCereHit*> data;

   if (tag == "") {
      const hddm_s::CereHitList &hits = record->getCereHits();
      hddm_s::CereHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         DCereHit *hit = new DCereHit;
         hit->sector = iter->getSector();
         hit->pe = iter->getPe();
         hit->t = iter->getT();
         data.push_back(hit);
      }
   }
   else if (tag == "TRUTH") {
      const hddm_s::CereTruthHitList &hits = record->getCereTruthHits();
      hddm_s::CereTruthHitList::iterator iter;
      for (iter = hits.begin(); iter != hits.end(); ++iter) {
         DCereHit *hit = new DCereHit;
         hit->sector = iter->getSector();
         hit->pe = iter->getPe();
         hit->t = iter->getT();
         data.push_back(hit);
      }
   }

   // copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

/*
 * Extract_DRichTruthHit
 * added by yqiang Oct 7, 2013
 */
jerror_t DEventSourceHDDM::Extract_DRichTruthHit(hddm_s::HDDM *record,
                                   JFactory<DRichTruthHit>* factory,
                                   string tag)
{
   if (factory == NULL)
      return OBJECT_NOT_AVAILABLE;
   if (tag != "")
      return OBJECT_NOT_AVAILABLE;

   vector<DRichTruthHit*> data;

   const hddm_s::RichTruthPointList &points = record->getRichTruthPoints();
   hddm_s::RichTruthPointList::iterator iter;
   for (iter = points.begin(); iter != points.end(); ++iter) {
      DRichTruthHit *hit = new DRichTruthHit;
      hit->x       = iter->getX();
      hit->y       = iter->getY();
      hit->z       = iter->getZ();
      hit->px      = iter->getPx();
      hit->py      = iter->getPy();
      hit->pz      = iter->getPz();
      hit->t       = iter->getT();
      hit->E       = iter->getE();
      hit->track   = iter->getTrack();
      hit->primary = iter->getPrimary();
      hit->ptype   = iter->getPtype();
      data.push_back(hit);
   }

   factory->CopyTo(data);
   return NOERROR;
}

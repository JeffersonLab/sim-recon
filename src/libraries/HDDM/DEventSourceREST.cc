//
// Author: Richard Jones  June 29, 2012
//
//
// DEventSourceREST methods
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <climits>

#include <JANA/JFactory_base.h>
#include <JANA/JEventLoop.h>
#include <JANA/JEvent.h>
#include <DANA/DStatusBits.h>

#include <DVector2.h>
#include <DEventSourceREST.h>

//----------------
// Constructor
//----------------
DEventSourceREST::DEventSourceREST(const char* source_name)
 : JEventSource(source_name)
{
   /// Constructor for DEventSourceREST object
   ifs = new std::ifstream(source_name);
   if (ifs && ifs->is_open()) {
      // hddm_r::istream constructor can throw a std::runtime_error
      // which is not being caught here -- policy question in JANA:
      // who catches the exceptions, top-level user code or here?
      fin = new hddm_r::istream(*ifs);
   }
   else {
      // One might want to throw an exception or report an error here.
      fin = NULL;
   }
}

//----------------
// Destructor
//----------------
DEventSourceREST::~DEventSourceREST()
{  
  if (fin) {
    delete fin;
  }
  if (ifs) {
    delete ifs;
  }
}

//----------------
// GetEvent
//----------------
jerror_t DEventSourceREST::GetEvent(JEvent &event)
{
   /// Implementation of JEventSource virtual function

   if (!fin) {
      return EVENT_SOURCE_NOT_OPEN;
   }

   // Each open hddm file takes up about 1M of memory so it's
   // worthwhile to close it as soon as we can.
   if (ifs->eof()) {
      delete fin;
      fin = NULL;
      delete ifs;
      ifs = NULL;

      return NO_MORE_EVENTS_IN_SOURCE;
   }

#if HDDM_SETPOSITION_EXAMPLE
   static std::ifstream fevlist("events.list");
   static int events_to_go = 0;
   if (events_to_go-- == 0 && fevlist.good()) {
      uint64_t start;
      uint32_t offset, status;
      fevlist >> start >> offset >> status >> events_to_go;
      if (fevlist.good())
         fin->setPosition(hddm_r::streamposition(start, offset, status));
   }
#endif

#if HDDM_GETPOSITION_EXAMPLE
   hddm_r::streamposition pos(fin->getPosition());
   // Later on below, if this event passes all of your selection cuts
   // then you might want to write the event position to output, as in
   // std::cout << "interesting event found at " 
   //           << pos.start << "," << pos.offset << "," << pos.status
   //           << std::endl;
#endif

   hddm_r::HDDM *record = new hddm_r::HDDM();
   try{
      if (! (*fin >> *record)) {
         delete fin;
         fin = NULL;
         delete ifs;
         ifs = NULL;
	     return NO_MORE_EVENTS_IN_SOURCE;
      }
   }catch(std::runtime_error &e){
      cerr << "Exception caught while trying to read REST file!" << endl;
	  cerr << e.what() << endl;
	  _DBG__;
	  return NO_MORE_EVENTS_IN_SOURCE;
   }

   // Copy the reference info into the JEvent object
   while (true) {
      hddm_r::ReconstructedPhysicsEvent &re
            = record->getReconstructedPhysicsEvent();
      int runno = re.getRunNo();
      int eventno = re.getEventNo();
      if (runno == 0 && eventno == 0) {
         // found a comment record, print comment strings and continue
         const hddm_r::CommentList &comments = re.getComments();
         hddm_r::CommentList::iterator iter;
         for (iter = comments.begin(); iter != comments.end(); ++iter) {
            std::cout << "   | " << iter->getText() << std::endl;
         }

         //set version string
         const hddm_r::DataVersionStringList& locVersionStrings = re.getDataVersionStrings();
         hddm_r::DataVersionStringList::iterator Versioniter;
         for (Versioniter = locVersionStrings.begin(); Versioniter != locVersionStrings.end(); ++Versioniter) {
        	 string HDDM_DATA_VERSION_STRING = Versioniter->getText();
             gPARMS->SetDefaultParameter("REST:DATAVERSIONSTRING", HDDM_DATA_VERSION_STRING);
         }

         //set REST calib context
         const hddm_r::CcdbContextList& locContextStrings = re.getCcdbContexts();
         hddm_r::CcdbContextList::iterator Contextiter;
         for (Contextiter = locContextStrings.begin(); Contextiter != locContextStrings.end(); ++Contextiter) {
        	 string REST_JANA_CALIB_CONTEXT = Contextiter->getText();
             gPARMS->SetDefaultParameter("REST:JANACALIBCONTEXT", REST_JANA_CALIB_CONTEXT);
         }

         if (! (*fin >> *record)) {
            delete fin;
            fin = NULL;
            delete ifs;
            ifs = NULL;
	        return NO_MORE_EVENTS_IN_SOURCE;
         }

         continue;
      }
      event.SetEventNumber(re.getEventNo());
      event.SetRunNumber(re.getRunNo());
      event.SetJEventSource(this);
      event.SetRef(record);
      event.SetStatusBit(kSTATUS_REST);
      event.SetStatusBit(kSTATUS_FROM_FILE);
	  event.SetStatusBit(kSTATUS_PHYSICS_EVENT);

	  ++Nevents_read;
      break;
   }
 
   return NOERROR;
}

//----------------
// FreeEvent
//----------------
void DEventSourceREST::FreeEvent(JEvent &event)
{
   hddm_r::HDDM *record = (hddm_r::HDDM*)event.GetRef();
   delete record;
}

//----------------
// GetObjects
//----------------
jerror_t DEventSourceREST::GetObjects(JEvent &event, JFactory_base *factory)
{
   /// This gets called through the virtual method of the
   /// JEventSource base class. It creates the objects of the type
   /// on which factory is based.  It uses the HDDM_Element object
   /// kept in the ref field of the JEvent object passed.

   // We must have a factory to hold the data
   if (!factory) {
      throw RESOURCE_UNAVAILABLE;
   }

   // The ref field of the JEvent is just the HDDM record pointer.
   hddm_r::HDDM *record = (hddm_r::HDDM*)event.GetRef();
   if (!record) {
      throw RESOURCE_UNAVAILABLE;
   }

   JEventLoop* locEventLoop = event.GetJEventLoop();
   string dataClassName = factory->GetDataClassName();
   
	//Get target center
		//multiple reader threads can access this object: need lock
	bool locNewRunNumber = false;
	unsigned int locRunNumber = event.GetRunNumber();
	LockRead();
	{
		locNewRunNumber = (dTargetCenterZMap.find(locRunNumber) == dTargetCenterZMap.end());
	}
	UnlockRead();
	if(locNewRunNumber)
	{
		DApplication* dapp = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
		DGeometry* locGeometry = dapp->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
		double locTargetCenterZ = 0.0;
		locGeometry->GetTargetZ(locTargetCenterZ);

		vector<double> locBeamPeriodVector;
        if(locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector))
            throw JException("Could not load CCDB table: PHOTON_BEAM/RF/beam_period");
		double locBeamBunchPeriod = locBeamPeriodVector[0];

		LockRead();
		{
			dTargetCenterZMap[locRunNumber] = locTargetCenterZ;
			dBeamBunchPeriodMap[locRunNumber] = locBeamBunchPeriod;
		}
		UnlockRead();
	}

   if (dataClassName =="DMCReaction") {
      return Extract_DMCReaction(record,
                     dynamic_cast<JFactory<DMCReaction>*>(factory), locEventLoop);
   }
   if (dataClassName =="DRFTime") {
      return Extract_DRFTime(record,
                     dynamic_cast<JFactory<DRFTime>*>(factory), locEventLoop);
   }
   if (dataClassName =="DBeamPhoton") {
      return Extract_DBeamPhoton(record,
                     dynamic_cast<JFactory<DBeamPhoton>*>(factory),
                     locEventLoop);
   }
   if (dataClassName =="DMCThrown") {
      return Extract_DMCThrown(record,
                     dynamic_cast<JFactory<DMCThrown>*>(factory));
   }
   if (dataClassName =="DTOFPoint") {
      return Extract_DTOFPoint(record,
                     dynamic_cast<JFactory<DTOFPoint>*>(factory));
   }
   if (dataClassName =="DSCHit") {
      return Extract_DSCHit(record,
                     dynamic_cast<JFactory<DSCHit>*>(factory));
   }
   if (dataClassName =="DFCALShower") {
      return Extract_DFCALShower(record,
                     dynamic_cast<JFactory<DFCALShower>*>(factory));
   }
   if (dataClassName =="DBCALShower") {
      return Extract_DBCALShower(record,
                     dynamic_cast<JFactory<DBCALShower>*>(factory));
   }
   if (dataClassName =="DTrackTimeBased") {
      return Extract_DTrackTimeBased(record,
                     dynamic_cast<JFactory<DTrackTimeBased>*>(factory), locEventLoop);
   }
   if (dataClassName =="DTrigger") {
      return Extract_DTrigger(record,
                     dynamic_cast<JFactory<DTrigger>*>(factory));
   }
   if (dataClassName =="DDetectorMatches") {
      return Extract_DDetectorMatches(locEventLoop, record,
                     dynamic_cast<JFactory<DDetectorMatches>*>(factory));
   }

   return OBJECT_NOT_AVAILABLE;
}

//------------------
// Extract_DMCReaction
//------------------
jerror_t DEventSourceREST::Extract_DMCReaction(hddm_r::HDDM *record,
                                   JFactory<DMCReaction> *factory, JEventLoop* locEventLoop)
{
   /// Copies the data from the Reaction hddm class. This is called
   /// from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
   
   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   std::string tag = (factory->Tag())? factory->Tag() : "";

	double locTargetCenterZ = 0.0;
	int locRunNumber = locEventLoop->GetJEvent().GetRunNumber();
	LockRead();
	{
		locTargetCenterZ = dTargetCenterZMap[locRunNumber];
	}
	UnlockRead();
	DVector3 locPosition(0.0, 0.0, locTargetCenterZ);

   vector<DMCReaction*> dmcreactions;

   // loop over reaction records
   const hddm_r::ReactionList &reactions = record->getReactions();
   hddm_r::ReactionList::iterator iter;
   for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      if (iter->getJtag() != tag) {
         continue;
      }
      DMCReaction *mcreaction = new DMCReaction;
      dmcreactions.push_back(mcreaction);
      mcreaction->type = iter->getType();
      mcreaction->weight = iter->getWeight();
      double Ebeam = iter->getEbeam();

      hddm_r::Origin &origin = iter->getVertex().getOrigin();
      double torig = origin.getT();
      double zorig = origin.getVz();

      mcreaction->beam.setPosition(locPosition);
      mcreaction->beam.setMomentum(DVector3(0.0, 0.0, Ebeam));
      mcreaction->beam.setMass(0.0);
      mcreaction->beam.setCharge(0.0);
      mcreaction->beam.clearErrorMatrix();
      mcreaction->beam.setT0(0.0, 0.0, SYS_NULL);
      mcreaction->beam.setT1(0.0, 0.0, SYS_NULL);
      mcreaction->beam.setTime(torig - (zorig - locTargetCenterZ)/29.9792458);
      mcreaction->beam.setPID(Gamma);

      mcreaction->target.setPosition(locPosition);
      mcreaction->target.setMomentum(DVector3(0.0, 0.0, 0.0));
      Particle_t ttype = iter->getTargetType();
      mcreaction->target.setPID((Particle_t)ttype);
      mcreaction->target.setMass(ParticleMass(ttype));
      mcreaction->target.setCharge(ParticleCharge(ttype));
      mcreaction->target.clearErrorMatrix();
      mcreaction->target.setT0(0.0, 0.0, SYS_NULL);
      mcreaction->target.setT1(0.0, 0.0, SYS_NULL);
      mcreaction->target.setTime(torig - (zorig - locTargetCenterZ)/29.9792458);
   }
   
   // Copy into factories
   factory->CopyTo(dmcreactions);

   return NOERROR;
}

//------------------
// Extract_DRFTime
//------------------
jerror_t DEventSourceREST::Extract_DRFTime(hddm_r::HDDM *record,
                                   JFactory<DRFTime> *factory, JEventLoop* locEventLoop)
{
   if (factory==NULL)
      return OBJECT_NOT_AVAILABLE;
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DRFTime*> locRFTimes;

   // loop over RF-time records
   const hddm_r::RFtimeList &rftimes = record->getRFtimes();
   hddm_r::RFtimeList::iterator iter;
   for (iter = rftimes.begin(); iter != rftimes.end(); ++iter)
   {
      if (iter->getJtag() != tag)
         continue;
      DRFTime *locRFTime = new DRFTime;
      locRFTime->dTime = iter->getTsync();
      locRFTime->dTimeVariance = 0.0; //SET ME!!
      locRFTimes.push_back(locRFTime);
   }

	if(!locRFTimes.empty())
	{
		//found in the file, copy into factory and return
		factory->CopyTo(locRFTimes);
		return NOERROR;
	}

	//Not found in the file, so either:
		//Experimental data & it's missing: bail
		//MC data: generate it

	vector<const DBeamPhoton*> locMCGENPhotons;
	locEventLoop->Get(locMCGENPhotons, "MCGEN");
	if(locMCGENPhotons.empty())
		return OBJECT_NOT_AVAILABLE; //Experimental data & it's missing: bail

	//Is MC data. Either:
		//No tag: return photon_time propagated by +/- n*locBeamBunchPeriod to get close to 0.0
		//TRUTH tag: get exact t from DBeamPhoton tag MCGEN

	if(tag == "TRUTH")
	{
		DRFTime *locRFTime = new DRFTime;
		locRFTime->dTime = locMCGENPhotons[0]->time();
		locRFTime->dTimeVariance = 0.0;
		locRFTimes.push_back(locRFTime);
	}
	else
	{
		double locBeamBunchPeriod = 0.0;
		int locRunNumber = locEventLoop->GetJEvent().GetRunNumber();
		LockRead();
		{
			locBeamBunchPeriod = dBeamBunchPeriodMap[locRunNumber];
		}
		UnlockRead();

		//start with true RF time, increment/decrement by multiples of locBeamBunchPeriod ns until closest to 0
		double locTime = locMCGENPhotons[0]->time();
		int locNumRFBuckets = int(locTime/locBeamBunchPeriod);
		locTime -= double(locNumRFBuckets)*locBeamBunchPeriod;
		while(locTime > 0.5*locBeamBunchPeriod)
			locTime -= locBeamBunchPeriod;
		while(locTime < -0.5*locBeamBunchPeriod)
			locTime += locBeamBunchPeriod;

		DRFTime *locRFTime = new DRFTime;
		locRFTime->dTime = locTime;
		locRFTime->dTimeVariance = 0.0;
		locRFTimes.push_back(locRFTime);
	}

   // Copy into factories
   factory->CopyTo(locRFTimes);

   return NOERROR;
}

//------------------
// Extract_DBeamPhoton
//------------------
jerror_t DEventSourceREST::Extract_DBeamPhoton(hddm_r::HDDM *record,
                                   JFactory<DBeamPhoton> *factory,
                                   JEventLoop *eventLoop)
{
   /// This is called from JEventSourceREST::GetObjects. If factory is NULL,
   /// return OBJECT_NOT_AVAILABLE immediately. If factory tag="MCGEN" then
   /// copy the beam photon data from the Reaction hddm class.

   if (factory==NULL)
      return OBJECT_NOT_AVAILABLE;
   string tag = (factory->Tag())? factory->Tag() : "";

	vector<DBeamPhoton*> dbeam_photons;

	// extract the TAGH geometry
   vector<const DTAGHGeometry*> taghGeomVect;
   eventLoop->Get(taghGeomVect);
   if (taghGeomVect.empty())
      return OBJECT_NOT_AVAILABLE;
   const DTAGHGeometry* taghGeom = taghGeomVect[0];

   // extract the TAGM geometry
   vector<const DTAGMGeometry*> tagmGeomVect;
   eventLoop->Get(tagmGeomVect);
   if (tagmGeomVect.empty())
      return OBJECT_NOT_AVAILABLE;
   const DTAGMGeometry* tagmGeom = tagmGeomVect[0];

   if(tag == "MCGEN")
	{
		vector<const DMCReaction*> dmcreactions;
		eventLoop->Get(dmcreactions);

		for(size_t loc_i = 0; loc_i < dmcreactions.size(); ++loc_i)
		{
		   DBeamPhoton *beamphoton = new DBeamPhoton;
		   *(DKinematicData*)beamphoton = dmcreactions[loc_i]->beam;
		   if(!tagmGeom->E_to_column(beamphoton->energy(), beamphoton->dCounter))
		   	  taghGeom->E_to_counter(beamphoton->energy(), beamphoton->dCounter);
		   dbeam_photons.push_back(beamphoton);
		}

		// Copy into factories
		factory->CopyTo(dbeam_photons);

	   return NOERROR;
	}

	double locTargetCenterZ = 0.0;
	int locRunNumber = eventLoop->GetJEvent().GetRunNumber();
	LockRead();
	{
		locTargetCenterZ = dTargetCenterZMap[locRunNumber];
	}
	UnlockRead();

	DVector3 pos(0.0, 0.0, locTargetCenterZ);

   const hddm_r::TagmBeamPhotonList &locTagmBeamPhotonList = record->getTagmBeamPhotons();
   hddm_r::TagmBeamPhotonList::iterator locTAGMiter;
   for(locTAGMiter = locTagmBeamPhotonList.begin(); locTAGMiter != locTagmBeamPhotonList.end(); ++locTAGMiter)
   {
      if (locTAGMiter->getJtag() != tag)
         continue;

      DBeamPhoton* gamma = new DBeamPhoton();

		DVector3 mom(0.0, 0.0, locTAGMiter->getE());
		gamma->setPID(Gamma);
		gamma->setMomentum(mom);
		gamma->setPosition(pos);
		gamma->setCharge(0);
		gamma->setMass(0);
		gamma->setTime(locTAGMiter->getT());
		gamma->setT0(locTAGMiter->getT(), 0.200, SYS_TAGM);

        tagmGeom->E_to_column(locTAGMiter->getE(), gamma->dCounter);
		dbeam_photons.push_back(gamma);
   }


   const hddm_r::TaghBeamPhotonList &locTaghBeamPhotonList = record->getTaghBeamPhotons();
   hddm_r::TaghBeamPhotonList::iterator locTAGHiter;
   for(locTAGHiter = locTaghBeamPhotonList.begin(); locTAGHiter != locTaghBeamPhotonList.end(); ++locTAGHiter)
   {
      if (locTAGHiter->getJtag() != tag)
         continue;

      DBeamPhoton* gamma = new DBeamPhoton();

		DVector3 mom(0.0, 0.0, locTAGHiter->getE());
		gamma->setPID(Gamma);
		gamma->setMomentum(mom);
		gamma->setPosition(pos);
		gamma->setCharge(0);
		gamma->setMass(0);
		gamma->setTime(locTAGHiter->getT());
		gamma->setT0(locTAGHiter->getT(), 0.350, SYS_TAGH);

		taghGeom->E_to_counter(locTAGHiter->getE(), gamma->dCounter);
		dbeam_photons.push_back(gamma);
   }


	// Copy into factories
	factory->CopyTo(dbeam_photons);

   return NOERROR;
}

//------------------
// Extract_DMCThrown
//------------------
jerror_t DEventSourceREST::Extract_DMCThrown(hddm_r::HDDM *record,
                                   JFactory<DMCThrown> *factory)
{
   /// Copies the data from the hddm vertex records. This is called
   /// from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DMCThrown*> data;

   // loop over vertex records
   hddm_r::VertexList vertices = record->getVertices();
   hddm_r::VertexList::iterator iter;
   for (iter = vertices.begin(); iter != vertices.end(); ++iter) {
      if (iter->getJtag() != tag) {
         continue;
      }
      const hddm_r::Origin &orig = iter->getOrigin();
      double vx = orig.getVx();
      double vy = orig.getVy();
      double vz = orig.getVz();
      double vt = orig.getT();
      const hddm_r::ProductList &products = iter->getProducts();
      hddm_r::ProductList::iterator piter;
      for (piter = products.begin(); piter != products.end(); ++piter) {
         double E  = piter->getMomentum().getE();
         double px = piter->getMomentum().getPx();
         double py = piter->getMomentum().getPy();
         double pz = piter->getMomentum().getPz();
         double mass = sqrt(E*E - (px*px + py*py + pz*pz));
         if (!isfinite(mass)) {
            mass = 0.0;
         }
         DMCThrown *mcthrown = new DMCThrown;
         int pdgtype = piter->getPdgtype();
         Particle_t ptype = PDGtoPType(pdgtype);
         mcthrown->type = ptype;
         mcthrown->pdgtype = pdgtype;
         mcthrown->myid = piter->getId();
         mcthrown->parentid = piter->getParentId();
         mcthrown->mech = 0;
         mcthrown->setPID(ptype);
         mcthrown->setMass(mass);
         mcthrown->setMomentum(DVector3(px, py, pz));
         mcthrown->setPosition(DVector3(vx, vy, vz));
         mcthrown->setCharge(ParticleCharge(ptype));
         mcthrown->setT0(vt,0,SYS_NULL);
         mcthrown->setTime(vt);
         data.push_back(mcthrown);
      }
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DTOFPoint
//------------------
jerror_t DEventSourceREST::Extract_DTOFPoint(hddm_r::HDDM *record,
                                   JFactory<DTOFPoint>* factory)
{
   /// Copies the data from the tofPoint hddm record. This is called
   /// from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DTOFPoint*> data;

   // loop over tofPoint records
   const hddm_r::TofPointList &tofs = record->getTofPoints();
   hddm_r::TofPointList::iterator iter;
   for (iter = tofs.begin(); iter != tofs.end(); ++iter) {
      if (iter->getJtag() != tag) {
         continue;
      }
      DTOFPoint *tofpoint = new DTOFPoint();
      tofpoint->pos = DVector3(iter->getX(),iter->getY(),iter->getZ());
      tofpoint->t = iter->getT();
      tofpoint->dE = iter->getDE();
      tofpoint->tErr = iter->getTerr();

      //Status
		const hddm_r::TofStatusList& locTofStatusList = iter->getTofStatuses();
	   hddm_r::TofStatusList::iterator locStatusIterator = locTofStatusList.begin();
		if(locStatusIterator == locTofStatusList.end())
		{
			tofpoint->dHorizontalBar = 0;
			tofpoint->dVerticalBar = 0;
			tofpoint->dHorizontalBarStatus = 3;
			tofpoint->dVerticalBarStatus = 3;
		}
		else //should only be 1
		{
			for(; locStatusIterator != locTofStatusList.end(); ++locStatusIterator)
			{
				int locStatus = locStatusIterator->getStatus(); //horizontal_bar + 45*vertical_bar + 45*45*horizontal_status + 45*45*4*vertical_status
				tofpoint->dVerticalBarStatus = locStatus/(45*45*4);
				locStatus %= 45*45*4; //Assume compiler optimizes multiplication
				tofpoint->dHorizontalBarStatus = locStatus/(45*45);
				locStatus %= 45*45;
				tofpoint->dVerticalBar = locStatus/45;
				tofpoint->dHorizontalBar = locStatus % 45;
			}
		}

      data.push_back(tofpoint);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//------------------
// Extract_DSCHit
//------------------
jerror_t DEventSourceREST::Extract_DSCHit(hddm_r::HDDM *record,
                                   JFactory<DSCHit>* factory)
{
   /// Copies the data from the startHit hddm record. This is called
   /// from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DSCHit*> data;

   // loop over startHit records
   const hddm_r::StartHitList &starts = record->getStartHits();
   hddm_r::StartHitList::iterator iter;
   for (iter = starts.begin(); iter != starts.end(); ++iter) {
      if (iter->getJtag() != tag) {
         continue;
      }
      DSCHit *start = new DSCHit();
      start->sector = iter->getSector();
      start->dE = iter->getDE();
      start->t = iter->getT();
      data.push_back(start);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//-----------------------
// Extract_DTrigger
//-----------------------
jerror_t DEventSourceREST::Extract_DTrigger(hddm_r::HDDM *record, JFactory<DTrigger>* factory)
{
	/// Copies the data from the trigger hddm record. This is
	/// call from JEventSourceREST::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.

	if (factory==NULL)
		return OBJECT_NOT_AVAILABLE;
	string tag = (factory->Tag())? factory->Tag() : "";

	vector<DTrigger*> data;

	// loop over trigger records
	const hddm_r::TriggerList &triggers = record->getTriggers();
	hddm_r::TriggerList::iterator iter;
	for (iter = triggers.begin(); iter != triggers.end(); ++iter)
	{
		if (iter->getJtag() != tag)
			continue;

		DTrigger *locTrigger = new DTrigger();
		locTrigger->Set_L1TriggerBits(Convert_SignedIntToUnsigned(iter->getL1_trig_bits()));
		locTrigger->Set_L1FrontPanelTriggerBits(Convert_SignedIntToUnsigned(iter->getL1_fp_trig_bits()));
		data.push_back(locTrigger);
	}

	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

//-----------------------
// Extract_DFCALShower
//-----------------------
jerror_t DEventSourceREST::Extract_DFCALShower(hddm_r::HDDM *record,
                                   JFactory<DFCALShower>* factory)
{
   /// Copies the data from the fcalShower hddm record. This is
   /// call from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DFCALShower*> data;

   // loop over fcal shower records
   const hddm_r::FcalShowerList &showers =
                 record->getFcalShowers();
   hddm_r::FcalShowerList::iterator iter;
   for (iter = showers.begin(); iter != showers.end(); ++iter) {
      if (iter->getJtag() != tag)
         continue;

      DFCALShower *shower = new DFCALShower();
      shower->setPosition(DVector3(iter->getX(),iter->getY(),iter->getZ()));
      shower->setEnergy(iter->getE());
      shower->setTime(iter->getT());

      TMatrixFSym covariance(5);
	  covariance(0,0) = iter->getEerr()*iter->getEerr();
	  covariance(1,1) = iter->getXerr()*iter->getXerr();
	  covariance(2,2) = iter->getYerr()*iter->getYerr();
	  covariance(3,3) = iter->getZerr()*iter->getZerr();
	  covariance(4,4) = iter->getTerr()*iter->getTerr();
	  covariance(1,2) = covariance(2,1) = iter->getXycorr()*iter->getXerr()*iter->getYerr();
	  covariance(1,3) = covariance(3,1) = iter->getXzcorr()*iter->getXerr()*iter->getZerr();
	  covariance(2,3) = covariance(3,2) = iter->getYzcorr()*iter->getYerr()*iter->getZerr();
	  covariance(0,3) = covariance(3,0) = iter->getEzcorr()*iter->getEerr()*iter->getZerr();
	  covariance(3,4) = covariance(4,3) = iter->getTzcorr()*iter->getTerr()*iter->getZerr();

	  // further correlations (an extension of REST format, so code is different.)
	  const hddm_r::FcalCorrelationsList& locFcalCorrelationsList = iter->getFcalCorrelationses();
	  hddm_r::FcalCorrelationsList::iterator locFcalCorrelationsIterator = locFcalCorrelationsList.begin();
	  if(locFcalCorrelationsIterator != locFcalCorrelationsList.end()) {
	  	  covariance(0,4) = covariance(4,0) = locFcalCorrelationsIterator->getEtcorr()*iter->getEerr()*iter->getTerr();
	  	  covariance(0,1) = covariance(1,0) = locFcalCorrelationsIterator->getExcorr()*iter->getEerr()*iter->getXerr();
	  	  covariance(0,2) = covariance(2,0) = locFcalCorrelationsIterator->getEycorr()*iter->getEerr()*iter->getYerr();
	  	  covariance(1,4) = covariance(4,1) = locFcalCorrelationsIterator->getTxcorr()*iter->getTerr()*iter->getXerr();
	  	  covariance(2,4) = covariance(4,2) = locFcalCorrelationsIterator->getTycorr()*iter->getTerr()*iter->getYerr();
	  }
	  shower->ExyztCovariance = covariance;

      data.push_back(shower);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//-----------------------
// Extract_DBCALShower
//-----------------------
jerror_t DEventSourceREST::Extract_DBCALShower(hddm_r::HDDM *record,
                                   JFactory<DBCALShower>* factory)
{
   /// Copies the data from the bcalShower hddm record. This is
   /// call from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DBCALShower*> data;

   // loop over bcal shower records
   const hddm_r::BcalShowerList &showers =
                 record->getBcalShowers();
   hddm_r::BcalShowerList::iterator iter;
   for (iter = showers.begin(); iter != showers.end(); ++iter) {
      if (iter->getJtag() != tag)
         continue;

      DBCALShower *shower = new DBCALShower();
      shower->E = iter->getE();
      shower->E_raw = -1;
      shower->x = iter->getX();
      shower->y = iter->getY();
      shower->z = iter->getZ();
      shower->t = iter->getT();
      TMatrixFSym covariance(5);
	  covariance(0,0) = iter->getEerr()*iter->getEerr();
	  covariance(1,1) = iter->getXerr()*iter->getXerr();
	  covariance(2,2) = iter->getYerr()*iter->getYerr();
	  covariance(3,3) = iter->getZerr()*iter->getZerr();
	  covariance(4,4) = iter->getTerr()*iter->getTerr();
	  covariance(1,2) = covariance(2,1) = iter->getXycorr()*iter->getXerr()*iter->getYerr();
	  covariance(1,3) = covariance(3,1) = iter->getXzcorr()*iter->getXerr()*iter->getZerr();
	  covariance(2,3) = covariance(3,2) = iter->getYzcorr()*iter->getYerr()*iter->getZerr();
	  covariance(0,3) = covariance(3,0) = iter->getEzcorr()*iter->getEerr()*iter->getZerr();
	  covariance(3,4) = covariance(4,3) = iter->getTzcorr()*iter->getTerr()*iter->getZerr();

	  // further correlations (an extension of REST format, so code is different.)
	  const hddm_r::BcalCorrelationsList& locBcalCorrelationsList = iter->getBcalCorrelationses();
	  hddm_r::BcalCorrelationsList::iterator locBcalCorrelationsIterator = locBcalCorrelationsList.begin();
	  if(locBcalCorrelationsIterator != locBcalCorrelationsList.end()) {
		  covariance(0,4) = covariance(4,0) = locBcalCorrelationsIterator->getEtcorr()*iter->getEerr()*iter->getTerr();
		  covariance(0,1) = covariance(1,0) = locBcalCorrelationsIterator->getExcorr()*iter->getEerr()*iter->getXerr();
		  covariance(0,2) = covariance(2,0) = locBcalCorrelationsIterator->getEycorr()*iter->getEerr()*iter->getYerr();
		  covariance(1,4) = covariance(4,1) = locBcalCorrelationsIterator->getTxcorr()*iter->getTerr()*iter->getXerr();
		  covariance(2,4) = covariance(4,2) = locBcalCorrelationsIterator->getTycorr()*iter->getTerr()*iter->getYerr();
	  }
	  shower->ExyztCovariance = covariance;

		// preshower
		const hddm_r::PreshowerList& locPreShowerList = iter->getPreshowers();
		hddm_r::PreshowerList::iterator locPreShowerIterator = locPreShowerList.begin();
		if(locPreShowerIterator == locPreShowerList.end())
			shower->E_preshower = 0.0;
		else //should only be 1
		{
			for(; locPreShowerIterator != locPreShowerList.end(); ++locPreShowerIterator)
				shower->E_preshower = locPreShowerIterator->getPreshowerE();
		}

		const hddm_r::BcalClusterList& locBcalClusterList = iter->getBcalClusters();
		hddm_r::BcalClusterList::iterator locBcalClusterIterator = locBcalClusterList.begin();
		if(locBcalClusterIterator == locBcalClusterList.end())
			shower->N_cell = -1;
		else //should only be 1
		{
			for(; locBcalClusterIterator != locBcalClusterList.end(); ++locBcalClusterIterator)
				shower->N_cell = locBcalClusterIterator->getNcell();
		}

      data.push_back(shower);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//--------------------------------
// Extract_DTrackTimeBased
//--------------------------------
jerror_t DEventSourceREST::Extract_DTrackTimeBased(hddm_r::HDDM *record,
                                   JFactory<DTrackTimeBased>* factory, JEventLoop* locEventLoop)
{
   /// Copies the data from the chargedTrack hddm record. This is
   /// call from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";
   uint64_t locEventNumber = locEventLoop->GetJEvent().GetEventNumber();

   vector<DTrackTimeBased*> data;

   // Allocate enough DReferenceTrajectory objects
   // for all DTrackTimeBased objects found in this event.

   const hddm_r::ChargedTrackList &tracks = record->getChargedTracks();

   // loop over chargedTrack records
   hddm_r::ChargedTrackList::iterator iter;
   for (iter = tracks.begin(); iter != tracks.end(); ++iter) {
      if (iter->getJtag() != tag) {
         continue;
      }
      DTrackTimeBased *tra = new DTrackTimeBased();
      tra->trackid = 0;
      tra->candidateid = iter->getCandidateId();
      Particle_t ptype = iter->getPtype();
      tra->setMass(ParticleMass(ptype));
      tra->setCharge(ParticleCharge(ptype));
      tra->setPID(ptype);

      const hddm_r::TrackFit &fit = iter->getTrackFit();
      tra->Ndof = fit.getNdof();
      tra->chisq = fit.getChisq();
      tra->FOM = TMath::Prob(tra->chisq, tra->Ndof);
      tra->setT0(fit.getT0(),fit.getT0err(),(DetectorSystem_t)fit.getT0det());
      tra->setTime(fit.getT0());
      DVector3 track_pos(fit.getX0(),fit.getY0(),fit.getZ0());
      DVector3 track_mom(fit.getPx(),fit.getPy(),fit.getPz());
      tra->setPosition(track_pos);
      tra->setMomentum(track_mom);

      TMatrixFSym* loc5x5ErrorMatrix = (dynamic_cast<DApplication*>(japp))->Get_CovarianceMatrixResource(5, locEventNumber);
      (*loc5x5ErrorMatrix)(0,0) = fit.getE11();
      (*loc5x5ErrorMatrix)(0,1) = (*loc5x5ErrorMatrix)(1,0) = fit.getE12();
      (*loc5x5ErrorMatrix)(0,2) = (*loc5x5ErrorMatrix)(2,0) = fit.getE13();
      (*loc5x5ErrorMatrix)(0,3) = (*loc5x5ErrorMatrix)(3,0) = fit.getE14();
      (*loc5x5ErrorMatrix)(0,4) = (*loc5x5ErrorMatrix)(4,0) = fit.getE15();
      (*loc5x5ErrorMatrix)(1,1) = fit.getE22();
      (*loc5x5ErrorMatrix)(1,2) = (*loc5x5ErrorMatrix)(2,1) = fit.getE23();
      (*loc5x5ErrorMatrix)(1,3) = (*loc5x5ErrorMatrix)(3,1) = fit.getE24();
      (*loc5x5ErrorMatrix)(1,4) = (*loc5x5ErrorMatrix)(4,1) = fit.getE25();
      (*loc5x5ErrorMatrix)(2,2) = fit.getE33();
      (*loc5x5ErrorMatrix)(2,3) = (*loc5x5ErrorMatrix)(3,2) = fit.getE34();
      (*loc5x5ErrorMatrix)(2,4) = (*loc5x5ErrorMatrix)(4,2) = fit.getE35();
      (*loc5x5ErrorMatrix)(3,3) = fit.getE44();
      (*loc5x5ErrorMatrix)(3,4) = (*loc5x5ErrorMatrix)(4,3) = fit.getE45();
      (*loc5x5ErrorMatrix)(4,4) = fit.getE55();
      tra->setTrackingErrorMatrix(loc5x5ErrorMatrix);

      // Convert from cartesian coordinates to the 5x1 state vector corresponding to the tracking error matrix.
      double vect[5];
      vect[2]=tan(M_PI_2 - track_mom.Theta());
      vect[1]=track_mom.Phi();
      double sinphi=sin(vect[1]);
      double cosphi=cos(vect[1]);
      vect[0]=tra->charge()/track_mom.Perp();
      vect[4]=track_pos.Z();
      vect[3]=track_pos.Perp();

      if ((track_pos.X() > 0 && sinphi>0) || (track_pos.Y() <0 && cosphi>0) || (track_pos.Y() >0 && cosphi<0) || (track_pos.X() <0 && sinphi<0))
        vect[3] *= -1.; 
      tra->setTrackingStateVector(vect[0], vect[1], vect[2], vect[3], vect[4]);

      // Set the 7x7 covariance matrix.
      TMatrixFSym* loc7x7ErrorMatrix = (dynamic_cast<DApplication*>(japp))->Get_CovarianceMatrixResource(7, locEventNumber);
	  Get7x7ErrorMatrix(tra->mass(), vect, loc5x5ErrorMatrix, loc7x7ErrorMatrix);
      tra->setErrorMatrix(loc7x7ErrorMatrix);

		// Hit layers
      const hddm_r::HitlayersList& locHitlayersList = iter->getHitlayerses();
	   hddm_r::HitlayersList::iterator locHitlayersIterator = locHitlayersList.begin();
		if(locHitlayersIterator == locHitlayersList.end())
		{
			tra->dCDCRings = 0;
			tra->dFDCPlanes = 0;
		}
		else //should only be 1
		{
			for(; locHitlayersIterator != locHitlayersList.end(); ++locHitlayersIterator)
			{
				tra->dCDCRings = locHitlayersIterator->getCDCrings();
				tra->dFDCPlanes = locHitlayersIterator->getFDCplanes();
			}
		}

		// MC match hit info
      const hddm_r::McmatchList& locMCMatchesList = iter->getMcmatchs();
	   hddm_r::McmatchList::iterator locMcmatchIterator = locMCMatchesList.begin();
		if(locMcmatchIterator == locMCMatchesList.end())
		{
			tra->dCDCRings = 0;
			tra->dFDCPlanes = 0;
		}
		else //should only be 1
		{
			for(; locMcmatchIterator != locMCMatchesList.end(); ++locMcmatchIterator)
			{
				tra->dMCThrownMatchMyID = locMcmatchIterator->getIthrown();
				tra->dNumHitsMatchedToThrown = locMcmatchIterator->getNumhitsmatch();
			}
		}

      // add the drift chamber dE/dx information
      const hddm_r::DEdxDCList &el = iter->getDEdxDCs();
      hddm_r::DEdxDCList::iterator diter = el.begin();
      if (diter != el.end()) {
         tra->dNumHitsUsedFordEdx_FDC = diter->getNsampleFDC();
         tra->dNumHitsUsedFordEdx_CDC = diter->getNsampleCDC();
         tra->ddEdx_FDC = diter->getDEdxFDC();
         tra->ddEdx_CDC = diter->getDEdxCDC();
         tra->ddx_FDC = diter->getDxFDC();
         tra->ddx_CDC = diter->getDxCDC();
         tra->setdEdx((tra->dNumHitsUsedFordEdx_CDC >= tra->dNumHitsUsedFordEdx_FDC) ? tra->ddEdx_CDC : tra->ddEdx_FDC);
      }
      else {
         tra->dNumHitsUsedFordEdx_FDC = 0;
         tra->dNumHitsUsedFordEdx_CDC = 0;
         tra->ddEdx_FDC = 0.0;
         tra->ddEdx_CDC = 0.0;
         tra->ddx_FDC = 0.0;
         tra->ddx_CDC = 0.0;
         tra->setdEdx(0.0);
      }

      data.push_back(tra);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//--------------------------------
// Extract_DDetectorMatches
//--------------------------------
jerror_t DEventSourceREST::Extract_DDetectorMatches(JEventLoop* locEventLoop, hddm_r::HDDM *record, JFactory<DDetectorMatches>* factory)
{
   /// Copies the data from the detectorMatches hddm record. This is
   /// called from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if(factory==NULL)
     return OBJECT_NOT_AVAILABLE;

   string tag = (factory->Tag())? factory->Tag() : "";
   vector<DDetectorMatches*> data;

   vector<const DTrackTimeBased*> locTrackTimeBasedVector;
   locEventLoop->Get(locTrackTimeBasedVector);

   vector<const DSCHit*> locSCHits;
   locEventLoop->Get(locSCHits);

   vector<const DTOFPoint*> locTOFPoints;
   locEventLoop->Get(locTOFPoints);

   vector<const DBCALShower*> locBCALShowers;
   locEventLoop->Get(locBCALShowers);

   vector<const DFCALShower*> locFCALShowers;
   locEventLoop->Get(locFCALShowers);

   const hddm_r::DetectorMatchesList &detectormatches = record->getDetectorMatcheses();

   // loop over chargedTrack records
   hddm_r::DetectorMatchesList::iterator iter;
   for(iter = detectormatches.begin(); iter != detectormatches.end(); ++iter)
	{
      if(iter->getJtag() != tag)
         continue;

      DDetectorMatches *locDetectorMatches = new DDetectorMatches();

      const hddm_r::BcalMatchParamsList &bcalList = iter->getBcalMatchParamses();
      hddm_r::BcalMatchParamsList::iterator bcalIter = bcalList.begin();
      for(; bcalIter != bcalList.end(); ++bcalIter)
      {
         size_t locShowerIndex = bcalIter->getShower();
         size_t locTrackIndex = bcalIter->getTrack();

         DBCALShowerMatchParams locShowerMatchParams;
         locShowerMatchParams.dTrack = locTrackTimeBasedVector[locTrackIndex];
         locShowerMatchParams.dBCALShower = locBCALShowers[locShowerIndex];
         locShowerMatchParams.dx = bcalIter->getDx();
         locShowerMatchParams.dFlightTime = bcalIter->getTflight();
         locShowerMatchParams.dFlightTimeVariance = bcalIter->getTflightvar();
         locShowerMatchParams.dPathLength = bcalIter->getPathlength();
         locShowerMatchParams.dDeltaPhiToShower = bcalIter->getDeltaphi();
         locShowerMatchParams.dDeltaZToShower = bcalIter->getDeltaz();

         locDetectorMatches->Add_Match(locTrackTimeBasedVector[locTrackIndex], locBCALShowers[locShowerIndex], locShowerMatchParams);
      }

      const hddm_r::FcalMatchParamsList &fcalList = iter->getFcalMatchParamses();
      hddm_r::FcalMatchParamsList::iterator fcalIter = fcalList.begin();
      for(; fcalIter != fcalList.end(); ++fcalIter)
      {
         size_t locShowerIndex = fcalIter->getShower();
         size_t locTrackIndex = fcalIter->getTrack();

         DFCALShowerMatchParams locShowerMatchParams;
         locShowerMatchParams.dTrack = locTrackTimeBasedVector[locTrackIndex];
         locShowerMatchParams.dFCALShower = locFCALShowers[locShowerIndex];

         locShowerMatchParams.dx = fcalIter->getDx();
         locShowerMatchParams.dFlightTime = fcalIter->getTflight();
         locShowerMatchParams.dFlightTimeVariance = fcalIter->getTflightvar();
         locShowerMatchParams.dPathLength = fcalIter->getPathlength();
         locShowerMatchParams.dDOCAToShower = fcalIter->getDoca();

         locDetectorMatches->Add_Match(locTrackTimeBasedVector[locTrackIndex], locFCALShowers[locShowerIndex], locShowerMatchParams);
      }

      const hddm_r::ScMatchParamsList &scList = iter->getScMatchParamses();
      hddm_r::ScMatchParamsList::iterator scIter = scList.begin();
      for(; scIter != scList.end(); ++scIter)
      {
         size_t locHitIndex = scIter->getHit();
         size_t locTrackIndex = scIter->getTrack();

         DSCHitMatchParams locSCHitMatchParams;
         locSCHitMatchParams.dTrack = locTrackTimeBasedVector[locTrackIndex];
         locSCHitMatchParams.dSCHit = locSCHits[locHitIndex];

         locSCHitMatchParams.dEdx = scIter->getDEdx();
         locSCHitMatchParams.dHitTime = scIter->getThit();
         locSCHitMatchParams.dHitTimeVariance = scIter->getThitvar();
         locSCHitMatchParams.dHitEnergy = scIter->getEhit();
         locSCHitMatchParams.dFlightTime = scIter->getTflight();
         locSCHitMatchParams.dFlightTimeVariance = scIter->getTflightvar();
         locSCHitMatchParams.dPathLength = scIter->getPathlength();
         locSCHitMatchParams.dDeltaPhiToHit = scIter->getDeltaphi();

         locDetectorMatches->Add_Match(locTrackTimeBasedVector[locTrackIndex], locSCHits[locHitIndex], locSCHitMatchParams);
      }

      const hddm_r::TofMatchParamsList &tofList = iter->getTofMatchParamses();
      hddm_r::TofMatchParamsList::iterator tofIter = tofList.begin();
      for(; tofIter != tofList.end(); ++tofIter)
      {
         size_t locHitIndex = tofIter->getHit();
         size_t locTrackIndex = tofIter->getTrack();

         DTOFHitMatchParams locTOFHitMatchParams;
         locTOFHitMatchParams.dTrack = locTrackTimeBasedVector[locTrackIndex];
         locTOFHitMatchParams.dTOFPoint = locTOFPoints[locHitIndex];

         locTOFHitMatchParams.dHitTime = tofIter->getThit();
         locTOFHitMatchParams.dHitTimeVariance = tofIter->getThitvar();
         locTOFHitMatchParams.dHitEnergy = tofIter->getEhit();

         locTOFHitMatchParams.dEdx = tofIter->getDEdx();
         locTOFHitMatchParams.dFlightTime = tofIter->getTflight();
         locTOFHitMatchParams.dFlightTimeVariance = tofIter->getTflightvar();
         locTOFHitMatchParams.dPathLength = tofIter->getPathlength();
         locTOFHitMatchParams.dDeltaXToHit = tofIter->getDeltax();
         locTOFHitMatchParams.dDeltaYToHit = tofIter->getDeltay();

         locDetectorMatches->Add_Match(locTrackTimeBasedVector[locTrackIndex], locTOFPoints[locHitIndex], locTOFHitMatchParams);
      }

      const hddm_r::BcalDOCAtoTrackList &bcaldocaList = iter->getBcalDOCAtoTracks();
      hddm_r::BcalDOCAtoTrackList::iterator bcaldocaIter = bcaldocaList.begin();
      for(; bcaldocaIter != bcaldocaList.end(); ++bcaldocaIter)
      {
         size_t locShowerIndex = bcaldocaIter->getShower();
         double locDeltaPhi = bcaldocaIter->getDeltaphi();
         double locDeltaZ = bcaldocaIter->getDeltaz();
         locDetectorMatches->Set_DistanceToNearestTrack(locBCALShowers[locShowerIndex], locDeltaPhi, locDeltaZ);
      }

      const hddm_r::FcalDOCAtoTrackList &fcaldocaList = iter->getFcalDOCAtoTracks();
      hddm_r::FcalDOCAtoTrackList::iterator fcaldocaIter = fcaldocaList.begin();
      for(; fcaldocaIter != fcaldocaList.end(); ++fcaldocaIter)
      {
         size_t locShowerIndex = fcaldocaIter->getShower();
         double locDOCA = fcaldocaIter->getDoca();
         locDetectorMatches->Set_DistanceToNearestTrack(locFCALShowers[locShowerIndex], locDOCA);
      }

      const hddm_r::TflightPCorrelationList &correlationList = iter->getTflightPCorrelations();
      hddm_r::TflightPCorrelationList::iterator correlationIter = correlationList.begin();
      for(; correlationIter != correlationList.end(); ++correlationIter)
      {
         size_t locTrackIndex = correlationIter->getTrack();
         DetectorSystem_t locDetectorSystem = (DetectorSystem_t)correlationIter->getSystem();
         double locCorrelation = correlationIter->getCorrelation();
         locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBasedVector[locTrackIndex], locDetectorSystem, locCorrelation);
      }

      data.push_back(locDetectorMatches);
   }

   // Copy data to factory
   factory->CopyTo(data);

   return NOERROR;
}

// Transform the 5x5 tracking error matrix into a 7x7 error matrix in cartesian
// coordinates.
// This was copied and transformed from DKinFit.cc
void DEventSourceREST::Get7x7ErrorMatrix(double mass, const double vec[5], const TMatrixFSym* C5x5, TMatrixFSym* loc7x7ErrorMatrix)
{
  DMatrix J(7,5);

  // State vector
  double q_over_pt=vec[0];
  double phi=vec[1];
  double tanl=vec[2];
  double D=vec[3];

  double pt=1./fabs(q_over_pt);
  double pt_sq=pt*pt;
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double q=(q_over_pt>0)?1.:-1.;

  J(0, 0)=-q*pt_sq*cosphi;
  J(0, 1)=-pt*sinphi;
  
  J(1, 0)=-q*pt_sq*sinphi;
  J(1, 1)=pt*cosphi;
  
  J(2, 0)=-q*pt_sq*tanl;
  J(2, 2)=pt;
  
  J(3, 1)=-D*cosphi;
  J(3, 3)=-sinphi;
  
  J(4, 1)=-D*sinphi;
  J(4, 3)=cosphi;
  
  J(5, 4)=1.;

  // C'= JCJ^T
  TMatrixFSym locTempMatrix = *C5x5;
  loc7x7ErrorMatrix=locTempMatrix.Similarity(J);
}

uint32_t DEventSourceREST::Convert_SignedIntToUnsigned(int32_t locSignedInt) const
{
	//Convert uint32_t to int32_t
	//Reverse scheme (from DEventWriterREST):
		//If is >= 0, then the int32_t is the same as the uint32_t (last bit not set)
		//If is the minimum int: bit 32 is 1, and all of the other bits are zero
		//Else, bit 32 is 1, then the uint32_t is -1 * int32_t, + add the last bit
	if(locSignedInt >= 0)
		return uint32_t(locSignedInt); //bit 32 is zero
	else if(locSignedInt == numeric_limits<int32_t>::min()) // -(2^31)
		return uint32_t(0x80000000); //bit 32 is 1, all others are 0
	return uint32_t(-1*locSignedInt) + uint32_t(0x80000000); //bit 32 is 1, all others are negative of signed int (which was negative)
}

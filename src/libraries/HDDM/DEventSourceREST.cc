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

#include <DVector2.h>
#include <DEventSourceREST.h>


//----------------
// Constructor
//----------------
DEventSourceREST::DEventSourceREST(const char* source_name)
 : JEventSource(source_name), dapp(0), 
   saved_bfield(0), saved_geom(0), 
   saved_runnumber(-1)
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
   pthread_mutex_init(&rt_mutex, NULL);
}

//----------------
// Destructor
//----------------
DEventSourceREST::~DEventSourceREST()
{
  // Check for DReferenceTrajectory objects we need to delete
  pthread_mutex_lock(&rt_mutex);
  while (rt_pool.size()>0){
    std::list<DReferenceTrajectory*>::iterator it=rt_pool.begin();
    delete *it;
    rt_pool.pop_front();
  }
  rt_pool.clear();
  pthread_mutex_unlock(&rt_mutex);
  
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

      // Check for DReferenceTrajectory objects we need to delete
      pthread_mutex_lock(&rt_mutex);
      while (rt_pool.size()>0){
	std::list<DReferenceTrajectory*>::iterator it=rt_pool.begin();
	delete *it;
	rt_pool.pop_front();
      }
      rt_pool.clear();
      pthread_mutex_unlock(&rt_mutex);

      return NO_MORE_EVENTS_IN_SOURCE;
   }

   hddm_r::HDDM *record = new hddm_r::HDDM();
   *fin >> *record;


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
         *fin >> *record;
         continue;
      }
      event.SetEventNumber(re.getEventNo());
      event.SetRunNumber(re.getRunNo());
      event.SetJEventSource(this);
      event.SetRef(record);
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

   // Check for DReferenceTrajectory objects we need to delete
   pthread_mutex_lock(&rt_mutex);
   std::map<hddm_r::HDDM*, std::vector<DReferenceTrajectory*> >::iterator
                           iter = rt_by_event.find(record);
   if (iter != rt_by_event.end()) {
     std::vector<DReferenceTrajectory*> &rtvector = iter->second;
     unsigned int i=0;
     while (i<rtvector.size()&&rt_pool.size()<max_rt_pool_size){
       rt_pool.push_back(rtvector[i]);
       i++;
     }
     /*
      for (unsigned int i=0; i < rtvector.size(); ++i) {
      rt_pool.push_back(rtvector[i]);
      }
     */
     for(; i < rtvector.size(); ++i) {
       delete rtvector[i];
     }
     
     rt_by_event.erase(iter);
   }
   pthread_mutex_unlock(&rt_mutex);
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

   if (dapp == 0) {
      JEventLoop *loop = event.GetJEventLoop();
      dapp = dynamic_cast<DApplication*>(loop->GetJApplication());    
      if (!dapp) {
         throw RESOURCE_UNAVAILABLE;
      }
      map<pthread_t,double> rates_by_thread;
      dapp->GetInstantaneousThreadRates(rates_by_thread);
      max_rt_pool_size=13*rates_by_thread.size();
   }

   string dataClassName = factory->GetDataClassName();
	
   if (dataClassName =="DMCReaction") {
      return Extract_DMCReaction(record,
                     dynamic_cast<JFactory<DMCReaction>*>(factory));
   }
   if (dataClassName =="DBeamPhoton") {
      return Extract_DBeamPhoton(record,
                     dynamic_cast<JFactory<DBeamPhoton>*>(factory));
   }
   if (dataClassName =="DMCThrown") {
      return Extract_DMCThrown(record,
                     dynamic_cast<JFactory<DMCThrown>*>(factory));
   }
   if (dataClassName =="DTagger") {
      return Extract_DTagger(record,
                     dynamic_cast<JFactory<DTagger>*>(factory));
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
                     dynamic_cast<JFactory<DTrackTimeBased>*>(factory));
   }
   if (dataClassName =="DMCTrigger") {
      return Extract_DMCTrigger(record,
                     dynamic_cast<JFactory<DMCTrigger>*>(factory));
   }
#if 0
   // one day we will need a class to hold RF timing information
   if (dataClassName =="DRFTime") {
      return Extract_DRFTime(record,
                     dynamic_cast<JFactory<DRFTime>*>(factory));
   }
#endif

   return OBJECT_NOT_AVAILABLE;
}

//------------------
// Extract_DMCReaction
//------------------
jerror_t DEventSourceREST::Extract_DMCReaction(hddm_r::HDDM *record,
                                   JFactory<DMCReaction> *factory)
{
   /// Copies the data from the Reaction hddm class. This is called
   /// from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.
	
   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   std::string tag = (factory->Tag())? factory->Tag() : "";

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
      mcreaction->beam.setPosition(DVector3(0.0, 0.0, 65.0));
      mcreaction->beam.setMomentum(DVector3(0.0, 0.0, Ebeam));
      mcreaction->beam.setMass(0.0);
      mcreaction->beam.setCharge(0.0);
      mcreaction->beam.clearErrorMatrix();
      mcreaction->beam.setT0(0.0, 0.0, SYS_NULL);
      mcreaction->beam.setT1(0.0, 0.0, SYS_NULL);
      mcreaction->beam.setTime(0.0);
      mcreaction->beam.setPID(Gamma);
      mcreaction->target.setPosition(DVector3(0.0, 0.0, 65.0));
      mcreaction->target.setMomentum(DVector3(0.0, 0.0, 0.0));
      Particle_t ttype = iter->getTargetType();
      mcreaction->target.setPID((Particle_t)ttype);
      mcreaction->target.setMass(ParticleMass(ttype));
      mcreaction->target.setCharge(ParticleCharge(ttype));
      mcreaction->target.clearErrorMatrix();
      mcreaction->target.setT0(0.0, 0.0, SYS_NULL);
      mcreaction->target.setT1(0.0, 0.0, SYS_NULL);
      mcreaction->target.setTime(0.0);
   }
	
   // Copy into factories
   factory->CopyTo(dmcreactions);

   return NOERROR;
}

//------------------
// Extract_DBeamPhoton
//------------------
jerror_t DEventSourceREST::Extract_DBeamPhoton(hddm_r::HDDM *record,
                                   JFactory<DBeamPhoton> *factory)
{
   /// Copies the data from the Reaction hddm class. This is called
   /// from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DBeamPhoton*> dbeam_photons;

   // loop over reaction records
   const hddm_r::ReactionList &reactions = record->getReactions();
   hddm_r::ReactionList::iterator iter;
   for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      if (iter->getJtag() != tag) {
         continue;
      }
      DBeamPhoton *beamphoton = new DBeamPhoton;
      double Ebeam = iter->getEbeam();
      beamphoton->setPID(Gamma);
      beamphoton->setPosition(DVector3(0.0, 0.0, 65.0));
      beamphoton->setMomentum(DVector3(0.0, 0.0, Ebeam));
      beamphoton->setMass(0.0);
      beamphoton->setCharge(0.0);
      beamphoton->clearErrorMatrix();
      beamphoton->setT0(0.0, 0.0, SYS_NULL);
//      double zint = iter->getVertex(0).getOrigin().getVz();
//      beamphoton->setTime((zint-65.0)/SPEED_OF_LIGHT);
      beamphoton->setTime(0.0); //0 because position is defined at 0, 0, 65.  Would be non-zero if position wasn't forced...
      dbeam_photons.push_back(beamphoton);
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
         if (!finite(mass)) {
            mass = 0.0;
         }
         DMCThrown *mcthrown = new DMCThrown;
         int pdgtype = piter->getPdgtype();
         Particle_t ptype = PDGtoPtype(pdgtype);
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
// Extract_DTagger
//------------------
jerror_t DEventSourceREST::Extract_DTagger(hddm_r::HDDM *record,
                                   JFactory<DTagger>* factory)
{
   /// Copies the data from the taggerHit hddm record. This is called
   /// from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DTagger*> data;

   // loop over taggerHit records
   const hddm_r::TaggerHitList &tags = record->getTaggerHits();
   hddm_r::TaggerHitList::iterator iter;
   for (iter = tags.begin(); iter != tags.end(); ++iter) {
      if (iter->getJtag() != tag) {
         continue;
      }
      DTagger *tagger = new DTagger();
      tagger->E = iter->getE();
      tagger->t = iter->getT();
      tagger->row = 0;
      tagger->column = 0;
      data.push_back(tagger);
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
// Extract_DFCALShower
//-----------------------
jerror_t DEventSourceREST::Extract_DFCALShower(hddm_r::HDDM *record,
                                   JFactory<DFCALShower>* factory)
{
   /// Copies the data from the calorimeterCluster hddm record. This is
   /// call from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DFCALShower*> data;

   // loop over calorimeterCluster records
   const hddm_r::CalorimeterClusterList &clusters =
                 record->getCalorimeterClusters();
   hddm_r::CalorimeterClusterList::iterator iter;
   for (iter = clusters.begin(); iter != clusters.end(); ++iter) {
      if (iter->getJtag() != tag || iter->getZ() < 600) {
         continue;
      }
      DFCALShower *shower = new DFCALShower();
      shower->setPosition(DVector3(iter->getX(),iter->getY(),iter->getZ()));
      shower->setPosError(iter->getXerr(),iter->getYerr(),iter->getZerr());
      shower->setEnergy(iter->getE());
      shower->setTime(iter->getT());
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
   /// Copies the data from the calorimeterCluster hddm record. This is
   /// call from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DBCALShower*> data;

   // loop over calorimeterCluster records
   const hddm_r::CalorimeterClusterList &clusters =
                 record->getCalorimeterClusters();
   hddm_r::CalorimeterClusterList::iterator iter;
   for (iter = clusters.begin(); iter != clusters.end(); ++iter) {
      if (iter->getJtag() != tag || iter->getZ() > 600) {
         continue;
      }
      DBCALShower *shower = new DBCALShower();
      shower->E = iter->getE();
      shower->E_raw = -1;
      shower->x = iter->getX();
      shower->y = iter->getY();
      shower->z = iter->getZ();
      shower->t = iter->getT();
      shower->xErr = iter->getXerr();
      shower->yErr = iter->getYerr();
      shower->zErr = iter->getZerr();
      shower->tErr = iter->getTerr();
      shower->N_cell = -1;
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
                                   JFactory<DTrackTimeBased>* factory)
{
   /// Copies the data from the chargedTrack hddm record. This is
   /// call from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }
   string tag = (factory->Tag())? factory->Tag() : "";

   vector<DTrackTimeBased*> data;

   // Allocate enough DReferenceTrajectory objects
   // for all DTrackTimeBased objects found in this event.

   const hddm_r::ChargedTrackList &tracks = record->getChargedTracks();
   std::vector<DReferenceTrajectory*> free_rts;
   std::vector<DReferenceTrajectory*> used_rts;
   pthread_mutex_lock(&rt_mutex);
   while ((int)free_rts.size() < tracks.size()) {
      if (rt_pool.size() > 0) {
         free_rts.push_back(rt_pool.back());
         rt_pool.pop_back();
      }
      else {
         if (saved_bfield == 0) {
            saved_bfield = dapp->GetBfield();
         }
         free_rts.push_back(new DReferenceTrajectory(saved_bfield));
      }
   }
   pthread_mutex_unlock(&rt_mutex);

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
      DMatrixDSym mat(5);
      mat(0,0) = fit.getE11();
      mat(0,1) = mat(1,0) = fit.getE12();
      mat(0,2) = mat(2,0) = fit.getE13();
      mat(0,3) = mat(3,0) = fit.getE14();
      mat(0,4) = mat(4,0) = fit.getE15();
      mat(1,1) = fit.getE22();
      mat(1,2) = mat(2,1) = fit.getE23();
      mat(1,3) = mat(3,1) = fit.getE24();
      mat(1,4) = mat(4,1) = fit.getE25();
      mat(2,2) = fit.getE33();
      mat(2,3) = mat(3,2) = fit.getE34();
      mat(2,4) = mat(4,2) = fit.getE35();
      mat(3,3) = fit.getE44();
      mat(3,4) = mat(4,3) = fit.getE45();
      mat(4,4) = fit.getE55();
      tra->setTrackingErrorMatrix(mat);

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
		tra->setErrorMatrix(Get7x7ErrorMatrix(tra->mass(), vect, mat));

      // configure the DReferenceTrajectory object
      DReferenceTrajectory *rt = free_rts.back();
      free_rts.pop_back();
      rt->Reset();
      rt->SetMass(tra->mass());
      int run = record->getReconstructedPhysicsEvent().getRunNo();
      if (run != saved_runnumber || saved_geom == 0) {
         saved_geom = dapp->GetDGeometry(run);
         saved_runnumber = run;
      }
      rt->SetDGeometry(saved_geom);
      rt->Swim(track_pos,track_mom,tra->charge());
      used_rts.push_back(rt);
      tra->rt = rt;

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

   // Post the DReferenceTrajectory objects that were allocated above
   // to a global list so that they can be reused when we are done with them.
   pthread_mutex_lock(&rt_mutex);
   rt_by_event[record] = used_rts;
   pthread_mutex_unlock(&rt_mutex);
   
   return NOERROR;
}

//--------------------------------
// Extract_DMCTrigger
//--------------------------------
jerror_t DEventSourceREST::Extract_DMCTrigger(hddm_r::HDDM *record,
                                   JFactory<DMCTrigger>* factory)
{
	/// Copies the data from the trigger hddm record. This is
	/// called from JEventSourceREST::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.

	if (factory==NULL) {
	  return OBJECT_NOT_AVAILABLE;
	}
	string tag = (factory->Tag())? factory->Tag() : "";

	vector<DMCTrigger*> data;

	const hddm_r::TriggerList &triggers = record->getTriggers();

	// loop over chargedTrack records
	hddm_r::TriggerList::iterator iter;
	for (iter = triggers.begin(); iter != triggers.end(); ++iter) {
		if (iter->getJtag() != tag) {
			continue;
		}
      DMCTrigger *trigger = new DMCTrigger();
      trigger->L1a_fired = iter->getL1a();
      trigger->L1b_fired = iter->getL1b();
      trigger->L1c_fired = iter->getL1c();
	}

	// Copy data to factory
	factory->CopyTo(data);

	return NOERROR;
}

Particle_t DEventSourceREST::PDGtoPtype(int pdgtype)
{
   for (int ptype=0; ptype < 99; ++ptype) {
      if (PDGtype((Particle_t)ptype) == pdgtype) {
         return (Particle_t)ptype;
      }
   }
   return (Particle_t)0;
}

// Transform the 5x5 tracking error matrix into a 7x7 error matrix in cartesian
// coordinates.
// This was copied and transformed from DKinFit.cc
DMatrixDSym DEventSourceREST::Get7x7ErrorMatrix(double mass, const double vec[5], const DMatrixDSym& C5x5)
{
  DMatrixDSym C7x7(7);
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
  DMatrixDSym locTempMatrix = C5x5;
  C7x7=locTempMatrix.Similarity(J);
  
  return C7x7;
}


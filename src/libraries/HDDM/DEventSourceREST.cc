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

#include "BCAL/DBCALShower.h"
#include "FCAL/DFCALShower.h"

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
      return NO_MORE_EVENTS_IN_SOURCE;
   }

   hddm_r::HDDM *record = new hddm_r::HDDM();
   *fin >> *record;
   ++Nevents_read;

   // Copy the reference info into the JEvent object
   hddm_r::ReconstructedPhysicsEvent &re
         = record->getReconstructedPhysicsEvent();
   event.SetEventNumber(re.getEventNo());
   event.SetRunNumber(re.getRunNo());
   event.SetJEventSource(this);
   event.SetRef(record);
 
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

   // HDDM doesn't exactly support tagged factories, but the tag
   // can be used to direct filling of the correct factory.
   string tag = (factory->Tag())? factory->Tag() : "";

   // Get name of data class we're trying to extract
   string dataClassName = factory->GetDataClassName();
	
   if (dataClassName =="DMCReaction" && tag=="") {
      return Extract_DMCReaction(record,
                     dynamic_cast<JFactory<DMCReaction>*>(factory));
   }
   if (dataClassName =="DBeamPhoton" && tag=="") {
      return Extract_DBeamPhoton(record,
                     dynamic_cast<JFactory<DBeamPhoton>*>(factory));
   }
   if (dataClassName =="DMCThrown" && tag=="") {
      return Extract_DMCThrown(record,
                     dynamic_cast<JFactory<DMCThrown>*>(factory));
   }
   if (dataClassName =="DTagger" && tag=="") {
      return Extract_DTagger(record,
                     dynamic_cast<JFactory<DTagger>*>(factory));
   }
   if (dataClassName =="DFCALShower" && tag=="") {
      return Extract_DFCALShower(record,
                     dynamic_cast<JFactory<DFCALShower>*>(factory));
   }
   if (dataClassName =="DBCALShower" && tag=="") {
      return Extract_DBCALShower(record,
                     dynamic_cast<JFactory<DBCALShower>*>(factory));
   }
   if (dataClassName =="DNeutralShower" && tag=="") {
      return Extract_DNeutralShower(record,
                     dynamic_cast<JFactory<DNeutralShower>*>(factory));
   }
   if (dataClassName =="DChargedTrackHypothesis" && tag=="") {
      return Extract_DChargedTrackHypothesis(record,
                     dynamic_cast<JFactory<DChargedTrackHypothesis>*>(factory));
   }

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

   vector<DMCReaction*> dmcreactions;

   // loop over reaction records
   const hddm_r::ReactionList &reactions = record->getReactions();
   hddm_r::ReactionList::iterator iter;
   for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
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
      mcreaction->target.setPosition(DVector3(0.0, 0.0, 65.0));
      mcreaction->target.setMomentum(DVector3(0.0, 0.0, 0.0));
      mcreaction->target.setMass(ParticleMass(Proton));
      mcreaction->target.setCharge(ParticleCharge(Proton));
      mcreaction->target.clearErrorMatrix();
      mcreaction->target.setT0(0.0, 0.0, SYS_NULL);
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

   vector<DBeamPhoton*> dbeam_photons;

   // loop over reaction records
   const hddm_r::ReactionList &reactions = record->getReactions();
   hddm_r::ReactionList::iterator iter;
   for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      DBeamPhoton *beamphoton = new DBeamPhoton;
      double Ebeam = iter->getEbeam();
      beamphoton->setPosition(DVector3(0.0, 0.0, 65.0));
      beamphoton->setMomentum(DVector3(0.0, 0.0, Ebeam));
      beamphoton->setMass(0.0);
      beamphoton->setCharge(0.0);
      beamphoton->clearErrorMatrix();
      beamphoton->setT0(0.0, 0.0, SYS_NULL);
      double zint = iter->getVertex(0).getOrigin().getVz();
      beamphoton->t = (zint-65.0)/SPEED_OF_LIGHT;
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

   vector<DMCThrown*> data;

   // loop over vertex records
   hddm_r::VertexList vertices = record->getVertices();
   hddm_r::VertexList::iterator iter;
   for (iter = vertices.begin(); iter != vertices.end(); ++iter) {
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
         mcthrown->setMass(mass);
         mcthrown->setMomentum(DVector3(px, py, pz));
         mcthrown->setPosition(DVector3(vx, vy, vz));
         mcthrown->setCharge(ParticleCharge(ptype));
         mcthrown->setT0(vt,0,SYS_NULL);
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
   /// Copies the data from the tagger hddm record. This is called
   /// from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }

   vector<DTagger*> data;

   // loop over taggerHit records
   const hddm_r::TaggerHitList &tags = record->getTaggerHits();
   hddm_r::TaggerHitList::iterator iter;
   for (iter = tags.begin(); iter != tags.end(); ++iter) {
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

   vector<DFCALShower*> data;

   // loop over calorimeterCluster records
   const hddm_r::CalorimeterClusterList &clusters =
                 record->getCalorimeterClusters();
   hddm_r::CalorimeterClusterList::iterator iter;
   for (iter = clusters.begin(); iter != clusters.end(); ++iter) {
      if (iter->getZ() < 600) {
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

   vector<DBCALShower*> data;

   // loop over calorimeterCluster records
   const hddm_r::CalorimeterClusterList &clusters =
                 record->getCalorimeterClusters();
   hddm_r::CalorimeterClusterList::iterator iter;
   for (iter = clusters.begin(); iter != clusters.end(); ++iter) {
      if (iter->getZ() > 600) {
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

//-----------------------
// Extract_DNeutralShower
//-----------------------
jerror_t DEventSourceREST::Extract_DNeutralShower(hddm_r::HDDM *record,
                                   JFactory<DNeutralShower>* factory)
{
   /// Copies the data from the calorimeterCluster hddm record. This is
   /// call from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }

   vector<DNeutralShower*> data;

   // loop over calorimeterCluster records
   const hddm_r::CalorimeterClusterList &clusters =
                 record->getCalorimeterClusters();
   hddm_r::CalorimeterClusterList::iterator iter;
   for (iter = clusters.begin(); iter != clusters.end(); ++iter) {
      if (!iter->getIsNeutral()) {
         continue;
      }
      DNeutralShower *shower;
      double x = iter->getX();
      double y = iter->getY();
      double z = iter->getZ();
      double t = iter->getT();
      double E = iter->getE();
      if (z > 600) {
         DFCALShower s;
         shower = new DNeutralShower(&s);
      }
      else {
         DBCALShower s;
         shower = new DNeutralShower(&s);
      }
      double xerr = iter->getXerr();
      double yerr = iter->getYerr();
      double zerr = iter->getZerr();
      double terr = iter->getTerr();
      double Eerr = iter->getEerr();
      shower->dSpacetimeVertex.SetXYZT(x,y,z,t);
      shower->dSpacetimeVertexUncertainties.SetXYZT(xerr,yerr,zerr,terr);
      shower->dEnergy = E;
      shower->dEnergyUncertainty = Eerr;
      data.push_back(shower);
   }

   // Copy into factory
   factory->CopyTo(data);

   return NOERROR;
}

//--------------------------------
// Extract_DChargedTrackHypothesis
//--------------------------------
jerror_t DEventSourceREST::Extract_DChargedTrackHypothesis(hddm_r::HDDM *record,
                                   JFactory<DChargedTrackHypothesis>* factory)
{
   /// Copies the data from the chargedTrack hddm record. This is
   /// call from JEventSourceREST::GetObjects. If factory is NULL, this
   /// returns OBJECT_NOT_AVAILABLE immediately.

   if (factory==NULL) {
      return OBJECT_NOT_AVAILABLE;
   }

   vector<DChargedTrackHypothesis*> data;

   // loop over chargedTrack records
   const hddm_r::ChargedTrackList &tracks = record->getChargedTracks();
   hddm_r::ChargedTrackList::iterator iter;
   for (iter = tracks.begin(); iter != tracks.end(); ++iter) {
      DChargedTrackHypothesis *hyp = new DChargedTrackHypothesis();
      hyp->candidateid = iter->getCandidateId();
      Particle_t ptype = iter->getPtype();
      hyp->dPID = ptype;

      const hddm_r::TrackFit &fit = iter->getTrackFit();
      hyp->dNDF_Track = fit.getNdof();
      hyp->dChiSq_Track = fit.getChisq();
      double pt0 = fit.getT0();
      double pt0err = fit.getT0err();
      hyp->dProjectedStartTime = pt0;
      hyp->dProjectedStartTimeUncertainty = pt0err;
      double x0 = fit.getX0();
      double y0 = fit.getY0();
      double z0 = fit.getZ0();
      hyp->setPosition(DVector3(x0,y0,z0));
      double px = fit.getPx();
      double py = fit.getPy();
      double pz = fit.getPz();
      hyp->setMomentum(DVector3(px,py,pz));
      hyp->setMass(ParticleMass(ptype));
      hyp->setCharge(ParticleCharge(ptype));
      hyp->dNDF_Timing = fit.getNdof_RFtime();
      hyp->dChiSq_Timing = fit.getChisq_RFtime();
      double tflight = fit.getFlightTime();
      double tferror = fit.getFlightTimeErr();
      double pathlen = fit.getPathLength();
      double patherr = fit.getPathLengthErr();
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
      hyp->setTrackingErrorMatrix(mat);

      hyp->dNDF_DCdEdx = 0;
      const hddm_r::DCdEdxList &el = iter->getDCdEdxs();
      hddm_r::DCdEdxList::iterator diter;
      for (diter = el.begin(); diter != el.end(); ++diter) {
         hyp->dNDF_DCdEdx = diter->getNdof();
         hyp->dChiSq_DCdEdx = diter->getChisq();
         hyp->setdEdx(diter->getDEdx());
      }

      hyp->dNDF = hyp->dNDF_Timing+hyp->dNDF_DCdEdx;
      hyp->dChiSq = hyp->dChiSq_Timing+hyp->dChiSq_DCdEdx;
      hyp->dFOM = (hyp->dNDF > 0)? TMath::Prob(hyp->dChiSq,hyp->dNDF) : NAN;

      hyp->setT0(NAN,0,SYS_NULL);
      const hddm_r::StartList &st = iter->getStarts();
      hddm_r::StartList::iterator siter;
      for (siter = st.begin(); siter != st.end(); ++siter) {
         double SCt0 = siter->getT0();
         double SCt0err = siter->getT0err();
         hyp->setT0(SCt0,SCt0err,SYS_START);
      }
      hyp->setT1(tflight+hyp->t0(),
                 sqrt(tferror*tferror-pow(hyp->t0_err(),2)),SYS_NULL);
      hyp->setPathLength(pathlen,patherr);

      data.push_back(hyp);
   }

   // Copy into factory
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

//
// JEventProcessor_danarest.cc
//
// JANA event processor plugin writes out rest events to a file
//
// Richard Jones, 1-July-2012

#include <JANA/JApplication.h>
#include <HDDM/DEventSourceREST.h>

#include "JEventProcessor_danarest.h"
#include "PID/DNeutralShower.h"
#include "TRIGGER/DMCTrigger.h"


// rest output file name, use rest:FILENAME configuration parameter to override
static string restFileName = "dana_rest.hddm";

// mutex for serializing writing to file
static pthread_mutex_t hddmMutex = PTHREAD_MUTEX_INITIALIZER;

// Make us a plugin
// for initializing plugins
extern "C" {
   void InitPlugin(JApplication *app) {
      InitJANAPlugin(app);
      app->AddProcessor(new JEventProcessor_danarest(), true);
   }
} // "extern C"

//-------------------------------
// Constructor
//-------------------------------
JEventProcessor_danarest::JEventProcessor_danarest()
{
   jout << std::endl
        << "  Default JEventProcessor_danarest invoked" 
        << std::endl << std::endl;

   // Check for rest:FILENAME output file name parameter
   gPARMS->SetDefaultParameter("rest:FILENAME",restFileName);
   jout << std::endl << "  rest output file name is " << restFileName
       << std::endl << std::endl;

   ofs = NULL;
   fout = NULL;
   Nevents_written = 0;
}  

//-------------------------------
// Destructor
//-------------------------------
JEventProcessor_danarest::~JEventProcessor_danarest()
{
   if (fout) {
      delete fout;
   }
   if (ofs) {
      delete ofs;
   }
}

//-------------------------------
// init
//-------------------------------
jerror_t JEventProcessor_danarest::init(void)
{
   return NOERROR;
}

//-------------------------------
// brun
//-------------------------------
jerror_t JEventProcessor_danarest::brun(JEventLoop *loop, int runnumber)
{
   // If file is already open, don't reopen it. Just keep adding to it.
   // Not sure what to do with the argument runnumber, just ignore it.

   if (fout) {
      return NOERROR;
   }

   ofs = new ofstream(restFileName.c_str());
   if (!ofs->is_open()) {
      return UNRECOVERABLE_ERROR;
   }
   fout = new hddm_r::ostream(*ofs);

   // enable on-the-fly bzip2 compression on output stream
   fout->setCompression(hddm_r::k_bz2_compression);

   // write a comment record at the head of the file
   hddm_r::HDDM record;
   hddm_r::ReconstructedPhysicsEventList res =
           record.addReconstructedPhysicsEvents(1);
   hddm_r::CommentList comment = res().addComments();
   comment().setText("this is a REST event stream, yadda yadda");
   *fout << record;
   record.clear(); // is this needed?

   Nevents_written = 0;
   return NOERROR;
}

//-------------------------------
// evnt
//-------------------------------
jerror_t JEventProcessor_danarest::evnt(JEventLoop *loop, int eventnumber)
{
   // Write this event to the rest output stream.

   hddm_r::HDDM record;
   hddm_r::ReconstructedPhysicsEventList res =
           record.addReconstructedPhysicsEvents(1);

   // load the run and event numbers
   JEvent& event = loop->GetJEvent();
   res().setRunNo(event.GetRunNumber());
   res().setEventNo(event.GetEventNumber());

   // push any DMCReaction objects to the output record
   std::vector<const DMCReaction*> reactions;
   loop->Get(reactions);
   for (unsigned int i=0; i < reactions.size(); i++) {
      hddm_r::ReactionList rea = res().addReactions(1);
      rea().setType(reactions[i]->type);
      rea().setWeight(reactions[i]->weight);
      rea().setEbeam(reactions[i]->beam.energy());
      rea().setTargetType(reactions[i]->target.PID());

      // Right now the DMCThrown object does not tell which of the listed
      // reactions gave rise to it, so associate them all to the first one.
      if (i == 0) {
         std::vector<const DMCThrown*> throwns;
         loop->Get(throwns);
         double vx,vy,vz,vt;
         vx = vy = vz = vt = -9999;
         hddm_r::VertexList ver = rea().getVertices();
         for (unsigned int it=0; it < throwns.size(); ++it) {
            DVector3 orig(throwns[it]->x(),throwns[it]->y(),throwns[it]->z());
            if (vx != throwns[it]->x() || vy != throwns[it]->y() ||
                vz != throwns[it]->z() || vt != throwns[it]->time() ) {
               ver = rea().addVertices(1);
               hddm_r::OriginList ori = ver().addOrigins(1);
               ori().setT(vt=throwns[it]->time());
               ori().setVx(vx=throwns[it]->x());
               ori().setVy(vy=throwns[it]->y());
               ori().setVz(vz=throwns[it]->z());
            }
            hddm_r::ProductList pro = ver().addProducts(1);
            pro().setId(throwns[it]->myid);
            pro().setParentId(throwns[it]->parentid);
            int pdgtype = throwns[it]->pdgtype;
            if (pdgtype == 0) {
               pdgtype = PDGtype((Particle_t)throwns[it]->type);
            }
            pro().setPdgtype(pdgtype);
            hddm_r::MomentumList mom = pro().addMomenta(1);
            mom().setE(throwns[it]->energy());
            mom().setPx(throwns[it]->px());
            mom().setPy(throwns[it]->py());
            mom().setPz(throwns[it]->pz());
         }
      }
   }

   // push any DTagger objects to the output record
   std::vector<const DTagger*> taggerhits;
   loop->Get(taggerhits);
   for (unsigned int i=0; i < taggerhits.size(); i++) {
      hddm_r::TaggerHitList hit = res().addTaggerHits(1);
      hit().setT(taggerhits[i]->t);
      hit().setE(taggerhits[i]->E);
   }

   // push any DFCALShower objects to the output record
   std::vector<const DFCALShower*> fcalshowers;
   loop->Get(fcalshowers);
   for (unsigned int i=0; i < fcalshowers.size(); i++) {
      hddm_r::CalorimeterClusterList cal = res().addCalorimeterClusters(1);
      DVector3 pos = fcalshowers[i]->getPosition();
      DVector3 poserr = fcalshowers[i]->getPositionError();
      cal().setX(pos(0));
      cal().setY(pos(1));
      cal().setZ(pos(2));
      cal().setT(fcalshowers[i]->getTime());
      cal().setE(fcalshowers[i]->getEnergy());
      cal().setXerr(poserr(0));
      cal().setYerr(poserr(1));
      cal().setZerr(poserr(2));
      cal().setTerr(0);
      cal().setEerr(0);
      cal().setXycorr(0);
      cal().setXzcorr(0);
      cal().setYzcorr(0);
      cal().setEzcorr(0);
      cal().setTzcorr(0);
   }

   // determine which subclass tag to use for DBCALShowers
   static std::string bcalClusterTag;
   static int bcalClusterTagInit=0;
   if (!bcalClusterTagInit) {
      std::vector<const DNeutralShower*> neutralshowers;
      loop->Get(neutralshowers);
      bcalClusterTagInit = 1;
      int useKloeClusters;
      gPARMS->GetParameter("BCALRECON:USE_KLOE", useKloeClusters);
      if (useKloeClusters) {
         bcalClusterTag = "KLOE";
      }
   }

   // push any DBCALShower objects to the output record
   std::vector<const DBCALShower*> bcalshowers;
   loop->Get(bcalshowers,bcalClusterTag.c_str());
   for (unsigned int i=0; i < bcalshowers.size(); i++) {
      hddm_r::CalorimeterClusterList cal = res().addCalorimeterClusters(1);
      cal().setJtag(bcalClusterTag);
      DVector3 pos(bcalshowers[i]->x,bcalshowers[i]->y,bcalshowers[i]->z);
      cal().setX(bcalshowers[i]->x);
      cal().setY(bcalshowers[i]->y);
      cal().setZ(bcalshowers[i]->z);
      cal().setT(bcalshowers[i]->t);
      cal().setE(bcalshowers[i]->E);
      cal().setXerr(bcalshowers[i]->xErr);
      cal().setYerr(bcalshowers[i]->yErr);
      cal().setZerr(bcalshowers[i]->zErr);
      cal().setTerr(bcalshowers[i]->tErr);
      cal().setEerr(0);
      cal().setXycorr(0);
      cal().setXzcorr(0);
      cal().setYzcorr(0);
      cal().setEzcorr(0);
      cal().setTzcorr(0);
   }

   // push any DTOFPoint objects to the output record
   std::vector<const DTOFPoint*> tofpoints;
   loop->Get(tofpoints);
   for (unsigned int i=0; i < tofpoints.size(); i++) {
      hddm_r::TofPointList tof = res().addTofPoints(1);
      tof().setX(tofpoints[i]->pos(0));
      tof().setY(tofpoints[i]->pos(1));
      tof().setZ(tofpoints[i]->pos(2));
      tof().setT(tofpoints[i]->t);
      tof().setDE(tofpoints[i]->dE);
   }

   // push any DSCHit objects to the output record
   std::vector<const DSCHit*> starthits;
   loop->Get(starthits);
   for (unsigned int i=0; i < starthits.size(); i++) {
      hddm_r::StartHitList hit = res().addStartHits(1);
      hit().setSector(starthits[i]->sector);
      hit().setT(starthits[i]->t);
      hit().setDE(starthits[i]->dE);
   }

   // push any DTrackTimeBased objects to the output record
   std::vector<const DTrackTimeBased*> tracks;
   loop->Get(tracks);
   for (unsigned int i=0; i < tracks.size(); ++i) {
      hddm_r::ChargedTrackList tra = res().addChargedTracks(1);
      tra().setCandidateId(tracks[i]->candidateid);
      tra().setPtype(tracks[i]->PID());
      hddm_r::TrackFitList fit = tra().addTrackFits(1);
      fit().setNdof(tracks[i]->Ndof);
      fit().setChisq(tracks[i]->chisq);
      fit().setX0(tracks[i]->x());
      fit().setY0(tracks[i]->y());
      fit().setZ0(tracks[i]->z());
      fit().setPx(tracks[i]->px());
      fit().setPy(tracks[i]->py());
      fit().setPz(tracks[i]->pz());
      fit().setT0(tracks[i]->t0());
      fit().setT0err(tracks[i]->t0_err());
      fit().setT0det(tracks[i]->t0_detector());
      DMatrixDSym errors = tracks[i]->TrackingErrorMatrix();
      fit().setE11(errors(0,0));
      fit().setE12(errors(0,1));
      fit().setE13(errors(0,2));
      fit().setE14(errors(0,3));
      fit().setE15(errors(0,4));
      fit().setE22(errors(1,1));
      fit().setE23(errors(1,2));
      fit().setE24(errors(1,3));
      fit().setE25(errors(1,4));
      fit().setE33(errors(2,2));
      fit().setE34(errors(2,3));
      fit().setE35(errors(2,4));
      fit().setE44(errors(3,3));
      fit().setE45(errors(3,4));
      fit().setE55(errors(4,4));
      if (tracks[i]->dNumHitsUsedFordEdx_FDC +
          tracks[i]->dNumHitsUsedFordEdx_CDC) {
         hddm_r::DEdxDCList elo = tra().addDEdxDCs(1);
         elo().setNsampleFDC(tracks[i]->dNumHitsUsedFordEdx_FDC);
         elo().setNsampleCDC(tracks[i]->dNumHitsUsedFordEdx_CDC);
         elo().setDxFDC(tracks[i]->ddx_FDC);
         elo().setDxCDC(tracks[i]->ddx_CDC);
         elo().setDEdxFDC(tracks[i]->ddEdx_FDC);
         elo().setDEdxCDC(tracks[i]->ddEdx_CDC);
      }

   }

   // push any DMCTrigger objects to the output record
   std::vector<const DMCTrigger*> triggers;
   loop->Get(triggers);
   for (unsigned int i=0; i < triggers.size(); i++) {
      hddm_r::TriggerList trigger = res().addTriggers(1);
      trigger().setL1a(triggers[i]->L1a_fired);
      trigger().setL1b(triggers[i]->L1b_fired);
      trigger().setL1c(triggers[i]->L1c_fired);
   }

   // write the resulting record to the output stream
   pthread_mutex_lock(&hddmMutex);
   *fout << record;
   Nevents_written++;
   pthread_mutex_unlock(&hddmMutex);
   record.clear();  // is this needed?

   return NOERROR;
}

//-------------------------------
// erun
//-------------------------------
jerror_t JEventProcessor_danarest::erun(void)
{
   return NOERROR;
}

//-------------------------------
// fini
//-------------------------------
jerror_t JEventProcessor_danarest::fini(void)
{
   if (fout) {
      pthread_mutex_lock(&hddmMutex);
      delete fout;
      fout = 0;
      pthread_mutex_unlock(&hddmMutex);
      if (ofs) {
         delete ofs;
         ofs = 0;
      }
      std::cout << std::endl <<"Closed REST file" << std::endl;
   }
   std::cout << " " << Nevents_written << " events written to "
             << restFileName << std::endl;
   return NOERROR;
}

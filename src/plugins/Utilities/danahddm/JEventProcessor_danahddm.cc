// JEventProcessor_danahddm.cc
//
//
// JANA event processor plugin writes out hddm event to file
//
//
//  David Lawrence, 7-May-2010

#include <JANA/JApplication.h>
#include <HDDM/DEventSourceHDDM.h>
#include <TRACKING/DTrackTimeBased.h>

#include "JEventProcessor_danahddm.h"


// hddm output file name, use hddm:FILENAME configuration parameter to override
static string hddmFileName = "dana_events.hddm";


// mutex for serializing writing to file
static pthread_mutex_t hddmMutex = PTHREAD_MUTEX_INITIALIZER;

// Make us a plugin
// for initializing plugins
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_danahddm(), true);
  }
} // "extern C"

//-------------------------------
// Constructor
//-------------------------------
JEventProcessor_danahddm::JEventProcessor_danahddm()
{

  jout << std::endl << "  Default JEventProcessor_danahddm invoked" 
       << std::endl << std::endl;

  // Check for hddm:FILENAME output file name parameter
  gPARMS->SetDefaultParameter("hddm:FILENAME",hddmFileName);
  jout << std::endl << "  hddm output file name is " << hddmFileName
       << std::endl << std::endl;
  
  HDDM_USE_COMPRESSION = true;
  gPARMS->SetDefaultParameter("HDDM:USE_COMPRESSION", HDDM_USE_COMPRESSION,
                         "Turn on/off compression of the output HDDM stream."
                         " Set to \"0\" to turn off (it's on by default)");
  HDDM_USE_INTEGRITY_CHECKS = true;
  gPARMS->SetDefaultParameter("HDDM:USE_INTEGRITY_CHECKS",
                               HDDM_USE_INTEGRITY_CHECKS,
                         "Turn on/off automatic integrity checking on the"
                         " output HDDM stream."
                         " Set to \"0\" to turn off (it's on by default)");
  file = NULL;
  fout = NULL;
  Nevents_written = 0;
}  

//-------------------------------
// Destructor
//-------------------------------
JEventProcessor_danahddm::~JEventProcessor_danahddm()
{
}

//-------------------------------
// init
//-------------------------------
jerror_t JEventProcessor_danahddm::init(void)
{
   return NOERROR;
}

//-------------------------------
// brun
//-------------------------------
jerror_t JEventProcessor_danahddm::brun(JEventLoop *loop, int runnumber)
{
   // If file is already open, don't reopen it. Just keep adding to it.
   if (file)
      return NOERROR;

   // We wait until here to open the output so that we can check if the 
   // input is hddm. If it's not, tell the user and exit immediately
   JEvent& event = loop->GetJEvent();
   JEventSource *source = event.GetJEventSource();
   DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
   if (! hddm_source) {
      std::cerr << " This program MUST be used with an HDDM file as input!" << std::endl;
      exit(-1);
   }

   // If we got here, it must be an HDDM source. Open a new file.
   file = new std::ofstream(hddmFileName.c_str());
   fout = new hddm_s::ostream(*file);
   Nevents_written = 0;
 
   // enable on-the-fly bzip2 compression on output stream
   if (HDDM_USE_COMPRESSION) {
      jout << " Enabling bz2 compression of output HDDM file stream" 
           << std::endl;
      fout->setCompression(hddm_s::k_bz2_compression);
   }
   else {
      jout << " HDDM compression disabled" << std::endl;
   }

   // enable a CRC data integrity check at the end of each event record
   if (HDDM_USE_INTEGRITY_CHECKS) {
      jout << " Enabling CRC data integrity check in output HDDM file stream"
           << std::endl;
      fout->setIntegrityChecks(hddm_s::k_crc32_integrity);
   }
   else {
      jout << " HDDM integrity checks disabled" << std::endl;
   }

   return NOERROR;
}

//-------------------------------
// evnt
//-------------------------------
jerror_t JEventProcessor_danahddm::evnt(JEventLoop *loop, int eventnumber)
{
   JEvent& event = loop->GetJEvent();
   JEventSource *source = event.GetJEventSource();
   DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
   if (! hddm_source) {
      std::cerr << " This program MUST be used only with HDDM files as inputs!"
                << std::endl;
      exit(-1);
   }
   hddm_s::HDDM *hddm = (hddm_s::HDDM*)event.GetRef();
   if (! hddm)
      return NOERROR;
   
   // Delete any data in the reconView branch of the event.
   hddm->getPhysicsEvent().deleteReconViews();

   // Fill in reconstructed banks, replacing any that are already there
   hddm_s::ReconViewList revs = hddm->getPhysicsEvent().addReconViews();
   Add_DTrackTimeBased(loop, revs.begin());

   // get write lock
   pthread_mutex_lock(&hddmMutex);

   // Write event to file and update counter
   *fout << *hddm;
   Nevents_written++;

   // unlock
   pthread_mutex_unlock(&hddmMutex);

   return NOERROR;
}

//-------------------------------
// erun
//-------------------------------
jerror_t JEventProcessor_danahddm::erun(void)
{
   return NOERROR;
}

//-------------------------------
// fini
//-------------------------------
jerror_t JEventProcessor_danahddm::fini(void)
{
   if (fout)
      delete fout;
   if (file) {
      delete file;
      std::cout << std::endl << "Closed HDDM file" << std::endl;
   }
   std::cout << " " << Nevents_written << " event written to "
             << hddmFileName << std::endl;

   return NOERROR;
}

//-------------------------------
// Add_DTrackTimeBased
//-------------------------------
void JEventProcessor_danahddm::Add_DTrackTimeBased(JEventLoop *loop, 
                               hddm_s::ReconViewList::iterator riter)
{
   // Get objects to write out
   vector<const DTrackTimeBased*> tracktimebaseds;
   loop->Get(tracktimebaseds);
   if (tracktimebaseds.size() == 0)
      return;

   // Allocate memory for all time based tracks
   unsigned int ntbts = tracktimebaseds.size();
   hddm_s::TracktimebasedList tbts = riter->addTracktimebaseds(ntbts);
   for (unsigned int i=0; i < ntbts; i++) {
      const DTrackTimeBased *tbt_dana = tracktimebaseds[i];
      DVector3 pos = tbt_dana->position();
      DVector3 mom = tbt_dana->momentum();
      
      tbts(i).setFOM(tbt_dana->FOM);
      tbts(i).setCandidateid(tbt_dana->candidateid);
      tbts(i).setTrackid(tbt_dana->trackid);
      tbts(i).setId(tbt_dana->id);
      tbts(i).setChisq(tbt_dana->chisq);
      tbts(i).setNdof(tbt_dana->Ndof);

      hddm_s::MomentumList tmoms = tbts(i).addMomenta();
      hddm_s::PropertiesList tpros = tbts(i).addPropertiesList();
      hddm_s::OriginList torig = tbts(i).addOrigins();
      hddm_s::ErrorMatrixList terrs = tbts(i).addErrorMatrixs();
      hddm_s::TrackingErrorMatrixList tters = 
                              tbts(i).addTrackingErrorMatrixs();

      tmoms().setE(tbt_dana->energy());
      tmoms().setPx(mom.x());
      tmoms().setPy(mom.y());
      tmoms().setPz(mom.z());

      tpros().setCharge((int)tbt_dana->charge());
      tpros().setMass(tbt_dana->mass());
      
      torig().setT(0.0);
      torig().setVx(pos.x());
      torig().setVy(pos.y());
      torig().setVz(pos.z());
      
      string vals = DMatrixDSymToString(tbt_dana->errorMatrix());
      terrs().setNcols(7);
      terrs().setNrows(7);
      terrs().setType("DMatrixDSym");
      terrs().setVals(vals.c_str());

      tters().setNcols(5);
      tters().setNrows(5);
      tters().setType("DMatrixDSym");
      tters().setVals(DMatrixDSymToString(tbt_dana->TrackingErrorMatrix()));
   }
}

//-------------------------------
// DMatrixDSymToString
//-------------------------------
string JEventProcessor_danahddm::DMatrixDSymToString(const DMatrixDSym &mat)
{
   // Convert the given symmetric matrix into a single string that
   // can be used in an HDDM file.

   stringstream ss;
   for (int irow=0; irow<mat.GetNrows(); irow++) {
      for (int icol=irow; icol<mat.GetNcols(); icol++) {
         ss << mat[irow][icol] << " ";
      }
   }
   
   return ss.str();
}

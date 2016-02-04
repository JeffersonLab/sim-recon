// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"

#include <JANA/JEvent.h>

#include <HDDM/DEventSourceHDDM.h>
#include <TRACKING/DMCThrown.h>

#include <TRIGGER/DMCTrigger.h>

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
   // open HDDM file
   filename = "bggen_filtered.hddm";
   ofs = new ofstream(filename.c_str());
   if (!ofs->is_open()) {
      fout = 0;
      return UNRECOVERABLE_ERROR;
   }
   fout = new hddm_s::ostream(*ofs);

   HDDM_USE_COMPRESSION = false;
   gPARMS->SetDefaultParameter("HDDM:USE_COMPRESSION", HDDM_USE_COMPRESSION,
                          "Turn on/off compression of the output HDDM stream."
                          " Set to \"1\" to turn on (it's off by default)");
   HDDM_USE_INTEGRITY_CHECKS = false;
   gPARMS->SetDefaultParameter("HDDM:USE_INTEGRITY_CHECKS",
                                HDDM_USE_INTEGRITY_CHECKS,
                          "Turn on/off automatic integrity checking on the"
                          " output HDDM stream."
                          " Set to \"1\" to turn on (it's off by default)");

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

   Nevents_written = 0;

   return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *loop, int eventnumber)
{
   JEvent& event = loop->GetJEvent();
   JEventSource *source = event.GetJEventSource();
   DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
   if (!hddm_source) {
      cerr << " This program MUST be used with an HDDM file as input!" << endl;
      exit(-1);
   }
   hddm_s::HDDM *hddm = (hddm_s::HDDM*)event.GetRef();
   if (!hddm)
      return NOERROR;
   
   
   // Initialize write_out flag. We set this to true to write the event
   // out and false to ignore it. Initialize it to false so that we can
   // pick the conditions needed to keep it.
   bool write_out=false;

   // Here we do whatever calculations are needed to determine if we keep
   // the event. Since this is basically a skeleton meant as an example
   // we do a trivial check on the momentum of the thrown particles.
   // In practice, one could request objects that require full reconstruction
   // as well so that filters could be built on those quantities as well.

   //---------------------- Filter Code Start ----------------------
   // Get data
   vector<const DMCThrown*> mcthrowns;
   vector<const DMCThrown*> mcthrownsFS;
   loop->Get(mcthrowns);
   loop->Get(mcthrownsFS,"FinalState");

/*   bool etaprime = false;
   bool piplus = false;
   bool piminus = false;
   bool eta = false;
   int myid;*/
   int etaid;
   int gammas = 0;
/*
   for (int i = 0; i < (int)mcthrowns.size(); i++) {
     if (mcthrowns[i]->pdgtype == 331) {
       etaprime = true;
       myid = mcthrowns[i]->myid;
       break;
     }
   }

   if (etaprime) {
     for (int i = 0; i < (int)mcthrowns.size(); i++) {
       if (mcthrowns[i]->pdgtype == 211 && mcthrowns[i]->parentid == myid) piplus = true;
     }
     for (int i = 0; i < (int)mcthrowns.size(); i++) {
       if (mcthrowns[i]->pdgtype == -211 && mcthrowns[i]->parentid == myid) piminus = true;
     }
     for (int i = 0; i < (int)mcthrowns.size(); i++) {
       if (mcthrowns[i]->pdgtype == 221 && mcthrowns[i]->parentid == myid) {
         eta = true;
         etaid = mcthrowns[i]->myid;
       }
     }

//     if (piplus && piminus && eta) write_out = true;
*/
//     if (piplus && piminus && eta) {
       etaid = 1;
       for (int i = 0; i < (int)mcthrowns.size(); i++) {
         if (mcthrowns[i]->pdgtype == 22 && mcthrowns[i]->parentid == etaid) gammas++;
       }
       if (gammas == 2) {
         if (mcthrownsFS.size() == 5) write_out = true;
       }
//     }
//   }

/*   if (mcthrowns.size() == 6) {
     int photoncount = 0;
     for (unsigned int i=0; i < mcthrowns.size(); i++) {
        const DMCThrown *mcthrown = mcthrowns[i];
        if (mcthrown->type == 1) photoncount++;
     }
     if (photoncount == 2) write_out = true;
   }*/
   
   //----------------------- Filter Code End -----------------------
   
   // If write_out flag is set, write this event to our output file
   // otherwise, just flush the memory.
   //
   // WARNING: If a plugin is used with this program that tries to
   // access objects in the s_HDDM_t structure (pretty much every
   // plugin out there), it will likely seg. fault due to the
   // memory already being freed!
   if (write_out) {
      *fout << *hddm;
      Nevents_written++;
   }
   else {
      hddm->clear();
   }
   
   return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
   if (fout) {
      delete fout;
      fout = 0;
   }
   if (ofs) {
      ofs->close();
      delete ofs;
      cout << endl << "Closed HDDM output file" << endl;
   }
   cout << " " << Nevents_written << " events written to " << filename << endl;
   
   return NOERROR;
}

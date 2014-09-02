// $Id$
//
// Created Oct 10, 2013  David Lawrence

#include "hddm_cull_events.h"

#include <HDDM/hddm_r.hpp>
using namespace hddm_r;


// NOTE: The REST format files claim they are in a compressed format that
// the C API cannot handle so we use the C++ API here


//-----------
// Process_r  --  HDDM REST format
//-----------
void Process_r(unsigned int &NEvents, unsigned int &NEvents_read)
{
   // Output file
   std::cout << " output file: " << OUTFILENAME << std::endl;
   ofstream ofs(OUTFILENAME);
   if (!ofs.is_open()) {
      std::cout << " Error opening output file \"" << OUTFILENAME << "\"!" << std::endl;
      exit(-1);
   }

   hddm_r::ostream ostr(ofs);
   if (HDDM_USE_COMPRESSION) {
      std::cout << " Enabling bz2 compression of output HDDM file stream" 
                << std::endl;
      ostr.setCompression(hddm_r::k_bz2_compression);
   }
   else {
      std::cout << " HDDM compression disabled on output" << std::endl;
   }

   if (HDDM_USE_INTEGRITY_CHECKS) {
      std::cout << " Enabling data integrity check on output HDDM file stream"
                << std::endl;
      ostr.setIntegrityChecks(hddm_r::k_crc32_integrity);
   }
   else {
      std::cout << " HDDM integrity checks disabled on output" << std::endl;
   }

   // Loop over input files
   time_t last_time = time(NULL);
   for (unsigned int i=0; i<INFILENAMES.size(); i++) {
      std::cout << " input file: " << INFILENAMES[i] << std::endl;

      // Open hddm file for reading
      ifstream ifs(INFILENAMES[i]);

      // Associate input file stream with HDDM record
      hddm_r::istream istr(ifs);

      // Loop over events
      while (!ifs.eof() && ifs.good()) {
         try{
            HDDM xrec;
            istr >> xrec;
            NEvents_read++;
            
            bool write_this_event = false;
            
            // Loop over physics events within this event and see if one
            // has the event number of interest
            class ReconstructedPhysicsEvent &reconstructedPhysicsEvent = xrec.getReconstructedPhysicsEvent();
            if (EVENT_TO_KEEP_MODE) { // need to check if reconstructedPhysicsEvent is valid!!
               int eventNo = reconstructedPhysicsEvent.getEventNo();
               if ((unsigned int)eventNo == SPECIFIC_EVENT_TO_KEEP) {
                  write_this_event = true;
                  QUIT = true;
               }
            }
            
            // Check if we're in the range of offsets to write out
            if (NEvents_read>EVENTS_TO_SKIP)write_this_event = true;
            
            // Write this output event to file and free its memory
            if (write_this_event) {
               ostr << xrec;
               NEvents++;
            }
         
            // Update ticker
            time_t now = time(NULL);
            if (now != last_time) {
               std::cout << "  " << NEvents_read << " events read     (" << NEvents << " event written) \r"; std::cout.flush();
               last_time = now;
            }
            
            // Quit as soon as we wrote all of the events we're going to
            if (NEvents_read>=(EVENTS_TO_SKIP+EVENTS_TO_KEEP))break;

            if (QUIT)break;
            

         }catch(...) {
            break;
         }
      }
   }
}


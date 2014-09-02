// $Id$
//
// Created Oct 10, 2013  David Lawrence

#include "hddm_cull_events.h"

#include <HDDM/hddm_s.hpp>
#include <iostream>
#include <fstream>

#include <stdlib.h>

//-----------
// Process_s  --  HDDM simulation format
//-----------
void Process_s(unsigned int &NEvents, unsigned int &NEvents_read)
{
   // Output file
   std::ofstream *ofs = new std::ofstream(OUTFILENAME);
   if (! ofs->is_open()) {
      std::cout << " Error opening output file \"" << OUTFILENAME << "\"!"
                << std::endl;
      exit(-1);
   }

   hddm_s::ostream *fout = new hddm_s::ostream(*ofs);
   if (HDDM_USE_COMPRESSION) {
      std::cout << " Enabling bz2 compression of output HDDM file stream" 
                << std::endl;
      fout->setCompression(hddm_s::k_bz2_compression);
   }
   else {
      std::cout << " HDDM compression disabled on output" << std::endl;
   }

   if (HDDM_USE_INTEGRITY_CHECKS) {
      std::cout << " Enabling data integrity check on output HDDM file stream"
                << std::endl;
      fout->setIntegrityChecks(hddm_s::k_crc32_integrity);
   }
   else {
      std::cout << " HDDM integrity checks disabled on output" << std::endl;
   }

   // Loop over input files
   time_t last_time = time(NULL);
   for (unsigned int i=0; i < INFILENAMES.size(); i++) {
      std::cout << " input file: " << INFILENAMES[i] << std::endl;
      std::ifstream *ifs = new std::ifstream(INFILENAMES[i]);
      if (! ifs->is_open()) {
         std::cout << " Error opening input file \"" << INFILENAMES[i]
                   << "\"!" << std::endl;
         exit(-1);
      }
      hddm_s::istream *fin = new hddm_s::istream(*ifs);
         
      // Loop over all events in input
      while (ifs->good()) {
         hddm_s::HDDM record;
         *fin >> record;
         NEvents_read++;
         
         bool write_this_event = false;
         
         // Loop over physics events within this event and see if one
         // has the event number of interest
         if (EVENT_TO_KEEP_MODE) {
            hddm_s::PhysicsEventList pes = record.getPhysicsEvents();
            hddm_s::PhysicsEventList::iterator eviter;
            for (eviter = pes.begin(); eviter != pes.end(); ++eviter) {
               int eventNo = eviter->getEventNo();
               if ((unsigned int)eventNo == SPECIFIC_EVENT_TO_KEEP) {
                  write_this_event = true;
                  QUIT = true;
               }
            }
         }

         // Check if we're in the range of offsets to write out
         if (NEvents_read > EVENTS_TO_SKIP)
            write_this_event = true;
         
         // Write this output event to file and free its memory
         if (write_this_event) {
            *fout << record;
            NEvents++;
         }
         record.clear();
      
         // Update ticker
         time_t now = time(NULL);
         if (now != last_time) {
            std::cout << "  " << NEvents_read << " events read     ("
                      << NEvents << " event written) \r";
            std::cout.flush();
            last_time = now;
         }
         
         // Quit as soon as we wrote all of the events we're going to
         if (NEvents_read >= (EVENTS_TO_SKIP+EVENTS_TO_KEEP))
            break;

         if (QUIT)break;
      }

      // Close input file
      delete fin;
      delete ifs;
   }
      
   // Close output file
   delete fout;
   delete ofs;
}

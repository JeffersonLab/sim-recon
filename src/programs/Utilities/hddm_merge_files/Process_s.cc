// $Id$
//
// Created Oct 25, 2013  Kei Moriya

#include "hddm_merge_files.h"

#include <HDDM/hddm_s.hpp>

//-----------
// Process_s  --  HDDM simulation format
//-----------
void Process_s(unsigned int &NEvents, unsigned int &NEvents_read)
{
   // Output file
   std::cout << " output file: " << OUTFILENAME << std::endl;
   std::ofstream ofs(OUTFILENAME);
   if (! ofs.is_open()) {
      std::cout << " Error opening output file \"" << OUTFILENAME 
                << "\"!" << std::endl;
      exit(-1);
   }
   hddm_s::ostream *fout = new hddm_s::ostream(ofs);
   if (HDDM_USE_COMPRESSION) {
      std::cout << " Enabling bz2 compression of output HDDM file stream" 
               << std::endl;
      fout->setCompression(hddm_s::k_bz2_compression);
   }
   else {
      std::cout << " HDDM compression disabled" << std::endl;
   }
   if (HDDM_USE_INTEGRITY_CHECKS) {
      std::cout << " Enabling CRC data integrity check in output HDDM"
                   " file stream" << std::endl;
      fout->setIntegrityChecks(hddm_s::k_crc32_integrity);
   }
   else {
      std::cout << " HDDM integrity checks disabled" << std::endl;
   }

   // Loop over input files
   time_t last_time = time(NULL);
   for (unsigned int i=0; i<INFILENAMES.size(); i++) {
      std::cout << " input file: " << INFILENAMES[i] << std::endl;
      std::ifstream ifs(INFILENAMES[i]);
      if (! ifs.is_open()) {
         std::cout << " Error opening input file \"" << INFILENAMES[i]
                   << "\"!" << std::endl;
         exit(-1);
      }
      hddm_s::istream *fin = new hddm_s::istream(ifs);

      // Loop over all events in input
      while (ifs.good()) {
         hddm_s::HDDM record;
         *fin >> record;
         NEvents_read++;
         
         // Write this output event to file
         *fout << record;
         NEvents++;
      
         // Update ticker
         time_t now = time(NULL);
         if (now != last_time) {
            std::cout << "  " << NEvents_read << " events read     (" 
                      << NEvents << " event written) \r";
            std::cout.flush();
            last_time = now;
         }

         if (QUIT)
            break;
      }

      // Close input file
      ifs.close();
      delete fin;
   }
      
   // Close output file
   ofs.close();
}

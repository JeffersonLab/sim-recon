// $Id$
//
// Created Oct 25, 2013  Kei Moriya

#include "hddm_merge_files.h"

#include <HDDM/hddm_r.hpp>
using namespace hddm_r;


//-----------
// Process_r  --  HDDM REST format
//-----------
void Process_r(unsigned int &NEvents, unsigned int &NEvents_read)
{
   // Output file
   std::cout << " output file: " << OUTFILENAME << std::endl;
   std::ofstream ofs(OUTFILENAME);
   if (! ofs.is_open()) {
      std::cout << " Error opening output file \"" << OUTFILENAME 
                << "\"!" << std::endl;
      exit(-1);
   }
   hddm_r::ostream *ostr = new hddm_r::ostream(ofs);
   if (HDDM_USE_COMPRESSION) {
      std::cout << " Enabling bz2 compression of output HDDM file stream" 
               << std::endl;
      ostr->setCompression(hddm_r::k_bz2_compression);
   }
   else {
      std::cout << " HDDM compression disabled" << std::endl;
   }
   if (HDDM_USE_INTEGRITY_CHECKS) {
      std::cout << " Enabling CRC data integrity check in output HDDM"
                   " file stream" << std::endl;
      ostr->setIntegrityChecks(hddm_r::k_crc32_integrity);
   }
   else {
      std::cout << " HDDM integrity checks disabled" << std::endl;
   }

   // Loop over input files
   time_t last_time = time(NULL);
   for (unsigned int i=0; i<INFILENAMES.size(); i++) {
      std::cout << " input file: " << INFILENAMES[i] << std::endl;

      // Open hddm file for reading
      std::ifstream ifs(INFILENAMES[i]);

      // Associate input file stream with HDDM record
      hddm_r::istream *istr = new hddm_r::istream(ifs);
      
      // Loop over events
      while (ifs.good()) {
         HDDM xrec;
         *istr >> xrec;
         NEvents_read++;

         *ostr << xrec;
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
   }
}

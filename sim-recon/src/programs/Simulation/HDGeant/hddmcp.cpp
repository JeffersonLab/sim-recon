/*
 * hddmcp - a utility program for copying from one hddm file to another,
 *          and in the process stripping out everything except a desired
 *          subset of the data fields.  The template for the output hddm
 *          file (a subset of what is found in the input hddm file) is
 *          assumed to have already been created, perhaps by grabbing the
 *          header from the input hddm file into your favorite text editor,
 *          stripping out the desired parts, and saving the resulting xml
 *          template into a new file called x.hddm, which then must be
 *          processed into c++ using the following command.
 *
 *            $ hddm-cpp x.hddm
 *
 * Richard Jones
 * GlueX collaboration
 * January 14, 2012
 *
 */

#include <fstream>
#include <stdlib.h>
#include <HDDM/hddm_s.hpp>

int main(int argc, char **argv)
{
   int input, output;
   if (argc < 2) {
      std::cerr << "Usage: hddmcp <inputfile1.hddm> ... <outputfile.hddm>"
                << std::endl;
      exit(1);
   }
   else {
      output = argc-1;
   }

   std::ofstream ofs(argv[output]);
   if (!ofs.is_open()) {
      std::cerr << "Usage: hddmcp <inputfile1.hddm> ... <outputfile.hddm>"
                << std::endl;
      exit(1);
   }
   hddm_s::ostream ostr(ofs);
   ostr.setCompression(hddm_s::k_bz2_compression);

   std::ifstream ifs;
   hddm_s::HDDM record;
   for (input=1; input<output; input++) {
      ifs.open(argv[input]);
      if (!ifs.is_open()) {
         std::cerr << "Error - could not open input file " << argv[input]
                   << std::endl;
         exit(1);
      }
      else {
         hddm_s::istream istr(ifs);
         //istr.skip(9999);
         int count=0;
         while (ifs.good()) {
            istr >> record;
            ostr << record;
            record.clear();
            ++count;
         }
      }
   }
}

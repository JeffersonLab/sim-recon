//
// nullgen - produces a file of null (empty) simulated events. This can
//           be used in combination with mcsmear if all you want to do
//           is to merge the hits from events drawn from multiple input
//           hddm streams into a a single stream of fatter events.
//
// author: richard.t.jones at uconn.edu
// version: march 17, 2017
//
// usage: nullgen -n <number> -o <output filename>
//
// example: 
// Suppose you have a set of three files x1.hddm, x2.hddm, and x3.hddm
// all containing events that have already been smeared, and you want
// to simply merge these such that the output events are merges of 2
// events from x1.hddm, 3 events from x2.hddm, and 1 from x3.hddm,
// then the following sequence of commands will generate 1000 events
// in xsum.hddm, the merged output file.
//
//    $ nullgen -n 1000 -o xnull.hddm
//    $ mcsmear x1.hddm:2 x2.hddm:3 x3.hddm:1 -o xsum.hddm

#include <HDDM/hddm_s.hpp>
#include <fstream>
#include <string>

void usage()
{
   std::cout << "usage: nullgen -n <count> "
             << "[-r <run number>] "
             << "-o <output filename>"
             << std::endl;
   exit(1);
}

int main(int argc, char *argv[])
{
   int event_count = 0;
   int run_number = 0;
   std::string output_file;

   for (int iarg=1; iarg < argc; ++iarg) {
      std::string arg(argv[iarg]);
      if (arg.substr(0,2) == "-n") {
         if (arg.size() > 2) {
            event_count = std::stoi(arg.substr(2));
         }
         else {
            event_count = std::stoi(argv[++iarg]);
         }
      }
      else if (arg.substr(0,2) == "-r") {
         if (arg.size() > 2) {
            run_number = std::stoi(arg.substr(2));
         }
         else {
            run_number = std::stoi(argv[++iarg]);
         }
      }
      else if (arg.substr(0,2) == "-o") {
         if (arg.size() > 2) {
            output_file = arg.substr(2);
         }
         else {
            output_file = argv[++iarg];
         }
      }
      else {
         usage();
      }
   }

   if (event_count == 0 || output_file.size() == 0) {
      usage();
   }
   std::ofstream outf(output_file);
   hddm_s::ostream fout(outf);
   hddm_s::HDDM record;
   hddm_s::PhysicsEventList elist = record.addPhysicsEvents();
   elist(0).setRunNo(run_number);
   for (int event=0; event < event_count; ++event) {
      elist(0).setEventNo(event);
      fout << record;
   }
   return 0;
}

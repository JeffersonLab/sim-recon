/*
 * cdccount - an example program for accessing the contents
 *            of events stored in a hddm file.
 *
 * Richard Jones
 * GlueX collaboration
 * January 14, 2012
 *
 */

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <HDDM/hddm_s.hpp>

int process_event(hddm_s::HDDM &event);

int main(int argc, char **argv)
{
   hddm_s::HDDM record;
   int events_with_hits = 0;
   int events = 0;
   for (int input=1; input<argc; input++) {
      std::ifstream ifs(argv[input]);
      if (!ifs.is_open()) {
         std::cerr << "Error - could not open input file "
                   << argv[input] << std::endl;
         exit(1);
      }
      hddm_s::istream istr(ifs);
      while (ifs.good()) {
         istr >> record;
         events_with_hits += process_event(record);
         record.clear();
         ++events;
      }
   }
   printf("Total events seen %d, with cdc hits %d.\n",events,events_with_hits);
}

int process_event(hddm_s::HDDM &event)
{
   hddm_s::CdcTruthPointList truthPoints = event.getCdcTruthPoints();
   return truthPoints.size();
}

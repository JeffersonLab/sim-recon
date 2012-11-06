/*
 * cdcdump - an example program for accessing the contents
 *           of events stored in a hddm file.
 *
 * Richard Jones
 * GlueX collaboration
 * January 10, 2005
 *
 */

#include <fstream>
#include <stdlib.h>
#include "hddm_s.hpp"

int process_event(hddm_s::HDDM &event);

int main(int argc, char **argv)
{
   hddm_s::HDDM record;
   int input;
   for (input=1; input<argc; input++) {
      std::ifstream ifs(argv[input]);
      if (!ifs.is_open()) {
         std::cerr << "Error - could not open input file "
                   << argv[input] << std::endl;
         exit(1);
      }
      hddm_s::istream istr(ifs);
      while (ifs.good()) {
         istr >> record;
         process_event(record);
         record.clear();
      }
   }
}

int process_event(hddm_s::HDDM &event)
{
   std::cout << "New event number " << event.getPhysicsEvent().getEventNo()
             << " run number " << event.getPhysicsEvent().getRunNo()
             << std::endl;

   hddm_s::CdcTruthPointList truths = event.getCdcTruthPoints();
   std::cout << " found " << truths.size() << " cdcTruthPoints!"
             << std::endl;

   hddm_s::CdcTruthPointList::iterator iter;
   for (iter = truths.begin(); iter != truths.end(); ++iter) {
      if (fabs(iter->getDradius()-19.5) < 0.5e5) {
         std::cout << "  dradius=" << iter->getDradius() << ","
                   << "  phi=" << iter->getPhi() << ","
                   << "  primary=" << ((iter->getPrimary())? "true,":"false,")
                   << "  r=" << iter->getR() << ","
                   << "  track=" << iter->getTrack() << ","
                   << "  z=" << iter->getZ() << ","
                   << "  dE/dx=" << iter->getDEdx()*1e6 << std::endl;
      }
   }
   return 1;
}

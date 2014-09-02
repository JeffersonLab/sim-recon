// $Id: mcsmear.cc 2388 2007-01-10 16:46:03Z davidl $
//
// Created August 24, 2007  David Lawrence

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include <signal.h>
#include <time.h>

#include "HDDM/hddm_s.hpp"

bool Filter(hddm_s::HDDM &record);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;


//-----------
// main
//-----------
int main(int narg,char* argv[])
{
   // Set up to catch SIGINTs for graceful exits
   signal(SIGINT,ctrlCHandle);

   ParseCommandLineArguments(narg, argv);
   
   std::cout << " input file: " << INFILENAME << std::endl;
   std::cout << " output file: " << OUTFILENAME << std::endl;
   
   // Open Input file
   std::ifstream *ifs = new ifstream(INFILENAME);
   if (! ifs->is_open()) {
      std::cout << " Error opening input file \"" << INFILENAME << "\"!"
                << std::endl;
      exit(-1);
   }
   hddm_s::istream *fin = new hddm_s::istream(*ifs);
   
   // Output file
   std::ofstream *ofs = new ofstream(OUTFILENAME);
   if (! ofs->is_open()) {
      std::cout << " Error opening output file \"" << OUTFILENAME << "\"!" 
                << std::endl;
      exit(-1);
   }
   hddm_s::ostream *fout = new hddm_s::ostream(*ofs);
   
   // Loop over events in input file
   int NEvents_read = 0;
   int NEvents_written = 0;
   time_t last_time = time(NULL);
   while (ifs->good()) {
      NEvents_read++;
      hddm_s::HDDM record;
      *fin >> record;
      time_t now = time(NULL);
      if (now != last_time) {
         std::cout << " " << NEvents_read << " events read -- " 
                   << NEvents_written << " events written      \r";
         std::cout.flush();
         last_time = now;
      }
      
      // Write or don't depending on return value of Filter()
      if (Filter(record)) {
         *fout << record;
         NEvents_written++;
      }
      
      if (QUIT)
         break;
   }
   
   // close input and output files
   delete fin;
   delete ifs;
   delete fout;
   delete ofs;

   std::cout << std::endl << "FINAL:" << endl;
   std::cout << " " << NEvents_read << " events read -- " 
             << NEvents_written << " events written" << std::endl;
   std::cout << "Output file has " 
             << 100.0*(double)NEvents_written/(double)NEvents_read 
             << "% of the events that were in the input file." << std::endl;
   return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{

   for (int i=1; i<narg; i++) {
      char *ptr = argv[i];
      
      if (ptr[0] == '-') {
         switch(ptr[1]) {
          case 'h': Usage();
            break;
         }
      }
      else {
         INFILENAME = argv[i];
      }
   }

   if (! INFILENAME) {
      std::cout << std::endl << "You must enter a filename!" 
                << std::endl << std::endl;
      Usage();
   }
   
   // Generate output filename based on input filename
   char *ptr, *path_stripped;
   path_stripped = ptr = strdup(INFILENAME);
   while ((ptr = strstr(ptr, "/")))
      path_stripped = ++ptr;
   ptr = strstr(path_stripped, ".hddm");
   if (ptr)
      *ptr=0;
   char str[256];
   sprintf(str, "%s_filtered.hddm", path_stripped);
   OUTFILENAME = strdup(str);
}

//-----------
// Usage
//-----------
void Usage(void)
{
   std::cout << std::endl << "Usage:" << endl;
   std::cout << "     filtergen [options] file.hddm" << std::endl;
   std::cout << std::endl;
   std::cout << " Read filter events from the given HDDM file while copying it"
             << std::endl;
   std::cout << "to another HDDM file. This was written to provide an easy way"
             << std::endl;
   std::cout << "to filter events from a file produced by pythia that one does"
             << std::endl;
   std::cout << "not want to waste time tracking through the whole detector."
             << std::endl;
   std::cout << std::endl;
   std::cout << "At this point, the filtering conditions are hardwired and"
             << std::endl;
   std::cout << "can only be changed by editing the source. In the future, we"
             << std::endl;
   std::cout << "will want command line options so the filter can be tuned" 
             << std::endl;
   std::cout << "without recompiling." << std::endl;
   std::cout << std::endl;
   std::cout << "  options:" << std::endl;
   std::cout << "    -h       Print this usage statement." << std::endl;
   std::cout << std::endl;

   exit(0);
}

//-----------------------------------------------------------------
// ctrlCHandle
//-----------------------------------------------------------------
void ctrlCHandle(int x)
{
   QUIT++;
   std::cerr << std::endl << "SIGINT received (" << QUIT << ")....." << endl;
}

// $Id: mcsmear.cc 2388 2007-01-10 16:46:03Z davidl $
//
// Created August 24, 2007  David Lawrence

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
using namespace std;

#include <signal.h>
#include <stdlib.h>
#include <time.h>

#include "HDDM/hddm_s.hpp"

bool Filter(hddm_s::HDDM &record);
bool Filter_eta_p(hddm_s::HDDM &record);
bool Filter_eta_p_pi0(hddm_s::HDDM &record);
bool Filter_eta_n_pip(hddm_s::HDDM &record);
bool Filter_eta_p_gamma(hddm_s::HDDM &record);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);
hddm_s::ostream *NewHDDMFile(std::string fname);
std::map<hddm_s::ostream*, std::ofstream*> ofslist;


char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;
bool CREATE_SINGLE_CHANNEL_FILES = false;


#ifndef _DBG_
#define _DBG_ std::cerr << __FILE__ << ":" << __LINE__ << " "
#define _DBG__ _DBG_ << std::endl
#endif

//-----------
// main
//-----------
int main(int narg,char* argv[])
{
   // Set up to catch SIGINTs for graceful exits
   signal(SIGINT,ctrlCHandle);

   ParseCommandLineArguments(narg, argv);
      
   // Open Input file
   std::cout << " input file: " << INFILENAME << std::endl;
   std::ifstream *ifs = new std::ifstream(INFILENAME);
   if (! ifs->is_open()) {
      std::cout << " Error opening input file \"" << INFILENAME << "\"!"
                << std::endl;
      exit(-1);
   }
   hddm_s::istream *fin = new hddm_s::istream(*ifs);
   
   // Output file
   hddm_s::ostream *fout = NewHDDMFile(OUTFILENAME);
   
   // Optionally create outputfiles for single channels
   hddm_s::ostream *fout_eta_p = NULL;
   hddm_s::ostream *fout_eta_p_pi0 = NULL;
   hddm_s::ostream *fout_eta_n_pip = NULL;
   hddm_s::ostream *fout_eta_p_gamma = NULL;
   int NEvents_written_eta_p = 0;
   int NEvents_written_eta_p_pi0 = 0;
   int NEvents_written_eta_n_pip = 0;
   int NEvents_written_eta_p_gamma = 0;
   if (CREATE_SINGLE_CHANNEL_FILES) {
      fout_eta_p = NewHDDMFile("eta_p.hddm");
      fout_eta_p_pi0 = NewHDDMFile("eta_p_pi0.hddm");
      fout_eta_n_pip = NewHDDMFile("eta_n_pip.hddm");
      fout_eta_p_gamma = NewHDDMFile("eta_p_gamma.hddm");
   }

   // Loop over events in input file
   int NEvents_read = 0;
   int NEvents_written = 0;
   time_t last_time = time(NULL);
   while (ifs->good()) {
      hddm_s::HDDM record;
      *fin >> record;
      NEvents_read++;
      time_t now = time(NULL);
      if (now != last_time) {
         std::cout << " " << NEvents_read << " events read -- "
                          << NEvents_written << " events written      \r";
         std::cout.flush();
         last_time = now;
      }
      
      // Write or don't depending on return value of Filter()
      if (Filter(record)) {
         if (fout_eta_p && Filter_eta_p(record)) {
            *fout_eta_p << record;
            NEvents_written_eta_p++;
         }
         if (fout_eta_p_pi0 && Filter_eta_p_pi0(record)) {
            *fout_eta_p_pi0 << record;
            NEvents_written_eta_p_pi0++;
         }
         if (fout_eta_n_pip && Filter_eta_n_pip(record)) {
            *fout_eta_n_pip << record;
            NEvents_written_eta_n_pip++;
         }
         if (fout_eta_p_gamma && Filter_eta_p_gamma(record)) {
            *fout_eta_p_gamma << record;
            NEvents_written_eta_p_gamma++;
         }

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
   delete ofslist[fout];

   std::cout << std::endl << "FINAL:" << std::endl;
   std::cout << " " << NEvents_read << " events read -- "
             << NEvents_written << " events written" << std::endl;
   double frac = (double)NEvents_written/(double)NEvents_read;
   std::cout << "Output eta file has " << 100.0*frac 
             << "%.  (about " << frac*122.0 << " microbarns)" << std::endl;

   if (fout_eta_p) {
      delete fout_eta_p;
      delete ofslist[fout_eta_p];
      double frac = (double)NEvents_written_eta_p/(double)NEvents_read;
      std::cout << "Output eta p file has       " << 100.0*frac 
                << "%.  (about " << frac*122000.0 << " nb)" << std::endl;
   }
   if (fout_eta_p_pi0) {
      delete fout_eta_p_pi0;
      delete ofslist[fout_eta_p_pi0];
      double frac = (double)NEvents_written_eta_p_pi0/(double)NEvents_read;
      std::cout << "Output eta p pi0 file has   " << 100.0*frac
                << "%.  (about " << frac*122000.0 << " nb)" << std::endl;
   }
   if (fout_eta_n_pip) {
      delete fout_eta_n_pip;
      delete ofslist[fout_eta_n_pip];
      double frac = (double)NEvents_written_eta_n_pip/(double)NEvents_read;
      std::cout << "Output eta n pi+ file has   " << 100.0*frac 
                << "%.  (about " << frac*122000.0 << " nb)" << std::endl;
   }
   if (fout_eta_p_gamma) {
      delete fout_eta_p_gamma;
      delete ofslist[fout_eta_p_gamma];
      double frac = (double)NEvents_written_eta_p_gamma/(double)NEvents_read;
      std::cout << "Output eta p gamma file has " << 100.0*frac 
                << "%.  (about " << frac*122000.0 << " nb)" << std::endl;
   }

   return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{

   for (int i=1; i < narg; i++) {
      char *ptr = argv[i];
      
      if (ptr[0] == '-') {
         switch(ptr[1]) {
          case 'h': 
            Usage();
            break;
          case 'c':
            CREATE_SINGLE_CHANNEL_FILES=true;
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
   std::cout << std::endl << "Usage:" << std::endl;
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
   std::cout << "    -c       Create single channel HDDM files as well." 
             << std::endl;
   std::cout << std::endl;

   exit(0);
}

//-----------------------------------------------------------------
// ctrlCHandle
//-----------------------------------------------------------------
void ctrlCHandle(int x)
{
   QUIT++;
   std::cerr << std::endl << "SIGINT received (" << QUIT << ")....." 
             << std::endl;
}

//-----------------------------------------------------------------
// NewHDDMFile
//-----------------------------------------------------------------
hddm_s::ostream* NewHDDMFile(std::string fname)
{
   // Open a new HDDM output file with the specified filename.
   // Print an error message and exit immediately if there is
   // a problem.

   std::ofstream *ofs = new std::ofstream(fname.c_str());
   if (! ofs->is_open()) {
      std::cout << " Error opening output file \"" << fname << "\"!" 
                << std::endl;
      exit(-1);
   }
   std::cout << " output file: " << fname << std::endl;
   hddm_s::ostream *fout = new hddm_s::ostream(*ofs);
   ofslist[fout] = ofs;
   
   return fout;
}

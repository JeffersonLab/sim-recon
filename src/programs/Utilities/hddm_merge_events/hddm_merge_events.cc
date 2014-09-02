// $Id: mcsmear.cc 2388 2007-01-10 16:46:03Z davidl $
//
// Created June 22, 2005  David Lawrence

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <stdlib.h>
#include <signal.h>
#include <time.h>

#include "HDDM/hddm_s.hpp"

void Smear(hddm_s::HDDM *record);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

std::vector<int> Nevents_to_merge;
std::vector<bool> loop_source;
std::vector<char*> INFILENAMES;
char *OUTFILENAME = NULL;
int QUIT = 0;
int MAXEVENTS=-1;
bool HDDM_USE_COMPRESSION=false;
bool HDDM_USE_INTEGRITY_CHECKS=false;


#define _DBG_ std::cout << __FILE__ << ":" << __LINE__ << " "
#define _DBG__ std::cout << __FILE__ << ":" << __LINE__ << std::endl


//-----------
// main
//-----------
int main(int narg,char* argv[])
{
   // Set up to catch SIGINTs for graceful exits
   signal(SIGINT,ctrlCHandle);

   ParseCommandLineArguments(narg, argv);
   
   // Dummy check
   if (Nevents_to_merge.size() != INFILENAMES.size()) {
      _DBG_ << "Size of Nevents_to_merge and INFILENAMES vectors"
               " not the same!" << std::endl;
      _DBG_ << "This indicates a bug in the program. Exiting ..."
            << std::endl;
      exit(-1);
   }
   if (Nevents_to_merge.size() != loop_source.size()) {
      _DBG_ << "Size of Nevents_to_merge and loop_source vectors"
               " not the same!" << std::endl;
      _DBG_ << "This indicates a bug in the program. Exiting ..."
            << std::endl;
      exit(-1);
   }
   
   // Open Input file(s)
   std::vector<std::ifstream*> infiles;
   std::vector<hddm_s::istream*> instreams;
   for (unsigned int i=0; i < INFILENAMES.size(); i++) {
      std::cout << " input file: " << INFILENAMES[i] 
                << "  (" << Nevents_to_merge[i] << " events)"
                << std::endl;
      std::ifstream *ifs = new std::ifstream(INFILENAMES[i]);
      if (! ifs->is_open()) {
         std::cout << " Error opening input file \"" << INFILENAMES[i]
                   << "\"!" << std::endl;
         exit(-1);
      }
      hddm_s::istream *fin = new hddm_s::istream(*ifs);
      instreams.push_back(fin);
      infiles.push_back(ifs);
   }
   
   // Output file
   std::cout << " output file: " << OUTFILENAME << std::endl;
   std::ofstream *ofs = new std::ofstream(OUTFILENAME);
   if (! ofs->is_open()) {
      std::cout << " Error opening output file \"" << OUTFILENAME
                << "\"!" << std::endl;
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
   
   // Loop over events from input files, interleaving them 
   // into the output file until one of the inputs does not
   // have enough events to make an output
   int NEvents = 0;
   int NEvents_read = 0;
   time_t last_time = time(NULL);
   while (true) {
      hddm_s::HDDM record;
      
      // Loop over inputs
      bool done=false;
      for (unsigned int i=0; i < instreams.size(); i++) {
         for (unsigned int j=0; j < (unsigned int)Nevents_to_merge[i]; j++) {
            if (! infiles[i]->good()) {
               // Looks like this source is out of events. Check if we need
               // to re-open this input to keep reading events from it.
               if (loop_source[i]) {
                  delete instreams[i];
                  delete infiles[i];
                  infiles[i] = new std::ifstream(INFILENAMES[i]);
                  instreams[i] = new hddm_s::istream(*infiles[i]);
                  std::cout << std::endl << "Reopened \"" << INFILENAMES[i] 
                            << "\" ..." << std::endl;
               }
               else {
                  done = true;
                  break;
               }
            }
            *instreams[i] >> record;
            *fout << record;
            NEvents_read++;
            NEvents++;
         }
         if (done)
            break;
      }
      if (done)
         break;
      
      // Update ticker
      time_t now = time(NULL);
      if (now != last_time) {
         std::cout << "  " << NEvents_read << " events read     ("
                   << NEvents << " event written) \r";
         std::cout.flush();
         last_time = now;
      }

      if (QUIT || (MAXEVENTS > 0 && NEvents >= MAXEVENTS))
         break;
   }
   
   // Close all inputs
   for (unsigned int i=0; i < instreams.size(); i++) {
      delete instreams[i];
      delete infiles[i];
   }
   delete fout;
   delete ofs;
   
   std::cout << " " << NEvents << " events read" << std::endl;

   return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
   int Nmerge = 1;
   bool loop = false;

   INFILENAMES.clear();
   loop_source.clear();
   Nevents_to_merge.clear();

   for (int i=1; i<narg; i++) {
      char *ptr = argv[i];
      
      if (ptr[0] == '-') {
         switch(ptr[1]) {
            case 'h':
               Usage();
               break;
            case 'N':
               Nmerge=atoi(&ptr[2]);
               break;
            case 'M':
               MAXEVENTS=atoi(&ptr[2]);
               break;
            case 'o':
               OUTFILENAME=&ptr[2];
               break;
            case 'l':
               loop=true;
               break;
            case 's':
               loop=false;
               break;
            case 'I':
               HDDM_USE_INTEGRITY_CHECKS=true;
               break;
            case 'C':
               HDDM_USE_COMPRESSION=true;
               break;
         }
      }
      else {
         INFILENAMES.push_back(argv[i]);
         loop_source.push_back(loop);
         Nevents_to_merge.push_back(Nmerge);
      }
   }

   if (INFILENAMES.size() == 0) {
      std::cout << std::endl << "You must enter a filename!" 
                << std::endl << std::endl;
      Usage();
   }
   
   if (OUTFILENAME == NULL) {
      char *default_outfile = (char*)malloc(20);
      OUTFILENAME = strcpy(default_outfile,"merged.hddm");
   }
}


//-----------
// Usage
//-----------
void Usage(void)
{
   std::cout << std::endl << "Usage:" << std::endl;
   std::cout << "     hddm_merge_events [-Mmaxevents] [-oOutputfile]"
                " [-l|s -Nnum] file1.hddm [-l|s -Nnum] file2.hddm ..."
             << std::endl;
   std::cout << std::endl;
   std::cout << "options:" << std::endl;
   std::cout << "    -oOutputfile  Set output filename (def. merged.hddm)"
             << std::endl;
   std::cout << "    -Mmaxevents   Do not write out more than maxevents"
                " events" << std::endl;
   std::cout << "    -Nnum         Merge together num events from each of"
                " the following input files" << std::endl;
   std::cout << "    -l            Continually loop over events from the"
                " following inputs" << std::endl;
   std::cout << "    -s            Stop when all events from the following"
                " inputs have been read" << std::endl;
   std::cout << "    -I            Enable data integrity checks on output"
                " hddm stream" << std::endl;
   std::cout << "    -C            Enable compression on output hddm stream"
             << std::endl;
   std::cout << std::endl;
   std::cout << " This will take events from 1 or more HDDM files"
                " and merge them" << std::endl;
   std::cout << " into events written to an output HDDM file." 
                "  This works by simply" << std::endl;
   std::cout << " copying the hddm objects into a single event that gets"
             << std::endl;
   std::cout << " written to the output file." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " Multiple input files can be specified and for each, the"
                " number of" << std::endl;
   std::cout << " events to merge. This allows one to do things such as merge"
                " 3 events" << std::endl;
   std::cout << " from a file of background events with a single event from"
                " a file of" << std::endl;
   std::cout << " signal events to get a sample of events with 3 times the "
                " background" << std::endl;
   std::cout << " rate overlayed on the signal." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " By specifying the -l flag, certain input sources can be"
                " looped over" << std::endl;
   std::cout << " to recycle thier events. This is useful if one has a"
                 " sample of" << std::endl;
   std::cout << " background events that is smaller than the sample of signal"
                " events." << std::endl;
   std::cout << " See the examples below." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " NOTE: Events are not merged such that double hits are" 
             << std::endl;
   std::cout << " merged. Hits are only copied so it is possible for" 
             << std::endl;
   std::cout << " the output event to have two hits on the same wire" 
             << std::endl;
   std::cout << " at the same time." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " Example 1:" << std::endl;
   std::cout << "     hddm_merge_events -N2 file.hddm" << std::endl;
   std::cout << " " << std::endl;
   std::cout << "  This will combine every 2 events from file.hddm into a"
                " single in the" << std::endl;
   std::cout << " output. Since no output file name is specified, the file"
                " \"merged.hddm\"" << std::endl;
   std::cout << " will be used." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " " << std::endl;
   std::cout << " Example 2:" << std::endl;
   std::cout << "     hddm_merge_events signal.hddm -N10 em_bkgnd.hddm"
                " -N1 hadronic_bkgnd.hddm" << std::endl;
   std::cout << " " << std::endl;
   std::cout << " This will combine 1 event from signal.hddm, 10 events from"
                " em_bkgnd.hddm," << std::endl;
   std::cout << " and 1 event from hadronic_bkgnd.hdm into a single output"
                " event." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " " << std::endl;
   std::cout << " Example 3:" << std::endl;
   std::cout << "     hddm_merge_events -N3 file1.hddm file2.hddm file3.hddm"
                " -N1 file4.hddm" << std::endl;
   std::cout << " " << std::endl;
   std::cout << " This will combine 3 events from each of file1.hddm,"
                " file2.hddm, and" << std::endl;
   std::cout << " file3.hddm with 1 event from file4.hddm." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " " << std::endl;
   std::cout << " Example 4:" << std::endl;
   std::cout << "     hddm_merge_events signal.hddm -l -N10 background.hddm" 
             << std::endl;
   std::cout << " " << std::endl;
   std::cout << " This will combine 1 signal event with 10 background"
                " events to create" << std::endl;
   std::cout << " an output event. The events in the background file"
                " will be looped" << std::endl;
   std::cout << " over continuosly as needed until all of the signal events"
                " have been" << std::endl;
   std::cout << " processed." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " " << std::endl;
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

// $Id$
//
// Created Dec 22, 2007  David Lawrence

#include "hddm_cull_events.h"

void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

string HDDM_CLASS = "s";
vector<char*> INFILENAMES;
char *OUTFILENAME = NULL;
int QUIT = 0;
unsigned int EVENTS_TO_SKIP = 0;
unsigned int EVENTS_TO_KEEP = 1;
unsigned int SPECIFIC_OFFSET_TO_KEEP = 0;
unsigned int SPECIFIC_EVENT_TO_KEEP = 0;
bool EVENT_TO_KEEP_MODE = false;
bool HDDM_USE_COMPRESSION = false;
bool HDDM_USE_INTEGRITY_CHECKS = false;


//-----------
// main
//-----------
int main(int narg,char* argv[])
{
   // Set up to catch SIGINTs for graceful exits
   signal(SIGINT,ctrlCHandle);

   ParseCommandLineArguments(narg, argv);
   
   std::cout << "Skipping " << EVENTS_TO_SKIP << std::endl;
   std::cout << "Keeping  " << EVENTS_TO_KEEP << std::endl;
   unsigned int NEvents = 0;
   unsigned int NEvents_read = 0;

   // Each HDDM class must have it's own cull routine
   if (HDDM_CLASS == "s") {
      Process_s(NEvents, NEvents_read);
   }
   else if (HDDM_CLASS == "r") {
      Process_r(NEvents, NEvents_read);
   }
   else {
      std::cout << "Don't know how to process HDDM class \"" << HDDM_CLASS << "\"!" << std::endl;
      return -1;
   }
   
   std::cout << std::endl;
   std::cout << " " << NEvents_read << " events read, " << NEvents << " events written" << std::endl;

   return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
  INFILENAMES.clear();
  
  for (int i=1; i < narg; i++) {
    char *ptr = argv[i];
    
    if (ptr[0] == '-') {
      switch(ptr[1]) {
       case 'h':
         Usage();
         break;
      case 'o': 
         OUTFILENAME=&ptr[2];
         break;
      case 's':
         EVENTS_TO_SKIP=atoi(&ptr[2]);
         break;
      case 'k':
         EVENTS_TO_KEEP=atoi(&ptr[2]);
         break;
      case 'e':
         SPECIFIC_OFFSET_TO_KEEP=atoi(&ptr[2]);
         break;
      case 'E':
         SPECIFIC_EVENT_TO_KEEP=atoi(&ptr[2]);
         EVENT_TO_KEEP_MODE=true;
         break;
      case 'r':
         HDDM_CLASS = "r";
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
    }
  }
  
  if (INFILENAMES.size() == 0) {
    std::cout << std::endl << "You must enter a filename!" 
              << std::endl << std::endl;
    Usage();
  }
  
  if (SPECIFIC_OFFSET_TO_KEEP>0) {
    EVENTS_TO_KEEP=1;
    EVENTS_TO_SKIP=SPECIFIC_OFFSET_TO_KEEP-1;
  }
  
  if (SPECIFIC_EVENT_TO_KEEP>0) {
    EVENTS_TO_KEEP=1;
    EVENTS_TO_SKIP=1000000000; // Large number of events to read in while looking for the specified event
  }
  
  if (OUTFILENAME == NULL) {
    if (SPECIFIC_OFFSET_TO_KEEP>0) {
      OUTFILENAME = new char[256];
      sprintf(OUTFILENAME,"evt%d.hddm", SPECIFIC_OFFSET_TO_KEEP);
    }
    else {
      OUTFILENAME = (char*)"culled.hddm";
    }
  }
}

//-----------
// Usage
//-----------
void Usage(void)
{
   std::cout << std::endl << "Usage:" << std::endl;
   std::cout << "     hddm_cull_events [-oOutputfile] [-sNeventsToSkip] [-kNeventsToKeep] file1.hddm file2.hddm ..." << std::endl;
   std::cout << std::endl;
   std::cout << "options:" << std::endl;
   std::cout << "    -oOutputfile     Set output filename (def. culled.hddm)" << std::endl;
   std::cout << "    -sNeventsToSkip  Set number of events to skip (def. 0)" << std::endl;
   std::cout << "    -kNeventsToKeep  Set number of events to keep (def. 1)" << std::endl;
   std::cout << "    -eSingleEvent    Keep only the single, specified event (file pos.)" << std::endl;
   std::cout << "    -ESingleEvent    Keep only the single, specified event (event number)" << std::endl;
   std::cout << "    -r               Input file is in REST format (def. hdgeant format)" << std::endl;
   std::cout << "    -I               Enable data integrity checks on output HDDM stream" << std::endl;
   std::cout << "    -C               Enable compression on output HDDM stream" << std::endl;
   std::cout << std::endl;
   std::cout << " This will copy a continguous set of events from the combined event streams" << std::endl;
   std::cout << " into a seperate output file. The primary use for this would be to copy" << std::endl;
   std::cout << " a single, problematic event into a seperate file for easier debugging." << std::endl;
   std::cout << std::endl;
   std::cout << " If the -eNNN option is used then only a single event is extracted" << std::endl;
   std::cout << " (the NNN-th event) and written to a file with the name evtNNN.hddm." << std::endl;
   std::cout << " Note that the event is counted from the begining of the file, starting" << std::endl;
   std::cout << " with \"1\". This does NOT look at the event number stored in the event itself." << std::endl;
   std::cout << " " << std::endl;
   std::cout << " If the -ENNN option is used then only a single event is extracted" << std::endl;
   std::cout << " (the specified event number) and written to a file with the name evtNNN.hddm." << std::endl;
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

// $Id$
//
// Created Dec 22, 2007  David Lawrence

#include "hddm_merge_files.h"

void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

string HDDM_CLASS = "s";
vector<char*> INFILENAMES;
char *OUTFILENAME = NULL;
int QUIT = 0;
bool HDDM_USE_COMPRESSION = false;
bool HDDM_USE_INTEGRITY_CHECKS = false;


//-----------
// main
//-----------
int main(int narg,char* argv[])
{
   // Set up to catch SIGINTs for graceful exits
   signal(SIGINT, ctrlCHandle);

   ParseCommandLineArguments(narg, argv);
   
   unsigned int NEvents = 0;
   unsigned int NEvents_read = 0;

   // Each HDDM class must have it's own cull routine
   if (HDDM_CLASS == "s")
      Process_s(NEvents, NEvents_read);
   else if (HDDM_CLASS == "r")
      Process_r(NEvents, NEvents_read);
   else {
      std::cout << "Don't know how to process HDDM class \"" << HDDM_CLASS 
                << "\"!" << std::endl;
      return -1;
   }

   std::cout << std::endl;
   std::cout << " " << NEvents_read << " events read, " 
             << NEvents << " events written" << std::endl;
   return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
   INFILENAMES.clear();

   for (int i=1; i<narg; i++) {
      char *ptr = argv[i];
      
      if (ptr[0] == '-') {
         switch(ptr[1]) {
            case 'h':
               Usage();
               break;
            case 'o':
               OUTFILENAME=&ptr[2];
               break;
            case 'r':
               HDDM_CLASS = "r";
               break;
            case 'C':
               HDDM_USE_COMPRESSION = true;
               break;
            case 'I':
               HDDM_USE_INTEGRITY_CHECKS = true;
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
   
   if (OUTFILENAME == NULL) {
      OUTFILENAME = new char[256];
      sprintf(OUTFILENAME,"merged_files.hddm");
   }
}


//-----------
// Usage
//-----------
void Usage(void)
{
   std::cout << std::endl << "Usage:" << std::endl;
   std::cout << "     hddm_merge_files [options] "
                "file1.hddm file2.hddm ..." << std::endl;
   std::cout << std::endl;
   std::cout << "options:" << std::endl;
   std::cout << "    -oOutputfile  Set output filename "
             << "(def. merged_files.hddm)" << std::endl;
   std::cout << "    -I            Enable data integrity checks on"
                " the output hddm stream" << std::endl;
   std::cout << "    -C            Enable data compression on"
                " the output hddm stream" << std::endl;
   std::cout << std::endl;
   std::cout << " This will merge 1 or more HDDM files "
                "into a single HDDM file." << std::endl;
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

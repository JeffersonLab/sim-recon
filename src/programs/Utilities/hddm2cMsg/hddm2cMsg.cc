// $Id$
//
// Created June 22, 2005  David Lawrence

#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

#include <stdlib.h>
#include <signal.h>
#include <time.h>

#if CMSG_PACKAGE_INSTALLED
#include <cMsg.hxx>
#endif

#include <HDDM/hddm_s.hpp>

void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

char *INFILENAME = NULL;
char UDL[] = "cMsg:cMsg//localhost:3456";
int QUIT = 0;

//-----------
// main
//-----------
int main(int narg,char* argv[])
{
   // Set up to catch SIGINTs for graceful exits
   signal(SIGINT, ctrlCHandle);

   ParseCommandLineArguments(narg, argv);
   
   std::cout << " input file: " << INFILENAME << std::endl;
   std::cout << " output UDL: " << UDL << std::endl;
   
#if CMSG_PACKAGE_INSTALLED
   // Open cMsg output
   cMsg cMsgSys(UDL, "hddm2cMsg", "Feed events from HDDM file to cMsg");
   try {
      std::cout << "Connecting to cMsg system ..." << std::endl;
      cMsgSys.connect();
   }
   catch (cMsgException e) {
      std::cout << e.toString() << std::endl;
   }
   
   cMsgSys.disconnect();
   exit(0);
#endif

   // Open Input file
   std::ifstream ifs(INFILENAME);
   if (! ifs.is_open()) {
      std::cout << " Error opening input file \"" << INFILENAME 
                << "\"!" << std::endl;
      exit(-1);
   }
   hddm_s::istream *fin = new hddm_s::istream(ifs);

   // Loop over events in input file
   int NEvents = 0;
   hddm_s::HDDM record;
   time_t last_time = time(NULL);
   while (ifs.good()) {
      *fin >> record;
      NEvents++;
      time_t now = time(NULL);
      if (now != last_time) {
         std::cout << "  " << NEvents << " events processed      \r";
         std::cout.flush();
         last_time = now;
      }
      if (QUIT)
         break;
   }
   
   // close input and output files
   delete fin;
   ifs.close();

   std::cout << " " << NEvents << " events read" << std::endl;

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
         }
      }
      else {
         INFILENAME = argv[i];
      }
   }

   if (!INFILENAME) {
      std::cout << std::endl << "You must enter a filename!" 
                << std::endl << std::endl;
      Usage();
   }
}


//-----------
// Usage
//-----------
void Usage(void)
{
   std::cout << std::endl << "Usage:" << std::endl;
   std::cout << "     hddm2cMsg [options] file.hddm" << std::endl;
   std::cout << std::endl;
   std::cout << " Read the given, Geant-produced HDDM file as input and send"
             << std::endl;
   std::cout << "the events to the specified cMsg UDL." << std::endl;
   std::cout << std::endl;
   std::cout << "  options:" << std::endl;
   std::cout << "    -h       Print this usage statement." << std::endl;
   std::cout << std::endl;
   std::cout << " Example:" << std::endl;
   std::cout << std::endl;
   std::cout << "     hddm2cMsg hdgeant.hddm" << std::endl;
   std::cout << std::endl;
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

// $Id$
//
// Created Oct 09, 2013  Kei Moriya

#include "hddm_select_events.h"
#include <fstream>

void Usage(void);
void ctrlCHandle(int x);

string INFILENAME  = "";
string OUTFILENAME = "";
bool saveRemainder = false;
string OUTFILENAME_REMAINDER = "";
int selectType = 0;
unsigned int MAX = 100 * 1000 * 1000;
bool debug = false;
string HDDM_CLASS = "s";
int QUIT = 0;
int seed = 0;

bool HDDM_USE_COMPRESSION = false;
bool HDDM_USE_INTEGRITY_CHECKS = false;

TRandom2 *rndm;

//-----------
// main
//-----------
int main(int argc,char* argv[]) {
  // Set up to catch SIGINTs for graceful exits
  signal(SIGINT,ctrlCHandle);

  extern char* optarg;
  // Check command line arguments
  int c;
  while ((c = getopt(argc,argv,"ho:i:ars:M:dR:")) != -1) {
    switch(c) {
    case 'h':
      Usage();
      exit(-1);
      break;
    case 'i':
      // Specify infile name.
      // This will nullify options -d and -c
      INFILENAME = optarg;
      std::cout << "infile name: " << INFILENAME << std::endl;
      break;
    case 'o':
      OUTFILENAME = optarg;
      std::cout << "outfile name: " << OUTFILENAME << std::endl;
      break;
    case 'r':
      HDDM_CLASS = "r";
      break;
    case 'a': 
      saveRemainder = true;
      if (OUTFILENAME.find(".hddm") != std::string::npos) {
        // if we found .hddm, then add "_remainder" in front
        OUTFILENAME_REMAINDER = OUTFILENAME;
        OUTFILENAME_REMAINDER.replace(OUTFILENAME.find(".hddm"),0,"_remainder");
      }
      else {
        OUTFILENAME_REMAINDER = OUTFILENAME + "_remainder";
      }
      std::cout << "remainder outfile name: " << OUTFILENAME_REMAINDER
                << std::endl;
      break;
    case 's':
      selectType = atoi(optarg);
      std::cout << "event selection type: " << selectType << std::endl;
      break;
    case 'M':
      MAX = atoi(optarg);
      std::cout << "maximum number of events: " << MAX << std::endl;
      break;
    case 'd':
      debug = true;
      std::cout << "debug mode" << std::endl;
      break;
    case 'R':
      if (selectType != 4) {
         std::cout << "Random seed is only needed for select type 4" 
                   << std::endl;
         abort();
      }
      seed = atoi(optarg);
      std::cout << "random seed: " << seed << std::endl;
      break;
    case 'C':
      HDDM_USE_COMPRESSION = true;
      break;
    case 'I':
      HDDM_USE_INTEGRITY_CHECKS = true;
      break;
    default:
      break;
    }
  }
  //___________________________________________________________________________________________

  if (INFILENAME == "" || OUTFILENAME == "" || 
      (selectType < 1 || 7 < selectType))
  {
    Usage();
  }

  // if selectType == 4, we need the random generator
  rndm = new TRandom2(seed);
   
  // standard hddm
  hddm_s::istream *istr_s;
  hddm_s::ostream *ostr_s;
  hddm_s::ostream *ostr_s_remainder;
  // REST
  hddm_r::istream *istr_r;
  hddm_r::ostream *ostr_r;
  hddm_r::ostream *ostr_r_remainder;

  ifstream ifs;
  ofstream ofs;
  ofstream ofs_remainder;

  if (HDDM_CLASS == "s") {
    ifs.open(INFILENAME.c_str());
    if (! ifs.is_open()) {
      std::cout << " Error opening input file \"" << INFILENAME 
                << "\"!" << std::endl;
      exit(-1);
    }
    istr_s = new hddm_s::istream(ifs);

    ofs.open(OUTFILENAME.c_str());
    if (! ofs.is_open()) {
      std::cout << " Error opening output file \"" << OUTFILENAME 
                << "\"!" << std::endl;
      exit(-1);
    }
    ostr_s = new hddm_s::ostream(ofs);
    if (HDDM_USE_COMPRESSION) {
       ostr_s->setCompression(hddm_r::k_bz2_compression);
    }
    if (HDDM_USE_INTEGRITY_CHECKS) {
       ostr_s->setIntegrityChecks(hddm_r::k_crc32_integrity);
    }

    if (saveRemainder) {
      ofs_remainder.open(OUTFILENAME_REMAINDER.c_str());
      if (! ofs_remainder.is_open()) {
         std::cout << " Error opening output file \"" << OUTFILENAME_REMAINDER
                   << "\"!" << std::endl;
         exit(-1);
      }
      ostr_s_remainder = new hddm_s::ostream(ofs_remainder);
      if (HDDM_USE_COMPRESSION) {
         ostr_s_remainder->setCompression(hddm_r::k_bz2_compression);
      }
      if (HDDM_USE_INTEGRITY_CHECKS) {
         ostr_s_remainder->setIntegrityChecks(hddm_r::k_crc32_integrity);
      }
    }
    else {
      ostr_s_remainder = NULL;
    }

    // Make sure compiler doesn't give warnings
    istr_r = NULL;
    ostr_r = NULL;
    ostr_r_remainder = NULL;
  }

  else {
    ifs.open(INFILENAME.c_str());
    if (! ifs.is_open()) {
      std::cout << " Error opening input file \"" << INFILENAME 
                << "\"!" << std::endl;
      exit(-1);
    }
    istr_r = new hddm_r::istream(ifs);

    ofs.open(OUTFILENAME.c_str());
    if (! ofs.is_open()) {
      std::cout << " Error opening output file \"" << OUTFILENAME 
                << "\"!" << std::endl;
      exit(-1);
    }
    ostr_r = new hddm_r::ostream(ofs);
    if (HDDM_USE_COMPRESSION) {
       ostr_r->setCompression(hddm_r::k_bz2_compression);
    }
    if (HDDM_USE_INTEGRITY_CHECKS) {
       ostr_r->setIntegrityChecks(hddm_r::k_crc32_integrity);
    }
    if (saveRemainder) {
      ofs_remainder.open(OUTFILENAME_REMAINDER.c_str());
      if (! ofs_remainder.is_open()) {
         std::cout << " Error opening output file \"" << OUTFILENAME_REMAINDER
                   << "\"!" << std::endl;
         exit(-1);
      }
      ostr_r_remainder = new hddm_r::ostream(ofs_remainder);
      if (HDDM_USE_COMPRESSION) {
         ostr_r_remainder->setCompression(hddm_r::k_bz2_compression);
      }
      if (HDDM_USE_INTEGRITY_CHECKS) {
         ostr_r_remainder->setIntegrityChecks(hddm_r::k_crc32_integrity);
      }
    }
    else
      ostr_r_remainder = NULL;

    // Make sure compiler doesn't give warnings
    istr_s = NULL;
    ostr_s = NULL;
    ostr_s_remainder = NULL;
  }

  // Loop over input files
  unsigned int NEvents = 0;
  unsigned int NEvents_read = 0;
  time_t last_time = time(NULL);
         
  // Loop over all events in input
  while (true && NEvents_read < MAX) {

    /////////////////////////////////////////////////////
    //                                                 //
    //  At this stage we have the current event in     //
    //  hddm_s, so we can choose our events with       //
    //  any information that is contained.             //
    //                                                 //
    /////////////////////////////////////////////////////

    if (HDDM_CLASS == "s") {
      if (! ifs.good())
         break;
      hddm_s::HDDM record;
      *istr_s >> record;
      NEvents_read++;
      if (debug)
        std::cout << NEvents_read << std::endl;
      if (selectEvent_s(selectType, record, NEvents_read, debug)) {
        // Write this output event to file
        *ostr_s << record;
        NEvents++;
      }
      else if (saveRemainder)
        *ostr_s_remainder << record;
    }
    else {
      if (! ifs.good())
        break;
      hddm_r::HDDM record;
      *istr_r >> record;
      NEvents_read++;
      if (debug)
        std::cout << NEvents_read << std::endl;
      hddm_r::ReconstructedPhysicsEvent &re =
              record.getReconstructedPhysicsEvent();
      int runno = re.getRunNo();
      int eventno = re.getEventNo();
      if (debug)
        std::cout << runno << "\t" << eventno << std::endl;
      if (selectEvent_r(selectType, record, NEvents_read, debug)) {
         // Write this output event to file 
        *ostr_r << record;
        NEvents++;
      }
      else if (saveRemainder)
        *ostr_r_remainder << record;
    }

    // Update ticker
    time_t now = time(NULL);
    if (now != last_time) {
      std::cout << "  " << NEvents_read << " events read     ("
                << NEvents << " event written) \r";
      std::cout.flush();
      last_time = now;
    }
         
    if (QUIT)
      break;
  }

  if (HDDM_CLASS == "s") {
    // Close input file
    delete istr_s;
    ifs.close();
    
    // Close output files
    delete ostr_s;
    ofs.close();
    if (ofs_remainder.is_open()) {
      delete ostr_s_remainder;
      ofs_remainder.close();
    }
  }
  else {
    delete istr_r;
    ifs.close();
    delete ostr_r;
    ofs.close();
    if (ofs_remainder.is_open()) {
      delete ostr_r_remainder;
      ofs_remainder.close();
    }
  }

  std::cout << std::endl;
  std::cout << " " << NEvents_read << " events read, " 
            << NEvents << " events written" << std::endl;

  return 0;
}

//-----------
// Usage
//-----------
void Usage(void)
{
  std::cout << std::endl << "Usage:" << std::endl;
  std::cout << "     hddm_select_events [-i Inputfile] [-o Outputfile]"
               " [-r REST format] [-a save remainder events too]"
               " [-s selection type] [-M maximum number of events]"
            << std::endl;
  std::cout << std::endl;
  std::cout << "options:" << std::endl;
  std::cout << "    -i Inputfile     Set input  filename" << std::endl;
  std::cout << "    -o Outputfile    Set output filename" << std::endl;
  std::cout << "    -r               Input is REST format" << std::endl;
  std::cout << "    -a               Save remainder events in a file too" 
            << std::endl;
  std::cout << "    -s selectType    Set which type of cut to do" 
            << std::endl;
  std::cout << "    -M MAX           Set maximum number of events"
            << std::endl;
  std::cout << "    -C               Enable compression in the output"
               " hddm streams" << std::endl;
  std::cout << "    -I               Enable data integrity checks in the"
               " output hddm streams" << std::endl;
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

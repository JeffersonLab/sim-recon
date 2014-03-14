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

TRandom2 *rndm;

//-----------
// main
//-----------
int main(int argc,char* argv[]){
  // Set up to catch SIGINTs for graceful exits
  signal(SIGINT,ctrlCHandle);

  extern char* optarg;
  // Check command line arguments
  int c;
  while((c = getopt(argc,argv,"ho:i:ars:M:dR:")) != -1){
    switch(c){
    case 'h':
      Usage();
      exit(-1);
      break;
    case 'i':
      // Specify infile name.
      // This will nullify options -d and -c
      INFILENAME = optarg;
      cout << "infile name: " << INFILENAME << endl;
      break;
    case 'o':
      OUTFILENAME = optarg;
      cout << "outfile name: " << OUTFILENAME << endl;
      break;
    case 'r':
      HDDM_CLASS = "r";
      break;
    case 'a':
      {
	saveRemainder = true;
	size_t suffix_pos = OUTFILENAME.find(".hddm");
	if(suffix_pos!=std::string::npos){
	  // if we found .hddm, then add "_remainder" in front
	  OUTFILENAME_REMAINDER = OUTFILENAME;
	  OUTFILENAME_REMAINDER.replace(suffix_pos,0,"_remainder");
	}else{
	  OUTFILENAME_REMAINDER = OUTFILENAME + "_remainder";
	}
	cout << "remainder outfile name: " << OUTFILENAME_REMAINDER << endl;
      }
      break;
    case 's':
      selectType = atoi(optarg);
      cout << "event selection type: " << selectType << endl;
      break;
    case 'M':
      MAX = atoi(optarg);
      cout << "maximum number of events: " << MAX << endl;
      break;
    case 'd':
      debug = true;
      cout << "debug mode" << endl;
      break;
    case 'R':
      if(selectType!=4){
	cout << "Random seed is only needed for select type 4" << endl;
	abort();
      }
      seed = atoi(optarg);
      cout << "random seed: " << seed << endl;
      break;
    default:
      break;
    }
  }
  //___________________________________________________________________________________________

  if(INFILENAME=="" || OUTFILENAME=="" || (selectType < 1 || 7 < selectType)){
    Usage();
  }

  // if selectType==4, we need the random generator
  rndm = new TRandom2(seed);
	
  // Input/Output file
  char filename[200];

  // standard hddm
  s_iostream_t *fin_s;
  s_iostream_t *fout_s;
  s_iostream_t *fout_s_remainder;
  // REST
  hddm_r::istream *istr_r;
  hddm_r::ostream *ostr_r;
  hddm_r::ostream *ostr_r_remainder;

  ifstream ifs;
  ofstream ofs;
  ofstream ofs_remainder;

  if(HDDM_CLASS=="s"){
    sprintf(filename,"%s",INFILENAME.c_str());
    fin_s = open_s_HDDM(filename);
    if(!fin_s){
      cout << " Error opening input file \"" << INFILENAME << "\"!" << endl;
      exit(-1);
    }

    sprintf(filename,"%s",OUTFILENAME.c_str());
    fout_s = init_s_HDDM(filename);
    if(!fout_s){
      cout << " Error opening output file \"" << OUTFILENAME << "\"!" << endl;
      exit(-1);
    }
  
    if(saveRemainder){
      sprintf(filename,"%s",OUTFILENAME_REMAINDER.c_str());
      fout_s_remainder = init_s_HDDM(filename);
      if(!fout_s_remainder){
	cout << " Error opening output file \"" << OUTFILENAME_REMAINDER << "\"!" << endl;
	exit(-1);
      }
    }else{
      fout_s_remainder = NULL;
    }

    // Make sure compiler doesn't give warnings
    istr_r = NULL;
    ostr_r = NULL;
    ostr_r_remainder = NULL;
  }else{

    ifs.open(INFILENAME.c_str());
    istr_r = new hddm_r::istream(ifs);
    ofs.open(OUTFILENAME.c_str());
    ostr_r = new hddm_r::ostream(ofs);
    ofs_remainder.open(OUTFILENAME_REMAINDER.c_str());
    if(saveRemainder) ostr_r_remainder = new hddm_r::ostream(ofs_remainder);
    else              ostr_r_remainder = NULL;

    // Make sure compiler doesn't give warnings
    fin_s = NULL;
    fout_s = NULL;
    fout_s_remainder = NULL;
  }

  // Loop over input files
  unsigned int NEvents = 0;
  unsigned int NEvents_read = 0;
  time_t last_time = time(NULL);
			
  // Loop over all events in input
  while(true && NEvents_read < MAX){

    bool write_this_event = false;

    /////////////////////////////////////////////////////
    //                                                 //
    //  At this stage we have the current event in     //
    //  hddm_s, so we can choose our events with       //
    //  any information that is contained.             //
    //                                                 //
    /////////////////////////////////////////////////////

    if(HDDM_CLASS=="s"){
      s_HDDM_t *hddm_s = read_s_HDDM(fin_s);
      if(!hddm_s) break;
      NEvents_read++;
      if(debug) cout << NEvents_read << endl;
      write_this_event = selectEvent_s(selectType, hddm_s, NEvents_read, debug);

      // Write this output event to file and free its memory
      if(write_this_event){
	flush_s_HDDM(hddm_s, fout_s);
	NEvents++;
      }else{
	if(saveRemainder) flush_s_HDDM(hddm_s, fout_s_remainder);
	else              flush_s_HDDM(hddm_s, NULL);
      }
    }else{
      hddm_r::HDDM *record = new hddm_r::HDDM();
      if(!ifs.eof() && ifs.good()){
	NEvents_read++;
	if(debug) cout << NEvents_read << endl;
	*istr_r >> *record;
	hddm_r::ReconstructedPhysicsEvent &re = record->getReconstructedPhysicsEvent();
	int runno = re.getRunNo();
	int eventno = re.getEventNo();
	if(debug) cout << runno << "\t" << eventno << endl;
	write_this_event = selectEvent_r(selectType, record, NEvents_read, debug);

	// Write this output event to file and free its memory
	if(write_this_event){
	  *ostr_r << *record;
	  NEvents++;
	}else{
	  if(saveRemainder) *ostr_r_remainder << *record;
	}
      }else{
	break;
      }
    }

    // Update ticker
    time_t now = time(NULL);
    if(now != last_time){
      cout << "  " << NEvents_read << " events read     (" << NEvents << " event written) \r";
      cout.flush();
      last_time = now;
    }
			
    if(QUIT)break;
  }

  if(HDDM_CLASS=="s"){
    // Close input file
    close_s_HDDM(fin_s);
    
    // Close output file
    close_s_HDDM(fout_s);
  }else{
    // Close REST
  }

  cout << endl;
  cout << " " << NEvents_read << " events read, " << NEvents << " events written" << endl;

  return 0;
}

//-----------
// Usage
//-----------
void Usage(void)
{
  cout<<endl<<"Usage:"<<endl;
  cout<<"     hddm_select_events [-i Inputfile] [-o Outputfile] [-r REST format] [-a save remainder events too] [-s selection type] [-M maximum number of events]"<<endl;
  cout<<endl;
  cout<<"options:"<<endl;
  cout<<"    -i Inputfile     Set input  filename"<<endl;
  cout<<"    -o Outputfile    Set output filename"<<endl;
  cout<<"    -r               Input is REST format"<<endl;
  cout<<"    -a               Save remainder events in a file too"<<endl;
  cout<<"    -s selectType    Set which type of cut to do"<<endl;
  cout<<"    -M MAX           Set maximum number of events"<<endl;
  cout<<endl;

  exit(0);
}

//-----------------------------------------------------------------
// ctrlCHandle
//-----------------------------------------------------------------
void ctrlCHandle(int x)
{
  QUIT++;
  cerr<<endl<<"SIGINT received ("<<QUIT<<")....."<<endl;
}

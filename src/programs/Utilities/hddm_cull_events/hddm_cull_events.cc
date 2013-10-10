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



//-----------
// main
//-----------
int main(int narg,char* argv[])
{
	// Set up to catch SIGINTs for graceful exits
	signal(SIGINT,ctrlCHandle);

	ParseCommandLineArguments(narg, argv);
	
	cout<<"Skipping "<<EVENTS_TO_SKIP<<endl;
	cout<<"Keeping  "<<EVENTS_TO_KEEP<<endl;
	unsigned int NEvents = 0;
	unsigned int NEvents_read = 0;

	// Each HDDM class must have it's own cull routine
	if(HDDM_CLASS == "s"){
		Process_s(NEvents, NEvents_read);
	}else if(HDDM_CLASS == "r"){
		Process_r(NEvents, NEvents_read);
	}else{
		cout << "Don't know how to process HDDM class \"" << HDDM_CLASS << "\"!" << endl;
		return -1;
	}
	
	cout<<endl;
	cout<<" "<<NEvents_read<<" events read, "<<NEvents<<" events written"<<endl;

	return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
  INFILENAMES.clear();
  
  for(int i=1; i<narg; i++){
    char *ptr = argv[i];
    
    if(ptr[0] == '-'){
      switch(ptr[1]){
      case 'h': Usage(); break;
      case 'o': OUTFILENAME=&ptr[2]; break;
      case 's': EVENTS_TO_SKIP=atoi(&ptr[2]); break;
      case 'k': EVENTS_TO_KEEP=atoi(&ptr[2]); break;
      case 'e': SPECIFIC_OFFSET_TO_KEEP=atoi(&ptr[2]); break;
      case 'E': SPECIFIC_EVENT_TO_KEEP=atoi(&ptr[2]); EVENT_TO_KEEP_MODE=true; break;
	  case 'r': HDDM_CLASS = "r";
      }
    }else{
      INFILENAMES.push_back(argv[i]);
    }
  }
  
  if(INFILENAMES.size()==0){
    cout<<endl<<"You must enter a filename!"<<endl<<endl;
    Usage();
  }
  
  if(SPECIFIC_OFFSET_TO_KEEP>0){
    EVENTS_TO_KEEP=1;
    EVENTS_TO_SKIP=SPECIFIC_OFFSET_TO_KEEP-1;
  }
  
  if(SPECIFIC_EVENT_TO_KEEP>0){
    EVENTS_TO_KEEP=1;
    EVENTS_TO_SKIP=1000000000; // Large number of events to read in while looking for the specified event
  }
  
  if(OUTFILENAME==NULL){
    if(SPECIFIC_OFFSET_TO_KEEP>0){
      OUTFILENAME = new char[256];
      sprintf(OUTFILENAME,"evt%d.hddm", SPECIFIC_OFFSET_TO_KEEP);
    }else{
      OUTFILENAME = (char*)"culled.hddm";
    }
  }
}

//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     hddm_cull_events [-oOutputfile] [-sNeventsToSkip] [-kNeventsToKeep] file1.hddm file2.hddm ..."<<endl;
	cout<<endl;
	cout<<"options:"<<endl;
	cout<<"    -oOutputfile     Set output filename (def. culled.hddm)"<<endl;
	cout<<"    -sNeventsToSkip  Set number of events to skip (def. 0)"<<endl;
	cout<<"    -kNeventsToKeep  Set number of events to keep (def. 1)"<<endl;
	cout<<"    -eSingleEvent    Keep only the single, specified event (file pos.)"<<endl;
	cout<<"    -ESingleEvent    Keep only the single, specified event (event number)"<<endl;
	cout<<"    -r               Input file is in REST format (def. hdgeant format)"<<endl;
	cout<<endl;
	cout<<" This will copy a continguous set of events from the combined event streams"<<endl;
	cout<<" into a seperate output file. The primary use for this would be to copy"<<endl;
	cout<<" a single, problematic event into a seperate file for easier debugging."<<endl;
	cout<<endl;
	cout<<" If the -eNNN option is used then only a single event is extracted"<<endl;
	cout<<" (the NNN-th event) and written to a file with the name evtNNN.hddm."<<endl;
	cout<<" Note that the event is counted from the begining of the file, starting"<<endl;
	cout<<" with \"1\". This does NOT look at the event number stored in the event itself."<<endl;
	cout<<" "<<endl;
	cout<<" If the -ENNN option is used then only a single event is extracted"<<endl;
	cout<<" (the specified event number) and written to a file with the name evtNNN.hddm."<<endl;
	cout<<" "<<endl;
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

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



//-----------
// main
//-----------
int main(int narg,char* argv[])
{
	// Set up to catch SIGINTs for graceful exits
	signal(SIGINT,ctrlCHandle);

	ParseCommandLineArguments(narg, argv);
	
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
				case 'h': Usage();						break;
				case 'o': OUTFILENAME=&ptr[2];		break;
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
	
	if(OUTFILENAME==NULL){
	        OUTFILENAME = new char[256];
		sprintf(OUTFILENAME,"merged_files.hddm");
	}
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     hddm_merge_files [-oOutputfile] file1.hddm file2.hddm ..."<<endl;
	cout<<endl;
	cout<<"options:"<<endl;
	cout<<"    -oOutputfile  Set output filename (def. merged_files.hddm)"<<endl;
	cout<<endl;
	cout<<" This will merge 1 or more HDDM files into a single HDDM file."<<endl;
	cout<<" "<<endl;
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

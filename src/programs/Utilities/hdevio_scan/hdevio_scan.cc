
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <thread>
using namespace std;

#include <TFile.h>

#include "DMapEVIOWords.h"

#include <DAQ/HDEVIO.h>


void Usage(string mess);
void ParseCommandLineArguments(int narg, char *argv[]);
void PrintSummary(void);
void MapEVIOWords(void);


vector<string> filenames;
bool   PRINT_SUMMARY = true;
bool   MAP_WORDS     = false;
string ROOT_FILENAME = "hdevio_scan.root";
uint64_t MAX_EVIO_EVENTS = 20000;
uint32_t BLOCK_SIZE = 20; // used for daq_block_size histogram

//----------------
// main
//----------------
int main(int narg, char *argv[])
{

	ParseCommandLineArguments(narg, argv);
	
	if(PRINT_SUMMARY) PrintSummary();
	
	if(MAP_WORDS    ) MapEVIOWords();

	return 0;
}

//----------------
// Usage
//----------------
void Usage(string mess="")
{
	cout << endl;
	cout << "Usage:" << endl;
	cout << endl;
	cout <<"    hdevio [options] file.evio [file2.evio ...]" << endl;
	cout << endl;
	cout << "options:" << endl;
	cout << "   -h, --help    Print this usage statement" << endl;
	cout << "   -w            Make histogram of population by word type" << endl;
	cout << "   -r file.root  Set name of ROOT file to save histo to. " << endl;
	cout << "                 (implies -w)" << endl;
	cout << "   -m max_events Max. EVIO events (not physics events) to process." << endl;
	cout << "   -b block_size EVIO events to add for daq_block_size histo." << endl;
	cout << endl;

	if(mess != "") cout << endl << mess << endl << endl;
	
	exit(0);
}

//----------------
// ParseCommandLineArguments
//----------------
void ParseCommandLineArguments(int narg, char *argv[])
{

	if(narg<2) Usage("You must supply a filename!");

	for(int i=1; i<narg; i++){
		string arg  = argv[i];
		string next = (i+1)<narg ? argv[i+1]:"";
		
		if(arg == "-h" || arg == "--help") Usage();
		else if(arg == "-w"){ MAP_WORDS = true; PRINT_SUMMARY = false; }
		else if(arg == "-r"){ MAP_WORDS = true; PRINT_SUMMARY = false; ROOT_FILENAME = next; i++;}
		else if(arg == "-m"){ MAX_EVIO_EVENTS = atoi(next.c_str()); i++;}
		else if(arg == "-b"){ BLOCK_SIZE = atoi(next.c_str()); i++;}
		else if(arg[0] == '-') {cout << "Unknown option \""<<arg<<"\" !" << endl; exit(-1);}
		else filenames.push_back(arg);
	}
}

//----------------
// PrintSummary
//----------------
void PrintSummary(void)
{
	// Loop over input files
	for(uint32_t i=0; i<filenames.size(); i++){
		string &filename = filenames[i];
		cout << "Processing file " << (i+1) << "/" << filenames.size() << " : " << filename << endl;

		HDEVIO *hdevio = new HDEVIO(filename);
		if(!hdevio->is_open){
			cout << hdevio->err_mess.str() << endl;
			continue;
		}
		
		time_t start_time = time(NULL);
		hdevio->PrintFileSummary();
		time_t end_time = time(NULL);
		
		delete hdevio;

		cout << (end_time - start_time) << " sec " << endl;
	}
}

//----------------
// MapEVIOWords
//----------------
void MapEVIOWords(void)
{
	cout << endl;
	cout << "Mapping will be limited to first " << MAX_EVIO_EVENTS << " events per input file" << endl;

	// Open ROOT file
	TFile *rootfile = new TFile(ROOT_FILENAME.c_str(), "RECREATE");

	DMapEVIOWords mapevio;

	// Loop over input files
	uint64_t Nevents = 0;
	for(uint32_t i=0; i<filenames.size(); i++){

		// Open EVIO file
		string &filename = filenames[i];
		cout << "Processing file " << (i+1) << "/" << filenames.size() << " : " << filename << endl;
		HDEVIO *hdevio = new HDEVIO(filename);
		if(!hdevio->is_open){
			cout << hdevio->err_mess.str() << endl;
			continue;
		}
		
		// Read all events in file
		uint32_t buff_len = 1000;
		uint32_t *buff = new uint32_t[buff_len];
		bool done = false;
		while(!done){
			hdevio->readNoFileBuff(buff, buff_len);
			switch(hdevio->err_code){
				case HDEVIO::HDEVIO_OK:
					mapevio.ParseEvent(buff);
					break;
				case HDEVIO::HDEVIO_USER_BUFFER_TOO_SMALL:
					buff_len = hdevio->last_event_len;
					delete buff;
					buff = new uint32_t[buff_len];
					break;
				case HDEVIO::HDEVIO_EOF:
					cout << endl << " end of file" << endl;
					done = true;
					break;
				default:
					cout << endl;
					cout << hdevio->err_mess.str() << endl;
					done = true;
					break;
				
			}
			
			if((++Nevents % 1000) == 0) {
				int percent_done = (100*Nevents)/MAX_EVIO_EVENTS;
				cout << " " << Nevents << "/" << MAX_EVIO_EVENTS << " (" << percent_done << "%) processed       \r"; cout.flush();
			}
			if(Nevents>MAX_EVIO_EVENTS) break;
		}
		cout << endl;

		// Close EVIO file
		delete hdevio;
	}
	
	// Flush and close ROOT file
	rootfile->Write();
	rootfile->Close();
	delete rootfile;
}




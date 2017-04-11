
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
bool   SAVE_FILE_MAP = false;
bool   SKIP_EVENT_MAPPING = false;
bool   MAP_WORDS     = false;
bool   GENERATE_ERROR_REPORT = false;
string ROOT_FILENAME = "hdevio_scan.root";
string MAP_FILENAME = "";
uint64_t MAX_EVIO_EVENTS = 20000;
uint32_t BLOCK_SIZE = 20; // used for daq_block_size histogram
uint32_t MAX_HISTORY_BUFF_SIZE = 400;

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
	cout << "   -e            Write details of bad event tag location to file" << endl;
	cout << "   -n max_buff   Max. events to keep timing info for (only valid)" << endl;
	cout << "                 with -w option)" << endl;
	cout << "   -s            Save file block/event map" << endl;
	cout << "   -blocksonly   Save only block map not events. (Only use with -s)" << endl;
	cout << "   -f file.map   Set name of file to save block/event to. " << endl;
	cout << "                 (implies -s)" << endl;
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
		else if(arg == "-e"){ GENERATE_ERROR_REPORT = true; }
		else if(arg == "-n"){ MAX_HISTORY_BUFF_SIZE = atoi(next.c_str()); i++;}		
		else if(arg == "-s"){ SAVE_FILE_MAP = true;}
		else if(arg == "-f"){ SAVE_FILE_MAP = true; MAP_FILENAME = next; i++;}
		else if(arg == "-blocksonly") { SKIP_EVENT_MAPPING = true;}
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
		
		if(SKIP_EVENT_MAPPING) hdevio->SKIP_EVENT_MAPPING = true;
		
		time_t start_time = time(NULL);
		hdevio->PrintFileSummary();
		time_t end_time = time(NULL);
		
		if(GENERATE_ERROR_REPORT){
			ofstream ofs("hdevio_scan.err");
			ofs << "#";
			ofs << "# hdevio_scan report for " << filename << endl;
			ofs << "#";
			ofs << "# The following list is for each EVIO event that contained an" << endl;
			ofs << "# unknown top-level bank tag. Since any corruption may have" << endl;
			ofs << "# started in the previous event, locations for that are also" << endl;
			ofs << "# provided (when available)." <<endl;
			ofs << "# columns are:" << endl;
			ofs << "#  1    block number in file (starting from 0)" << endl;
			ofs << "#  2    block offset in file (hex)" <<endl;
			ofs << "#  3    block number of previous event" << endl;
			ofs << "#  4    block offset of previous event" <<endl;
			ofs << "#  5    event number in block (starting from 1)" << endl;
			ofs << "#  6    number of events in block" << endl;
			ofs << "#  7    event offset in file (hex)" <<endl;
			ofs << "#  8    event length (inclusive)" <<endl;
			ofs << "#  9    event header (hex) <-- This value is what indicates a corrupt event" <<endl;			
			ofs << "# 10    event number of previous event" << endl;
			ofs << "# 11    event offset of previous event" <<endl;
			ofs << "# 12    event length of previous event" <<endl;
			ofs << "# 13    event header of previous event" <<endl;
			ofs << "#" << endl;
			ofs << "# 1      2        3      4        5   6      7          8          9      10     11         12       13" <<endl;

			vector<HDEVIO::EVIOBlockRecord> brs = hdevio->GetEVIOBlockRecords();
			
			HDEVIO::EVIOBlockRecord *br_prev = NULL;
			HDEVIO::EVIOEventRecord *er_prev = NULL;
			int32_t ibr = 0;
			for(HDEVIO::EVIOBlockRecord &br : brs){
				int32_t ier = 1;
				for(HDEVIO::EVIOEventRecord &er : br.evio_events){
					if(er.event_type == HDEVIO::kBT_UNKNOWN){
						char str[512];
						//              1    2      3    4     5   6    7     8    9     10    11    12   13
						sprintf(str, "%04u 0x%08x %04d 0x%08x %3u/%3u 0x%08x %10u 0x%08x %03d 0x%08x %10u 0x%08x"
							/*  1 */ , ibr
							/*  2 */ , (unsigned int)br.pos
							/*  3 */ , ibr-1
							/*  4 */ , br_prev!=NULL ? (unsigned int)br_prev->pos:0
							/*  5 */ , ier
							/*  6 */ , (unsigned int)br.evio_events.size()
							/*  7 */ , (unsigned int)er.pos
							/*  8 */ , (unsigned int)er.event_len
							/*  9 */ , (unsigned int)er.event_header
							/* 10 */ , ier-1
							/* 11 */ , er_prev!=NULL ? (unsigned int)er_prev->pos:0
							/* 12 */ , er_prev!=NULL ? (unsigned int)er_prev->event_len:0
							/* 13 */ , er_prev!=NULL ? (unsigned int)er_prev->event_header:0);
						ofs << str << endl;							
					}
					
					er_prev = &er;
					ier++;
				}
				
				br_prev = &br;
				ibr++;
			}
			
			ofs.close();
		}
		
		if(SAVE_FILE_MAP) hdevio->SaveFileMap(MAP_FILENAME);
		
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
	mapevio.max_history_buff_size = MAX_HISTORY_BUFF_SIZE;

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
					delete[] buff;
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




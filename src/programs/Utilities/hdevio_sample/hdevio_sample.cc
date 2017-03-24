
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

#include <DAQ/HDEVIO.h>


void Usage(string mess);
void ParseCommandLineArguments(int narg, char *argv[]);


vector<string> filenames;

bool   KEEP_BOR   = true;
bool   KEEP_EPICS = false;
bool   KEEP_CODA  = false;
string   OUTPUT_FILENAME = "hdevio_pare.evio";
 uint64_t PRESCALE = 100;

//----------------
// main
//----------------
int main(int narg, char *argv[])
{

	ParseCommandLineArguments(narg, argv);
	
	// Open output file
	ofstream ofs(OUTPUT_FILENAME);
	if(!ofs.is_open()){
		_DBG_ << "Unable to open " << OUTPUT_FILENAME << " for writing!" << endl;
		exit(-1);
	}
	
	// Create buffer for reading blocks into.
	uint32_t *buff = new uint32_t[1000000];
	
	// Loop over input files
	for(auto filename : filenames){

		// Open input file
		HDEVIO *hdevio = new HDEVIO(filename);
		if(!hdevio->is_open){
			cout << hdevio->err_mess.str() << endl;
			continue;
		}
		
		// Get file map
		vector<HDEVIO::EVIOBlockRecord> brs = hdevio->GetEVIOBlockRecords();
		
		// Have HDEVIO close file and open it ourselves
		delete hdevio;
		ifstream ifs(filename);

		// Loop over blocks
		uint64_t idx = 0;
		for(auto &br : brs){

			bool write_block = false;
			switch(br.block_type){
				case HDEVIO::kBT_BOR  : if(KEEP_BOR  ) write_block = true; break;
				case HDEVIO::kBT_EPICS: if(KEEP_EPICS) write_block = true; break;
				case HDEVIO::kBT_PRESTART:
				case HDEVIO::kBT_GO:
				case HDEVIO::kBT_PAUSE:
				case HDEVIO::kBT_END:
					if(KEEP_CODA) write_block = true;
					break;
				case HDEVIO::kBT_PHYSICS:
				default:
					break;
			}
			
			// Any events we're not keeping, prescale
			if(!write_block){
				write_block = ((idx++)%PRESCALE) == 0;
			}
			
			if(write_block){
				ifs.seekg(br.pos, ios::beg);
				int nbytes = br.block_len*4;
				ifs.read((char*)buff, nbytes);
				if( ifs.gcount() != nbytes ){
					_DBG_<<"Unable to read entire block (is file truncated?)" << endl;
				}
				
				ofs.write((char*)buff, nbytes);
			}
		}
		
		ifs.close();
	}
	
	ofs.close();
	if(buff) delete buff;

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
	cout <<"    hdevio_pare [options] file.evio [file2.evio ...]" << endl;
	cout << endl;
	cout << "options:" << endl;
	cout << "   -h, --help    Print this usage statement" << endl;
	cout << "   -o file.evio  Set name of EVIO ouput file" << endl;
	cout << "   -p prescale   Prescale factor for EVIO events (not L1 trigger events)" << endl;
	cout << "   -bor          Don't save BOR events" << endl;
	cout << "   -epics        Don't save EPICS events" << endl;
	cout << "   -coda         Don't save CODA control events" << endl;
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
		else if(arg == "-o"    ){ OUTPUT_FILENAME = next.c_str(); i++;}
		else if(arg == "-p"    ){ PRESCALE   = atoi(next.c_str()); i++;}
		else if(arg == "-bor"  ){ KEEP_BOR   = false;}
		else if(arg == "-epics"){ KEEP_EPICS = false;}		
		else if(arg == "-coda" ){ KEEP_CODA  = false;}		
		else if(arg[0] == '-') {cout << "Unknown option \""<<arg<<"\" !" << endl; exit(-1);}
		else filenames.push_back(arg);
	}
}


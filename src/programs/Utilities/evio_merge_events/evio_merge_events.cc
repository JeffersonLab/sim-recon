// $Id$
//
// Created June 10, 2014  David Lawrence

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

#include <signal.h>
#include <time.h>
#include <stdlib.h>

#include <evioFileChannel.hxx>
#include <evioUtil.hxx>
using namespace evio;


#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);
void Process(unsigned int &NEvents, unsigned int &NEvents_read);

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

	// Process all events
	Process(NEvents, NEvents_read);
	
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
		sprintf(OUTFILENAME,"merged_files.evio");
	}
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     evio_merge_events [-oOutputfile] file1.evio file2.evio ..."<<endl;
	cout<<endl;
	cout<<"options:"<<endl;
	cout<<"    -oOutputfile  Set output filename (def. merged_files.evio)"<<endl;
	cout<<endl;
	cout<<" This will merge events from 1 or more EVIO files into a single EVIO file."<<endl;
	cout<<"This is done at the event level by copying all EVIO banks from the top-level" << endl;
	cout<<"bank into the top-level bank of the output file." << endl;
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

//-----------
// Process
//-----------
void Process(unsigned int &NEvents, unsigned int &NEvents_read)
{
	// Output file
	cout<<" output file: "<<OUTFILENAME<<endl;
	evioFileChannel ochan(OUTFILENAME, "w");
	ochan.open();
	
	// Open all input files
	vector<evioFileChannel*> ichan;
	for(unsigned int i=0; i<INFILENAMES.size(); i++){
		cout << "Opening input file : \"" << INFILENAMES[i] << "\"" << endl;
		evioFileChannel *chan = new evioFileChannel(INFILENAMES[i], "r");
		chan->open();
		ichan.push_back(chan);
	}
	
	// Loop until an input file runs out of events
	time_t last_time = time(NULL);
	while(true){
	
		// Read in event from each input file, creating a DOM tree for each
		vector<evioDOMTree*> doms;
		for(unsigned int i=0; i<ichan.size(); i++){
			try{
				if(! ichan[i]->read() ) {
					cout << endl << "No more events in " << INFILENAMES[i] << endl;
					break;
				}
			}catch(evioException e){
				cerr << e.what() << endl;
				QUIT=true;
				break;
			}
			evioDOMTree *dom = new evioDOMTree(ichan[i]);
			doms.push_back(dom);
		}
		if(QUIT) break;
		if(doms.size() != ichan.size()) break;
		NEvents_read++;
		
		// Merge DOM Trees into single DOM tree
		evioDOMTreeP tree=NULL;
		evioDOMNodeP root=NULL;
		for(unsigned int i=0; i<doms.size(); i++){
			evioDOMNodeListP nodes = doms[i]->getNodeList();
			evioDOMNodeList::const_iterator ulIter;
			for(ulIter=nodes->begin(); ulIter!=nodes->end(); ulIter++) {
				evioDOMNodeP node = *ulIter;
				if(node->getParent() != NULL ) continue;  // only interested in top-level nodes
				if(tree==NULL){
					// Use the first node of the first file as the root node
					root = node;
					root->cut(); // remove from input file's DOM
					tree = new evioDOMTree(root);
				}else{
					// Move all daughter nodes of this input file's top-level node
					// into the root node of the merged event
					evioDOMNodeListP mynodes = node->getChildren();
					if(mynodes.get() == NULL) continue;
					evioDOMNodeList::const_iterator myiter;
					for(myiter=mynodes->begin(); myiter!=mynodes->end(); myiter++) {
						try{
							(*myiter)->move(root);
						}catch(...){
							QUIT = true;
							break;
						}
					}
					if(QUIT) break;
				}
				if(QUIT) break;
			}
			if(QUIT) break;
		}
		
		// Write event to output file
		if(!QUIT){
			ochan.write(tree);
			NEvents++;
		}

		// Update ticker
		time_t now = time(NULL);
		if(now != last_time){
			cout<<"  "<<NEvents_read<<" events read     ("<<NEvents<<" event written) \r";cout.flush();
			last_time = now;
		}
		
		if(QUIT)break;
	}
	
	// Close all input files
	for(unsigned int i=0; i<ichan.size(); i++){
		try{
			ichan[i]->close();
			delete ichan[i];
		}catch(...){}
	}
	
	// Close output file
	ochan.close();

}


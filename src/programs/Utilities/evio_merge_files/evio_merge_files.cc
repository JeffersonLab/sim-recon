// $Id:$
//
// Created Feb. 24, 2015  David Lawrence

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

#include <signal.h>
#include <time.h>
#include <stdlib.h>



#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);
void Process(unsigned int &NEvents, unsigned int &NEvents_read);

static vector<char*> INFILENAMES;
static char *OUTFILENAME = NULL;
static int QUIT = 0;
static int BLOCKSIZE = 10485760;   // 10 M
static bool VERBOSE = false;

//-----------
// GetFilesize
//-----------
ifstream::pos_type GetFilesize(const char* filename)
{
    std::ifstream in(filename, ifstream::ate | ifstream::binary);
    return in.tellg(); 
}

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
	
	//cout<<endl;
	//cout<<" "<<NEvents_read<<" events read, "<<NEvents<<" events written"<<endl;

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
		sprintf(OUTFILENAME,"merged.evio");
	}
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     evio_merge_files [-oOutputfile] file1.evio file2.evio ..."<<endl;
	cout<<endl;
	cout<<"options:"<<endl;
	cout<<"    -oOutputfile  Set output filename (def. merged.evio)"<<endl;
	cout<<endl;
	cout<<" This will merge together multiple EVIO files into a single EVIO file."<<endl;
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
    // An EVIO file is made of a series of independent blocks, so one should be
    // able to concatenate EVIO files together, except that the end-of-file
    // is indicated by an empty EVIO block, i.e., an 8-word EVIO header 
    // with no data payload.  Therefore, we concatenate the files leaving 
    // out the end-of-file headers, except for the final files.
    // This is different than the previous algorithm, which used the standard
    // EVIO library and required parsing each EVIO event.  This new algorithm
    // should be faster and reduce our dependency on the EVIO library.
    // Note that the merged files are currently unreliable.  We leave this
    // code in for future development, since the old code currently crashes.
    // sdobbs, 6/15/2016

    // Set up a buffer
    char *buffer = new char[BLOCKSIZE];

	// Output file
	cout<<" output file: "<<OUTFILENAME<<endl;
    ofstream outfile (OUTFILENAME, ios::out|ios::binary);
    if(!outfile.is_open()) {
        cerr<<"Could not open output file " << OUTFILENAME << endl;
        return;
    }
	
	// Open all input files
	for(unsigned int i=0; i<INFILENAMES.size(); i++){
		cout << "Opening input file : \"" << INFILENAMES[i] << "\"" << endl;
        ifstream infile (INFILENAMES[i] , ios::in|ios::binary);
        if(!infile.is_open()) {
            cerr<<"Could not open input file " << INFILENAMES[i] << endl;
            continue;
        }

        // figure out the size of the file
        ifstream::pos_type length_to_read = GetFilesize(INFILENAMES[i]);

        if(VERBOSE) {
            cout << "first block" << endl;
            infile.read(buffer, 32);
            unsigned int *int_buffer = (unsigned int*)buffer;
            for(int i=0; i<8; i++)
                cout << "0x" << hex << int_buffer[i] << endl;
            infile.seekg(-32, ios_base::cur);
        }
        // copy files in blocks to avoid reading in gigs of memory at once
        while(length_to_read > BLOCKSIZE) {
            infile.read(buffer, BLOCKSIZE);
            outfile.write(buffer, BLOCKSIZE);
            length_to_read -= BLOCKSIZE;
        }
        //strip out EOF blocks for all files but the last one
        if(i+1 != INFILENAMES.size())
            length_to_read -= 8*4L;
        infile.read(buffer, length_to_read);
        outfile.write(buffer, length_to_read);

        if(VERBOSE) {
            // dump end
            cout << "last block" << endl;
            infile.read(buffer, 32);
            //unsigned int *int_buffer = (unsigned int*)buffer;
            unsigned int *int_buffer = (unsigned int*)buffer;
            for(int i=0; i<8; i++)
                cout << "0x" << hex << int_buffer[i] << endl;
        }

        infile.close();
    }

    outfile.close();
}

#if 0

// old code disabled, see description above (sdobbs, 6/15/2016)

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
		evioFileChannel *ichan = new evioFileChannel(INFILENAMES[i], "r");
		ichan->open();
	
		// Loop until input file runs out of events
		time_t last_time = time(NULL);
		while(true){
		
			// Read in event from input file, creating a DOM tree
			try{
				if(! ichan->read() ) {
					cout << endl << "No more events in " << INFILENAMES[i] << endl;
					break;
				}
			}catch(evioException e){
				cerr << e.what() << endl;
				QUIT=true;
				break;
			}
			evioDOMTree *dom = new evioDOMTree(ichan);

			if(QUIT) break;
			NEvents_read++;
			
			// Merge DOM Trees into single DOM tree
			evioDOMTreeP tree=NULL;
			evioDOMNodeP root=NULL;
			evioDOMNodeListP nodes = dom->getNodeList();
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
	
		// Close input file
		try{
			ichan->close();
			delete ichan;
		}catch(...){}
	}
	
	// Close output file
	ochan.close();

}

#endif

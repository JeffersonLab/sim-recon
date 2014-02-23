// $Id$
//
// Created Oct 25, 2013  Kei Moriya

#include "hddm_merge_files.h"

#include <HDDM/hddm_r.hpp>
using namespace hddm_r;


// NOTE: The REST format files claim they are in a compressed format that
// the C API cannot handle so we use the C++ API here


//-----------
// Process_r  --  HDDM REST format
//-----------
void Process_r(unsigned int &NEvents, unsigned int &NEvents_read)
{
	// Output file
	cout<<" output file: "<<OUTFILENAME<<endl;
	ofstream ofs(OUTFILENAME);
	hddm_r::ostream ostr(ofs);

	if(!ofs.is_open()){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}

	// Loop over input files
	time_t last_time = time(NULL);
	for(unsigned int i=0; i<INFILENAMES.size(); i++){
		cout<<" input file: "<<INFILENAMES[i]<<endl;

		// Open hddm file for reading
		ifstream ifs(INFILENAMES[i]);

		// Associate input file stream with HDDM record
		hddm_r::istream istr(ifs);
		
		// Loop over events
		while(!ifs.eof() && ifs.good()){
			try{
				HDDM xrec;
				istr >> xrec;
				NEvents_read++;
				
				ostr << xrec;
				NEvents++;
			
				// Update ticker
				time_t now = time(NULL);
				if(now != last_time){
					cout<<"  "<<NEvents_read<<" events read     ("<<NEvents<<" event written) \r";cout.flush();
					last_time = now;
				}
				
				if(QUIT)break;
			}catch(...){
				break;
			}
		}
	}
}

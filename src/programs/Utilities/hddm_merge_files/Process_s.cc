// $Id$
//
// Created Oct 25, 2013  Kei Moriya

#include "hddm_merge_files.h"

#include <HDDM/hddm_s.h>



//-----------
// Process_s  --  HDDM simulation format
//-----------
void Process_s(unsigned int &NEvents, unsigned int &NEvents_read)
{
	// Output file
	cout<<" output file: "<<OUTFILENAME<<endl;
	s_iostream_t *fout = init_s_HDDM(OUTFILENAME);
	if(!fout){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}

	// Loop over input files
	time_t last_time = time(NULL);
	for(unsigned int i=0; i<INFILENAMES.size(); i++){
		cout<<" input file: "<<INFILENAMES[i]<<endl;
		s_iostream_t *fin = open_s_HDDM(INFILENAMES[i]);
		if(!fin){
			cout<<" Error opening input file \""<<INFILENAMES[i]<<"\"!"<<endl;
			exit(-1);
		}
			
		// Loop over all events in input
		while(true){
			s_HDDM_t *hddm_s = read_s_HDDM(fin);
			if(!hddm_s)break;
			NEvents_read++;
			
			// Write this output event to file and free its memory
			flush_s_HDDM(hddm_s, fout);
			NEvents++;
		
			// Update ticker
			time_t now = time(NULL);
			if(now != last_time){
				cout<<"  "<<NEvents_read<<" events read     ("<<NEvents<<" event written) \r";cout.flush();
				last_time = now;
			}

			if(QUIT)break;
		}

		// Close input file
		close_s_HDDM(fin);
	}
		
	// Close output file
	close_s_HDDM(fout);
}

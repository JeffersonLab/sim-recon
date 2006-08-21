/*----------------------------------------------------------------------------*
	GlueX file to ET  producer

	This program reads events from an EVIO(CODA) file(s) and then inserts
	them into an ET system.
 *----------------------------------------------------------------------------*/


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

#include <et.h>
#include <evio.h>

#define max_event_size 100000

// Globals
et_sys_id id;
et_att_id  attach;
string et_fname("");
useconds_t DELAY;
vector<string> source_files;
bool LOOP=false;

// Routines
int InsertEventIntoET(const char *buff, int nwords);
void Usage(void);

//------------------------
// main
//------------------------
int main(int narg,char **argv)
{
	// Loop over command line arguments.
	for(int i=1; i<narg; i++){
		string s(argv[i]);
		if(s=="-h" || s=="--help")Usage();
		if(s=="-loop"){
			LOOP = true;
			continue;
		}
		if(s=="-d"){
			DELAY = atoi(argv[++i]);
			continue;
		}
		if(s=="-f"){
			et_fname = argv[++i];
			continue;
		}
		source_files.push_back(s);
	}
	if(source_files.size()==0 || et_fname.size()==0)Usage();
  
	// Open the ET system for inserting events into
	et_openconfig openconfig;
	et_open_config_init(&openconfig);
	if(et_open(&id, et_fname.c_str(), openconfig) != ET_OK){
		cerr<<__FILE__<<":"<<__LINE__<<" Error opening ET system "<<et_fname<<endl;
		return -1;
	}
	et_open_config_destroy(openconfig);
  
	// set the debug level (everything)
	et_system_setdebug(id, ET_DEBUG_INFO);
   
	// attach to GRANDCENTRAL station since we are producing events
	if (et_station_attach(id, ET_GRANDCENTRAL, &attach) < 0) {
		cerr<<__FILE__<<":"<<__LINE__<<" Error attaching to GRAND CENTRAL on "<<et_fname<<endl;
		return -1;
	}

	// Outer loop for when LOOP is true
	if(LOOP){
		cout<<"Entering infinite loop over events. Hit ctl-C to quit ..."<<endl;
	}
	do{

		// Loop over source files
		for(unsigned int i=0; i<source_files.size(); i++){
			// Open EVIO file
			int handle;
			if(evOpen((char*)source_files[i].c_str(),"r", &handle)!=S_SUCCESS){
				cerr<<"Error opening \""<<source_files[i]<<"\""<<endl;
				continue;
			}
		
			// Loop over events in file
			char buff[max_event_size];
			while(evRead(handle, (unsigned long*)buff, max_event_size)==S_SUCCESS){
				// Event size in words should be first word in buffer
				int event_size = *(int*)buff + 1;
				if(event_size>max_event_size)event_size = max_event_size;
				InsertEventIntoET(buff, event_size);
			}

			// Close EVIO file
			evClose(handle);
		}
	}while(LOOP);

	return 0;
}

//------------------------
// InsertEventIntoET
//------------------------
int InsertEventIntoET(const char *buff, int nwords)
{
	et_event   *pe;
	int status;
	char *pdata;

	// Make sure ET system is still alive. Wait for it if not.
	if (!et_alive(id)) {
		et_wait_for_alive(id);
	}

	/* get new/unused event */  
	status = et_event_new(id, attach, &pe, ET_SLEEP, NULL, nwords*sizeof(int));
	if (status != ET_OK) {
		printf("et_producer: error in et_event_new\n");
		exit(0);
	}

	// Get pointer to data of new event
    et_event_getdata(pe,(void**)&pdata);

    // Copy contents into new event
    memcpy( (char*)pdata, (char*)buff, nwords*sizeof(int));
    
	// put event back into the ET system
	status = et_event_put(id, attach, pe);
	if (status != ET_OK) {
		printf("et_producer: put error\n");
			exit(0);
	}
	
	// Delay a little (if specified) to lower the event rate
	if(DELAY != 0)usleep(DELAY);
	
	return 0;
}

//------------------------
// Usage
//------------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"      file2et [-d delay] [-loop] -f et_system_file file [file] ..."<<endl;
	cout<<endl;
	cout<<"Read events from an EVIO (CODA formatted file) and write "<<endl;
	cout<<"them into the specified ET system."<<endl;
	cout<<endl;
	cout<<"The value of et_system_file is usually something like:"<<endl;
	cout<<endl;
	cout<<"/tmp/et_sys_mysession"<<endl;
	cout<<endl;
	cout<<"If the -d delay option is given, then a delay of \"delay\""<<endl;
	cout<<"microseconds will be inserted between events so that the"<<endl;
	cout<<"rate can be controlled. Note that most operating systems"<<endl;
	cout<<"have a minimum delay so if you specify a value below it,"<<endl;
	cout<<"the minimum delay will still occur. No additional delay is"<<endl;
	cout<<"incurred if the value of delay is 0 or the -d argument is"<<endl;
	cout<<"ommitted."<<endl;
	cout<<endl;
	cout<<"If the -loop option is specified, then the source files are"<<endl;
	cout<<"repeatedly read over and over until the program is killed."<<endl;
	cout<<endl;
	
	exit(0);
}




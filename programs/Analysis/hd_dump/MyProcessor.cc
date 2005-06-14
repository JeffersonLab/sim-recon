// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <stdio.h>
#include <unistd.h>

#include "MyProcessor.h"
#include "hddm_s.h"


int PAUSE_BETWEEN_EVENTS = 1;
int SKIP_BORING_EVENTS = 0;
int PRINT_ALL=0;

vector<string> toprint;

#define ansi_escape		((char)0x1b)
#define ansi_bold 		ansi_escape<<"[1m"
#define ansi_black		ansi_escape<<"[30m"
#define ansi_red			ansi_escape<<"[31m"
#define ansi_green		ansi_escape<<"[32m"
#define ansi_blue			ansi_escape<<"[34m"
#define ansi_normal		ansi_escape<<"[0m"
#define ansi_up(A)		ansi_escape<<"["<<(A)<<"A"
#define ansi_down(A)		ansi_escape<<"["<<(A)<<"B"
#define ansi_forward(A)	ansi_escape<<"["<<(A)<<"C"
#define ansi_back(A)		ansi_escape<<"["<<(A)<<"D"

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
derror_t MyProcessor::brun(DEventLoop *eventLoop, int runnumber)
{
	vector<string> factory_names = eventLoop->GetFactoryNames();

	usleep(100000); //this just gives the Main thread a chance to finish printing the "Launching threads" message
	cout<<endl;

	// If int PRINT_ALL is set then add EVERYTHING.
	if(PRINT_ALL){
		toprint = factory_names;
		SKIP_BORING_EVENTS = 0; // with PRINT_ALL, nothing is boring!
	}else{
		// make sure factories exist for all requested data types
		// If a factory isn't found, but one with a "D" prefixed
		// is, go ahead and correct the name.
		vector<string> really_toprint;
		for(unsigned int i=0; i<toprint.size();i++){
			int found = 0;
			int dfound = 0;
			for(unsigned int j=0;j<factory_names.size();j++){
				if(factory_names[j] == toprint[i])found = 1;
				if(factory_names[j] == "D" + toprint[i])dfound = 1;
			}
			if(found)
				really_toprint.push_back(toprint[i]);
			else if(dfound)
				really_toprint.push_back("D" + toprint[i]);
			else
				cout<<ansi_red<<"WARNING:"<<ansi_normal
					<<" Couldn't find factory for \""
					<<ansi_bold<<toprint[i]<<ansi_normal
					<<"\"!"<<endl;
		}
		
		toprint = really_toprint;
	}
	
	cout<<endl;

	return NOERROR;
}

//------------------------------------------------------------------
// evnt
//------------------------------------------------------------------
derror_t MyProcessor::evnt(DEventLoop *eventLoop, int eventnumber)
{

	// If we're skipping boring events (events with no rows for any of
	// the types we're printing) we must find out first if the event is
	// "boring".
	int event_is_boring = 1;
	if(SKIP_BORING_EVENTS){
		for(unsigned int i=0;i<toprint.size();i++){
			DFactory_base *factory = eventLoop->GetFactory(toprint[i]);
			if(!factory)factory = eventLoop->GetFactory("D" + toprint[i]);
			if(factory){
				if(factory->GetNrows()>0){
					event_is_boring=0;
					break;
				}
			}
		}
		if(event_is_boring)return NOERROR;
	}else{
		event_is_boring = 0;
	}
	
	// Print event separator
	cout<<"================================================================"<<endl;
	cout<<"Event: "<<eventnumber<<endl;

	// We want to print info about all factories results, even if we aren't
	// printing the actual data.
	eventLoop->PrintFactories(1);
	
	// Print data for all specified factories
	for(unsigned int i=0;i<toprint.size();i++){
		eventLoop->Print(toprint[i]);
	}
	
	// Wait for user input if pausing
	if(PAUSE_BETWEEN_EVENTS && !event_is_boring){
		cerr.flush();
		cout<<endl<<"< Hit return for the next event (P=prev. Q=quit) >";
		cout.flush();
		char c = getchar(); // see note in hd_dump.cc:main()
		if(c=='\n')cout<<ansi_up(1);
		cout<<endl;
		switch(toupper(c)){
			case 'Q':
				eventLoop->Quit();
				break;
			case 'P':
				//eventLoop->GotoEvent(eventnumber-1);
				break;
			case 'N':
				break;
		}

		cout<<ansi_up(1)<<"\r                                                     \r";
		cout.flush();
	}
	
	return NOERROR;
}


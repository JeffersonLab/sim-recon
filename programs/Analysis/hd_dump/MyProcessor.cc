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
// evnt   -Fill histograms here
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	// If int PRINT_ALL is set then add EVERYTHING.
	if(PRINT_ALL){
		toprint = eventLoop->GetFactoryNames();
		for(int i=0; i<toprint.size(); i++){
			string name = toprint[i];
			cout<<"Adding to print list : "<<name<<endl;
		}

		PRINT_ALL = 0; // clear PRINT_ALL flag so we don't re-add all factories
		SKIP_BORING_EVENTS = 0; // with PRINT_ALL, nothing is boring!
	}

	// If we're skipping boring events (events with no rows for any of
	// the types we're printing) we must find out first if the event is
	// "boring".
	int event_is_boring = 1;
	if(SKIP_BORING_EVENTS){
		for(int i=0;i<toprint.size();i++){
			DFactory_base *factory = eventLoop->GetFactory(toprint[i]);
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
	for(int i=0;i<toprint.size();i++){
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
				eventLoop->GotoEvent(eventnumber-1);
				break;
			case 'N':
				break;
		}

		cout<<ansi_up(1)<<"\r                                                     \r";
		cout.flush();
	}
	
	return NOERROR;
}


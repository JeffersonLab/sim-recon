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

#include "DFactory_CDCHits.h"
#include "DFactory_CDCClusters.h"
#include "DFactory_FCALHits.h"

int PAUSE_BETWEEN_EVENTS = 1;
int SKIP_BORING_EVENTS = 0;
int PRINT_ALL=0;

int Ntoprint = 0;
char *toprint[1024];

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
		DContainer *factoryNames = event_loop->GetFactoryNames();
		char **name = (char**)factoryNames->first();
		Ntoprint = 0; // Clear any data specified by -D so it doesn't print twice
		for(int i=0; i<factoryNames->nrows; i++, name++){
			cout<<"Adding to print list : "<<*name<<endl;
			toprint[Ntoprint++] = *name;
		}

		PRINT_ALL = 0; // clear PRINT_ALL flag so we don't re-add all factories
		SKIP_BORING_EVENTS = 0; // with PRINT_ALL, nothing is boring!
		delete factoryNames;
	}

	// If we're skipping boring events (events with no rows for any of
	// the types we're printing) we must find out first if the event is
	// "boring".
	int event_is_boring = 1;
	if(SKIP_BORING_EVENTS){
		for(int i=0;i<Ntoprint;i++){
			DContainer *c = event_loop->Get(toprint[i]);
			if(c){
				if(c->nrows>0){
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

	// We want to print info about all factories results, even if we aren't
	// printing the actual data. Hence we must call every factory's Get() method.
	DContainer *factoryNames = event_loop->GetFactoryNames();
	char **name = (char**)factoryNames->first();
	for(int i=0; i<factoryNames->nrows; i++, name++)event_loop->Get(*name);
	delete factoryNames;	
	event_loop->PrintFactories(1); // print sparsified factory info
	
	// Print data for all specified factories
	for(int i=0;i<Ntoprint;i++){
		event_loop->Print(toprint[i]);
	}
	
	// Wait for user input if pausing
	if(PAUSE_BETWEEN_EVENTS && !event_is_boring){
		char c;
		cerr.flush();
		cout<<endl<<"< Hit return for the next event (type Q to quit) >";
		cout.flush();
		while(!scanf("%c",&c)){cout<<"c="<<(int)c<<endl;usleep(100000);}
		if(c=='q' || c=='Q')event_loop->Quit();
		cout<<ansi_up(1)<<"\r                                                     \r";
		cout.flush();
	}
	
	return NOERROR;
}


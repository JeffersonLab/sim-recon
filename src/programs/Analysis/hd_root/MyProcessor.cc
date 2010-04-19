// Author: Edward Brash February 15, 2005
// revised severely 2006-2007  David Lawrence
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <TROOT.h>
#include <TTree.h>

#include "MyProcessor.h"

vector<string> toprint;
bool ACTIVATE_ALL=false;

extern string OUTPUT_FILENAME;

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
// MyProcessor
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
	ROOTfile = NULL;
}

//------------------------------------------------------------------
// ~MyProcessor
//------------------------------------------------------------------
MyProcessor::~MyProcessor()
{

}

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// open ROOT file
	ROOTfile = new TFile(OUTPUT_FILENAME.c_str(),"RECREATE","Produced by hd_root");
	cout<<"Opened ROOT file \""<<OUTPUT_FILENAME<<"\" ..."<<endl;

	return NOERROR;
}

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventLoop, int runnumber)
{
	vector<string> factory_names;
	eventLoop->GetFactoryNames(factory_names);

	usleep(100000); //this just gives the Main thread a chance to finish printing the "Launching threads" message
	cout<<endl;

	// If ACTIVATE_ALL is set then add EVERYTHING.
	if(ACTIVATE_ALL){
		toprint = factory_names;
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
	
	// At this point, toprint should contain a list of all factories
	// in dataClassName:tag format, that both exist and were requested.
	// Seperate the tag from the name and fill the fac_info vector.
	fac_info.clear();
	for(unsigned int i=0;i<toprint.size();i++){
		string name = toprint[i];
		string tag = "";
		unsigned int pos = name.rfind(":",name.size()-1);
		if(pos != (unsigned int)string::npos){
			tag = name.substr(pos+1,name.size());
			name.erase(pos);
		}
		factory_info_t f;
		f.dataClassName = name;
		f.tag = tag;
		fac_info.push_back(f);
	}
	
	cout<<endl;

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill tree here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{
	// Loop over factories explicitly mentioned on command line
	for(unsigned int i=0;i<toprint.size();i++){
		string name =fac_info[i].dataClassName;
		string tag = fac_info[i].tag;
		JFactory_base *factory = eventLoop->GetFactory(name,tag.c_str());
		if(!factory)factory = eventLoop->GetFactory("D" + name,tag.c_str());
		if(factory){
			try{
				factory->GetNrows();
			}catch(...){
				// someone threw an exception
			}
		}
	}
	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
	if(ROOTfile!=NULL){
		ROOTfile->Write();
		//ROOTfile->Close();
		//delete ROOTfile;
		ROOTfile=NULL;
		cout<<endl<<"Closed ROOT file"<<endl;
	}
	
	return NOERROR;
}


// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <TFile.h>
#include <TTree.h>

extern vector<string> toprint;
extern bool ACTIVATE_ALL;

class MyProcessor:public JEventProcessor
{
	public:
		MyProcessor();
		~MyProcessor();
	
		jerror_t init(void);				///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);						///< Called every event.
		jerror_t erun(void){return NOERROR;};				///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);				///< Called after last event of last event source has been processed.

		typedef struct{
			string dataClassName;
			string tag;
		}factory_info_t;
		vector<factory_info_t> fac_info;

		TFile *ROOTfile;
};

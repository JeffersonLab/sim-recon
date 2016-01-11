#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JApplication.h>
#include "DL1Trig_factory.h"
#include "DAQ/DCODAROCInfo.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DL1Trig_factory::init(void)
{

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DL1Trig_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DL1Trig_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{


	DL1Trig *l1trig = new DL1Trig;


	l1trig->timestamp     =  0;
	l1trig->event_type    =  0;	
	l1trig->trig_mask     =  0;
	l1trig->fp_trig_mask  =  0;


	vector<const DCODAROCInfo*> coda_info;
	loop->Get(coda_info);
	
	
	if(coda_info.size() == 0){
	  //	  cout << " CODA_INFO does not exist" << endl;
	  return NOERROR;
	}

	l1trig->event_type  =  1;

	//	cout << " Event = " << eventnumber << endl;


	for(unsigned int ii = 0; ii < coda_info.size(); ii++){
	  if(coda_info[ii]->rocid == 1){

	    l1trig->timestamp = coda_info[ii]->timestamp;
	    const vector<uint32_t> &info = coda_info[ii]->misc;
	    
	    for(unsigned int jj = 0; jj < info.size(); jj++){
	      if(jj == 0) l1trig->trig_mask     =  info[jj];
	      if(jj == 1) l1trig->fp_trig_mask  =  info[jj];
	      //	      cout << " jj =   "  << jj <<  "   val = " <<  info[jj] << endl;
	    }
	    
	    
	  }
	}

	_data.push_back(l1trig);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DL1Trig_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DL1Trig_factory::fini(void)
{
	return NOERROR;
}


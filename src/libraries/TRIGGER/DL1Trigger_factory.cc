#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JApplication.h>
#include <DAQ/DCODAROCInfo.h>
#include <DAQ/DL1Info.h>

using namespace jana;

#include "DL1Trigger_factory.h"

//------------------
// init
//------------------
jerror_t DL1Trigger_factory::init(void)
{

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DL1Trigger_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DL1Trigger_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// The L1 trigger latch words are added as "misc" words to the 
	// Physic's Event's Built Trigger Bank. These show up in the
	// "misc" member of the DCODAROCInfo object corresponding to the
	// Trigger Supervisor. The TS has rocid of "1" (which will hopefully
	// never change!)


        int l1_found = 0;  

	vector<const DCODAROCInfo*> codarocinfos;
	loop->Get(codarocinfos);
	
	vector<const DL1Info*> l1_info;
	loop->Get(l1_info);

	DL1Trigger *l1trigger = new DL1Trigger;


	for(uint32_t i = 0; i < codarocinfos.size(); i++){
		const DCODAROCInfo *cri = codarocinfos[i];
		
		// Only interested in Trigger Supervisor 
		if(cri->rocid != 1) continue;
		
		if(cri->misc.size()==2){
			// Potentially, the TS can be configured to only 
			// add the GTP latch word or only the FP latch
			// word (or neither or both). We have no way of 
			// knowing how it was configured, but if both are
			// read out, then the GTP word comes first. (See
			// e-mail sent by Bryan Moffit Dec. 15, 2015).
			l1trigger->trig_mask     = cri->misc[0]; // Global Trigger Processor latch word
			l1trigger->fp_trig_mask  = cri->misc[1]; // Front Panel latch word
			l1_found++;
		}

		l1trigger->timestamp =  cri->timestamp;

		l1trigger->AddAssociatedObject(cri);
			
		break; // Should only be 1 with rocid=1
	}

	if(l1_info.size() == 1){
	  const DL1Info *l1_sc   = l1_info[0];
	  l1trigger->nsync       = l1_sc->nsync;
	  l1trigger->trig_number = l1_sc->trig_number;
	  l1trigger->live        = l1_sc->live_time;
	  l1trigger->busy        = l1_sc->busy_time;
	  l1trigger->live_inst   = l1_sc->live_inst;
	  l1trigger->unix_time   = l1_sc->unix_time;
	  
	  l1trigger->gtp_sc    =  l1_sc->gtp_sc;
	  l1trigger->fp_sc     =  l1_sc->fp_sc;
	  l1trigger->gtp_rate  =  l1_sc->gtp_rate;
	  l1trigger->fp_rate   =  l1_sc->fp_rate;

	  l1_found++;	  
	}


	if(l1_found){
	  _data.push_back(l1trigger);	 
	}



	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DL1Trigger_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DL1Trigger_factory::fini(void)
{
	return NOERROR;
}


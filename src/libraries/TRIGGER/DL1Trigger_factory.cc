#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JApplication.h>
#include <DAQ/DCODAROCInfo.h>
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

	vector<const DCODAROCInfo*> codarocinfos;
	loop->Get(codarocinfos);
	
	for(uint32_t i=0; i<codarocinfos.size(); i++){
		const DCODAROCInfo *cri = codarocinfos[i];
		
		// Only interested in Trigger Supervisor 
		if(cri->rocid != 1) continue;
		
		DL1Trigger *l1trigger = new DL1Trigger;
		if(cri->misc.size()==2){
			// Potentially, the TS can be configured to only 
			// add the GTP latch word or only the FP latch
			// word (or neither or both). We have no way of 
			// knowing how it was configured, but if both are
			// read out, then the GTP word comes first. (See
			// e-mail sent by Bryan Moffit Dec. 15, 2015).
			l1trigger->trig_mask     = cri->misc[0]; // Global Trigger Processor latch word
			l1trigger->fp_trig_mask  = cri->misc[1]; // Front Panel latch word
		}
		l1trigger->AddAssociatedObject(cri);
		
		_data.push_back(l1trigger);
		break; // Should only be 1 with rocid=1
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


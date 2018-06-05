#include "DTrigger_factory.h"
#include <bitset>

//------------------
// evnt
//------------------
jerror_t DTrigger_factory::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	vector<const DL1Trigger*> locL1Triggers;
	locEventLoop->Get(locL1Triggers);
	const DL1Trigger* locL1Trigger = locL1Triggers.empty() ? NULL : locL1Triggers[0];

	DTrigger *locTrigger = new DTrigger;
	
	//SET LEVEL-1 TRIGGER INFO
	if(locL1Trigger != NULL)
	{
		locTrigger->Set_L1TriggerBits(locL1Trigger->trig_mask);
		locTrigger->Set_L1FrontPanelTriggerBits(locL1Trigger->fp_trig_mask);
	}
	else 
    {
        // IF TRIGGER INFO IS NOT IN THE DATA STREAM, LOAD TRIGGER SIMULATIONS LAZILY
        
      // cerr << " event status = " << std::bitset<32>(locEventLoop->GetJEvent().GetStatus()) << endl;

        // don't bother simulating the trigger for non-physics events
        // for now, just don't run this for EVIO (raw data) events
        if( locEventLoop->GetJEvent().GetStatusBit(kSTATUS_EVIO) ){
            _data.push_back(locTrigger);
            return NOERROR;
        }

        // realistic trigger simulation
        vector<const DL1MCTrigger*> locMCTriggers;
        locEventLoop->Get(locMCTriggers);
        const DL1MCTrigger* locMCTrigger = locMCTriggers.empty() ? NULL : locMCTriggers[0];

        if(locMCTrigger != NULL)
        {
            //IS MC DATA: USE SIMULATED TRIGGER INFORMATION IF AVAILABLE
            locTrigger->Set_L1TriggerBits(locMCTrigger->trig_mask);
            locTrigger->Set_L1FrontPanelTriggerBits(0);
        }
        else 
        {
	  locTrigger->Set_L1TriggerBits(0);
	  locTrigger->Set_L1FrontPanelTriggerBits(0); 
	}
    }

	//SET LEVEL-3 TRIGGER INFO HERE

	_data.push_back(locTrigger);

	return NOERROR;
}

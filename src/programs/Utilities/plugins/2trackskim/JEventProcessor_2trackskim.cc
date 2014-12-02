//
// JEventProcessor_2trackskim.cc
//
// JANA event processor plugin to skim 2-track events to an EVIO file
//
// Paul Mattione, 19-November-2014

#include "JEventProcessor_2trackskim.h"

// for initializing plugins
extern "C" {
   void InitPlugin(JApplication *app)
	{
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_2trackskim(), true);
   }
} // "extern C"

//-------------------------------
// init
//-------------------------------
jerror_t JEventProcessor_2trackskim::init(void)
{
	dEventWriterEVIO = NULL;
   return NOERROR;
}

//-------------------------------
// brun
//-------------------------------
jerror_t JEventProcessor_2trackskim::brun(JEventLoop *locEventLoop, int runnumber)
{
	locEventLoop->GetSingle(dEventWriterEVIO);
   return NOERROR;
}

//-------------------------------
// evnt
//-------------------------------
jerror_t JEventProcessor_2trackskim::evnt(JEventLoop *locEventLoop, int eventnumber)
{
	//Save EPICS events
	vector<const DEPICSvalue*> locEPICSValues;
	locEventLoop->Get(locEPICSValues);
	if(!locEPICSValues.empty())
	{
		dEventWriterEVIO->Write_EVIOEvent(locEventLoop, "2tracks");
		return NOERROR;
	}

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "PreSelect");
	if(locChargedTracks.size() >= 2)
	{
		dEventWriterEVIO->Write_EVIOEvent(locEventLoop, "2tracks");
		return NOERROR;
	}

   return NOERROR;
}

//-------------------------------
// erun
//-------------------------------
jerror_t JEventProcessor_2trackskim::erun(void)
{
   return NOERROR;
}

//-------------------------------
// fini
//-------------------------------
jerror_t JEventProcessor_2trackskim::fini(void)
{
   return NOERROR;
}


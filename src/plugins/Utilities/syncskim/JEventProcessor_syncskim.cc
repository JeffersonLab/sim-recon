// $Id$
//
//    File: JEventProcessor_syncskim.cc
// Created: Wed Feb 22 20:04:25 EST 2017
// Creator: davidl (on Linux gluon48.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//

#include "JEventProcessor_syncskim.h"
using namespace jana;

#include <DAQ/DL1Info.h>
#include <DAQ/DCODAEventInfo.h>


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_syncskim());
}
} // "C"


//------------------
// JEventProcessor_syncskim (Constructor)
//------------------
JEventProcessor_syncskim::JEventProcessor_syncskim()
{
	sum_n = 0.0;
	sum_x = 0.0;
	sum_y = 0.0;
	sum_xy = 0.0;
	sum_x2 = 0.0;
}

//------------------
// ~JEventProcessor_syncskim (Destructor)
//------------------
JEventProcessor_syncskim::~JEventProcessor_syncskim()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_syncskim::init(void)
{
	
	gPARMS->SetParameter("EVIO:LINK",              false); 
	gPARMS->SetParameter("EVIO:LINK_BORCONFIG",    false); 
	gPARMS->SetParameter("EVIO:PARSE_F250",        false); 
	gPARMS->SetParameter("EVIO:PARSE_F125",        false); 
	gPARMS->SetParameter("EVIO:PARSE_F1TDC",       false); 
	gPARMS->SetParameter("EVIO:PARSE_CAEN1290TDC", false); 
	gPARMS->SetParameter("EVIO:PARSE_CONFIG",      false); 
	gPARMS->SetParameter("EVIO:PARSE_BOR",         false); 
	gPARMS->SetParameter("EVIO:PARSE_EPICS",       false); 
	gPARMS->SetParameter("EVIO:PARSE_EVENTTAG",    false); 
	//gPARMS->SetParameter("EVIO:PARSE_TRIGGER",     false);
	gPARMS->SetParameter("EVIO:APPLY_TRANSLATION_TABLE", false);
	gPARMS->SetParameter("EVIO:F250_EMULATION_MODE", 0);
	gPARMS->SetParameter("EVIO:F125_EMULATION_MODE", 0);
	
	tree = new TTree("T", "Sync Events Tree");
	tree->Branch("run_number",    &synevt.run_number,    "run_number/i"    );
	tree->Branch("run_type",      &synevt.run_type,      "run_type/i"      );
	tree->Branch("event_number",  &synevt.event_number,  "event_number/l"  );
	tree->Branch("event_type",    &synevt.event_type,    "event_type/s"    );
	tree->Branch("avg_timestamp", &synevt.avg_timestamp, "avg_timestamp/l" );

	tree->Branch("nsync",         &synevt.nsync,         "nsync/i"         );
	tree->Branch("trig_number",   &synevt.trig_number,   "trig_number/i"   );
	tree->Branch("live_time",     &synevt.live_time,     "live_time/i"     );
	tree->Branch("busy_time",     &synevt.nsync,         "busy_time/i"     );
	tree->Branch("live_inst",     &synevt.live_inst,     "live_inst/i"     );
	tree->Branch("unix_time",     &synevt.unix_time,     "unix_time/i"     );
 
	tree->Branch("gtp_sc",        &synevt.gtp_sc,        "gtp_sc[32]/i"    );
	tree->Branch("fp_sc",         &synevt.fp_sc,         "fp_sc[16]/i"     );
	tree->Branch("gtp_rate",      &synevt.gtp_rate,      "gtp_rate[32]/i"  );
	tree->Branch("fp_rate",       &synevt.fp_rate,       "fp_rate[16]/i"   );
	tree->Print();

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_syncskim::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_syncskim::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	vector<const DL1Info*> l1infos;
	loop->Get(l1infos);
	if(l1infos.empty()) return NOERROR;
	
	vector<const DCODAEventInfo*> codainfos;
	loop->Get(codainfos);
	if(codainfos.empty()) return NOERROR;
	const DCODAEventInfo *codainfo = codainfos[0];

	japp->RootFillLock(this);
	
	for(auto l1info : l1infos){
		
		time_t t = l1info->unix_time;
		cout << "sync event at " << ctime(&t);
		
		synevt.run_number    = codainfo->run_number;
		synevt.run_type      = codainfo->run_type;
		synevt.event_number  = codainfo->event_number;
		synevt.event_type    = codainfo->event_type;
		synevt.avg_timestamp = codainfo->avg_timestamp;
	
		synevt.nsync       = l1info->nsync;
		synevt.trig_number = l1info->trig_number;
		synevt.live_time   = l1info->live_time;
		synevt.busy_time   = l1info->busy_time;
		synevt.live_inst   = l1info->live_inst;
		synevt.unix_time   = l1info->unix_time;
		for(int i=0; i<32; i++){
			synevt.gtp_sc[i]   = l1info->gtp_sc[i];
			synevt.gtp_rate[i] = l1info->gtp_rate[i];
			if(i<16){
				synevt.fp_sc[i]   = l1info->fp_sc[i];
				synevt.fp_rate[i] = l1info->fp_rate[i];
			}
		}		
		
		tree->Fill();
		
		// scale and shift x/y values to make sure range of sum doesn't cut them off
		double x = (double)synevt.avg_timestamp / 250.0E6;
		double y = (double)synevt.unix_time - 13.0E8;
		sum_n  += 1.0;
		sum_x  += x;
		sum_y  += y;
		sum_xy += x * y;
		sum_x2 += x * x;
	}
	
	japp->RootFillUnLock(this);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_syncskim::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_syncskim::fini(void)
{

	double m = (sum_n*sum_xy - sum_x*sum_y)/(sum_n*sum_x2 - sum_x*sum_x);
	double b = (sum_y*sum_x2 - sum_x*sum_xy)/(sum_n*sum_x2 - sum_x*sum_x);
	
	// scale shift back
	b += 13.0E8;
	double one_over_m = 250.0E6/m;
	
	cout << endl << "timestamp to unix time conversion: 1/m=" << one_over_m << " b=" << (uint64_t)b << endl << endl;

	return NOERROR;
}


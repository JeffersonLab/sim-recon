// $Id$
//
//    File: JEventSource_ETEVIO.h
// Created: Mon Nov 26 10:48:42 EST 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2  x86_64)
//

#include "JEventSource_ETEVIO.h"
using namespace jana;

#include <evioETChannel.hxx>
using namespace evio;


//----------------
// Constructor
//----------------
JEventSource_ETEVIO::JEventSource_ETEVIO(const char* source_name):JEventSource_EVIO(source_name)
{
	// open event source here
	et_sys_id sys_id;
	et_att_id att_id;
	et_stat_id sta_id;
	
	// Split source name into session, station, etc...
	vector<string> fields;
	size_t cutAt;
	string str = source_name;
	while( (cutAt = str.find(":")) != str.npos ){
		if(cutAt > 0)fields.push_back(str.substr(0,cutAt));
		str = str.substr(cutAt+1);
	}
	if(str.length() > 0)fields.push_back(str);
	string session = fields.size()>1 ? fields[1]:"none";
	string station = fields.size()>2 ? fields[2]:"DANA";
	int Nevents    = fields.size()>3 ? atoi(fields[3].c_str()):1;
	
	cout<<"Opening ET session:"<<session<<"  station:"<<station<<endl;

	et_openconfig openconfig;
	et_open_config_init(&openconfig);

	// create station config in case no station exists
	et_statconfig et_station_config;
	et_station_config_init(&et_station_config);
	et_station_config_setblock(et_station_config,ET_STATION_BLOCKING);
	et_station_config_setselect(et_station_config,ET_STATION_SELECT_ALL);
	et_station_config_setuser(et_station_config,ET_STATION_USER_MULTI);
	et_station_config_setrestore(et_station_config,ET_STATION_RESTORE_OUT);
	et_station_config_setcue(et_station_config,Nevents);
	et_station_config_setprescale(et_station_config,0);
	cout<<"ET station configured\n";

	// connect to the ET system
	string fname = string("/tmp/et_sys_") + session;
	et_open_config_init(&openconfig);
	if(et_open(&sys_id,fname.c_str(),openconfig)!=ET_OK){
		cout<<__FILE__<<":"<<__LINE__<<" Problem opening ET system"<<endl;
		return;
	}
	
	// create station if not already created
	int status=et_station_create(sys_id,&sta_id,station.c_str(),et_station_config);
	if((status!=ET_OK)&&(status!=ET_ERROR_EXISTS)) { 
		et_close(sys_id);
		cerr << "Unable to create station " << station << endl;
		return;
	}
	cout<<"ET station created\n";

	status=et_station_attach(sys_id,sta_id,&att_id);
	if(status!=ET_OK) {
		et_close(sys_id);
		cerr << "Unable to attach to station " << station << endl;
		return;
	}

	cout << "...now connected to ET system: " << fname 
		<< ",   station: " << station << " (station id=" << sta_id << ", attach id=" << att_id <<")" << endl;
		
	chan = new evioETChannel(sys_id, att_id);

	if(chan)chan->open();
}

//----------------
// ReadEVIOEvent
//----------------
jerror_t JEventSource_ETEVIO::ReadEVIOEvent(void)
{
	if(!chan->read())return NO_MORE_EVENTS_IN_SOURCE;

	return NOERROR;
}

//----------------
// Destructor
//----------------
JEventSource_ETEVIO::~JEventSource_ETEVIO()
{
	// chan is closed and deleted by base class destructor
}

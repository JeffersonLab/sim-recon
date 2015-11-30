// $Id$
//
//    File: JEventProcessor_janaded.cc
// Created: Fri 20 Jul 2012 10:03:57 AM EDT 
// Creator: garmon
//

#include <iostream>
#include <fstream>
#include <queue>
using namespace std;

#include <JANA/JApplication.h>
#include "JEventProcessor_janaded.h"

#include <pthread.h>

using namespace jana;

static pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t cond;
static queue<string> commands;



// Routine used to allow us to register our JEventSourceGenerator
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_janaded());
}
} // "C"

// Tokenize a string
static inline void Tokenize(string str, vector<string> &tokens, const char delim=' ')
{
	tokens.clear();
	unsigned int cutAt;
	while( (cutAt = str.find(delim)) != (unsigned int)str.npos ){
		if(cutAt > 0)tokens.push_back(str.substr(0,cutAt));
		str = str.substr(cutAt+1);
	}
	if(str.length() > 0)tokens.push_back(str);
}

//-----------------------------------------
// JEventProcessor_janaded (constructor)
//-----------------------------------------
JEventProcessor_janaded::JEventProcessor_janaded()
{
	string myName = "JEventProcessor_janaded";
  	string myDescription = "janaded";
  	string UDL = "cMsg:cMsg://localhost/cMsg/test";

	cMsg* conn = new cMsg(UDL, myName, myDescription);
	conn->connect();
	conn->subscribe("janaded", "test", this, NULL);
	conn->start();

	japp->monitor_heartbeat=false;

}

//------------------
// init
//------------------
jerror_t JEventProcessor_janaded::init(void)
{
	JANADED_VERBOSE=0;
	app->GetJParameterManager()->SetDefaultParameter("JANADED_VERBOSE", JANADED_VERBOSE);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_janaded::brun(JEventLoop *eventLoop, int32_t runnumber)
{

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_janaded::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	while(true) {
		pthread_mutex_lock(&mutex1);
		pthread_cond_wait( &cond, &mutex1 );
		string command = commands.front();
		commands.pop();
		if(command=="NEXT_EVENT") {
			cout << "Next Event Command" <<endl;
			return NOERROR;
		}
		else if (command=="Whatever") {
			//do whatever
		}
		else {
			cout << "Command was" << command <<endl;
		}

		pthread_mutex_unlock(&mutex1);
	}
	
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_janaded::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_janaded::fini(void)
{

	return NOERROR;
}

void JEventProcessor_janaded::callback(cMsgMessage *msg, void *arg) {  
	pthread_mutex_lock(&mutex1);
	commands.push(msg->getText());
	pthread_mutex_unlock(&mutex1);	
	//pthread_cond_broadcast( &cond );
	pthread_cond_signal( &cond ); 
	
	delete msg;
}




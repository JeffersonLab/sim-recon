// $Id$
//
//    File: JEventProcessor_CODA_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)

// Note...user must set eventCount and dataCount in statistics_thread for monitoring purposes

// ejw, 8-Nov-2013


// still to do:

//   need unique name
//   set event count and data count
//   implement Vardan's new status reporting scheme
//   check for missing thread.join(), likely in codaObject library
//   turn off debug
//   is singleton the correct pattern?
//   is unique_ptr needed since pointer to object is passed to JANA?


#include "JEventProcessor_CODA_online.h"
#include <JANA/JApplication.h>


using namespace std;
using namespace codaObject;
using namespace jana;


// for testing
static bool debug = true;
int CodaObject::debug = 1;


//----------------------------------------------------------------------------------


/**
 * Initializes the plugin.
 *
 * Creates JEventProcessor/RunObject.  Only called once by JANA framework.
 */
extern "C" {
  void InitPlugin(JApplication *app) {

    string UDL        = "cMsg://localhost/cMsg";
    string name       = "JEventProcessor_codaObject";  // ???
    string descr      = "JEventProcessor coda object"; // ???
    string theSession = "halldsession";


    // initialize the plugin
    InitJANAPlugin(app);


    // get parameters from command line
    gPARMS->SetDefaultParameter("CODAOBJECT:UDL",            UDL, "UDL to use for connecting to cMsg server");
    gPARMS->SetDefaultParameter("CODAOBJECT:NAME",          name, "Name to use for connecting to cMsg server");
    gPARMS->SetDefaultParameter("CODAOBJECT:DESCR",        descr, "Description to use for connecting to cMsg server");
    gPARMS->SetDefaultParameter("CODAOBJECT:SESSION", theSession, "Session name to use");


    // create object, add to JANA event processor list
    app->AddProcessor(new JEventProcessor_CODA_online(UDL,name,descr,theSession));


    if(debug)jout << "leaving InitPlugin" << endl;
  }
}


//----------------------------------------------------------------------------------


/**
 * Constructor sets session, launches statistics thread, automatically starts processing and reporting.
 */
JEventProcessor_CODA_online::JEventProcessor_CODA_online(const string &UDL, const string &name, const string &descr, 
                                                         const string &theSession) 
  : RunObject(UDL,name,descr), done(false) {


  // set session if specified
  if(!theSession.empty())handleSetSession(theSession);


  // launch statistics thread
  stat_thread.reset(new thread(&JEventProcessor_CODA_online::statistics_thread,this));


  // set state and status
  changeState("active");
  changeStatus("ok");


  // start processing
  startProcessing();


  // start reporting
  handleStartReporting(nullptr);

  
  if(debug)jout << "leaving JEventProcessor_CODA_online constructor" << endl;
}


//----------------------------------------------------------------------------------


/**
 * Destructor stops processing.
 */
JEventProcessor_CODA_online::~JEventProcessor_CODA_online() throw() {
  if(debug)jout << "in destructor" << endl;
  done=true;
  stopProcessing();
}


//-----------------------------------------------------------------------------


/**
 * Called during configure transition.
 */
bool JEventProcessor_CODA_online::userConfigure(const string& s) throw(CodaException) {
  if(debug)jout << "in userConfigure" << endl;
  return(true);
}



//-----------------------------------------------------------------------------


/**
 * Called during download transition.
 */
bool JEventProcessor_CODA_online::userDownload(const string& s) throw(CodaException)   {
  if(debug)jout << "in userDownload" << endl;
  return(true);
}


//-----------------------------------------------------------------------------


/**
 * Called during prestart transition.
 */
bool JEventProcessor_CODA_online::userPrestart(const string& s) throw(CodaException)   {
  if(debug)jout << "in userPrestart" << endl;
  return(true);
}


//-----------------------------------------------------------------------------


/**
 * Called during go transition.
 */
bool JEventProcessor_CODA_online::userGo(const string& s) throw(CodaException)   {
  if(debug)jout << "in userGo" << endl;
  return(true);
}


//-----------------------------------------------------------------------------


/**
 * Called during pause transition.
 */
bool JEventProcessor_CODA_online::userPause(const string& s) throw(CodaException)   {
  if(debug)jout << "in userPause" << endl;
  japp->Pause();
  return(true);
}


//-----------------------------------------------------------------------------


/**
 * Called during resume transition.
 */
bool JEventProcessor_CODA_online::userResume(const string& s) throw(CodaException)   {
  if(debug)jout << "in userResum" << endl;
  japp->Resume();
  return(true);
}


//-----------------------------------------------------------------------------


/**
 * Called during end transition.
 */
bool JEventProcessor_CODA_online::userEnd(const string& s) throw(CodaException)   {
  if(debug)jout << "in userEnd" << endl;
  return(true);
}


//-----------------------------------------------------------------------------


/**
 * Called during reset transition.
 */
bool JEventProcessor_CODA_online::userReset(const string& s) throw(CodaException) {
  if(debug)jout << "in userReset" << endl;
  return(true);
}


//-----------------------------------------------------------------------------


/**
 * Called when exit command received.
 */
void JEventProcessor_CODA_online::exit(const string& s) throw(CodaException) {
  if(debug)jout << "in exit" << endl;
  done=true;
  stopProcessing();
  stat_thread->join();
  japp->Quit();
}


//-----------------------------------------------------------------------------


/**
 * Called when unknown message type received.
 */
void JEventProcessor_CODA_online::userMsgHandler(cMsgMessage *msgp, void *userArg) throw(CodaException) {
  unique_ptr<cMsgMessage> msg(msgp);
  jerr << "?JEventProcessor_CODA_online...received unknown message subject,type: "
       << msg->getSubject() << "," << msg->getType() << endl << endl;
}


//-----------------------------------------------------------------------------


/**
 * fills rc/report/status message
 *
 * @param m message
 */
void JEventProcessor_CODA_online::fillReport(cMsgMessage *m) throw() {
  RunObject::fillReport(m);
  m->add("test1",0);
  m->add("test2","fred");
}


//-----------------------------------------------------------------------------


/**
 * Thread in which user must set event number and data count manually, used for monitoring purposes.
 */
const void JEventProcessor_CODA_online::statistics_thread(void) throw() {

  if(debug)jout << "entering JEventProcessor_CODA_online statistics_thread" << endl;

  while(!done) {
    sleep(1);
    eventNumber=0;
    dataCount=0;
    rcConn->setMonitoringString("      <test> This is a test monitoring string </test>\n");
  }
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CODA_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//-------------------------------------------------------------------------------
// end class definition
//----------------------------------------------------------------------------------

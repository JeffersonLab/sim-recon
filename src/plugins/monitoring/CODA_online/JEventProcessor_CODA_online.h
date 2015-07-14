// $Id$
//
//    File: JEventProcessor_CODA_online.h
// Created: Fri Nov  9 10:56:00 EST 2013
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_CODA_online_
#define _JEventProcessor_CODA_online_

#include<memory>
#include<thread>


// for JANA
#include <JANA/JApplication.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JParameterManager.h>


// for coda object package
#include <RunObject.hxx>


using namespace std;
using namespace jana;
using namespace codaObject;


/**
 * Plugin extends RunObject to interface between coda object and JANA event processor for farm management.
 *
 * Uses gPARMS flags CODAOBJECT:UDL, NAME, DESCR, and SESSION.
 */
class JEventProcessor_CODA_online : public JEventProcessor, public RunObject {

 public:
  const char* className(void){return "JEventProcessor_CODA_online";}

  JEventProcessor_CODA_online(void) = delete;
  JEventProcessor_CODA_online(const JEventProcessor_CODA_online&) = delete;
  JEventProcessor_CODA_online& operator=(const JEventProcessor_CODA_online&) = delete;
  JEventProcessor_CODA_online(const string &UDL, const string &name, const string &descr, const string &theSession);
  ~JEventProcessor_CODA_online() throw();


  // for JANA
 private:
  jerror_t fini(void);	/**<Called after last event of last event source has been processed.*/


 public:
  virtual bool userConfigure(const string& s) throw(CodaException) override;
  virtual bool userDownload(const string& s) throw(CodaException) override;
  virtual bool userPrestart(const string& s) throw(CodaException) override;
  virtual bool userGo(const string& s) throw(CodaException) override;
  virtual bool userPause(const string& s) throw(CodaException) override;
  virtual bool userResume(const string& s) throw(CodaException) override;
  virtual bool userEnd(const string& s) throw(CodaException) override;
  virtual bool userReset(const string& s) throw(CodaException) override;
  virtual void exit(const string& s) throw(CodaException) override;

  virtual void userMsgHandler(cMsgMessage *msgp, void *userArg) throw(CodaException) override;


 public:
  virtual const void statistics_thread(void) throw();
  void fillReport(cMsgMessage *m) throw();


 protected:
  unique_ptr<thread> stat_thread;  /**<Sets event and data counts for monitoring.*/
  
 public:
  bool done;            /**<Tells coda object it is done.*/


};

#endif // _JEventProcessor_CODA_online_


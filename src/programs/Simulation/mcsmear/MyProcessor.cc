// Author: David Lawrence  Sat Jan 29 09:37:37 EST 2011
//
//
// MyProcessor.cc
//

#include <iostream>
#include <cmath>
using namespace std;

#include <strings.h>

#include "MyProcessor.h"

#include <JANA/JEvent.h>

#include <HDDM/DEventSourceHDDM.h>
#include <TRACKING/DMCThrown.h>

extern void Smear(hddm_s::HDDM *record);
extern char *OUTFILENAME;

static pthread_mutex_t output_file_mutex;
static pthread_t output_file_mutex_last_owner;

void mcsmear_thread_HUP_sighandler(int sig)
{
   jerr<<" Caught HUP signal for thread 0x"<<hex<<pthread_self()<<dec<<" thread exiting..."<<endl;

   // We use output_file_mutex_owner to keep track (sort of)
   // of which thread has the mutex locked. This mutex is only
   // locked at the end of the evnt method. Once the lock is
   // obtained, this value is set to hold the id of the thread
   // that locked it. It may help in debugging to know the last
   // known thread to have locked the mutex when the signal
   // handler was called
   jerr<<endl;
   jerr<<" Last thread to lock output file mutex: 0x"<<hex<<pthread_self()<<dec<<endl;
   jerr<<" Attempting to unlock mutex to avoid deadlock." <<endl;
   jerr<<" However, the output file is likely corrupt at" <<endl;
   jerr<<" this point and the process should be restarted ..." <<endl;
   jerr<<endl;
   pthread_mutex_unlock(&output_file_mutex);
   pthread_exit(NULL);
}


//------------------------------------------------------------------
// init   -Open output file 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
   // open HDDM file
   ofs = new ofstream(OUTFILENAME);
   if (!ofs->is_open()){
      cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
      exit(-1);
   }
   fout = new hddm_s::ostream(*ofs);
   Nevents_written = 0;

   HDDM_USE_COMPRESSION = true;
   gPARMS->SetDefaultParameter("HDDM:USE_COMPRESSION", HDDM_USE_COMPRESSION,
                          "Turn on/off compression of the output HDDM stream."
                          " Set to \"0\" to turn off (it's on by default)");
   HDDM_USE_INTEGRITY_CHECKS = true;
   gPARMS->SetDefaultParameter("HDDM:USE_INTEGRITY_CHECKS",
                                HDDM_USE_INTEGRITY_CHECKS,
                          "Turn on/off automatic integrity checking on the"
                          " output HDDM stream."
                          " Set to \"0\" to turn off (it's on by default)");

   // enable on-the-fly bzip2 compression on output stream
   if (HDDM_USE_COMPRESSION) {
      jout << " Enabling bz2 compression of output HDDM file stream" 
           << std::endl;
      fout->setCompression(hddm_s::k_bz2_compression);
   }
   else {
      jout << " HDDM compression disabled" << std::endl;
   }

   // enable a CRC data integrity check at the end of each event record
   if (HDDM_USE_INTEGRITY_CHECKS) {
      jout << " Enabling CRC data integrity check in output HDDM file stream"
           << std::endl;
      fout->setIntegrityChecks(hddm_s::k_crc32_integrity);
   }
   else {
      jout << " HDDM integrity checks disabled" << std::endl;
   }

   // We set the mutex type to "ERRORCHECK" so that if the
   // signal handler is called, we can unlock the mutex
   // safely whether we have it locked or not.
   pthread_mutexattr_t attr;
   pthread_mutexattr_init(&attr);
   pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK);
   pthread_mutex_init(&output_file_mutex, NULL);
   
   // pthreads does not provide an "invalid" value for 
   // a pthread_t that we can initialize with. Furthermore,
   // the pthread_t may be a simple as an integer or as
   // a complicated structure. Hence, to make this portable
   // we clear it with bzero.
   bzero(&output_file_mutex_last_owner, sizeof(pthread_t));
   
   return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *loop, int eventnumber)
{
   JEvent& event = loop->GetJEvent();
   JEventSource *source = event.GetJEventSource();
   DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
   if (!hddm_source) {
      cerr << " This program MUST be used with an HDDM file as input!" << endl;
      exit(-1);
   }
   hddm_s::HDDM *record = (hddm_s::HDDM*)event.GetRef();
   if (!record)
      return NOERROR;
   
   // Smear values and add noise hits
   Smear(record);
   
   // Write event to output file
   pthread_mutex_lock(&output_file_mutex);
   output_file_mutex_last_owner = pthread_self();
   *fout << *record;
   Nevents_written++;
   pthread_mutex_unlock(&output_file_mutex);

   return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
   if (fout)
      delete fout;
   if (ofs) {
      ofs->close();
      cout << endl << "Closed HDDM file" << endl;
   }
   cout << " " << Nevents_written << " event written to " << OUTFILENAME
        << endl;
   
   return NOERROR;
}

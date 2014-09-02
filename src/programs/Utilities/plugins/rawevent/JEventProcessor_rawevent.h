// $Id$
//
//    File: JEventProcessor_rawevent.h
// Created: Fri Jun 24 12:05:19 EDT 2011
// Creator: wolin (on Linux stan.jlab.org 2.6.18-194.11.1.el5 x86_64)
//

#ifndef _JEventProcessor_rawevent_
#define _JEventProcessor_rawevent_


// temporary root stuff
#include <TH1.h>


#include <vector>
#include <map>


#include <JANA/JApplication.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>


#include <evioFileChannel.hxx>
#include <evioUtil.hxx>

#include <BCAL/DBCALHit.h>
#include <BCAL/DBCALTDCHit.h>
#include <CDC/DCDCHit.h>
#include <FCAL/DFCALHit.h>
#include <FDC/DFDCHit.h>
#include <START_COUNTER/DSCHit.h>
#include <TOF/DTOFHit.h>
#include <TAGGER/DTAGMHit.h>
#include <TAGGER/DTAGHHit.h>


using namespace std;
using namespace jana;
using namespace evio;


// holds crate/slot/channel
typedef struct {
  int crate;
  int slot;
  int channel;
} cscVal;
typedef const cscVal &cscRef;



//----------------------------------------------------------------------------


class JEventProcessor_rawevent : public jana::JEventProcessor {

	public:
		JEventProcessor_rawevent();
		~JEventProcessor_rawevent();
		const char* className(void){return "JEventProcessor_rawevent";}


	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);


                // these routines read and fill the translation tables
                static void readTranslationTable(void);
                static void StartElement(void *userData, const char *xmlname, const char **atts);
                static void EndElement(void *userData, const char *xmlname);


                // these routines access the translation tables
                cscRef DTOFHitTranslationADC(const DTOFHit* hit) const;
                cscRef DTOFHitTranslationTDC(const DTOFHit* hit) const;

                cscRef DBCALHitTranslationADC(const DBCALHit* hit) const;
                //cscRef DBCALHitTranslationTDC(const DBCALHit* hit) const;
                cscRef DBCALHitTranslationTDC(const DBCALTDCHit* hit) const;

                cscRef DFCALHitTranslationADC(const DFCALHit* hit) const;

                cscRef DFDCAnodeHitTranslation(const DFDCHit* hit) const;
                cscRef DFDCCathodeHitTranslation(const DFDCHit* hit) const;

                cscRef DCDCHitTranslationADC(const DCDCHit* hit) const;

                cscRef DSTHitTranslationADC(const DSCHit* hit) const;
                cscRef DSTHitTranslationTDC(const DSCHit* hit) const;

                cscRef DTAGMHitTranslationADC(const DTAGMHit* hit) const;
                cscRef DTAGMHitTranslationTDC(const DTAGMHit* hit) const;

                cscRef DTAGHHitTranslationADC(const DTAGHHit* hit) const;
                cscRef DTAGHHitTranslationTDC(const DTAGHHit* hit) const;
};

#endif // _JEventProcessor_rawevent_


//----------------------------------------------------------------------------

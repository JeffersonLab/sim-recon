// $Id$
//
//    File: JEventProcessor_rawevent.h
// Created: Fri Jun 24 12:05:19 EDT 2011
// Creator: wolin (on Linux stan.jlab.org 2.6.18-194.11.1.el5 x86_64)
//

#ifndef _JEventProcessor_rawevent_
#define _JEventProcessor_rawevent_


#include <vector>
#include <map>


#include <JANA/JApplication.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>


#include <evioFileChannel.hxx>
#include <evioUtil.hxx>


#include "FCAL/DFCALHit.h"
#include "BCAL/DBCALHit.h"
#include "TOF/DTOFRawHit.h"
#include "CDC/DCDCHit.h"
#include "FDC/DFDCHit.h"
#include "START_COUNTER/DSCHit.h"
#include "TAGGER/DTagger.h"


#include<boost/tuple/tuple.hpp>


using namespace std;
using namespace jana;
using namespace evio;
using namespace boost;


typedef tuple<int,int,int>  cscVal;
typedef const tuple<int,int,int> &cscRef;



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
                static void startElement(void *userData, const char *xmlname, const char **atts);


                // these routines access the translation tables
                cscRef DTOFRawHitTranslationADC(const DTOFRawHit* hit) const;
                cscRef DTOFRawHitTranslationTDC(const DTOFRawHit* hit) const;

                cscRef DBCALHitTranslationADC(const DBCALHit* hit) const;
                cscRef DBCALHitTranslationTDC(const DBCALHit* hit) const;

                cscRef DFCALHitTranslationADC(const DFCALHit* hit) const;

                cscRef DFDCAnodeHitTranslation(const DFDCHit* hit) const;
                cscRef DFDCCathodeHitTranslation(const DFDCHit* hit) const;

                cscRef DCDCHitTranslationADC(const DCDCHit* hit) const;

                cscRef DSCHitTranslationADC(const DSCHit* hit) const;
                cscRef DSCHitTranslationTDC(const DSCHit* hit) const;

                cscRef DTaggerTranslationADC(const DTagger* hit) const;
                cscRef DTaggerTranslationTDC(const DTagger* hit) const;
};

#endif // _JEventProcessor_rawevent_


//----------------------------------------------------------------------------

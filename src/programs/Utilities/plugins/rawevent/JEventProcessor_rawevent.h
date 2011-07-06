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


typedef tuple<int,int,int> cscVal;



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


                // these routines access the translation table
                void readTranslationTable(void);
                static void startElement(void *userData, const char *xmlname, const char **atts);

                cscVal DTOFRawHitTranslationADC(const DTOFRawHit* hit);
                cscVal DTOFRawHitTranslationTDC(const DTOFRawHit* hit);

                cscVal DBCALHitTranslationADC(const DBCALHit* hit);
                cscVal DBCALHitTranslationTDC(const DBCALHit* hit);

                cscVal DFCALHitTranslationADC(const DFCALHit* hit);

                cscVal DFDCHitTranslation(const DFDCHit* hit);

                cscVal DCDCHitTranslationADC(const DCDCHit* hit);

                cscVal DSCHitTranslationADC(const DSCHit* hit);
                cscVal DSCHitTranslationTDC(const DSCHit* hit);

                cscVal DTaggerTranslationADC(const DTagger* hit);
                cscVal DTaggerTranslationTDC(const DTagger* hit);


                // maps
                static map< string, tuple<int,int,int> >   fcalMap;
                static map< string, tuple<int,int,int> >   bcalMap;
                static map< string, tuple<int,int,int> >   cdcMap;
                static map< string, tuple<int,int,int> >   scMap;
                static map< string, tuple<int,int,int> >   hodoscopeMap;
                static map< string, tuple<int,int,int> >   fdcCathodeMap;
                static map< string, tuple<int,int,int> >   fdcAnodeMap;
                static map< string, tuple<int,int,int> >   tofMap;
                static map< string, tuple<int,int,int> >   microscopeMap;
};

#endif // _JEventProcessor_rawevent_


//----------------------------------------------------------------------------

// $Id$
//
//    File: DFDCPseudo_factory_WIRESONLY.h
// Created: Fri Nov  9 09:57:12 EST 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DFDCPseudo_factory_WIRESONLY_
#define _DFDCPseudo_factory_WIRESONLY_

#include <JANA/JFactory.h>
using namespace jana;

#include "DFDCPseudo.h"

class DFDCPseudo_factory_WIRESONLY:public JFactory<DFDCPseudo>{
	public:
		DFDCPseudo_factory_WIRESONLY(){};
		~DFDCPseudo_factory_WIRESONLY(){};
		const char* Tag(void){return "WIRESONLY";}

	private:
		//jerror_t init(void);						///< Called once at program start.
		//jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void MakePseudo(const DFDCHit *hit, const DFDCWire *wire, const DVector3 &pos);

};

#endif // _DFDCPseudo_factory_WIRESONLY_


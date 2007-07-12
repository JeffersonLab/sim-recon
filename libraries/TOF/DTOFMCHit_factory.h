// $Id: DTOFMCHit_factory.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DTOFMCHit_factory.h
// Created: Mon Jul  9 16:21:12 EDT 2007
// Creator: B.Zihlmann 
//

#ifndef _DTOFMCHit_factory_
#define _DTOFMCHit_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DTOFMCHit.h"

class DTOFMCHit_factory:public JFactory<DTOFMCHit>{
        public:
                DTOFMCHit_factory(){};
                ~DTOFMCHit_factory(){};
                const string toString(void);

        protected:
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);  ///< Called every event.
};

#endif // _DTOFMCHit_factory_


// $Id$
//
//    File: DTOFMCResponse_factory.h
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DTOFMCResponse_factory_
#define _DTOFMCResponse_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DTOFMCResponse.h"

class DTOFMCResponse_factory:public JFactory<DTOFMCResponse>{
	public:
		DTOFMCResponse_factory(){};
		~DTOFMCResponse_factory(){};
		const string toString(void);


	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DTOFMCResponse_factory_


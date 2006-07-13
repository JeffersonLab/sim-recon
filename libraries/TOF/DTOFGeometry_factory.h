// $Id$
//
//    File: DTOFGeometry_factory.h
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DTOFGeometry_factory_
#define _DTOFGeometry_factory_

#include "JANA/JFactory.h"
#include "DTOFGeometry.h"

class DTOFGeometry_factory:public JFactory<DTOFGeometry>{
	public:
		DTOFGeometry_factory(){};
		~DTOFGeometry_factory(){};
		const string toString(void);


	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		jerror_t init(void);
};

#endif // _DTOFGeometry_factory_


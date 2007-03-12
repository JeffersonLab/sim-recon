// $Id: DFCALPhoton_factory.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALPhoton_factory.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALPhoton_factory_
#define _DFCALPhoton_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DFCALPhoton.h"
#include "DFCALCluster.h"


class DFCALPhoton_factory:public JFactory<DFCALPhoton>{
	public:
		DFCALPhoton_factory();
		~DFCALPhoton_factory(){};
		const string toString(void);
	
	private:
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

};


#endif // _DFCALPhoton_factory_


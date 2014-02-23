// $Id$
//
//    File: DTrackHitSelector_factory.h
// Created: Thu Feb  5 13:34:58 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackHitSelector_factory_
#define _DTrackHitSelector_factory_

#include <JANA/JFactory.h>
#include "DTrackHitSelector.h"

class DTrackHitSelector_factory:public jana::JFactory<DTrackHitSelector>{
	public:
		DTrackHitSelector_factory(){};
		~DTrackHitSelector_factory(){};


	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// This is a trivial class that simply implements the 
			// ALT1 tagged factory as the default. It is here so
			// that the default can be changed easily by simply
			// changing the tag here or on the command line.
			vector<const DTrackHitSelector*> selectors;
			loop->Get(selectors, "ALT1");
			for(unsigned int i=0; i< selectors.size(); i++){
				_data.push_back(const_cast<DTrackHitSelector*>(selectors[i]));
			}
			SetFactoryFlag(NOT_OBJECT_OWNER);

			return NOERROR;
		}
};

#endif // _DTrackHitSelector_factory_


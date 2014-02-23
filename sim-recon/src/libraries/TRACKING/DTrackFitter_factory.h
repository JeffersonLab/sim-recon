// $Id$
//
//    File: DTrackFitter_factory_ALT1.h
// Created: Mon Sep  1 10:29:51 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#ifndef _DTrackFitter_factory_
#define _DTrackFitter_factory1_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitter.h>

class DTrackFitter_factory:public jana::JFactory<DTrackFitter>{
	public:
		DTrackFitter_factory(){};
		~DTrackFitter_factory(){};

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// This is a trivial class that simply implements the
			// ALT1 tagged factory as the default. It is here so 
			// that the default can be changed easily by simply
			// changing the tag here or on the command line.
			vector<const DTrackFitter*> fitters;
			loop->Get(fitters, "KalmanSIMD_ALT1");
			for(unsigned int i=0; i< fitters.size(); i++){
				_data.push_back(const_cast<DTrackFitter*>(fitters[i]));
			}
			SetFactoryFlag(NOT_OBJECT_OWNER);

			return NOERROR;
		}
};

#endif // _DTrackFitter_factory_


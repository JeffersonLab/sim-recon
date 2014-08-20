// $Id$
//
//    File: DTrackFinder_factory.h
// Created: Mon Aug 18 11:00:33 EDT 2014
// Creator: staylor (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DTrackFinder_factory_
#define _DTrackFinder_factory_

#include <JANA/JFactory.h>
#include "DTrackFinder.h"

class DTrackFinder_factory:public jana::JFactory<DTrackFinder>{
 public:
  DTrackFinder_factory(){};
  ~DTrackFinder_factory(){};

 private:
  jerror_t evnt(jana::JEventLoop *loop, int eventnumber){
    
    SetFactoryFlag(PERSISTANT);
    ClearFactoryFlag(WRITE_TO_OUTPUT);
   
    _data.push_back(new DTrackFinder(loop));   

    return NOERROR;
  }
};

#endif // _DTrackFinder_factory_


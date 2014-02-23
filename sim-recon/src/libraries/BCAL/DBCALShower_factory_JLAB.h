// $Id: DBCALShower_factory_JLAB.h 
//
//    File: DBCALShower_factory_JLAB.h
// Created: Mon Mar 18 09:42:29 EDT 2013
// Creator: Benedikt Zihlmann version 0.1
//


#ifndef _DBCALShower_factory_JLAB_
#define _DBCALShower_factory_JLAB_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <BCAL/DBCALShower.h>
#include <DBCALClump_factory.h>


/// The showers produced here are based on the output of the DBCALClump_factory

class DBCALShower_factory_JLAB:public JFactory<DBCALShower>{
  
 public:
  
  DBCALShower_factory_JLAB();
  ~DBCALShower_factory_JLAB(){};
  
  const char* Tag(void){return "JLAB";}
  
 private:
  jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
  
};

#endif // _DBCALShower_factory_JLAB_


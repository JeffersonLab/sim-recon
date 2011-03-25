/*
 *  DBCALPhoton_factory.h
 *  Hall D
 *
 *  Created by Matthew Shepherd on 7/23/07.
 *
 */

#ifndef _DBCALPhoton_factory_
#define _DBCALPhoton_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
using namespace jana;


#include <BCAL/DBCALShower.h>
#include <BCAL/DBCALPhoton.h>

class DVertex;

class DBCALPhoton_factory : public JFactory< DBCALPhoton > { 
    
public:
    
  DBCALPhoton_factory();
  ~DBCALPhoton_factory(){}
    
private:

  jerror_t init();
  
  jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt( JEventLoop *loop, int eventnumber );
  DBCALPhoton* MakeDBCALPhoton(const DBCALShower* shower, const DVertex *vertex);

  double m_zTarget;
  double m_bcalIR;
  
  int USE_KLOE;
};

#endif

/*
 *  DFCALPhoton_factory.h
 *  Hall D
 *
 */

#ifndef _DFCALPhoton_factory_
#define _DFCALPhoton_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
using namespace jana;


#include <FCAL/DFCALShower.h>
#include <FCAL/DFCALPhoton.h>

class DVertex;

class DFCALPhoton_factory : public JFactory< DFCALPhoton > { 
    
public:
    
  DFCALPhoton_factory();
  ~DFCALPhoton_factory(){}
    
private:

  
 // jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt( JEventLoop *loop, int eventnumber );
  DFCALPhoton* MakeDFCALPhoton(const DFCALShower* shower, const DVertex *vertex);

  double m_zTarget;

};

#endif

#ifndef _DBCALShower_factory_IU_
#define _DBCALShower_factory_IU_

/*
 *  DBCALShower_factory_IU.h
 *  (formerly DBCALShower_factory.h)
 *
 *  Created by Matthew Shepherd on 3/24/11.
 *
 */

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALShower.h"

class DBCALShower_factory_IU : public JFactory< DBCALShower > {
  
public:
  
  DBCALShower_factory_IU();
  ~DBCALShower_factory_IU(){}

  const char* Tag(void){return "IU";}
  
private:
  
  jerror_t evnt(JEventLoop *loop, uint64_t eventnumber);
  jerror_t brun(JEventLoop *loop, int32_t runnumber);

  double m_zTarget;

// energy calibration parameters
  
  float m_scaleZ_p0;
  float m_scaleZ_p1;
  float m_scaleZ_p2;
  float m_scaleZ_p3;
  
  float m_nonlinZ_p0;
  float m_nonlinZ_p1;
  float m_nonlinZ_p2;
  float m_nonlinZ_p3;
  
};

#endif

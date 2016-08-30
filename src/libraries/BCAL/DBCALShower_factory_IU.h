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

  double LOAD_CCDB_CONSTANTS;
  double energy_cutoff;
  double linear_intercept;
  double linear_slope;
  double exponential_parm0;
  double exponential_parm1;
  double exponential_parm2;

  double m_zTarget;

};

#endif

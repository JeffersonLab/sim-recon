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

#include "TH2F.h"
#include <DMatrixDSym.h>

class DBCALShower_factory_IU : public JFactory< DBCALShower > {
  
public:
  
  DBCALShower_factory_IU();
  ~DBCALShower_factory_IU(){}

  const char* Tag(void){return "IU";}
  jerror_t LoadCovarianceLookupTables();
  jerror_t FillCovarianceMatrix(DBCALShower* shower);

private:
  
  jerror_t evnt(JEventLoop *loop, uint64_t eventnumber);
  jerror_t brun(JEventLoop *loop, int32_t runnumber);
  //jerror_t CreateCovarianceMatrix();

  int VERBOSE;
  string COVARIANCEFILENAME;
  TH2F *CovarianceLookupTable[5][5];

  double LOAD_CCDB_CONSTANTS;
  double energy_cutoff;
  double linear_intercept;
  double linear_slope;
  double exponential_param0;
  double exponential_param1;
  double exponential_param2;

  double m_zTarget;

};

#endif

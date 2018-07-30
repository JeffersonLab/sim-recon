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
#include "BCAL/DBCALGeometry.h"

#include "TH2F.h"
#include <DMatrixDSym.h>

class DBCALShower_factory_IU : public JFactory< DBCALShower > {
  
public:
  
  DBCALShower_factory_IU();
  ~DBCALShower_factory_IU(){}

  const char* Tag(void){return "IU";}
  jerror_t LoadCovarianceLookupTables(JEventLoop *eventLoop);
  jerror_t FillCovarianceMatrix(DBCALShower* shower);

private:
  
  jerror_t evnt(JEventLoop *loop, uint64_t eventnumber);
  jerror_t brun(JEventLoop *loop, int32_t runnumber);
  jerror_t erun(void);
  //jerror_t CreateCovarianceMatrix();

  int VERBOSE;
  string COVARIANCEFILENAME;
  TH2F *CovarianceLookupTable[5][5];

  double LOAD_CCDB_CONSTANTS;
  double const_term;
  double first_term_scale_factor;
  double first_exp_param0;
  double first_exp_param1;
  double second_term_scale_factor;
  double second_exp_const_term;
  double second_exp_scale_factor;
  double second_exp_param0;
  double second_exp_param1;

  double m_zTarget;

  const DBCALGeometry *dBCALGeom;

};

#endif

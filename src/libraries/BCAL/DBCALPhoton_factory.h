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
    ~DBCALPhoton_factory(){};

    
private:

    jerror_t brun(JEventLoop *loop, int runnumber);
    jerror_t evnt( JEventLoop *loop, int eventnumber );
	 DBCALPhoton* MakeDBCALPhoton(const DBCALShower* shower, const DVertex *vertex);
	 DBCALPhoton* MergeDBCALPhotons(vector<DBCALPhoton*> &showers);
    
	 double MIN_CLUSTER_SEPARATION_XY;
	 double MIN_CLUSTER_SEPARATION_Z;
	 
    double m_scaleZ_p0LT;
    double m_scaleZ_p1LT;
    double m_scaleZ_p2LT;
   
    double m_scaleZ_p0GE;
    double m_scaleZ_p1GE;
    double m_scaleZ_p2GE;
    double m_scaleZ_p3GE;

    double m_scaleZ_p0;
    double m_scaleZ_p1;
    double m_scaleZ_p2;
    double m_scaleZ_p3;
    double m_scaleZ_p4;
   
    double m_nonlinZ_p0LT;
    double m_nonlinZ_p1LT;
    double m_nonlinZ_p2LT;
    double m_nonlinZ_p3LT;

    double m_nonlinZ_p0GE;
    double m_nonlinZ_p1GE;
    double m_nonlinZ_p2GE;
    double m_nonlinZ_p3GE;

    double m_nonlinZ_p0;
    double m_nonlinZ_p1;
    double m_nonlinZ_p2;
    double m_nonlinZ_p3;
    double m_nonlinZ_p4;
 
    double m_linZ_p0LT;
    double m_linZ_p1LT;
    double m_linZ_p2LT;
    double m_linZ_p3LT;
  
    double m_linZ_p0GE;
    double m_linZ_p1GE;
    double m_linZ_p2GE;
    double m_linZ_p3GE;
 
    double m_linZ_p0;
    double m_linZ_p1;
    double m_linZ_p2;
    double m_linZ_p3;   

    double m_bcalIR;
    double m_zTarget;

    double nonlin;
    double scale;
    double lin;
};

#endif

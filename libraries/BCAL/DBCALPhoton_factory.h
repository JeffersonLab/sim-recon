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


#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALPhoton.h"

class DBCALPhoton_factory : public JFactory< DBCALPhoton > { 
    
public:
    
    DBCALPhoton_factory();
    ~DBCALPhoton_factory(){};

    
private:

    jerror_t brun(JEventLoop *loop, int runnumber);
    jerror_t evnt( JEventLoop *loop, int eventnumber );
    
    double m_scaleZ_p0LT;
    double m_scaleZ_p1LT;
    double m_scaleZ_p2LT;
   
    double m_scaleZ_p0GE;
    double m_scaleZ_p1GE;
    double m_scaleZ_p2GE;
    
    double m_nonlinZ_p0LT;
    double m_nonlinZ_p1LT;
    double m_nonlinZ_p2LT;
    double m_nonlinZ_p3LT;

  
    double m_nonlinZ_p0GE;
    double m_nonlinZ_p1GE;
    double m_nonlinZ_p2GE;
    double m_nonlinZ_p3GE;
    
    double m_bcalIR;
    double m_zTarget;

    double nonlin;
    double scale;

};

#endif

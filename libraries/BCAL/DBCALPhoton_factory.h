/*
 *  DBCALPhoton_factory.h
 *  Hall D
 *
 *  Created by Matthew Shepherd on 7/23/07.
 *
 */

#ifndef _DBCALPhoton_factory_
#define _DBCALPhoton_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALPhoton.h"

class DBCALPhoton_factory : public JFactory< DBCALPhoton > { 
    
public:
    
    DBCALPhoton_factory();
    ~DBCALPhoton_factory(){};

    const string toString( void );
    
private:

    jerror_t DBCALPhoton_factory::brun(JEventLoop *loop, int runnumber);
    jerror_t evnt( JEventLoop *loop, int eventnumber );
    
    double m_scaleZ_p0;
    double m_scaleZ_p1;
    double m_scaleZ_p2;
    double m_scaleZ_p3;
    
    double m_nonlinZ_p0;
    double m_nonlinZ_p1;
    double m_nonlinZ_p2;
    
    double m_bcalIR;
    double m_zTarget;
};

#endif

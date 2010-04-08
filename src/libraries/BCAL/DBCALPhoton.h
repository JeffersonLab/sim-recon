/*
 *  DBCALPhoton.h
 *  Hall D
 *
 *  Created by Matthew Shepherd on 7/23/07.
 *
 *  This class takes final BCALShowers and calibrates them
 *  to create photons.
 */

#ifndef _DBCALPhoton_
#define _DBCALPhoton_

#include "DVector3.h"
#include "DLorentzVector.h"

#include "BCAL/DBCALShower.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DBCALPhoton : public JObject {

public:

    JOBJECT_PUBLIC( DBCALPhoton );
    
    DBCALPhoton(){}
    ~DBCALPhoton(){}
        
    inline DVector3 showerPosition() const { return m_position; }
    inline double showerTime() const { return m_time; }
    inline DLorentzVector lorentzMomentum() const { return m_p4; }

    inline DVector3 fitLayPoint()    const { return m_fitLayPoint;    }
    inline DVector3 fitLayPointErr() const { return m_fitLayPointErr; }
    inline DVector3 fitLaySlope()    const { return m_fitLayPoint;    }
    inline DVector3 fitLaySlopeErr() const { return m_fitLaySlopeErr; }
    
    void setShowerPosition( const DVector3& vec3 ) { m_position = vec3; }
    void setShowerTime( const double t ) { m_time = t; }
    void setLorentzMomentum( const DLorentzVector& p4 ) { m_p4 = p4; }
    
    void setFitLayPoint( const DVector3& vec3 )    { m_fitLayPoint = vec3;    }
    void setFitLayPointErr( const DVector3& vec3 ) { m_fitLayPointErr = vec3; }
    void setFitLaySlope( const DVector3& vec3 )    { m_fitLaySlope = vec3;    }
    void setFitLaySlopeErr( const DVector3& vec3 ) { m_fitLaySlopeErr = vec3; }

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "x", "%3.2f", m_position.X());
			AddString(items, "y", "%3.2f", m_position.Y());
			AddString(items, "z", "%3.2f", m_position.Z());
			AddString(items, "px", "%5.2f", m_p4.Px());
			AddString(items, "py", "%5.2f", m_p4.Py());
			AddString(items, "pz", "%5.2f", m_p4.Pz());
			AddString(items, "E" , "%5.2f", m_p4.E());
			AddString(items, "t" , "%5.2f", m_time);
		}

private:
      
    DVector3 m_position;
    double m_time;
    DLorentzVector m_p4;
    
    DVector3 m_fitLayPoint;
    DVector3 m_fitLayPointErr;
    DVector3 m_fitLaySlope;
    DVector3 m_fitLaySlopeErr;
};


#endif // _DBCALPhoton_


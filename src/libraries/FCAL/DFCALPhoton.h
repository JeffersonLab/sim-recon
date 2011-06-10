/*
 *  DFCALPhoton.h
 *  Hall D
 *
 */

#ifndef _DFCALPhoton_
#define _DFCALPhoton_

#include "DVector3.h"
#include "DLorentzVector.h"

#include "FCAL/DFCALShower.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DFCALPhoton : public JObject {

public:

    JOBJECT_PUBLIC( DFCALPhoton );
    
    DFCALPhoton(){}
    ~DFCALPhoton(){}
        
    inline DVector3 showerPosition() const { return m_position; }
    inline DVector3 showerPositionErr() const { return m_positionErr; }
    inline double showerTime() const { return m_time; }
    inline DLorentzVector lorentzMomentum() const { return m_p4; }
    
    void setShowerPosition( const DVector3& vec3 ) { m_position = vec3; }
    void setShowerPositionErr( const DVector3& vec3 ) { m_positionErr = vec3; }
    void setShowerTime( const double t ) { m_time = t; }
    void setLorentzMomentum( const DLorentzVector& p4 ) { m_p4 = p4; }
    
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
    DVector3 m_positionErr;
    double m_time;
    DLorentzVector m_p4;
};


#endif // _DFCALPhoton_


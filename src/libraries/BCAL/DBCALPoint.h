#ifndef _DBCALPoint_
#define _DBCALPoint_

/*
 *  DBCALPoint.h
 *
 *  Created by Matthew Shepherd on 3/13/11.
 *
 */

#include "BCAL/DBCALHit.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

using namespace jana;

class DBCALPoint : public JObject {

public:

  JOBJECT_PUBLIC(DBCALPoint);
  
  DBCALPoint(){}
  
  // this constructor uses two hits to obtain a local z position
  DBCALPoint( const DBCALHit& hit1, const DBCALHit& hit2 );
  
  // this constructor is helpful for single-ended hits when
  // z is known -- z measured with respect to target
  DBCALPoint( const DBCALHit& hit, float zTarget );
  
  float E() const { return m_E; }
  float t() const { return m_t; }

  // assuming a photon, this gives time at the inner radius of BCAL
  // by extrapolating back on path from center of cell to target
  float tInnerRadius() const;
    
  // spherical coordinates, origin at the target
  float rho()   const { return m_rho; }
  float sigRho() const { return m_sig_rho; }
  
  float theta() const { return m_theta; }
  float sigTheta() const { return m_sig_theta; }

  float phi() const { return m_phi; }
  float sigPhi() const { return m_sig_phi; }

  // cylindrical coordinates, origin at the target
  float z()   const { return m_z; }
  float r()   const { return m_r; }

  int module() const {return m_module;}
  int layer() const {return m_layer;}
  int sector() const {return m_sector;}
                      
  // these rotate point in phi by 2 pi -- effectively leaves point unchanged
  void add2Pi() const;
  void sub2Pi() const;

  void toStrings(vector<pair<string,string> > &items) const {
    AddString(items, "E(GeV)", "%5.3f", m_E);
    AddString(items, "t(ns)", "%5.1f", m_t);
    AddString(items, "z(cm)", "%5.1f", m_z);
    AddString(items, "r(cm)", "%5.1f", m_r);
    AddString(items, "phi", "%5.3f", m_phi);
    AddString(items, "module", "%i", m_module);
    AddString(items, "layer", "%i", m_layer);
    AddString(items, "sector", "%i", m_sector);
    
  }
  
private:
  
  void convertCylindricalToSpherical();
  
  float m_E;
  float m_t;

  int m_module, m_layer, m_sector;
  
  // cylindrical coordinate locations
  float m_zLocal;
  float m_z, m_sig_z;
  float m_r, m_sig_r;
  float m_phi, m_sig_phi;
  
  // spherical coordinate locations
  float m_rho, m_sig_rho;
  float m_theta, m_sig_theta;

};


#endif

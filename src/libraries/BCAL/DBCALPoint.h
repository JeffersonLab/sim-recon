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

  // To take into account TDC hits, timewalk corrections must be done before
  // determining z-position. You should pass in corrected times to this
  // function. If no TDC hits are present, tUp and tUp_ADC should both be
  // specified with the same value
  DBCALPoint(int module, int layer, int sector, float EUp, float EDown, float tUp, float tDown, float tUp_ADC, float tDown_ADC);
  
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
  float sigZ() const { return m_sig_z; }
  float r()   const { return m_r; }

  int module() const {return m_module;}
  int layer() const {return m_layer;}
  int sector() const {return m_sector;}
                      
  // these rotate point in phi by 2 pi -- effectively leaves point unchanged
  void add2Pi() const;
  void sub2Pi() const;

  //Energies mesaured at the two BCAL ends
  //These energies are in units of "GeV", but have not been corrected for
  //attenuation.
  //They are saved so that any calibration is done only once
  //(in DBCALPoint_factory) and any class upstream can just grab the corrected
  //values from here
  float EUp() const {return m_EUp;}
  float EDown() const {return m_EDown;}

  //Hit times (after all corrections/calibration) using only information
  //from fADC hits (i.e. DBCALHit objects).
  //Again, these are saved so that upstream objects don't have to redo any
  //calibrations.
  float tUp_ADC() const {return m_tUp_ADC;}
  float tDown_ADC() const {return m_tDown_ADC;}

  //Hit times (after corrections/calibration) using the most precise
  //information available
  //If a TDC hit exists, the (corrected) TDC time will be returned.
  //If no TDC hit exists, this function should return the same values
  //as tUp_ADC()/tDown_ADC()
  float tUp() const {return m_tUp;}
  float tDown() const {return m_tDown;}

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

  float m_EUp, m_EDown;
  float m_tUp, m_tDown;
  float m_tUp_ADC, m_tDown_ADC;

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

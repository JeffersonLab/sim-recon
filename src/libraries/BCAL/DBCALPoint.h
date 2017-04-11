#ifndef _DBCALPoint_
#define _DBCALPoint_

/*
 *  DBCALPoint.h
 *
 *  Created by Matthew Shepherd on 3/13/11.
 *
 */

#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALUnifiedHit.h"
#include "BCAL/DBCALGeometry.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

using namespace jana;

/** 
   This object gives a z position to BCAL hits using timing information.
   The z position is relative to the target center because higher objects use 
   the position of hits relative to the target in spherical coordinates.

*/

class DBCALPoint : public JObject {

public:

  JOBJECT_PUBLIC(DBCALPoint);
  
  // this constructor uses two hits to obtain a local z position
  DBCALPoint(const DBCALUnifiedHit& hit1, const DBCALUnifiedHit& hit2, double z_target_center, 
  			 double attenutation_length, double c_effective, double track_p0, double track_p1, double track_p2, 
  			 const DBCALGeometry *locGeom);
  
  float E() const { return m_E; }
  float E_US() const { return m_E_US; }  ///< Return the attenuation corrected Energy of US Hit
  float E_DS() const { return m_E_DS; }  ///< Return the attenuation corrected Energy of DS Hit
  float t() const { return m_t; }
  float t_US() const { return m_t_US; }  ///< Return the time of US Hit
  float t_DS() const { return m_t_DS; }  ///< Return the time of DS Hit

  // assuming a photon, this gives time at the inner radius of BCAL
  // by extrapolating back on path from center of cell to target
  float tInnerRadius() const;
    
  // spherical coordinates, origin at the center of the target
  float rho()   const { return m_rho; }
  float sigRho() const { return m_sig_rho; }
  
  float theta() const { return m_theta; }
  float sigTheta() const { return m_sig_theta; }

  float phi() const { return m_phi; }
  float sigPhi() const { return m_sig_phi; }

  // cylindrical coordinates, origin at the center of the target
  float z()   const { return m_z; }
  float sigZ() const { return m_sig_z; }
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
  
  float m_E;                     ///< Energy of the Point used in higher objects
  float m_E_US;                  ///< Attenuation corrected Energy of US Hit that contributed to the Point
  float m_E_DS;                  ///< Attenuation corrected Energy of DS Hit that contributed to the Point
  float m_t;                     ///< Arrival time
  float m_t_US;                  ///< Time of DS Hit that contributed to the Point
  float m_t_DS;                  ///< Time of DS Hit that contributed to the Point

  int m_module, m_layer, m_sector;

  // cylindrical coordinate locations
  float m_zGlobal;                ///< z-coordinate relative to the beginning of the BCAL
  float m_zLocal;                ///< z-coordinate relative to the center of BCAL
  float m_z, m_sig_z;            ///< z-coordinate relative to the center of the target
  float m_r, m_sig_r;            ///< distance from beam axis
  float m_phi, m_sig_phi;        ///< azimuthal angle
  
  // spherical coordinate locations
  float m_rho, m_sig_rho;        ///< spherical distance wrt target center
  float m_theta, m_sig_theta;    ///< polar angle wrt target center
  
  const DBCALGeometry *m_BCALGeom;

};


#endif

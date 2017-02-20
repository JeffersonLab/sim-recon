//***************************************************
// DFDCGeometry - temporary geometry class for FDC.
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
//***************************************************

#ifndef DFDCGEOMETRY_H
#define DFDCGEOMETRY_H

#include "HDGEOMETRY/DGeometry.h"

#define FDC_NUM_LAYERS			24
#define FDC_NUM_PLANES                  72
#define FDC_ACTIVE_RADIUS		48.5

//----- These were cut from HDGeant/hitFDC.c -----
#define DRIFT_SPEED           .0055
#define CATHODE_ROT_ANGLE     1.309  // 75 degrees
//#define CATHODE_ROT_ANGLE     0.7854 // 45 degrees
#define WIRE_DEAD_ZONE_RADIUS 3.2
#define ANODE_CATHODE_SPACING 0.5
#define TWO_HIT_RESOL         250.
#define WIRES_PER_PLANE       96
#define WIRE_SPACING          1.0
#define U_OF_WIRE_ZERO        (-((WIRES_PER_PLANE-1)*WIRE_SPACING)/2)
#define STRIPS_PER_PLANE      192 
#define STRIP_SPACING         0.500
#define U_OF_STRIP_ZERO		   (-((STRIPS_PER_PLANE-1)*STRIP_SPACING)/2)
#define STRIP_ZERO_OFFSET     -95.5
#define STRIP_GAP             0.1
#define MAX_HITS             100
//#define K2                  1.15
#define STRIP_NODES           3
#define THRESH_KEV           1.
#define THRESH_STRIPS        5.   /* mV */
#define ELECTRON_CHARGE 1.6022e-4 /* fC */
//------------------------------------------------

#include <math.h>

#include <JANA/JObject.h>
#include "DFDCHit.h"
#include "DFDCWire.h"
#include "DFDCCathode.h"

///
/// class DFDCGeometry: definition of a class providing basic geometry
/// methods for FDC reconstruction.
///
class DFDCGeometry : public JObject {
 public:
  DFDCGeometry(){}
  ~DFDCGeometry(){}
  
  ///
  /// DFDCGeometry::gLayer():
  /// Compute the global layer (detection layer 1-24) number based on module and layer
  ///
  static inline int gLayer(const DFDCHit* h) {
    return 3*(h->module - 1) + (h->layer - 1) + 1;
  }		
  
  ///
  /// DFDCGeometry::gPlane():
  /// Compute the global plane (1-74) number based on module, layer, and plane
  ///
  static inline int gPlane(const DFDCHit* h) {
    return 9*(h->module-1) + 3*(h->layer-1) + (h->plane-1) + 1;
  }
  
  ///		
  /// DFDCGeometry::getWireR():
  /// Return X coordinate of a wire
  /// 
  static inline float getWireR(const DFDCHit* h) {
    return WIRE_SPACING*(h->element-1)+U_OF_WIRE_ZERO;
    
  }
  
  ///
  /// DFDCGeometry::getStripR():
  /// Return coordinate in U or V space of a strip
  ///
  static inline float getStripR(const DFDCHit* h) {
    return STRIP_SPACING*(h->element-1)+U_OF_STRIP_ZERO;
  }

  ///
  /// DFDCGeometry::getXLocalStrips()
  ///  Compute the x-coordinate in the layer coordinate system
  /// from the strip data.
  ///
  static inline float getXLocalStrips(float u, float v){
    return STRIP_SPACING*(u-1 + v-1 - (STRIPS_PER_PLANE-1))
      /2./cos(CATHODE_ROT_ANGLE);
  }  
  static inline float getXLocalStrips(float u, float phi_u,float v, float phi_v){
    float cosPhiU=cos(phi_u);
    float cosPhiV=cos(phi_v);
    float sinPhiU=sin(phi_u);
    float sinPhiV=sin(phi_v);
    return (u*sinPhiV-v*sinPhiU)/(cosPhiV*sinPhiU-cosPhiU*sinPhiV);
  }
 

  ///
  /// DFDCGeometry::getYLocalStrips()
  ///  Compute the y-coordinate in the layer coordinate system
  /// from the strip data
  ///
  static inline float getYLocalStrips(float u, float v){
    return  STRIP_SPACING*(u-v)/2./sin(CATHODE_ROT_ANGLE);
  } 
  static inline float getYLocalStrips(float u, float phi_u,float v, float phi_v){
    float cosPhiU=cos(phi_u);
    float cosPhiV=cos(phi_v);
    float sinPhiU=sin(phi_u);
    float sinPhiV=sin(phi_v);
    return (u*cosPhiV-v*cosPhiU)/(cosPhiU*sinPhiV-cosPhiV*sinPhiU);
  }
 


};

#endif // DFDCGEOMETRY_H

//***************************************************
// DFDCGeometry - temporary geometry class for FDC.
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
//***************************************************

#ifndef DFDCGEOMETRY_H
#define DFDCGEOMETRY_H

#define FDC_NUM_LAYERS			24
#define FDC_ACTIVE_RADIUS		53.6

//----- These were cut from HDGeant/hitFDC.c -----
#define DRIFT_SPEED           .0055
#define WIRE_DEAD_ZONE_RADIUS 3.5
#define ANODE_CATHODE_SPACING 0.5
#define TWO_HIT_RESOL         250.
#define WIRES_PER_PLANE       96
#define WIRE_SPACING          1.116
#define U_OF_WIRE_ZERO        (-((WIRES_PER_PLANE-1)*WIRE_SPACING)/2)
#define STRIPS_PER_PLANE      216
#define STRIP_SPACING         0.5
#define U_OF_STRIP_ZERO		   (-((STRIPS_PER_PLANE-1)*STRIP_SPACING)/2)
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

///
/// class DFDCGeometry: definition of a class providing basic geometry
/// methods for FDC reconstruction.
///
class DFDCGeometry : public JObject {
	public:
	
		DFDCGeometry(void);
				

		// Get z-coordinate for first wire plane in a package
		static inline float GetZpackage(int layer){
		  int package_number=(layer-1)/6;
		  float z=212.0+95.5;
		  switch (package_number){
		  case 0:
		    z+=-92.5-2.0;
		    break;
		  case 1:
		    z+= -32.5-2.0;
		    break;
		  case 2:
		    z+=+26.5-2.0;
		    break;
		  case 3:
		    z+=+86.5-2.0;
		    break;		    
		  }
		  return z;
		}


		static const DFDCWire* GetDFDCWire(int layer, int wire);
	
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
		/// DFDCGeometry::getLayerZ():
		/// Get the Z position of a layer
		///
		static float getLayerZ(const DFDCHit* h) {
			if (h->gPlane <= 18)  
				return 227.5 + h->gPlane*0.5;			// Thing 1
			if (h->gPlane <= 36)	   
				return 283.5 + (h->gPlane - 18)*0.5;	// Thing 2
			if (h->gPlane <= 54)	   
				return 358.5 + (h->gPlane - 36)*0.5;	// Thing 3
			return 394.5 + (h->gPlane - 540*0.5); 		// Thing 4
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
		  return STRIP_SPACING*(u-1 + v-1 - (STRIPS_PER_PLANE-1))/sqrt(2.);
		}

		///
		/// DFDCGeometry::getYLocalStrips()
		///  Compute the y-coordinate in the layer coordinate system
		/// from the strip data
		///
		static inline float getYLocalStrips(float u, float v){
		  return  STRIP_SPACING*(u-v)/sqrt(2.);
		} 
	
		///
		/// DFDCGeometry::getLayerRotation():
		/// Compute the rotation of a detection layer (0, 60, -60)
		///
		static inline float getLayerRotation(const int gLayer) {
			switch ((gLayer-1) % 3) {
				case 0:
					return 0.0;
				case 1:
					return +M_PI / 3; // +60 degrees
				case 2:
					return -M_PI / 3; // -60 degrees
				return 0.0;
			}
		}
};

#endif // DFDCGEOMETRY_H

//**************************************************
// DFDCGeometry - temporary geometry class for FDC.
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
//**************************************************

#ifndef DFDCGEOMETRY_H
#define DFDCGEOMETRY_H

#include "DObject.h"
#include "DFDCHit.h"

class DFDCGeometry : public DObject {
	public:
		inline int globalLayer(const DFDCHit* h) {
			return 3*(h->module - 1) + (h->layer - 1) + 1;
		}		
		
		inline int globalPlane(const DFDCHit* h) {
			return 9*(h->module-1) + 3*(h->layer-1) + (h->plane-1) + 1;
		}
				
		float getLayerZ(const DFDCHit* h) {
			if (h->globalPlane <= 18)  
				return 227.5 + h->globalPlane*0.5;			// Thing 1
			if (h->globalPlane <= 36)	   
				return 283.5 + (h->globalPlane - 18)*0.5;	// Thing 2
			if (h->globalPlane <= 54)	   
				return 358.5 + (h->globalPlane - 36)*0.5;	// Thing 3
			return 394.5 + (h->globalPlane - 540*0.5); 		// Thing 4
		}
				
		inline float getWireR(const DFDCHit* h) {
			return h->element - 60.0;
		}
		
		inline float getStripR(const DFDCHit* h) {
			return (h->element - 119.5)*0.5;
		}
};

#endif // DFDCGEOMETRY_H

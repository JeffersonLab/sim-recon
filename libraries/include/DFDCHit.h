///
/// DFDChit.h - base class for all FDC hit types. 
///
/// Author: Craig Bookwalter, David Lawrence
/// Date:	March 2006
///

#ifndef _DFDCHit_
#define _DFDCHit_

#include "DObject.h"
#include "DFactory.h"

class DFDCHit : public DObject{
	public:
		HDCLASSDEF(DFDCHit);		
		int layer;				// 1, 2, or 3
		int module;				// 1 through 8
		int element;			// wire or strip number
	    int plane;				// for cathodes only: u or v plane
	    int globalPlane;		// 1 through 72
	    int globalLayer;		// 1 through 24
	    float dE;				// charge deposited
	    float t;				// drift time
	    float r;				// perpendicular distance from center of chamber to wire/strip center
	    int type;				// cathode=1, anode=0
};

bool DFDCHit_dE_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->dE < b->dE;
}

bool DFDCHit_layer_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->globalLayer < b->globalLayer;
}

bool DFDCHit_plane_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->globalPlane < b->globalPlane;
}


#endif // _DFDCHit_


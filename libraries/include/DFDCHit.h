// $Id$
//
//    File: DFDCHit.h
// Created: Sun Apr  3 10:39:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFDCHit_
#define _DFDCHit_

#include "DFactory.h"

class DFDCHit{
	public:
		HDCLASSDEF(DFDCHit);
		
		int layer;
		int module;
		float tau;
		float z;
		float u;
		float dE;
		float t;
		int type; ///< 0=anode, 1=cathode
};

#endif // _DFDCHit_


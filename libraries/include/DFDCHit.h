// $Id$
//
//    File: DFDCHit.h
// Created: Thu Jun  9 10:25:22 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFDCHit_
#define _DFDCHit_

#include "DObject.h"
#include "DFactory.h"

class DFDCHit:public DObject{
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


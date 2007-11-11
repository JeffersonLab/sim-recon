// $Id$
//
//    File: DFDCIntersection.h
// Created: Tue Oct 30 11:24:53 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DFDCIntersection_
#define _DFDCIntersection_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <DVector3.h>
#include "FDC/DFDCHit.h"
#include "FDC/DFDCWire.h"

class DFDCIntersection:public JObject{
	public:
		HDCLASSDEF(DFDCIntersection);
		
		const DFDCHit *hit1;
		const DFDCHit *hit2;
		const DFDCWire *wire1;
		const DFDCWire *wire2;
		DVector3 pos;
};

#endif // _DFDCIntersection_


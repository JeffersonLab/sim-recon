// $Id$
//
//    File: FDC_cathode.h
// Created: Fri Jun 24 15:26:54 EDT 2011
// Creator: davidl (on Linux ifarm1101 2.6.18-128.7.1.el5 #1 SMP x86_64)
//

#ifndef _FDC_cathode_
#define _FDC_cathode_

#include <TObject.h>


class FDC_cathode:public TObject{
	public:
		
		// ctor
		FDC_cathode(void){}
		
		// dtor
		~FDC_cathode(void){}
		
		// Clear
		void Clear(void){
			L1a_fired = false;
			L1b_fired = false;

			gPlane = 0;
			element = 0;
		}
		
		bool L1a_fired;
		bool L1b_fired;
		
		unsigned int gPlane;
		unsigned int element;

	private:
		ClassDef(FDC_cathode,1);
};



#endif // _FDC_cathode_


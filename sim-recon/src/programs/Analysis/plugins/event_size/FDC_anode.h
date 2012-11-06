// $Id$
//
//    File: FDC_anode.h
// Created: Fri Jun 24 15:26:54 EDT 2011
// Creator: davidl (on Linux ifarm1101 2.6.18-128.7.1.el5 #1 SMP x86_64)
//

#ifndef _FDC_anode_
#define _FDC_anode_

#include <TObject.h>


class FDC_anode:public TObject{
	public:
		
		// ctor
		FDC_anode(void){}
		
		// dtor
		~FDC_anode(void){}
		
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
		ClassDef(FDC_anode,1);
};



#endif // _FDC_anode_


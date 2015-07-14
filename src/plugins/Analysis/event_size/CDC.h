// $Id$
//
//    File: CDC.h
// Created: Fri Jun 24 15:26:54 EDT 2011
// Creator: davidl (on Linux ifarm1101 2.6.18-128.7.1.el5 #1 SMP x86_64)
//

#ifndef _CDC_
#define _CDC_

#include <TObject.h>


class CDC:public TObject{
	public:
		
		// ctor
		CDC(void){}
		
		// dtor
		~CDC(void){}
		
		// Clear
		void Clear(void){
			L1a_fired = false;
			L1b_fired = false;

			ring = 0;
			straw = 0;
		}
		
		bool L1a_fired;
		bool L1b_fired;
		
		unsigned int ring;
		unsigned int straw;

	private:
		ClassDef(CDC,1);
};



#endif // _CDC_


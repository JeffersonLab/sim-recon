// $Id$
//
//    File: TOF.h
// Created: Fri Jun 24 15:26:54 EDT 2011
// Creator: davidl (on Linux ifarm1101 2.6.18-128.7.1.el5 #1 SMP x86_64)
//

#ifndef _TOF_
#define _TOF_

#include <TObject.h>


class TOF:public TObject{
	public:
		
		// ctor
		TOF(void){}
		
		// dtor
		~TOF(void){}
		
		// Clear
		void Clear(void){
			L1a_fired = false;
			L1b_fired = false;

			plane = 0;
			bar = 0;
			lr = 0;
		}
		
		bool L1a_fired;
		bool L1b_fired;
		
		unsigned int plane;
		unsigned int bar;
		unsigned int lr;

	private:
		ClassDef(TOF,1);
};



#endif // _TOF_


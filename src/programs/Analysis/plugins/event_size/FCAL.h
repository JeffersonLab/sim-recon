// $Id$
//
//    File: FCAL.h
// Created: Fri Jun 24 15:26:54 EDT 2011
// Creator: davidl (on Linux ifarm1101 2.6.18-128.7.1.el5 #1 SMP x86_64)
//

#ifndef _FCAL_
#define _FCAL_

#include <TObject.h>


class FCAL:public TObject{
	public:
		
		// ctor
		FCAL(void){}
		
		// dtor
		~FCAL(void){}
		
		// Clear
		void Clear(void){
			L1a_fired = false;
			L1b_fired = false;

			row = 0;
			column = 0;
		}
		
		bool L1a_fired;
		bool L1b_fired;
		
		unsigned int row;
		unsigned int column;

	private:
		ClassDef(FCAL,1);
};



#endif // _FCAL_


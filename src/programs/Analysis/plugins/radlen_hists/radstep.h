// $Id$
//
//    File: radstep.h
// Created: Fri Apr 17 13:27:38 EDT 2009
// Creator: davidl
//

#ifndef _radstep_
#define _radstep_

#include <TObject.h>
#include <TVector3.h>


class radstep:public TObject{

	public:
		TVector3 pthrown;
		TVector3 pos;
		double s;
		double stot;
		double Xo;
		double ix_over_Xo;

	private:
		ClassDef(radstep,1);

};

#endif // _radstep_


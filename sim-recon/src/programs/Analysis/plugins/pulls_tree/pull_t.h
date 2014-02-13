// $Id$
//
//    File: pull_t.h
// Created: Fri. Feb. 19, 2010
// Creator: davidl
//

#ifndef _pull_t_
#define _pull_t_

#include <TObject.h>
#include <TVector3.h>


class pull_t:public TObject{

	public:
		int eventnumber;
		float resi;
		float err;
		float s;
		float pull;
		float trk_chisq;
		float trk_Ndof;
		TVector3 pthrown;

	private:
		ClassDef(pull_t,1);

};

#endif // _pull_t_


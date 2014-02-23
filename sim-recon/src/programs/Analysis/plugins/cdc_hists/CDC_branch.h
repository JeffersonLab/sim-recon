// $Id$
//
//    File: FDC_branch.h.h
// Created: Mon Dec 17 2007
// Creator: davidl
//

#ifndef _CDC_branch_
#define _CDC_branch_

#include <TObject.h>
#include <TVector3.h>


class CDC_branch:public TObject{

	public:

		TVector3 pos_truth;

	private:
		ClassDef(CDC_branch,1);

};

#endif // _CDC_branch_


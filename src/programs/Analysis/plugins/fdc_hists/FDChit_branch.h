// $Id$
//
//    File: FDChit_branch.h.h
// Created: Mon Dec 17 2007
// Creator: davidl
//

#ifndef _FDChit_branch_
#define _FDChit_branch_

#include <TObject.h>
#include <TVector3.h>


class FDChit_branch:public TObject{

	public:

		int layer;
		int module;
		int element;
		int plane;
		int gPlane;
		int gLayer;
		float q;
		float t;
		float r;
		int type;

	private:
		ClassDef(FDChit_branch,1);

};

#endif // _FDChit_branch_


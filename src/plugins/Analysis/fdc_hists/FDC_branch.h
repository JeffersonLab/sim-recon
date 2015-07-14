// $Id$
//
//    File: FDC_branch.h.h
// Created: Mon Dec 17 2007
// Creator: davidl
//

#ifndef _FDC_branch_
#define _FDC_branch_

#include <TObject.h>
#include <TVector3.h>


class FDC_branch:public TObject{

	public:

  float dX;
  float dY;
  float dPhi;
  int layer;
  int N;

	private:
		ClassDef(FDC_branch,1);

};

#endif // _FDC_branch_


// $Id$
//
//    File: FDC_c_branch.h.h
// Created: Tue March 25 2014
// Creator: staylor
//

#ifndef _FDC_c_branch_
#define _FDC_c_branch_

#include <TObject.h>
#include <TVector3.h>


class FDC_c_branch:public TObject{

	public:

  float dU;
  float dPhiU;
  float dV;
  float dPhiV;
  int layer;
  int N;

	private:
		ClassDef(FDC_c_branch,1);

};

#endif // _FDC_c_branch_


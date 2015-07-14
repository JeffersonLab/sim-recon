// $Id$
//
//    File: dchit.h
// Created: Tues. Feb. 5, 2008
// Creator: davidl
//

#ifndef _dchit_
#define _dchit_

#include <TObject.h>
#include <TVector3.h>


class dchit:public TObject{

	public:
		int eventnumber;
		int wire;
		int layer;
		float t;
		float tof;
		float doca;
		float resi;
		float u;
		float u_pseudo;
		float u_lorentz;
		float resic;
		float trk_chisq;
		float trk_Ndof;
		int LRis_correct;
		int LRfit;
		TVector3 pos_doca;
		TVector3 pos_wire;

	private:
		ClassDef(dchit,1);

};

#endif // _dchit_


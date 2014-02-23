// $Id$
//
//    File: DEventProcessor_trackeff_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _track_info_
#define _track_info_

#include <iostream>

#include <TObject.h>
#include <TVector3.h>

class track_info:public TObject{
	public:
		
		track_info(void):trk_chisq(1.0E6),trk_Ndof(-1),Ncdc(0),Nfdc(0){}
		
		TVector3 p;
		float trk_chisq;
		int trk_Ndof;
		int Ncdc;
		int Nfdc;
		
		track_info& operator=(const track_info &ti){
			this->p = ti.p;
			this->trk_chisq = ti.trk_chisq;
			this->trk_Ndof = ti.trk_Ndof;
			this->Ncdc = ti.Ncdc;
			this->Nfdc = ti.Nfdc;
			return *this;
		}
		
	private:
		ClassDef(track_info,1);
};



#endif // _track_info_


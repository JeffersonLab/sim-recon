// $Id$
//
//    File: DEventProcessor_trackeff_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _trackpar_
#define _trackpar_

#include <TObject.h>
#include <TVector3.h>


class trackpar:public TObject{

	public:

		int eventnumber;
		TVector3 pthrown;
		TVector3 pfit;
		float z_thrown;
		float z_fit;
		float z_can;
		float r_fit;
		int NLRcorrect;
		int NLRincorrect;

	private:
		ClassDef(trackpar,1);

};

#endif // _trackpar_


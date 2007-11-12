// $Id$
//
//    File: DEventProcessor_trackeff_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _TrkEff_Leaf_
#define _TrkEff_Leaf_

#include <TObject.h>
#include <TVector3.h>


class TrkEff_Leaf:public TObject{

	public:

		TVector3 pthrown;
		TVector3 pcan;
		float z_thrown;
		float z_can;
		int ncdc_hits_thrown;
		int ncdc_hits_can;
		int ncdc_hits_thrown_and_can;
		int nfdc_hits_thrown;
		int nfdc_hits_can;
		int nfdc_hits_thrown_and_can;
		int ncdc_hits;
		int nfdc_hits;
		float cdc_chisq_can;
		float fdc_chisq_can;
		int status;

	private:
		ClassDef(TrkEff_Leaf,1);

};

#endif // _TrkEff_Leaf_


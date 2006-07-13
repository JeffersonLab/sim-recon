// $Id$
//
//    File: DHDDMBCALHit_factory.cc
// Created: Thu Jun  9 10:25:22 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//
// adapted by Chuncheng Xu from FDC facotries On Dec.22,2005
// 

#include "DBCALHit.h"
#include "DHDDMBCALHit_factory.h"

//------------------
// evnt
//------------------
jerror_t DHDDMBCALHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
  
  vector<const DBCALHit*> bcalhits;
  eventLoop->Get(bcalhits);
   int end=0; 
  for (unsigned int i = 0; i < bcalhits.size(); i++) {
     const DBCALHit *bcalhit = bcalhits[i];
     DHDDMBCALHit  *hddmhit    = new  DHDDMBCALHit;
     hddmhit->module = bcalhit->module;
     hddmhit->layer  = bcalhit->layer;
     hddmhit->sector = bcalhit->sector;
     if(bcalhit->end==DBCALHit::UPSTREAM)  end=0;
     if(bcalhit->end==DBCALHit::DOWNSTREAM)end=1;
     hddmhit->end    = end;
     hddmhit->E      = bcalhit->E;
     hddmhit->t      = bcalhit->t;     
     _data.push_back(hddmhit);
  }



	return NOERROR;
}




//------------------
// toString
//------------------
const string DHDDMBCALHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("mod:   lay:   sec:   end:   E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
          	DHDDMBCALHit *s = _data[i];
	        printf("%d  %d  %d  %d  %f  %f\n",s->module,s->layer,s->sector,s->end,s->E,s->t);

	  //	  		DHDDMBCALHit *hddmhit = _data[i];			  //	printnewrow();
	  //	printcol("%d",	i);
	  //	printcol("%d",	hddmhit->module);
	  //	printcol("%d",	hddmhit->layer);
	  //	printcol("%d",	hddmhit->sector);
	  //	printcol(hddmhit->end==0 ? "upstream":"downstream");
	  //	printcol("%2.3f",	hddmhit->E);
	  //	printcol("%1.3f",	hddmhit->t);
	  //	printrow();
	}

	//  printnewrow();
	//  printrow();

	return _table;

}

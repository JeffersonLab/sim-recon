// $Id$
//
//    File: DFactory_DBCALMCResponse.cc
// Created: Thu Nov 17 09:56:05 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#include <cassert>	
#include "DFactory_DBCALMCResponse.h"
#include "DFactory_DHDDMBCALHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DBCALMCResponse::evnt(DEventLoop *loop, int eventnumber)
{

  vector<const DHDDMBCALHit*> hddmhits;
  eventLoop->Get(hddmhits);
    
  for (unsigned int i = 0; i < hddmhits.size(); i++) {
    const DHDDMBCALHit *hddmhit = hddmhits[i];
    DBCALMCResponse *response = new DBCALMCResponse;
   
     response->module =hddmhit->module;
     response->layer = hddmhit->layer;
     response->sector = hddmhit->sector;
     response->end = hddmhit->end;
     response->E = hddmhit->E;
     response->t = hddmhit->t;
     
    _data.push_back(response);
  }

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DBCALMCResponse::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	//		printheader("row:    x:     y:");
	//
	// 	for(int i=0; i<_data.size(); i++){
	//			DBCALMCResponse *myDBCALMCResponse = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDBCALMCResponse->x);
	//			printcol("%3.2f",	myDBCALMCResponse->y);
	//			printrow();
	//		}
	//
	printheader("row:   module:  layer:  sector:         end:     E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DBCALMCResponse *s = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%d",	s->module);
		printcol("%d",	s->layer);
		printcol("%d",	s->sector);
		printcol(s->end==0? "upstream":"downstream");
		printcol("%2.3f",	s->E);
		printcol("%2.3f",	s->t);
		printrow();
	}


	return _table;

}

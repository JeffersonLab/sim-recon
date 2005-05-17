// $Id$
//
//    File: DFactory_DFCALShower.cc
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include "DFactory_DFCALShower.h"
#include "DFactory_DFCALHit.h"
#include "DEvent.h"

//------------------
// evnt
//    Setting up the structure of the FCAL reconstruction code.
//    Nothing is done at this point except a copy of the DFCALHits
//      object into the DFCALShower object.
//------------------
derror_t DFactory_DFCALShower::evnt(int eventnumber)
{

        vector<const DFCALHit*> fcalhits;
        event->Get(fcalhits);

        for (unsigned int i = 0; i < fcalhits.size(); i++){

          const DFCALHit *fcalhit = fcalhits[i];
          DFCALShower *fcalshower = new DFCALShower;

          fcalshower->x = fcalhit->x;
          fcalshower->y = fcalhit->y;
          fcalshower->E = fcalhit->E;
          fcalshower->t = fcalhit->t;

          _data.push_back(fcalshower);

        }

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DFCALShower::toString(void)
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
	//			DFCALShower *myDFCALShower = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDFCALShower->x);
	//			printcol("%3.2f",	myDFCALShower->y);
	//			printrow();
	//		}
	//
	return _table;

}

// $Id$
//
//    File: DBeamPhoton_factory.h
// Created: Thu Feb 14 10:11:52 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.10.1 i386)
//

#ifndef _DBeamPhoton_factory_
#define _DBeamPhoton_factory_

#include <string>
#include <cmath>
using std::string;

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <PID/DBeamPhoton.h>

class DBeamPhoton_factory: public JFactory<DBeamPhoton>{
	public:
		DBeamPhoton_factory();
		virtual ~DBeamPhoton_factory();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DBeamPhoton";}
		
		double t; ///< Time at which photon arrives at interaction vertex location
		
		inline const string toString(void);
};

//------------------
// toString
//------------------
inline const string DBeamPhoton_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("  row:    E(GeV):  theta(degrees):   phi(degrees):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DBeamPhoton * photon = _data[i];
		DVector3 mom = photon->momentum();

		printnewrow();
		
		printcol("%d", i);
		printcol("%3.3f", mom.Mag());
		printcol("%3.3f", mom.Theta()*180.0/M_PI);
		printcol("%3.1f", mom.Phi()*180.0/M_PI);
		printcol("%3.1f", photon->t);

		printrow();
		
		//mom.Print();
	}
	
	return _table;
}


#endif // _DBeamPhoton_factory_


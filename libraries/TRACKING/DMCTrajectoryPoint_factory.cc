// $Id$
//
//    File: DMCTrajectoryPoint_factory.cc
// Created: Mon Jun 12 09:29:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.6.0 powerpc)
//

#include <cassert>	

#include "HDDM/hddm_s.h"
#include "DMCTrajectoryPoint_factory.h"

//------------------
// evnt
//------------------
jerror_t DMCTrajectoryPoint_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// This doesn't do anything. All of the work is done in  Extract_HDDM()

	return NOERROR;
}

//------------------
// toString
//------------------
const string DMCTrajectoryPoint_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The JFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	printheader("pri-trk: trk: part:     x:     y:     z:   t(ns):     px:     py:     pz: E(GeV):     step: radlen: dE(MeV): mech:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DMCTrajectoryPoint *p = _data[i];

		printnewrow();
		printcol("%d"   ,	p->primary_track);
		printcol("%d"   ,	p->track);
		printcol("%d"   ,	p->part);
		printcol("%3.1f",	p->x);
		printcol("%3.1f",	p->y);
		printcol("%3.1f",	p->z);
		printcol("%3.2f",	p->t*1.0E9);
		printcol("%3.3f",	p->px);
		printcol("%3.3f",	p->py);
		printcol("%3.3f",	p->pz);
		printcol("%3.3f",	p->E);
		printcol("%3.3g",	p->step);
		printcol("%3.1f",	p->radlen);
		printcol("%3.3f",	p->dE*1.0E3);
		printcol("%d"   ,	p->mech);
		printrow();
	}

	return _table;

}

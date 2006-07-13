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
// Extract_HDDM
//------------------
jerror_t DMCTrajectoryPoint_factory::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects.
	
	v.clear();

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->mcTrajectory == HDDM_NULL ||
			hits->mcTrajectory->mcTrajectoryPoints == HDDM_NULL)continue;

		s_McTrajectoryPoints_t *points = hits->mcTrajectory->mcTrajectoryPoints;
		for(unsigned int i=0; i<points->mult; i++){
			DMCTrajectoryPoint *p = new DMCTrajectoryPoint;
			
			p->x = points->in[i].x;
			p->y = points->in[i].y;
			p->z = points->in[i].z;
			p->t = points->in[i].t;
			p->px = points->in[i].px;
			p->py = points->in[i].py;
			p->pz = points->in[i].pz;
			p->E = points->in[i].E;

			p->dE = points->in[i].dE;
			p->track = points->in[i].track;
			p->part = points->in[i].part;
			
			
			v.push_back(p);
		}		
	}

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
	printheader("track: part:     x:     y:     z:   t(ns):     px:     py:     pz: E(GeV):    dE(MeV):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DMCTrajectoryPoint *p = _data[i];

		printnewrow();
		printcol("%d",	p->track);
		printcol("%d",	p->part);
		printcol("%3.1f",	p->x);
		printcol("%3.1f",	p->y);
		printcol("%3.1f",	p->z);
		printcol("%3.2f",	p->t*1.0E9);
		printcol("%3.3f",	p->px);
		printcol("%3.3f",	p->py);
		printcol("%3.3f",	p->pz);
		printcol("%3.3f",	p->E);
		printcol("%3.3f",	p->dE*1.0E3);
		printrow();
	}

	return _table;

}

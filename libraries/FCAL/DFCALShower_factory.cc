// $Id$
//
//    File: DFCALShower_factory.cc
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <math.h>

#include "DFCALShower_factory.h"
#include "DFCALHit.h"
#include "JANA/JEvent.h"

// Used to sort hits by Energy
bool FCALHitSort_C(const DFCALHit* const &thit1, const DFCALHit* const &thit2) {
	return thit1->E > thit2->E;
}

//----------------
// Constructor
//----------------
DFCALShower_factory::DFCALShower_factory()
{
	// Set defaults
	MAX_SHOWER_DIST = 9.0; // cm
	
	gPARMS->SetDefaultParameter("FCAL:MAX_SHOWER_DIST",		MAX_SHOWER_DIST);

}


//------------------
// evnt
//    Trivial calorimeter reconstruction. (D. Lawrence)
//------------------
jerror_t DFCALShower_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	vector<const DFCALHit*> fcalhits;
	eventLoop->Get(fcalhits);
	
	// Sort hits by energy
	sort(fcalhits.begin(), fcalhits.end(), FCALHitSort_C);
	
	// Make boolean array to keep track of used hits
	vector<bool> used;
	for(unsigned int i=0; i<fcalhits.size(); i++)used.push_back(false);

	// Look for showers. Take first, unused hit and assume he is
	// center of shower. Consider all hits within 9 cm(about 2.5
	// Moliere radii) part of the same shower.
	do{
		// Unused hit with largest energy is center of shower
		unsigned int center_index = 0;
		for (center_index = 0; center_index < fcalhits.size(); center_index++){
			if(!used[center_index])break;
		}
		if(center_index>=fcalhits.size())break;

		// All subsequent hits within 9cm are part of this shower
		used[center_index] = true;
		double E = fcalhits[center_index]->E;
		double logE = log(E);
		double logEsum = logE;
		double X = fcalhits[center_index]->x*logE;
		double Y = fcalhits[center_index]->y*logE;
		for(unsigned int i=center_index+1; i<fcalhits.size(); i++){
			if(used[i])continue;
			double deltaX = fcalhits[center_index]->x - fcalhits[i]->x;
			double deltaY = fcalhits[center_index]->y - fcalhits[i]->y;
			double r2 = deltaX*deltaX + deltaY*deltaY;
			if(r2 > MAX_SHOWER_DIST*MAX_SHOWER_DIST)continue;
			
			used[i] = true;
			E += fcalhits[i]->E;
			logE = log(E);
			logEsum += logE;
			X += fcalhits[i]->x*logE;
			Y += fcalhits[i]->y*logE;
		}
		
		X/=logEsum;
		Y/=logEsum;

		DFCALShower *fcalshower = new DFCALShower;

		fcalshower->x = X;
		fcalshower->y = Y;
		fcalshower->E = E*1.90; // this number is empirical and should be a calibration constant
		fcalshower->t = fcalhits[center_index]->t;

		_data.push_back(fcalshower);

	}while(true);

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFCALShower_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   x(cm):   y(cm):   E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALShower *fcalhit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%3.1f", fcalhit->x);
		printcol("%3.1f", fcalhit->y);
		printcol("%2.3f", fcalhit->E);
		printcol("%4.0f", fcalhit->t);
		printrow();
	}

	return _table;
}

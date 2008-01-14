// $Id: DFCALPhoton_factory.cc 2280 2006-12-07 17:10:06Z davidl $
//
//    File: DFCALPhoton_factory.cc
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <math.h>
#include <DVector3.h>
#include <TLorentzVector.h>

#include "DFCALPhoton_factory.h"
//#include "DFCALPhoton.h"
#include "DFCALCluster.h"
#include "DFCALHit.h"
#include "JANA/JEvent.h"


//----------------
// Constructor
//----------------
DFCALPhoton_factory::DFCALPhoton_factory()
{
	// Set defaults
	
//	gPARMS->SetDefaultParameter("FCAL:EGAMMA_NORM", EGAMMA_EPSILON);
//	gPARMS->SetDefaultParameter("FCAL:EGAMMA_EPSILON", EGAMMA_EPSILON);

}


//------------------
// evnt
//    Trivial calorimeter reconstruction. (D. Lawrence)
//------------------
jerror_t DFCALPhoton_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	vector<const DFCALCluster*> fcalClusters;
	eventLoop->Get(fcalClusters);
	
       for ( unsigned int i = 0; i < fcalClusters.size(); i++ ) {

                DFCALPhoton* fcalPhoton = makePhoton( fcalClusters[i] );

		_data.push_back(fcalPhoton);

       } 

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFCALPhoton_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   E(GeV):   Px(GeV):   Py(GeV):   Pz(GeV):   X(cm):   Y(cm):   Z(cm):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALPhoton *fcalphot = _data[i];
               
		printnewrow();
		printcol("%d",	i);
		printcol("%6.2f", fcalphot->getEnergy());
		printcol("%6.2f", fcalphot->getMom3().X());
		printcol("%6.2f", fcalphot->getMom3().Y());
		printcol("%6.2f", fcalphot->getMom3().Z());
		printcol("%7.2f", fcalphot->getPosition().X());
		printcol("%7.2f", fcalphot->getPosition().Y());
		printcol("%7.2f", fcalphot->getPosition().Z());
		printrow();
	}

	return _table;
}

// Non-linear and depth corrections should be fixed within DFCALPhoton member functions
DFCALPhoton* DFCALPhoton_factory::makePhoton(const DFCALCluster* cluster) 
{

        DFCALPhoton* photon = new DFCALPhoton;
       
// Do non-linar energy correction first,
	photon->fixEnergy( cluster->getEnergy() );

// than depth. Idealy, the two previous steps should be done simultaneously.
	photon->fixDepth( photon->getEnergy(), cluster->getCentroid() );   

        photon->setErrorXY(cluster->getRMS_x(), cluster->getRMS_y());

// Than set momentum to units of GeV.
	photon->setMom3( photon->getEnergy(), photon->getPosition() );   
        photon->setMom4();

        return photon;
}


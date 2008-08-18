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
#include <JANA/JEvent.h>
using namespace jana;

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


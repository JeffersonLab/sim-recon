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
	
    
       for (vector<const DFCALCluster*>::const_iterator cluster  = fcalClusters.begin(); 
                                                        cluster != fcalClusters.end(); 
							cluster++) {

//		DFCALPhoton *fcalPhoton = new DFCALPhoton;

	        // Apply simple non-linear and depth corrections to clusters
               
/*                const double Ein = (**cluster).getEnergy();
                const DVector3 mom( (**cluster).getCentroid().x(), 
                                    (**cluster).getCentroid().y(),
                                    (**cluster).getCentroid().z());
                    
                fcalPhoton->setEnergy( (**cluster).getEnergy() );
                fcalPhoton->set3Mom( fcalPhoton->getEnergy(), (**cluster).getCentroid() );
                */

                TLorentzVector gamma( (**cluster).getCentroid(), (**cluster).getEnergy() );
                DFCALPhoton* fcalPhoton = makePhoton( gamma );

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

	printheader("row:   E(GeV):   Px(GeV):	  Py(GeV):    Pz(GeV):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALPhoton *fcalphot = _data[i];
               
		printnewrow();
		printcol("%d",	i);
		printcol("%3.1f", fcalphot->getEnergy());
		printcol("%3.1f", fcalphot->getMom4().X());
		printcol("%3.1f", fcalphot->getMom4().Y());
		printcol("%3.1f", fcalphot->getMom4().Z());
		printrow();
	}

	return _table;
}

// Non-linear and depth corrections should be fixed within DFCALPhoton member functions
DFCALPhoton* DFCALPhoton_factory::makePhoton(const TLorentzVector gamma) 
{

        DFCALPhoton* photon = new DFCALPhoton;
        
	photon->setEnergy( gamma.T() );   
	photon->setMom3(photon->getEnergy(), gamma.Vect());   
        photon->setMom4();

        return photon;
}


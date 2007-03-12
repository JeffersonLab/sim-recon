// $Id: DFCALPhoton_factory.cc 2280 2006-12-07 17:10:06Z davidl $
//
//    File: DFCALPhoton_factory.cc
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <math.h>
#include <TVector3.h>

#include "DFCALPhoton_factory.h"
#include "DFCALPhoton.h"
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


		DFCALPhoton *fcalPhoton = new DFCALPhoton;

	        // Apply simple non-linear and depth corrections to clusters
               
                const double Ein = (**cluster).getEnergy();
                const TVector3 pos( (**cluster).getCentroid().x, 
                                    (**cluster).getCentroid().y,
                                    (**cluster).getCentroid().z);

                fcalPhoton->setEnergy(Ein);
                fcalPhoton->set3Mom(Ein,pos);

//                DFCALPhoton::makePhoton(); // Apply corrections 

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
		printcol("%3.1f", fcalphot->get3Mom().x());
		printcol("%3.1f", fcalphot->get3Mom().y());
		printcol("%3.1f", fcalphot->get3Mom().z());
		printrow();
	}

	return _table;
}


DFCALPhoton::DFCALPhoton()
{
   fEnergy = 0;
   fMom.SetY(0) ;
   fMom.SetX(0) ;
   fMom.SetZ(0) ;
}

DFCALPhoton::~DFCALPhoton()
{
}

#define	FCAL_RADIATION_LENGTH  3.1
#define	FCAL_CRITICAL_ENERGY  0.01455
#define	FCAL_SHOWER_OFFSET   3.0 

// These change with the FCAL attenuation length (L=166cm)
#define	EGAMMA_NORM 0.6348
#define	EGAMMA_EPSILON   0.029
 
// Simple non-linear correction
void DFCALPhoton::setEnergy(const double energy) 
{

	double const A=1/EGAMMA_NORM; 				// Normalization factor 
        double const power = 1/(1+EGAMMA_EPSILON); 	        // Non-linear factor
 	fEnergy = pow(A*energy,power);
}

// Simple depth correction: parameters imported from Radphi. 
void DFCALPhoton::set3Mom(const double energy, const TVector3 pos) 
{

	double x = pos.X();
	double y = pos.Y();

//  estimate shower depth based on shower maximum
        double zMax = (FCAL_RADIATION_LENGTH*( 
			FCAL_SHOWER_OFFSET + log(energy/FCAL_CRITICAL_ENERGY)));

	double z = FCAL_Zmin - Shower_Vertex_Z + zMax;

// normalization factor to momenta [GeV]
        double f = energy/sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
        fMom.SetX(x*f) ;
    	fMom.SetY(y*f) ;
	fMom.SetZ(z*f) ;
}

// Non-linear and depth corrections should be fixed together
void DFCALPhoton::makePhoton() {

// Apply final non-linear and depth corrections;

}

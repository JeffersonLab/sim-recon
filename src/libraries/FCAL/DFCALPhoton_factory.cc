// $Id: DFCALPhoton_factory.cc 2280 2006-12-07 17:10:06Z davidl $
//
//    File: DFCALPhoton_factory.cc
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <math.h>
#include <DVector3.h>
#include <TLorentzVector.h>
using namespace std;

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

// Set of coefficients for non-linear energy corrections 
        NON_LIN_COEF_A = 0.5334;
        NON_LIN_COEF_B = 2.9598; 
        NON_LIN_COEF_C = 2.6635;
        NON_LIN_COEF_alfa = 1.0073;

// Parameters to make shower-depth correction taken from Radphi, 
// slightly modifed to match photon-polar angle
        FCAL_RADIATION_LENGTH = 3.1;
        FCAL_CRITICAL_ENERGY = 0.035;
        FCAL_SHOWER_OFFSET = 1.0;

	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_A", NON_LIN_COEF_A);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_B", NON_LIN_COEF_B);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_C", NON_LIN_COEF_C);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_alfa", NON_LIN_COEF_alfa);

	gPARMS->SetDefaultParameter("FCAL:FCAL_RADIATION_LENGTH", FCAL_RADIATION_LENGTH);
	gPARMS->SetDefaultParameter("FCAL:FCAL_CRITICAL_ENERGY", FCAL_CRITICAL_ENERGY);
	gPARMS->SetDefaultParameter("FCAL:FCAL_SHOWER_OFFSET", FCAL_SHOWER_OFFSET);
	
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
       
// Non-linar energy correction are done here
       int MAXITER = 1000;

       double A  = NON_LIN_COEF_A;
       double B  = NON_LIN_COEF_B;
       double C  = NON_LIN_COEF_C;
       double alfa  = NON_LIN_COEF_alfa;
       double Eclust = cluster->getEnergy();
       double Egamma = Eclust/A;

       for (int niter=0; 1; niter++) {

           double energy = Egamma;
           double non_lin_part = pow(Egamma,1+alfa)/(B+C*Egamma);
           Egamma = Eclust/A - non_lin_part;
           if ( fabs( (Egamma-energy)/energy ) < 0.001 ) {

               break;

           }
           else if ( niter > MAXITER ) {

               Egamma  = 0;
               cout << " Iteration failed for cluster energy " << Eclust << endl;
               break;

           }
          
       }

	photon->setEnergy( Egamma );

// than depth corrections 

        DVector3  pos = cluster->getCentroid();
        double z0 = FCAL_Zmin - Shower_Vertex_Z;
        double zMax = (FCAL_RADIATION_LENGTH*(
                       FCAL_SHOWER_OFFSET + log(Egamma/FCAL_CRITICAL_ENERGY)));
        double zed = z0;
        double zed1 = z0 + zMax;

        if ( Egamma > 0 ) { 

            double r0 = sqrt( pos.X()*pos.X() + pos.Y()*pos.Y() );

            int niter;
            for ( niter=0; niter<100; niter++) {

                double tt = r0/zed1;
                zed = z0 + zMax/sqrt( 1 + tt*tt );
                if ( fabs( (zed-zed1) ) < 0.001) {
                   break;
                }
                zed1 = zed;
            }

        }
    
        pos.SetZ( zed );

	photon->setPosition( pos );   

        photon->setPosError(cluster->getRMS_x(), cluster->getRMS_y(), fabs( zed-zed1 ) );

// Than set momentum to GeV units.
	photon->setMom3( Egamma, photon->getPosition() );   
        photon->setMom4();

        return photon;

}




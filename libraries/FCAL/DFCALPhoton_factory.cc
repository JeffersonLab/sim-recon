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
  
  //Regular Lead Glass Fit values for 30cm RHG radius
  NON_LIN_COEF_A1 = 0.53109;
  NON_LIN_COEF_B1 = 2.66426; 
  NON_LIN_COEF_C1 = 2.70763;
  NON_LIN_COEF_alfa1 = 1+0.0191858;

  //Radiation Hard Lead Glass for 30cm radius
  NON_LIN_COEF_A2 = 0.463044;
  NON_LIN_COEF_B2 = 2.4628; 
  NON_LIN_COEF_C2 = 2.39377;
  NON_LIN_COEF_alfa2 = 1+0.03614;


  
  BUFFER_RADIUS = 8.0;   //transition region buffer
  RHG_RADIUS = 30.0; //RHG radius

// Parameters to make shower-depth correction taken from Radphi, 
// slightly modifed to match photon-polar angle
        FCAL_RADIATION_LENGTH = 3.1;
        FCAL_CRITICAL_ENERGY = 0.035;
        FCAL_SHOWER_OFFSET = 1.0;

	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_A1", NON_LIN_COEF_A1);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_B1", NON_LIN_COEF_B1);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_C1", NON_LIN_COEF_C1);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_alfa1", NON_LIN_COEF_alfa1);

	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_A2", NON_LIN_COEF_A2);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_B2", NON_LIN_COEF_B2);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_C2", NON_LIN_COEF_C2);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_alfa2", NON_LIN_COEF_alfa2);

	gPARMS->SetDefaultParameter("FCAL:BUFFER_RADIUS",BUFFER_RADIUS);
	gPARMS->SetDefaultParameter("FCAL:RHG_RADIUS",RHG_RADIUS);
	
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
       DVector3  pos = cluster->getCentroid();
       float x0 = pos.Px();
       float y0 = pos.Py();
       float hrad = sqrt(x0*x0+y0*y0);
       float x;
       float y;
   
       
       double A;
       double B;
       double C;
       double alfa;
       
      

       int blocks = cluster->getHits();
       DVector3 Ipos[blocks];
       double efrac[blocks];
       float Ein=0;
       float Eout=0;
       float Etot=0;
       float erad;
       
       //transition region
       if(fabs(RHG_RADIUS-hrad) < BUFFER_RADIUS ){
	 cluster->getHitsEf(efrac,blocks);
	 cluster->getHitsPos(Ipos,blocks);
	 for (int h=0;h<blocks;h++){
	   x=Ipos[h].Px();
	   y=Ipos[h].Py();
	   erad = sqrt(x*x+y*y);
	   if(erad<RHG_RADIUS){
	     
	     Ein=Ein+efrac[h];
   
	   }
	   else{
	     
	     Eout = Eout+efrac[h];

	   }
	 }
	 
	 Etot=Eout+Ein;
	 A  = Eout/Etot*NON_LIN_COEF_A1+Ein/Etot*NON_LIN_COEF_A2;
	 B  = Eout/Etot*NON_LIN_COEF_B1+Ein/Etot*NON_LIN_COEF_B2;
	 C  = Eout/Etot*NON_LIN_COEF_C1+Ein/Etot*NON_LIN_COEF_C2;
	 alfa  = Eout/Etot*NON_LIN_COEF_alfa1+Ein/Etot*NON_LIN_COEF_alfa2;
	 
       }
       
       //Inner region
       else if(hrad<RHG_RADIUS){
	 
	 A  = NON_LIN_COEF_A2;
	 B  = NON_LIN_COEF_B2;
	 C  = NON_LIN_COEF_C2;
	 alfa  = NON_LIN_COEF_alfa2;
	 
       }
       
       //Outer Region
       else{
	 
	 A  = NON_LIN_COEF_A1;
	 B  = NON_LIN_COEF_B1;
	 C  = NON_LIN_COEF_C1;
	 alfa  = NON_LIN_COEF_alfa1;
	 
       }

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




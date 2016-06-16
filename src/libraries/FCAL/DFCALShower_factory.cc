//
//    File: DFCALShower_factory.cc
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)

#include <math.h>
#include <DVector3.h>
using namespace std;

#include "FCAL/DFCALShower_factory.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALCluster.h"
#include "FCAL/DFCALHit.h"
#include <JANA/JEvent.h>
#include <JANA/JApplication.h>
using namespace jana;

//----------------
// Constructor
//----------------
DFCALShower_factory::DFCALShower_factory()
{
  // should we use CCDB constants?
  LOAD_CCDB_CONSTANTS = 1.;
  gPARMS->SetDefaultParameter("FCAL:LOAD_NONLIN_CCDB", LOAD_CCDB_CONSTANTS);

  SHOWER_ENERGY_THRESHOLD = 50*k_MeV;
  gPARMS->SetDefaultParameter("FCAL:SHOWER_ENERGY_THRESHOLD", SHOWER_ENERGY_THRESHOLD);

  // these need to come from database to ensure accuracy
  // remove default value which might be close to the right solution,
  // but not quite correct -- allow command line tuning

  NON_LIN_COEF_A = 0; 
  NON_LIN_COEF_B = 0;
  NON_LIN_COEF_C = 0;
  NON_LIN_COEF_D = 0;
  NON_LIN_COEF_E = 0;

  linfit_slope = 0;
  linfit_intercept = 0;
  expfit_param1 = 0;
  expfit_param2 = 0;
  expfit_param3 = 0;


/*
  gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_A", NON_LIN_COEF_A);
  gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_B", NON_LIN_COEF_B);
  gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_C", NON_LIN_COEF_C);
  gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_D", NON_LIN_COEF_D);
  gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_E", NON_LIN_COEF_E);
*/

  gPARMS->SetDefaultParameter("FCAL:linfit_slope", linfit_slope);
  gPARMS->SetDefaultParameter("FCAL:linfit_intercept", linfit_intercept);
  gPARMS->SetDefaultParameter("FCAL:expfit_param1", expfit_param1);
  gPARMS->SetDefaultParameter("FCAL:expfit_param2", expfit_param2);
  gPARMS->SetDefaultParameter("FCAL:expfit_param3", expfit_param3);


  // Parameters to make shower-depth correction taken from Radphi, 
  // slightly modifed to match photon-polar angle
  FCAL_RADIATION_LENGTH = 3.1;
  FCAL_CRITICAL_ENERGY = 0.035;
  FCAL_SHOWER_OFFSET = 1.0;
	
  gPARMS->SetDefaultParameter("FCAL:FCAL_RADIATION_LENGTH", FCAL_RADIATION_LENGTH);
  gPARMS->SetDefaultParameter("FCAL:FCAL_CRITICAL_ENERGY", FCAL_CRITICAL_ENERGY);
  gPARMS->SetDefaultParameter("FCAL:FCAL_SHOWER_OFFSET", FCAL_SHOWER_OFFSET);

}

//------------------
// brun
//------------------
jerror_t DFCALShower_factory::brun(JEventLoop *loop, int32_t runnumber)
{
 
    // Get calibration constants
    map<string, double> fcal_parms;
    loop->GetCalib("FCAL/fcal_parms", fcal_parms);
    if (fcal_parms.find("FCAL_C_EFFECTIVE")!=fcal_parms.end()){
	FCAL_C_EFFECTIVE = fcal_parms["FCAL_C_EFFECTIVE"];
	if(debug_level>0)jout<<"FCAL_C_EFFECTIVE = "<<FCAL_C_EFFECTIVE<<endl;
    } else {
	jerr<<"Unable to get FCAL_C_EFFECTIVE from FCAL/fcal_parms in Calib database!"<<endl;
    }
  
    DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
    const DGeometry *geom = dapp->GetDGeometry(runnumber);
    
    if (geom) {
	geom->GetTargetZ(m_zTarget);
	geom->GetFCALZ(m_FCALfront);
    }
    else{
      
      cerr << "No geometry accessbile." << endl;
      return RESOURCE_UNAVAILABLE;
    }

    // by default, load non-linear shower corrections from the CCDB
    // but allow these to be overridden by command line parameters
    if(LOAD_CCDB_CONSTANTS > 0.1) {
	map<string, double> shower_calib_piecewise;
	loop->GetCalib("FCAL/shower_calib_piecewise", shower_calib_piecewise);
	linfit_slope = shower_calib_piecewise["linfit_slope"];
	linfit_intercept = shower_calib_piecewise["linfit_intercept"];
	expfit_param1 = shower_calib_piecewise["expfit_param1"];
	expfit_param2 = shower_calib_piecewise["expfit_param2"];
	expfit_param3 = shower_calib_piecewise["expfit_param3"];

	if(debug_level>0) {
	    jout << "linfit_slope = " << linfit_slope << endl;
	    jout << "linfit_intercept = " << linfit_intercept << endl;
	    jout << "expfit_param1 = " << expfit_param1 << endl;
	    jout << "expfit_param2 = " << expfit_param2<< endl;
	    jout << "expfit_param3 = " << expfit_param3 << endl;

	}
    }
    
    return NOERROR;
}


//------------------
// evnt
//------------------
jerror_t DFCALShower_factory::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
  vector<const DFCALCluster*> fcalClusters;
  eventLoop->Get(fcalClusters);
  if(fcalClusters.size()<1)return NOERROR;
 
  // Use the center of the target as an approximation for the vertex position
  DVector3 vertex(0.0, 0.0, m_zTarget);

  // Loop over list of DFCALCluster objects and calculate the "Non-linear" corrected
  // energy and position for each. We'll use a logarithmic energy-weighting to 
  // find the final position and error. 
  for( vector< const DFCALCluster* >::const_iterator clItr = fcalClusters.begin();
       clItr != fcalClusters.end();  ++clItr ){
    const DFCALCluster* cluster=*clItr;

    double cTime = cluster->getTime();
 		
    double errX = cluster->getRMS_x();
    double errY = cluster->getRMS_y();
    double errZ;  // will be filled by call to GetCorrectedEnergyAndPosition()
		
    // Get corrected energy, position, and errZ
    double Ecorrected;
    DVector3 pos_corrected;
    GetCorrectedEnergyAndPosition( cluster , Ecorrected, pos_corrected, errZ, &vertex);

    if (Ecorrected>0.){		
      //up to this point, all times have been times at which light reaches
      //the back of the detector. Here we correct for the time that it 
      //takes the Cherenkov light to reach the back of the detector
      //so that the t reported is roughly the time of the shower at the
      //position pos_corrected	
      cTime -= ( m_FCALfront + DFCALGeometry::blockLength() - pos_corrected.Z() )/FCAL_C_EFFECTIVE;

      // Make the DFCALShower object
      DFCALShower* shower = new DFCALShower;
      
      shower->setEnergy( Ecorrected );
      shower->setPosition( pos_corrected );   
      shower->setPosError( errX, errY, errZ );
      shower->setTime ( cTime );
      
      shower->AddAssociatedObject(cluster);

      _data.push_back(shower);
    }
  }

  return NOERROR;
}

//--------------------------------
// GetCorrectedEnergyAndPosition
//
// Non-linear and depth corrections should be fixed within DFCALShower member functions
//--------------------------------
void DFCALShower_factory::GetCorrectedEnergyAndPosition(const DFCALCluster* cluster, double &Ecorrected, DVector3 &pos_corrected, double &errZ, const DVector3 *vertex)
{
  // Non-linear energy correction are done here
  //int MAXITER = 1000;

  DVector3  posInCal = cluster->getCentroid();
  float x0 = posInCal.Px();
  float y0 = posInCal.Py();

  double Eclust = cluster->getEnergy();
  
  double A  = linfit_slope;
  double B  = linfit_intercept;
  double C  = expfit_param1;
  double D  = expfit_param2;
  double E  = expfit_param3;


 //double alfa  = NON_LIN_COEF_alfa;
	 
  double Egamma = 0.;
  
  // 06/02/2016 Shower Non-linearity Correction by Adesh. Uses new shower_calib_piecewise table in ccdb to include constants from linear fit.
  
  if ( Eclust <= 1.5 ) { 
  
  Egamma = Eclust/(A*Eclust + B); // Linear part
  
  }
  
  if ( Eclust > 1.5 ) { 
  
  Egamma = Eclust/(C - exp(-D*Eclust+ E)); // Non-linear part
  
  }
  
  // End Adesh's Correction  
  
 /*
 // Original Correction
  if ( A > 0 ) { 
    
    Egamma = Eclust/A;

    for ( int niter=0; 1; niter++) {

      double energy = Egamma;
      double non_lin_part = pow(Egamma,1+alfa)/(B+C*Egamma);
      Egamma = Eclust/A - non_lin_part;
      if ( fabs( (Egamma-energy)/energy ) < 0.001 ) {
	break;
	
      }
      else if ( niter > MAXITER ) {
	
	cout << " Iteration failed for cluster energy " << Eclust << endl;
	Egamma  = 0;
        
	break;
	
      }
      
    }
    
  }
  else {
    cout  << "Warning: DFCALShower : parameter A=" <<  NON_LIN_COEF_A 
	  << " is not valid!" << endl; 
  }
  */

  // then depth corrections 
  if ( Egamma > 0 ) { 
    float dxV = x0-vertex->X();
    float dyV = y0-vertex->Y();
    float zV = vertex->Z();
   
    double z0 = m_FCALfront - zV;
    double zMax = FCAL_RADIATION_LENGTH*(FCAL_SHOWER_OFFSET 
					 + log(Egamma/FCAL_CRITICAL_ENERGY));
    double zed = z0;
    double zed1 = z0 + zMax;

    double r0 = sqrt(dxV*dxV + dyV*dyV );

    int niter;
    for ( niter=0; niter<100; niter++) {
      double tt = r0/zed1;
      zed = z0 + zMax/sqrt( 1 + tt*tt );
      if ( fabs( (zed-zed1) ) < 0.001) {
	break;
      }
      zed1 = zed;
    }
    
    posInCal.SetZ( zed + zV );
    errZ = zed - zed1;
  }
  
  Ecorrected = Egamma;
  pos_corrected = posInCal;
}




//
//    File: DFCALShower_factory.cc
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
// Edited: B. Schaefer 3/23/2012 removed radiation hard insert functionality

#include <math.h>
#include <DVector3.h>
//#include <DLorentzVector.h>
using namespace std;

#include "DFCALShower_factory.h"
#include "DFCALGeometry.h"
//#include "DFCALCluster.h"
#include "DFCALHit.h"
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

// Set of coefficients for non-linear energy corrections 

    SHOWER_ENERGY_THRESHOLD = 50*k_MeV;

    NON_LIN_COEF_A = 0.439287;
    NON_LIN_COEF_B = 0.503378; 
    NON_LIN_COEF_C = 1.80842;
    NON_LIN_COEF_alfa = 1+0.0724789;

    //Loose timing cut
    SHOWER_TIMING_WINDOW = 5.0;
    
// Parameters to make shower-depth correction taken from Radphi, 
// slightly modifed to match photon-polar angle
        FCAL_RADIATION_LENGTH = 3.1;
        FCAL_CRITICAL_ENERGY = 0.035;
        FCAL_SHOWER_OFFSET = 1.0;

	gPARMS->SetDefaultParameter("FCAL:SHOWER_ENERGY_THRESHOLD", SHOWER_ENERGY_THRESHOLD);

	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_A", NON_LIN_COEF_A);
	if (NON_LIN_COEF_A<=0.) {
	  cout  << "Warning: DFCALShower : parameter A=" <<  NON_LIN_COEF_A 
		<< " is not valid!" << endl; 
	}


	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_B", NON_LIN_COEF_B);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_C", NON_LIN_COEF_C);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_alfa", NON_LIN_COEF_alfa);

	gPARMS->SetDefaultParameter("FCAL:SHOWER_TIMING_WINDOW", SHOWER_TIMING_WINDOW);
	
	gPARMS->SetDefaultParameter("FCAL:FCAL_RADIATION_LENGTH", FCAL_RADIATION_LENGTH);
	gPARMS->SetDefaultParameter("FCAL:FCAL_CRITICAL_ENERGY", FCAL_CRITICAL_ENERGY);
	gPARMS->SetDefaultParameter("FCAL:FCAL_SHOWER_OFFSET", FCAL_SHOWER_OFFSET);

	gPARMS->SetDefaultParameter("FCAL:LOAD_NONLIN_CCDB", LOAD_CCDB_CONSTANTS);

}

//------------------
// brun
//------------------
// take merging out
jerror_t DFCALShower_factory::brun(JEventLoop *loop, int runnumber)
{
  /*// Get calibration constants
	map<string, double> cluster_merging;
	loop->GetCalib("FCAL/cluster_merging", cluster_merging);
	if(cluster_merging.find("MIN_CLUSTER_SEPARATION")!=cluster_merging.end()){
		MIN_CLUSTER_SEPARATION = cluster_merging["MIN_CLUSTER_SEPARATION"];
		if(debug_level>0)jout<<"MIN_CLUSTER_SEPARATION = "<<MIN_CLUSTER_SEPARATION<<endl;
	}else{
		jerr<<"Unable to get from MIN_CLUSTER_SEPARATION FCAL/cluster_merging in Calib database!"<<endl;
		loop->GetJApplication()->Quit();
        }*/

    // Get calibration constants
    FCAL_C_EFFECTIVE=15.;
    map<string, double> fcal_parms;
    loop->GetCalib("FCAL/fcal_parms", fcal_parms);
    if (fcal_parms.find("FCAL_C_EFFECTIVE")!=fcal_parms.end()){
	FCAL_C_EFFECTIVE = fcal_parms["FCAL_C_EFFECTIVE"];
	if(debug_level>0)jout<<"FCAL_C_EFFECTIVE = "<<FCAL_C_EFFECTIVE<<endl;
    } else {
	jerr<<"Unable to get FCAL_C_EFFECTIVE from FCAL/fcal_parms in Calib database!"<<endl;
    }
  
    m_zTarget=65.0*k_cm;
    m_FCALfront= DFCALGeometry::fcalFaceZ();
    m_FCALback= m_FCALfront+DFCALGeometry::blockLength();
    DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
    if (dapp) {
	const DGeometry *geom = dapp->GetDGeometry(runnumber);
	geom->GetTargetZ(m_zTarget);
	geom->GetFCALZ(m_FCALfront);
	
	vector<double>block;
	geom->Get("//box[@name='LGBL']/@X_Y_Z",block);
	m_FCALback=m_FCALfront+block[2];
    }

    // by default, load non-linear shower corrections from the CCDB
    // but allow these to be overridden by command line parameters
    if(LOAD_CCDB_CONSTANTS > 0.1) {
	map<string, double> shower_calib;
	loop->GetCalib("FCAL/shower_calib", shower_calib);
	NON_LIN_COEF_A = shower_calib["FCAL_SHOWER_CALIB_A"];
	NON_LIN_COEF_B = shower_calib["FCAL_SHOWER_CALIB_B"];
	NON_LIN_COEF_C = shower_calib["FCAL_SHOWER_CALIB_C"];
	NON_LIN_COEF_alfa = shower_calib["FCAL_SHOWER_CALIB_D"];
	if(debug_level>0) {
	    jout << "NON_LIN_COEF_A = " << NON_LIN_COEF_A << endl;
	    jout << "NON_LIN_COEF_B = " << NON_LIN_COEF_B << endl;
	    jout << "NON_LIN_COEF_C = " << NON_LIN_COEF_C << endl;
	    jout << "NON_LIN_COEF_alfa = " << NON_LIN_COEF_alfa << endl;
	}
    }
    
    return NOERROR;
}


//------------------
// evnt
//------------------
jerror_t DFCALShower_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
  vector<const DFCALCluster*> fcalClusters;
  eventLoop->Get(fcalClusters);
  if(fcalClusters.size()<1)return NOERROR;
 
  // Use the center of the target as an approximation for the vertex position
  DVector3 target(0.0, 0.0, m_zTarget);

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
    GetCorrectedEnergyAndPosition( cluster , Ecorrected, pos_corrected, errZ, &target);

    if (Ecorrected>0.){		
      //up to this point, all times have been times at which light reaches
      //the back of the detector. Here we correct for the time that it 
      //takes the Cherenkov light to reach the back of the detector
      //so that the t reported is roughly the time of the shower at the
      //position pos_corrected	
      cTime -= ( m_FCALback - pos_corrected.Z() )/FCAL_C_EFFECTIVE;

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
  int MAXITER = 1000;

  DVector3  posInCal = cluster->getCentroid();
  float x0 = posInCal.Px();
  float y0 = posInCal.Py();

  double Eclust = cluster->getEnergy();
  
  double A  = NON_LIN_COEF_A;
  double B  = NON_LIN_COEF_B;
  double C  = NON_LIN_COEF_C;
  double alfa  = NON_LIN_COEF_alfa;
	 
  double Egamma = 0.;
  
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




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
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_B", NON_LIN_COEF_B);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_C", NON_LIN_COEF_C);
	gPARMS->SetDefaultParameter("FCAL:NON_LIN_COEF_alfa", NON_LIN_COEF_alfa);

	gPARMS->SetDefaultParameter("FCAL:SHOWER_TIMING_WINDOW", SHOWER_TIMING_WINDOW);
	
	gPARMS->SetDefaultParameter("FCAL:FCAL_RADIATION_LENGTH", FCAL_RADIATION_LENGTH);
	gPARMS->SetDefaultParameter("FCAL:FCAL_CRITICAL_ENERGY", FCAL_CRITICAL_ENERGY);
	gPARMS->SetDefaultParameter("FCAL:FCAL_SHOWER_OFFSET", FCAL_SHOWER_OFFSET);

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
  DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
  if (dapp) {
    const DGeometry *geom = dapp->GetDGeometry(runnumber);
    geom->GetTargetZ(m_zTarget);
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

	//------------- the following is a temporary kludge...
//	const DVertex *vertex=NULL;

	//loop->GetSingle(vertex);
//	vector<const DVertex *>vertices;
//	eventLoop->Get(vertices);
//	if (vertices.size())   vertex=vertices[0];

	// Return immediately if there isn't even one cluster
	if(fcalClusters.size()<1)return NOERROR;

// ------------------------------------------------
/*	
	// Here we try and merge FCAL clusters simply by looking to see if they are
	// within a certain distance of one another. We do this by first looping
	// through the list of clusters as many times as is needed until we get
	// a set of lists for which all clusters are at least MIN_CLUSTER_SEPARATION
	// from all clusters in all other lists.
	
	// Start by filling out out lists with one cluster each
	vector<vector<const DFCALCluster*> > merge_lists;
	for ( unsigned int i = 0; i < fcalClusters.size(); i++ ) {
		vector<const DFCALCluster*> merge_list;
		merge_list.push_back(fcalClusters[i]);
		merge_lists.push_back(merge_list);
	}
	
	// Loop until we find no more lists to merge
	bool clusters_were_merged;
	do{
		clusters_were_merged = false;
		
		// Loop over all pairs of cluster lists
		for(unsigned int i=0; i<merge_lists.size(); i++){
			vector<const DFCALCluster*> &merge_list1 = merge_lists[i];
			for(unsigned int j=i+1; j<merge_lists.size(); j++){
				vector<const DFCALCluster*> &merge_list2 = merge_lists[j];
				
				// Loop over all elements of both lists to see if any are
				// within MIN_CLUSTER_SEPARATION.
				for(unsigned int k=0; k<merge_list1.size(); k++){
					DVector3 pos1 = merge_list1[k]->getCentroid();
					for(unsigned int m=0; m<merge_list2.size(); m++){
						DVector3 pos2 = merge_list2[m]->getCentroid();
						
						double separation_xy = (pos1-pos2).Perp();
						if(separation_xy<MIN_CLUSTER_SEPARATION){
							// Phew! if we got here then we need to merge the 2 clusters
							// The easiest way to do this is to just add all clusters
							// from merge_list2 to merge_list1 and clear merge_list2.
							// Then, we ignore empty lists below.
							merge_list1.insert(merge_list1.end(), merge_list2.begin(), merge_list2.end());
							merge_list2.clear();
							clusters_were_merged = true;
						}
					}
					if(clusters_were_merged)break; // we'll need to do the outer "do" loop again anyway so bail now
				}
				if(clusters_were_merged)break; // we'll need to do the outer "do" loop again anyway so bail now
			}
			if(clusters_were_merged)break; // we'll need to do the outer "do" loop again anyway so bail now
		}
	
	}while(clusters_were_merged);

	// Now we loop over the lists of clusters to merge and make a
	// DFCALShower from the list. Note that it may well be that each
	// list is still only 1 element long!
	for ( unsigned int i = 0; i < merge_lists.size(); i++ ) {
		vector<const DFCALCluster*> &merge_list = merge_lists[i];
		if(merge_list.size()<1)continue; // ignore empty lists (see comment above)
	
		DFCALShower* fcalPhoton = makePhoton( merge_list, vertex );

		if ( fcalPhoton->getEnergy() <= 0  ) {
		  cout << "Deleting fcalPhoton " << endl;
		  delete fcalPhoton; 
		  continue;
		}else {
			_data.push_back(fcalPhoton);
		}
	} 
*/
// --------------------------
       for( vector< const DFCALCluster* >::const_iterator clItr = fcalClusters.begin();
	       clItr != fcalClusters.end();  ++clItr ){

		DFCALShower* fcalShower = makeFcalShower( *clItr );

		//Make a very loose timing cut here to eliminate out-of-time EM bkgd.
		//Photon travels dist1 (from center of target to shower maximum
		//position) at the speed of light and the remainder of the distance,
		//dist2, to the back of the FCAL (where timing is measured) at c_eff.
		//This travel time is the expected_time. Since this time is only used
		//to make a very loose timing cut, we can ignore other factors such as
		//the fact that not all particles originate exactly at the center of
		//the target.
		//Later analysis should use a tighter timing cut incoporating vertex
		//position and time.
		DVector3 pos = fcalShower->getPosition();
		DVector3 vertex(0.0, 0.0, m_zTarget);
		double dist1 = (pos-vertex).Mag();
		double dist2 = DFCALGeometry::fcalFaceZ() + DFCALGeometry::blockLength() - pos.Z();
		double expected_time = dist1/SPEED_OF_LIGHT + dist2/FCAL_C_EFFECTIVE;

		if ( fcalShower->getEnergy() <= 0  ) {
		  cout << "Deleting fcalShower " << endl;
		  delete fcalShower; 
		  continue;
		} else if ( fabs(fcalShower->getTime() - expected_time) > SHOWER_TIMING_WINDOW ){
		  delete fcalShower; 
		  continue;          
		} else {
			_data.push_back(fcalShower);
		}

        }

	return NOERROR;
}


//--------------------------------
// makePhoton
//--------------------------------
DFCALShower* DFCALShower_factory::makeFcalShower( const DFCALCluster* cluster ) 
{

	// Loop over list of DFCALCluster objects and calculate the "Non-linear" corrected
	// energy and position for each. We'll use a logarithmic energy-weighting to 
	// find the final position and error.

// Use target center if vertex does not exist 
        DVector3 target(0.0, 0.0, m_zTarget);

        double cTime = cluster->getTime();
 		
	double errX = cluster->getRMS_x();
	double errY = cluster->getRMS_y();
	double errZ;  // will be filled by call to GetCorrectedEnergyAndPosition()
		
		// Get corrected energy, position, and errZ
	double Ecorrected;
	DVector3 pos_corrected;
	GetCorrectedEnergyAndPosition( cluster , Ecorrected, pos_corrected, errZ, &target);
		
	
	// Make the DFCALShower object
	DFCALShower* shower = new DFCALShower;

	shower->setEnergy( Ecorrected );
	shower->setPosition( pos_corrected );   
	shower->setPosError( errX, errY, errZ );
	shower->setTime ( cTime );

	shower->AddAssociatedObject(cluster);

	return shower;
}

//--------------------------------
// GetCorrectedEnergyAndPosition
//
// Non-linear and depth corrections should be fixed within DFCALShower member functions
//--------------------------------
void DFCALShower_factory::GetCorrectedEnergyAndPosition(const DFCALCluster* cluster, double &Ecorrected, DVector3 &pos_corrected, double &errZ, const DVector3 *vertex)
{
// Non-linar energy correction are done here
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

          cout << "Warning: DFCALShower : parameter A" << A << " is not valid" << endl; 
       }


// then depth corrections 

       if ( Egamma > 0 ) { 
               float xV = vertex->X();
               float yV = vertex->Y();
               float zV = vertex->Z();
            

           double z0 = DFCALGeometry::fcalFaceZ() - zV;
           double zMax = (FCAL_RADIATION_LENGTH*(
                       FCAL_SHOWER_OFFSET + log(Egamma/FCAL_CRITICAL_ENERGY)));
           double zed = z0;
           double zed1 = z0 + zMax;

           double r0 = sqrt( (x0-xV)*(x0-xV) + (y0-yV)*(y0-yV) );

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




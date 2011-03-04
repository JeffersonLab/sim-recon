/*
 *  DBCALPhoton_factory.cc
 *  Hall D
 *
 *  Created by Matthew Shepherd on 7/23/07.
 *
 */

#include <iostream>
using namespace std;

#include <JANA/JApplication.h>

#include "BCAL/DBCALPhoton_factory.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALGeometry.h"
#include <PID/DVertex.h>

#include "DLorentzVector.h"

DBCALPhoton_factory::DBCALPhoton_factory()
{

    m_scaleZ_p0GE =  0.660299;
    m_scaleZ_p1GE =  0.00229035;
    m_scaleZ_p2GE =  -7.8725e-06;
    m_scaleZ_p3GE =  8.05729e-09;
     
    m_nonlinZ_p0GE =  0.117;
    m_nonlinZ_p1GE =  -3.79638e-04;
    m_nonlinZ_p2GE =  3.770e-07;    
    m_nonlinZ_p3GE =  1.9274e-10;

  /*    
    
    m_linZ_p0GE = -5.26475e-03;
    m_linZ_p1GE = -2.47419e-02;
    m_linZ_p2GE =  4.19082e01;    
    m_linZ_p3GE =  6.69810e01;
    
   
   //scaling parameter set for Z>370 (end of module)

    m_scaleZ_p0 =  0.8284;
    m_scaleZ_p1 =  3.3;
    m_scaleZ_p2 =  422.5;
    m_scaleZ_p3 =  12.04;
    m_scaleZ_p4 = 0.0;
     
    m_nonlinZ_p0 =  0.05136;
    m_nonlinZ_p1 = 1000.0;
    m_nonlinZ_p2 =  453.6;    
    m_nonlinZ_p3 = 17.21;

    m_linZ_p0 = -7.166e-03;
    m_linZ_p1 = -1000.0;
    m_linZ_p2 = 482.0;
    m_linZ_p3 = 24.19;
  


 // parameters for correcting with dark noise (different correction function)

    m_scaleZ_p0GE =  0.94795;
    m_scaleZ_p1GE =  -1.747;
    m_scaleZ_p2GE =  1000;
    m_scaleZ_p3GE = 251.5;
     
    m_nonlinZ_p0GE =  0.02597;
    m_nonlinZ_p1GE = -0.0347;
    m_nonlinZ_p2GE =  54.89;    
    m_nonlinZ_p3GE =  55.53;
    
    m_linZ_p0GE = -3.0533e-03;
    m_linZ_p1GE = -0.01671;
    m_linZ_p2GE =  14.363;    
    m_linZ_p3GE =  69.542;
    
   
   //scaling parameter set for Z>370 (end of module)

    m_scaleZ_p0 =  0.8776;
    m_scaleZ_p1 =  -10.0;
    m_scaleZ_p2 =  428.2;
    m_scaleZ_p3 =  12.22;
    m_scaleZ_p4 = 0.0;
     
    m_nonlinZ_p0 = 0.02604;
    m_nonlinZ_p1 = -10.1;
    m_nonlinZ_p2 =  433;    
    m_nonlinZ_p3 = 13.38;

    m_linZ_p0 = -3.3587e-03;
    m_linZ_p1 = -10.0;
    m_linZ_p2 = 456.7;
    m_linZ_p3 = 17.69;
   */

}

//------------------
// brun
//------------------
jerror_t DBCALPhoton_factory::brun(JEventLoop *loop, int runnumber)
{
    
    vector<const DBCALGeometry*> bcalGeomVect;
    loop->Get( bcalGeomVect );
    const DBCALGeometry& bcalGeom = *(bcalGeomVect[0]);

    m_bcalIR = bcalGeom.BCALINNERRAD;
    m_zTarget = 65;                    // global target position -- should come from database!

	// Get calibration constants
	map<string, double> cluster_merging;
	loop->GetCalib("BCAL/cluster_merging", cluster_merging);
	
	// Check that our constants are in the returned map. If not, print error and exit
	if(	cluster_merging.find("MIN_CLUSTER_SEPARATION_XY")==cluster_merging.end()
		||	cluster_merging.find("MIN_CLUSTER_SEPARATION_Z" )==cluster_merging.end()){
		jerr<<"Unable to get MIN_CLUSTER_SEPARATION from BCAL/cluster_merging in Calib database!"<<endl;
		exit(-1);
	}
	
	// Copy from container into local data members and optionally report the values
	MIN_CLUSTER_SEPARATION_XY = cluster_merging["MIN_CLUSTER_SEPARATION_XY"];
	MIN_CLUSTER_SEPARATION_Z = cluster_merging["MIN_CLUSTER_SEPARATION_Z"];
	if(debug_level>0)jout<<"MIN_CLUSTER_SEPARATION_XY = "<<MIN_CLUSTER_SEPARATION_XY<<endl;
	if(debug_level>0)jout<<"MIN_CLUSTER_SEPARATION_Z = "<<MIN_CLUSTER_SEPARATION_Z<<endl;
   
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALPhoton_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// This factory has been touched by several authors which may explain a little
	// why it is a bit of a mosaic. Here is how it works based on an incomplete
	// knowledge of the prior steps:
	//
	// Initial reconstruction involving clustering of hits and showers
	// is done in DBCALShower_factory.cc. This file applies corrections for dark noise
	// and non-linear scaling factors for calibration.
	//
	// At this point (where I come in) there are frequently multiple BCAL clusters
	// that are near to one another and should be merged into a single DBCALPhoton
	// object. This works by first making a DBCALPhoton object for every DBCALShower
	// and then merging the ones close in proximity into a single DBCALPhoton
	// (deleting the rest).
	// Calibration parameters controlling when to merge clusters are kept in
	// BCAL/cluster_merging.
	
	// This follows a similar algorithm used in DFCALPhoto_factory.cc. Specifically,
	// we generate a list of lists of DBCALPhoton objects that should be merged to
	// form a final DBCALPhoton object. Each list is initially just one element long
	// but grows or shrinks as lists are merged. When two lists are merged, one list
	// gets all of the elements and the other is cleared to have zero elements. Doing
	// it this way means the final merging is independent of the order of the 
	// shower objects (one can imagine the centroid of 2 merged showers is now out of
	// range for merging with a 3rd while one of the original 2 was within range.) 


	// Get DBCALShower objects
	vector< const DBCALShower* > showerVect;
	loop->Get( showerVect );

	//------------- the following is a temporary kludge...
	const DVertex *vertex=NULL;
	//loop->GetSingle(vertex);
	vector<const DVertex *>vertices;
	loop->Get(vertices);
	if (vertices.size()) vertex=vertices[0];

	// Make a list with a single DBCALPhoton for each DBCALShower object
	vector<vector<DBCALPhoton*> > merge_lists;
	for(unsigned int i=0; i<showerVect.size(); i++){
		vector<DBCALPhoton*> merge_list;
		merge_list.push_back(MakeDBCALPhoton(showerVect[i], vertex));
		merge_lists.push_back(merge_list);
	}
	
	// Loop until we find no more lists to merge
	bool photons_were_merged;
	do{
		photons_were_merged = false;
		
		// Loop over all pairs of photon lists
		for(unsigned int i=0; i<merge_lists.size(); i++){
			vector<DBCALPhoton*> &merge_list1 = merge_lists[i];
			for(unsigned int j=i+1; j<merge_lists.size(); j++){
				vector<DBCALPhoton*> &merge_list2 = merge_lists[j];
				
				// Loop over all elements of both lists to see if any are
				// within MIN_CLUSTER_SEPARATION.
				for(unsigned int k=0; k<merge_list1.size(); k++){
					DVector3 pos1 = merge_list1[k]->showerPosition();
					for(unsigned int m=0; m<merge_list2.size(); m++){
						DVector3 pos2 = merge_list2[m]->showerPosition();
						
						double separation_xy = (pos1-pos2).Perp();
						double separation_z = (pos1-pos2).Z();
						if(separation_xy<MIN_CLUSTER_SEPARATION_XY && separation_z<MIN_CLUSTER_SEPARATION_Z){
							// Phew! if we got here then we need to merge the 2 photons
							// The easiest way to do this is to just add all photons
							// from merge_list2 to merge_list1 and clear merge_list2.
							// Then, we ignore empty lists below.
							merge_list1.insert(merge_list1.end(), merge_list2.begin(), merge_list2.end());
							merge_list2.clear();
							photons_were_merged = true;
						}
					}
					if(photons_were_merged)break; // we'll need to do the outer "do" loop again anyway so bail now
				}
				if(photons_were_merged)break; // we'll need to do the outer "do" loop again anyway so bail now
			}
			if(photons_were_merged)break; // we'll need to do the outer "do" loop again anyway so bail now
		}
	
	}while(photons_were_merged);
	
	// Now we loop over the lists of clusters to merge and make a
	// DBCALPhoton from the list. Note that it may well be that each
	// list is still only 1 element long!
	for ( unsigned int i = 0; i < merge_lists.size(); i++ ) {
		vector<DBCALPhoton*> &merge_list = merge_lists[i];
		if(merge_list.size()<1)continue; // ignore empty lists (see comment above)
	
		DBCALPhoton* photon = MergeDBCALPhotons( merge_list );

		if(photon!=NULL)_data.push_back(photon);
	}
	
	return NOERROR;
}

//------------------
// MergeDBCALPhotons
//------------------
DBCALPhoton* DBCALPhoton_factory::MergeDBCALPhotons(vector<DBCALPhoton*> &photons)
{
	// Loop over all photons in the list of photons to be merged and calculate
	// the weighted sums of certain values to be used for the final "merged"
	// photon parameters below. Position of the merged photon is calculated
	// using a logarithmic weight of the energy. Photons with less than
	// 1 MeV are dropped.

	double sum_weight = 0.0;
	double Etotal = 0.0;
	double earliest_time = 1.0E6;
	DVector3 pos_total(0.0, 0.0, 0.0);
	double Ebiggest = 0.0;
	DVector3 fitLayPoint;
	DVector3 fitLayPointErr;
	DVector3 fitLaySlope;
	DVector3 fitLaySlopeErr;
	vector<const JObject*> all_associated_objects;
	for(unsigned int i=0; i<photons.size(); i++){
		DBCALPhoton *photon = photons[i];
		
		// Weight using log of energy
		double E = photon->lorentzMomentum().E();
		double weight = log10(E/0.001);
		if(weight<0.0) continue; // ignore clusters with less than 1MeV
		
		sum_weight += weight;
		Etotal += E;
		pos_total += weight*photon->showerPosition();
		if(photon->showerTime()<earliest_time) earliest_time = photon->showerTime();
		
		// Add list of associated objects of this photon to temporary
		// list that will be added back to the one, merged photon below.
		vector<const JObject*> myobjects;
		photon->GetT(myobjects);
		all_associated_objects.insert(all_associated_objects.end(), myobjects.begin(), myobjects.end());
		
		// ################### WARNING #################
		// WE NEED TO RECALCULATE THE VALUES FOR THE POINT, POINTERR, AND
		// SLOPE, SLOPEERR MEMBERS. WITHOUT KNOWING BETTER HOW TO DO IT,
		// WE JUST USE THE VALUES FROM THE MOST ENERGETIC CLUSTER.
		if(E>Ebiggest){
			Ebiggest = E;
			fitLayPoint = photon->fitLayPoint();
			fitLayPointErr = photon->fitLayPointErr();
			fitLaySlope = photon->fitLaySlope();
			fitLaySlopeErr = photon->fitLaySlopeErr();
		}
	}
	
	// Make sure at least one photon with >1MeV of energy was in the list. If so,
	// copy the "merged" values into the first photon object.
	if(sum_weight>0.0){

		pos_total *= 1.0/sum_weight;
		
		DBCALPhoton *photon = photons[0];
		DVector3 mom = pos_total - DVector3(0.0, 0.0, m_zTarget);
		mom.SetMag(Etotal);
		
		photon->setShowerPosition(pos_total);
		photon->setShowerTime(earliest_time);
		photon->setLorentzMomentum(DLorentzVector(mom, Etotal));
		
		photon->setFitLayPoint(fitLayPoint);
		photon->setFitLayPointErr(fitLayPointErr);
		photon->setFitLaySlope(fitLaySlope);
		photon->setFitLaySlopeErr(fitLaySlopeErr);
		
		// To ensure there are no duplicates, we need to remove any existing associated
		// objects in photon and then add back in the complete list from associated_objects_total
		vector<const JObject*> myobjects;
		photon->GetT(myobjects);
		for(unsigned int i=0; i<myobjects.size(); i++)photon->RemoveAssociatedObject(myobjects[i]);
		
		// Add back in all associated objects from all showers being merged
		for(unsigned int i=0; i<all_associated_objects.size(); i++){
			photon->AddAssociatedObject(all_associated_objects[i]);
		}
	}else{
		// If we get here then there were no photons in the list with
		// more than 1MeV of energy. Delete the first element (the
		// rest will be deleted below.)
		if(photons.size()>0)delete photons[0];
		photons[0] = NULL;
	}
	
	// Delete all but the first DBCALPhoton since they were all merged into it
	for(unsigned int i=1; i<photons.size(); i++)delete photons[i];
	
	return photons.size()>0 ? photons[0]:NULL;
}

//------------------
// MakeDBCALPhoton
//------------------
DBCALPhoton* DBCALPhoton_factory::MakeDBCALPhoton(const DBCALShower* shower, const DVertex *vertex)
{ 
                
	double xSh = shower->x;
	double ySh = shower->y;
	double zSh = shower->z;
	//     int nCell = shower->N_cell;        

	// Get vertex position as DVector3
	DVector3 my_vertex(0.0, 0.0, m_zTarget);
	if(vertex)my_vertex = vertex->x.Vect();
	
	// Momentum direction is vector pointing from vertex to shower center
	DVector3 pdir = DVector3(xSh, ySh, zSh) - my_vertex;
	pdir.SetMag(1.0);
	
	// We want to find the entry point of the photon into the BCAL. This is
	// determined by writing the vector pointing from the vertex to the
	// BCAL as:
	//
	//    v + alpha*pdir
	//
	// where v is the vertex position vector, pdir is the unit vector
	// pointing in the momentum direction, and alpha is a scalar value
	// to be solved for.
	//
	// The equation to solve is determined by setting the "R" value of the 
	// above equation to the inner BCAL radius
	//
	//   (v_x + alpha*p_x)^2 + (v_y + alpha*p_y)^2 = R^2
	//
	// This is quadratic in alpha with the coefficients
	//
	//   a = p_x^2 + p_y^2           (which is just pdir.Perp2())
	//   b = 2*(v_x*p_x + v_y*p_y)   
	//   c = v_x^2 + v_y^2 - R^2     (which is v.Perp2() - R^2)
	//
	// Note that there are 2 roots of the quadratic equation which should
	// correspond to a positive alpha solution and a negative alpha 
	// solution (assuming the vertex is inside the BCAL inner radius R).
	// We therefore just take the "+" solution of the Q.E. since that
	// has to be the more positive of the 2.
	double a = pdir.Perp2();
	double b = 2.0*(my_vertex.X()*pdir.X() + my_vertex.Y()*pdir.Y());
	double c = my_vertex.Perp2() - m_bcalIR*m_bcalIR;
	double alpha = (-b + sqrt(b*b - 4.0*a*c))/2.0/a;
	DVector3 x_surface = my_vertex + alpha*pdir;

	// get z where shower enters BCAL (this corresponds to generated z in tuning MC)
	//double zEntry = zSh - ( ( zSh - m_zTarget ) * 
	//							  ( 1 - m_bcalIR / ( sqrt( xSh * xSh + ySh * ySh ) ) ) ); 
	double zEntry = x_surface.Z();

	// calibrate energy:
	// Energy calibration has a z dependence -- the
	// calibration comes from fitting E_rec / E_gen to scale * E_gen^nonlin
	// for slices of z.  These fit parameters (scale and nonlin) are then plotted 
	// as a function of z and fit.
    
	/*
	if( zEntry < 370.0 ) {
		scale = (m_scaleZ_p0GE  + m_scaleZ_p1GE *(exp( -0.5 *(zEntry - m_scaleZ_p2GE )* (zEntry - m_scaleZ_p2GE ) / (m_scaleZ_p3GE * m_scaleZ_p3GE)   ) ) );
		nonlin =( m_nonlinZ_p0GE  + m_nonlinZ_p1GE *(exp( -0.5 *(zEntry - m_nonlinZ_p2GE )* (zEntry - m_nonlinZ_p2GE ) / (m_nonlinZ_p3GE * m_nonlinZ_p3GE)   ) ) ) ;
		lin = ( m_linZ_p0GE  + m_linZ_p1GE *(exp( -0.5 *(zEntry - m_linZ_p2GE )* (zEntry - m_linZ_p2GE ) / (m_linZ_p3GE * m_linZ_p3GE)   ) ) ) ;

		//	nonlin = 0.0; // fixed value for debug
		//	   lin = 0.0; // fixed value for debug
		//       scale = 1.0; // fixed value for debug
	}

	if( zEntry >= 370.0 ) {
		scale = m_scaleZ_p0 +  m_scaleZ_p1 *(exp( -0.5 *(zEntry - m_scaleZ_p2 )* (zEntry - m_scaleZ_p2 ) / (m_scaleZ_p3 * m_scaleZ_p3)   ) ) ;
		nonlin = m_nonlinZ_p0  + m_nonlinZ_p1 *(exp( -0.5 *(zEntry - m_nonlinZ_p2 )* (zEntry - m_nonlinZ_p2 ) / (m_nonlinZ_p3 * m_nonlinZ_p3)   ) )  ;
		lin = m_linZ_p0 + m_linZ_p1 *(exp( -0.5 *(zEntry - m_linZ_p2 )* (zEntry - m_linZ_p2 ) / (m_linZ_p3 * m_linZ_p3)   ) )  ;

		//  cout << scale << ' ' << nonlin << ' ' << lin << endl;    
	}*/

	scale = m_scaleZ_p0GE + m_scaleZ_p1GE*zEntry + m_scaleZ_p2GE*(zEntry*zEntry) + m_scaleZ_p3GE*(zEntry*zEntry*zEntry);
	nonlin = m_nonlinZ_p0GE + m_nonlinZ_p1GE*zEntry + m_nonlinZ_p2GE*(zEntry*zEntry)+ m_nonlinZ_p3GE*(zEntry*zEntry*zEntry);

	// if( zEntry < m_zTarget ) zEntry = m_zTarget;

		  //end of BCAL calibration        


		  
	// now turn E_rec into E_gen -->> E_gen = ( E_rec / scale ) ^ ( 1 / ( 1 + nonlin ) )
	double energy = pow( (shower->E ) / scale, 1 / ( 1 + nonlin ) );


	DBCALPhoton* photon = new DBCALPhoton();

	DVector3 mom = energy*pdir;
	photon->setLorentzMomentum( DLorentzVector( mom, energy ) );

	photon->setShowerPosition( DVector3( xSh, ySh, zSh ) );

	photon->setShowerTime( shower->t );

	photon->setFitLayPoint( DVector3( shower->Apx_x,
												shower->Apx_y,
												shower->Apx_z ) );

	photon->setFitLayPointErr( DVector3( shower->error_Apx_x,
													shower->error_Apx_y,
													shower->error_Apx_z ) );

	photon->setFitLaySlope( DVector3( shower->Cx,
												shower->Cy,
												shower->Cz ) );

	photon->setFitLaySlopeErr( DVector3( shower->error_Cx,
													shower->error_Cy,
													shower->error_Cz ) );

	photon->AddAssociatedObject(shower);
	photon->AddAssociatedObject(vertex);

	return photon;
}


/*
 *  DFCALPhoton_factory.cc
 *  Hall D
 *
 *  Created by Mihajlo Kornicer 5/31/11.
 *
 */

#include <iostream>
using namespace std;

#include <JANA/JApplication.h>

#include "FCAL/DFCALPhoton_factory.h"
#include "FCAL/DFCALShower.h"
//#include "FCAL/DFCALGeometry.h"
#include "PID/DVertex.h"
#include "DLorentzVector.h"

#include "units.h"

DFCALPhoton_factory::DFCALPhoton_factory(){
 
  m_zTarget = 65*k_cm;
  
}

//------------------
// brun
//------------------
//jerror_t DFCALPhoton_factory::brun(JEventLoop *loop, int runnumber)
//{
//    
//  m_zTarget = 65;                    // global target position -- should come from database!
//   
//  return NOERROR;
//}

//------------------
// evnt
//------------------
jerror_t DFCALPhoton_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Get DFCALShower objects
	vector< const DFCALShower* > showerVect;
        loop->Get( showerVect );
	
	const DVertex* vertex = NULL;
	vector< const DVertex* >vertices;
	loop->Get(vertices);
	if ( vertices.size() ) vertex = vertices[0];

  for( vector< const DFCALShower* >::const_iterator shItr = showerVect.begin();
      shItr != showerVect.end();
      ++shItr ){
    
    _data.push_back( MakeDFCALPhoton( *shItr, vertex ) );
    
  }
  
	
	return NOERROR;
}


//------------------
// MakeDFCALPhoton
//------------------
DFCALPhoton* DFCALPhoton_factory::MakeDFCALPhoton(const DFCALShower* shower, const DVertex *vertex)
{ 
                

	// Get vertex position as DVector3
	DVector3 my_vertex(0.0, 0.0, m_zTarget);
	if(vertex)my_vertex = vertex->x.Vect();
	
	// Momentum direction is vector pointing from vertex to shower center
	DVector3 pdir = shower->getPosition() - my_vertex;
	pdir.SetMag(1.0);
	
	double energy = shower->getEnergy();

	DFCALPhoton* photon = new DFCALPhoton();

	DVector3 mom = energy * pdir;

	photon->setLorentzMomentum( DLorentzVector( mom, energy ) );
	photon->setShowerPosition( shower->getPosition() );
        photon->setShowerPositionErr( shower->getPositionError() );
	photon->setShowerTime( shower->getTime() );
  

	if(shower)photon->AddAssociatedObject(shower);
	if(vertex)photon->AddAssociatedObject(vertex);

	return photon;
}


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
#include "PID/DVertex.h"
#include "DLorentzVector.h"

#include "units.h"

DBCALPhoton_factory::DBCALPhoton_factory(){
 
  m_zTarget = 65*k_cm;
  m_bcalIR = DBCALGeometry::BCALINNERRAD;
  
  USE_KLOE = 1;
  
  gPARMS->SetDefaultParameter( "BCALRECON:USE_KLOE", USE_KLOE );
}

jerror_t DBCALPhoton_factory::init(){
 
  if( USE_KLOE ){
    
    cout << "Using KLOE BCAL clustering." << endl;
  }
  else{
   
    cout << "Using alternative (experimental) BCAL clustering." << endl;
  }
  
  return NOERROR;
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
   
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALPhoton_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Get DBCALShower objects
	vector< const DBCALShower* > showerVect;

	if( USE_KLOE ){
    
    loop->Get( showerVect, "KLOE" );
  }
  else{
    
    loop->Get( showerVect );
  }
	
	const DVertex* vertex = NULL;
	vector< const DVertex* >vertices;
	loop->Get(vertices);
	if ( vertices.size() ) vertex = vertices[0];

  for( vector< const DBCALShower* >::const_iterator shItr = showerVect.begin();
      shItr != showerVect.end();
      ++shItr ){
    
    _data.push_back( MakeDBCALPhoton( *shItr, vertex ) );
    
  }
  
	
	return NOERROR;
}


//------------------
// MakeDBCALPhoton
//------------------
DBCALPhoton* DBCALPhoton_factory::MakeDBCALPhoton(const DBCALShower* shower, const DVertex *vertex)
{ 
                
	double xSh = shower->x;
	double ySh = shower->y;
	double zSh = shower->z;

	// Get vertex position as DVector3
	DVector3 my_vertex(0.0, 0.0, m_zTarget);
	if(vertex)my_vertex = vertex->x.Vect();
	
	// Momentum direction is vector pointing from vertex to shower center
	DVector3 pdir = DVector3(xSh, ySh, zSh) - my_vertex;
	pdir.SetMag(1.0);
	
	double energy = shower->E;

	DBCALPhoton* photon = new DBCALPhoton();

	DVector3 mom = energy * pdir;

	photon->setLorentzMomentum( DLorentzVector( mom, energy ) );
	photon->setShowerPosition( DVector3( xSh, ySh, zSh ) );
  photon->setShowerPositionErr( DVector3( shower->xErr, shower->yErr,
                                          shower->zErr ) );
	photon->setShowerTime( shower->t );
  

	if(shower)photon->AddAssociatedObject(shower);
	if(vertex)photon->AddAssociatedObject(vertex);

	return photon;
}


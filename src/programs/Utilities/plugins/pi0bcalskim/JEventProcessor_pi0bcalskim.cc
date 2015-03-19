// $Id$
//
//    File: JEventProcessor_pi0bcalskim.cc
// Created: Mon Dec  1 14:57:11 EST 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include <math.h>
#include <TLorentzVector.h>
#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include "JEventProcessor_pi0bcalskim.h"


#include "TRACKING/DMCThrown.h"
// Routine used to create our JEventProcessor
#include "PID/DVertex.h"
#include "DANA/DApplication.h"
#include "JANA/JApplication.h"
#include "JANA/JFactory.h"
#include "BCAL/DBCALShower.h"
#include "DLorentzVector.h"
#include "TTree.h"
#include "units.h"
#include "ANALYSIS/DAnalysisUtilities.h"

extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_pi0bcalskim());
  }
} // "C"


//------------------
// JEventProcessor_pi0bcalskim (Constructor)
//------------------
JEventProcessor_pi0bcalskim::JEventProcessor_pi0bcalskim()
{

  MIN_SH1_E = 0.2;
  MIN_SH2_E = 0.2;

  WRITE_EVIO = 1;

  gPARMS->SetDefaultParameter( "PI0BCALSKIM:WRITE_EVIO", WRITE_EVIO );
  gPARMS->SetDefaultParameter("PI0BCALSKIM:MIN_SH1_E" , MIN_SH1_E );
  gPARMS->SetDefaultParameter("PI0BCALSKIM:MIN_SH2_E" , MIN_SH2_E );



}

//------------------
// ~JEventProcessor_pi0bcalskim (Destructor)
//------------------
JEventProcessor_pi0bcalskim::~JEventProcessor_pi0bcalskim()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_pi0bcalskim::init(void)
{
  dEventWriterEVIO = NULL;

  //if( ! WRITE_EVIO) cerr << " output isnt working " << endl;

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_pi0bcalskim::brun(JEventLoop *eventLoop, int runnumber)
{
  eventLoop->GetSingle(dEventWriterEVIO);

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_pi0bcalskim::evnt(JEventLoop *loop, int eventnumber)
{
  vector< const DBCALShower* > locBCALShowers;
  loop->Get(locBCALShowers);
  vector< const DTrackTimeBased*> locTrackTimeBased;
  loop->Get(locTrackTimeBased);
  vector<const DVertex*> kinfitVertex;
  loop->Get(kinfitVertex);

  if(locBCALShowers.size() < 2 ) return NOERROR;

	bool Candidate = false;
	double sh1_E,sh2_E,inv_mass,kinfitVertexZ;
	vector <const DBCALShower *> matchedShowers;
	DVector3 mypos(0.0,0.0,0.0);
	for(unsigned int i = 0 ; i < locTrackTimeBased.size() ; ++i)
	{
		for(unsigned int j = 0 ; j < locBCALShowers.size() ; ++j)
		{
			double x = locBCALShowers[j]->x;
			double y = locBCALShowers[j]->y;
			double z = locBCALShowers[j]->z;
			double E = locBCALShowers[j]->E;
			DVector3 pos_bcal(x,y,z);
			double R = pos_bcal.Perp();
			double phi = pos_bcal.Phi();
			double L2 = 0.81*2.54+65.0;
			double L3 = L2 + 0.81*2.54*2;
			double L4 = L3 + 0.81*2.54*3;
			double L5 = L4 + 0.97*2.54*4;
			locTrackTimeBased[i]->rt->GetIntersectionWithRadius(R, mypos);
			double dPhi = TMath::Abs(mypos.Phi()-pos_bcal.Phi());
			double dZ = TMath::Abs(mypos.Z() - z);	
			if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R) {
				 matchedShowers.push_back(locBCALShowers[j]);
			}
		}
	}


	//japp->RootWriteLock();
	
	double kinfitVertexX, kinfitVertexY, kinfitVertexT;
	for (int i = 0 ; i < kinfitVertex.size(); i++)
	{
		kinfitVertexX = kinfitVertex[i]->dSpacetimeVertex.X();
		kinfitVertexY = kinfitVertex[i]->dSpacetimeVertex.Y();
		kinfitVertexZ = kinfitVertex[i]->dSpacetimeVertex.Z();
		kinfitVertexT = kinfitVertex[i]->dSpacetimeVertex.T();
		//		goodVertexZ->Fill(kinfitVertexZ);
	}
	

  for(unsigned int i=0; i<locBCALShowers.size() ; i++)	
  {
	if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[i]) != matchedShowers.end()) continue;
	 sh1_E = locBCALShowers[i]->E_raw;
	const DBCALShower *s1 = locBCALShowers[i];
	double sh1_x = s1->x - kinfitVertexX ;
	double sh1_y = s1->y - kinfitVertexY ;
	double sh1_z = s1->z - kinfitVertexZ ;
	double sh1_R = sqrt(sh1_x*sh1_x+sh1_y*sh1_y+sh1_z*sh1_z);
	TLorentzVector sh1_p(sh1_E*sh1_x/sh1_R,sh1_E*sh1_y/sh1_R,sh1_E*sh1_z/sh1_R,sh1_E);
	for(unsigned int j = i+1 ; j < locBCALShowers.size() ; j++)
	{
		if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[j]) != matchedShowers.end()) continue;
		const DBCALShower *s2 = locBCALShowers[j];
		 sh2_E = locBCALShowers[j]->E_raw;
		double sh2_x = s2->x - kinfitVertexX;
		double sh2_y = s2->y - kinfitVertexY;
		double sh2_z = s2->z - kinfitVertexZ;
		double sh2_R = sqrt(sh2_x*sh2_x + sh2_y*sh2_y + sh2_z*sh2_z);
		TLorentzVector sh2_p(sh2_E*sh2_x/sh2_R,sh2_E*sh2_y/sh2_R,sh2_E*sh2_z/sh2_R,sh2_E);
		TLorentzVector ptot = sh1_p+sh2_p;
		inv_mass = ptot.M();
		Candidate |= ( (sh2_E>0.67) && (kinfitVertexZ<68.0) &&(kinfitVertexZ>62.0)  && (inv_mass<0.30) );
	}
  }
  	if(Candidate){
	        if( WRITE_EVIO ) {
		  //	cout << " inv mass = " << inv_mass << " sh1 E = " << sh1_E << " sh2 E = " << sh2_E << " event num = " << eventnumber << endl;
     			 dEventWriterEVIO->Write_EVIOEvent( loop, "pi0bcalskim" );
  		  }
			}


  //japp->RootUnLock();
  
   



  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_pi0bcalskim::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// Fin
//------------------
jerror_t JEventProcessor_pi0bcalskim::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}


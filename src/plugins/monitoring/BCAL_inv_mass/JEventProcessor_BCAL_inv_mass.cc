// $Id$
//
//    File: DEventProcessor_BCAL_Shower.cc
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_BCAL_inv_mass.h"

#include <TLorentzVector.h>
#include "TMath.h"

#include "JANA/JApplication.h"
#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALHit.h"
#include "FCAL/DFCALCluster.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "PID/DVertex.h"
//#include "TRACKING/DTrackFinder.h"

#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

// Routine used to create our DEventProcessor


static TH1I* bcal_diphoton_mass_300 = NULL;
static TH1I* bcal_diphoton_mass_500 = NULL;
static TH1I* bcal_diphoton_mass_700 = NULL;
static TH1I* bcal_diphoton_mass_900 = NULL;
static TH1I* bcal_fcal_diphoton_mass_300 = NULL;
static TH1I* bcal_fcal_diphoton_mass_500 = NULL;
static TH1I* bcal_fcal_diphoton_mass_700 = NULL;
static TH1I* bcal_fcal_diphoton_mass_900 = NULL;


extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new JEventProcessor_BCAL_inv_mass()); //register this plugin
	}
} // "C"

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_inv_mass::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:

	if(bcal_diphoton_mass_500 != NULL){
	  return NOERROR;
	}

	TDirectory *main = gDirectory;
	gDirectory->mkdir("bcal_inv_mass")->cd();

        bcal_diphoton_mass_300 = new TH1I("bcal_diphoton_mass_300","bcal diphoton mass (Cluster E > 300 MeV)",100,0.0,1.0);
        bcal_diphoton_mass_300->GetXaxis()->SetTitle("invariant mass [GeV]");
        bcal_diphoton_mass_300->GetYaxis()->SetTitle("counts / 10 MeV");

	bcal_diphoton_mass_500 = new TH1I("bcal_diphoton_mass_500","bcal diphoton mass (Cluster E > 500 MeV)",100,0.0,1.0);
	bcal_diphoton_mass_500->GetXaxis()->SetTitle("invariant mass [GeV]");
	bcal_diphoton_mass_500->GetYaxis()->SetTitle("counts / 10 MeV");	

        bcal_diphoton_mass_700 = new TH1I("bcal_diphoton_mass_700","bcal diphoton mass (Cluster E > 700 MeV)",100,0.0,1.0);
        bcal_diphoton_mass_700->GetXaxis()->SetTitle("invariant mass [GeV]");
        bcal_diphoton_mass_700->GetYaxis()->SetTitle("counts / 10 MeV");

        bcal_diphoton_mass_900 = new TH1I("bcal_diphoton_mass_900","bcal diphoton mass (Cluster E > 900 MeV)",100,0.0,1.0);
        bcal_diphoton_mass_900->GetXaxis()->SetTitle("invariant mass [GeV]");
        bcal_diphoton_mass_900->GetYaxis()->SetTitle("counts / 10 MeV");

        bcal_fcal_diphoton_mass_300 = new TH1I("bcal_fcal_diphoton_mass_300","bcal and fcal diphoton mass (Cluster E > 300 MeV)",100,0.0,1.0);
        bcal_fcal_diphoton_mass_300->GetXaxis()->SetTitle("invariant mass [GeV]");
        bcal_fcal_diphoton_mass_300->GetYaxis()->SetTitle("counts / 10 MeV");
	
        bcal_fcal_diphoton_mass_500 = new TH1I("bcal_fcal_diphoton_mass_500","bcal and fcal diphoton mass (Cluster E > 500 MeV)",100,0.0,1.0);
        bcal_fcal_diphoton_mass_500->GetXaxis()->SetTitle("invariant mass [GeV]");
        bcal_fcal_diphoton_mass_500->GetYaxis()->SetTitle("counts / 10 MeV");

        bcal_fcal_diphoton_mass_700 = new TH1I("bcal_fcal_diphoton_mass_700","bcal and fcal diphoton mass (Cluster E > 700 MeV)",100,0.0,1.0);
        bcal_fcal_diphoton_mass_700->GetXaxis()->SetTitle("invariant mass [GeV]");
        bcal_fcal_diphoton_mass_700->GetYaxis()->SetTitle("counts / 10 MeV");

        bcal_fcal_diphoton_mass_900 = new TH1I("bcal_fcal_diphoton_mass_900","bcal and fcal diphoton mass (Cluster E > 900 MeV)",100,0.0,1.0);
        bcal_fcal_diphoton_mass_900->GetXaxis()->SetTitle("invariant mass [GeV]");
        bcal_fcal_diphoton_mass_900->GetYaxis()->SetTitle("counts / 10 MeV");

	

	//  ... create historgrams or trees ...

	 //	TDirectory *dir = new TDirectoryFile("BCAL","BCAL");
	 //	dir->cd();

	main->cd();

	return NOERROR;
}


//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_inv_mass::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{

	return NOERROR;
}

//------------------
// evnt
//------------------




jerror_t JEventProcessor_BCAL_inv_mass::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
{

	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// locEventLoop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	//
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// locEventLoop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software

	vector<const DBCALShower*> locBCALShowers;
	vector<const DFCALCluster*> locFCALClusters;
	vector<const DVertex*> kinfitVertex;
	//const DDetectorMatches* locDetectorMatches = NULL;
	//locEventLoop->GetSingle(locDetectorMatches);
	locEventLoop->Get(locBCALShowers);
	locEventLoop->Get(locFCALClusters);
	locEventLoop->Get(kinfitVertex);

	vector<const DTrackTimeBased*> locTrackTimeBased;
	locEventLoop->Get(locTrackTimeBased);

	vector <const DBCALShower*> matchedShowers;
	vector <const DFCALCluster*> matchedFCALClusters;
	vector <const DTrackTimeBased*> matchedTracks;
	DVector3 mypos(0.0,0.0,0.0);

	for (unsigned int i=0; i < locTrackTimeBased.size() ; ++i){
	  for (unsigned int j=0; j< locBCALShowers.size(); ++j){
	
	  	double x = locBCALShowers[j]->x;
		double y = locBCALShowers[j]->y;
		double z = locBCALShowers[j]->z;
		DVector3 pos_bcal(x,y,z);
		double R = pos_bcal.Perp();
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(R, mypos);
		 locTrackTimeBased[i]->momentum().Mag();
		double dPhi = TMath::Abs(mypos.Phi()-pos_bcal.Phi());
		double dZ = TMath::Abs(mypos.Z() - z);
		
		if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R) {
		  matchedShowers.push_back(locBCALShowers[j]);
	          matchedTracks.push_back(locTrackTimeBased[i]);

		}

	  }
	}

   for(unsigned int i = 0 ; i < locTrackTimeBased.size(); ++i)
        {
                for(unsigned int j = 0 ; j < locFCALClusters.size(); ++j)
                {
                        const DFCALCluster *c1 = locFCALClusters[j];
                        double x = c1->getCentroid().X();
                        double y = c1->getCentroid().Y();
                        double z = c1->getCentroid().Z();
                        DVector3 fcalpos(x,y,z);
                        //cout << " x = " << x << " y = " << y << endl;
                        DVector3 norm(0.0,0.0,-1);
                        DVector3 pos;
                        locTrackTimeBased[i]->rt->GetIntersectionWithPlane(fcalpos,norm,pos);
                        double diffX = TMath::Abs(x - pos.X());
                        double diffY = TMath::Abs(y - pos.Y());
                        if(diffX < 3.0 && diffY < 3.0)
                        {
                             matchedFCALClusters.push_back(locFCALClusters[j]);
                                                                                                                                                                                                                             }
                  }
	}                        

 	vector <const DChargedTrackHypothesis*> locParticles;
	double kinfitVertexX = 0.0, kinfitVertexY = 0.0, kinfitVertexZ = 0.0;
	for (unsigned int i = 0 ; i < kinfitVertex.size(); i++)
	{
		kinfitVertexX = kinfitVertex[i]->dSpacetimeVertex.X();
		kinfitVertexY = kinfitVertex[i]->dSpacetimeVertex.Y();
		kinfitVertexZ = kinfitVertex[i]->dSpacetimeVertex.Z();
		//p
		//kinfitVertexT = kinfitVertex[i]->dSpacetimeVertex.T();
	}
	
	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	for(unsigned int i=0; i<locBCALShowers.size(); i++)
	{
	     //   if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[i]))
	  // continue;
		if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[i]) != matchedShowers.end()) continue;
		const DBCALShower *s1 = locBCALShowers[i];
		double dx1 = s1->x - kinfitVertexX;
		double dy1 = s1->y - kinfitVertexY;
		double dz1 = s1->z - kinfitVertexZ;
		double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
		double  E1 = s1->E;
		double  E1_raw = s1->E_raw;
		TLorentzVector sh1_p(E1*dx1/R1, E1*dy1/R1, E1*dz1/R1, E1);
		TLorentzVector sh1_p_raw(E1_raw*dx1/R1, E1_raw*dy1/R1, E1_raw*dz1/R1, E1_raw);
			for(unsigned int j=i+1; j<locBCALShowers.size(); j++){
				const DBCALShower *s2 = locBCALShowers[j];
		     		if (find(matchedShowers.begin(), matchedShowers.end(),s2) != matchedShowers.end()) continue;
				double dx2 = s2->x - kinfitVertexX;
				double dy2 = s2->y - kinfitVertexY;
				double dz2 = s2->z - kinfitVertexZ; // shift to coordinate relative to center of target
				double R2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
				double E2 = s2->E;
				double E2_raw = s2->E_raw;
				TLorentzVector sh2_p(E2*dx2/R2, E2*dy2/R2, E2*dz2/R2, E2);	
				TLorentzVector sh2_p_raw(E2_raw*dx2/R2, E2_raw*dy2/R2, E2_raw*dz2/R2, E2_raw);		
				TLorentzVector ptot = sh1_p+sh2_p;
				TLorentzVector ptot_raw = sh1_p_raw + sh2_p_raw ;
				double inv_mass_raw = ptot_raw.M();
				if(E1_raw>.3&&E2_raw>.3) bcal_diphoton_mass_300->Fill(inv_mass_raw);
				if(E1_raw>.5&&E2_raw>.5) bcal_diphoton_mass_500->Fill(inv_mass_raw);
               	                if(E1_raw>.9&&E2_raw>.9) bcal_diphoton_mass_900->Fill(inv_mass_raw);
               		        if(E1_raw>.7&&E2_raw>.7) bcal_diphoton_mass_700->Fill(inv_mass_raw);
			}		
			for(unsigned int j=0; j<locFCALClusters.size(); j++){
				if (find(matchedFCALClusters.begin(), matchedFCALClusters.end(),locFCALClusters[j]) != matchedFCALClusters.end()) continue;
				const DFCALCluster *cl2 = locFCALClusters[j];
				double dx2 = cl2->getCentroid().X()-kinfitVertexX;
	                        double dy2 = cl2->getCentroid().Y()-kinfitVertexY;
                                double dz2 = cl2->getCentroid().Z()-kinfitVertexZ;
				double fcal_E = cl2->getEnergy();
				double R2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
				TLorentzVector cl2_p(fcal_E*dx2/R2, fcal_E*dy2/R2, fcal_E*dz2/R2, fcal_E);
				TLorentzVector ptot_fcal_bcal = sh1_p_raw + cl2_p;
				double inv_mass = ptot_fcal_bcal.M();
                                if(E1_raw>.3&&fcal_E>.3) bcal_fcal_diphoton_mass_300->Fill(inv_mass);
				if(E1_raw>.5&&fcal_E>.5) bcal_fcal_diphoton_mass_500->Fill(inv_mass);
			        if(E1_raw>.7&&fcal_E>.7) bcal_fcal_diphoton_mass_700->Fill(inv_mass);
                                if(E1_raw>.9&&fcal_E>.9) bcal_fcal_diphoton_mass_900->Fill(inv_mass);
			}
	}   


	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK


	/*
	//Optional: Save event to output REST file. Use this to create skims.
	dEventWriterREST->Write_RESTEvent(locEventLoop, "BCAL_Shower"); //string is part of output file name
	*/

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_inv_mass::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_inv_mass::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}


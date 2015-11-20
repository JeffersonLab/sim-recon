// $Id$
//
//    File: DEventProcessor_BCAL_gainmatrix.cc
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "DEventProcessor_BCAL_gainmatrix.h"

#include <TLorentzVector.h>
#include "TMath.h"

#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "PID/DVertex.h"

#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

// Routine used to create our DEventProcessor

extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new DEventProcessor_BCAL_gainmatrix()); //register this plugin
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_BCAL_gainmatrix::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	
	//  ... create historgrams or trees ...

	int n_channels = 768;
	int m_nElements = 768;

	h2D_mC = new TH2F("h2D_mC", "C Matrix as TH2F ", n_channels, 0.0, n_channels, n_channels, 0.0, n_channels);
	h1D_mD = new TH1F("h1D_mD", "D Matrix as TH1F", n_channels, 0.0, n_channels);
	h1D_mL = new TH1F("h1D_mL", " L Matrix as TH1F", n_channels, 0.0, n_channels);
	h1D_massbias = new TH1F( "h1D_massbias", "Mass Bias Value (in bin 2)",5,0.0,5.0);
	h1D_nhits = new TH1F("h1D_nhits", "Number of hits as function of channel number", n_channels,0,n_channels);
	h1D_nhits->SetXTitle("channel number");
	h1D_nhits->SetYTitle("counts");

	m_mD.ResizeTo(m_nElements,1);
	m_mC.ResizeTo(m_nElements,m_nElements);
	m_mL.ResizeTo(m_nElements,1);
	m_nhits.ResizeTo(m_nElements,1);
	m_mD.Zero();
	m_mL.Zero();
	m_mC.Zero();
	m_nhits.Zero();

	m_massbias = 0.0;
	BCAL_Neutrals = new TTree("BCALNeutrals","BCALNeutrals");
	BCAL_Neutrals->Branch("eventnum",&eventnum,"eventnum/i");
//	BCAL_Neutrals->Branch("sh1_E",&E1);
	BCAL_Neutrals->Branch("sh1_E_raw",&E1_raw);
//	BCAL_Neutrals->Branch("sh2_E",&E2);
       	BCAL_Neutrals->Branch("sh2_E_raw",&E2_raw);
//	BCAL_Neutrals->Branch("sh1_t",&t1);
//	BCAL_Neutrals->Branch("sh2_t",&t2);	
//	BCAL_Neutrals->Branch("inv_mass",&inv_mass);
	BCAL_Neutrals->Branch("inv_mass_raw",&inv_mass_raw);
	BCAL_Neutrals->Branch("sh1_z",&z1);
	BCAL_Neutrals->Branch("sh2_z",&z2);
//	BCAL_Neutrals->Branch("sh1_x",&x1);
//	BCAL_Neutrals->Branch("sh2_x",&x2);
//	BCAL_Neutrals->Branch("sh1_y",&y1);
//	BCAL_Neutrals->Branch("sh2_y",&y2);
	BCAL_Neutrals->Branch("psi",&psi);
	BCAL_Neutrals->Branch("cl1_contains_layer4",&logical1);
	BCAL_Neutrals->Branch("cl2_contains_layer4",&logical2);
	BCAL_Neutrals->Branch("vertexZ",&vertexz);
	BCAL_Neutrals->Branch("Run_Number", &Run_Number);
	
	dEventWriterROOT = NULL;
	dEventWriterREST = NULL;

	return NOERROR;
}


//------------------
// brun
//------------------
jerror_t DEventProcessor_BCAL_gainmatrix::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes
	Run_Number = locRunNumber;
	/*
	//Optional: Retrieve REST writer for writing out skims
	locEventLoop->GetSingle(dEventWriterREST);
	
	//Recommeded: Create output ROOT TTrees (nothing is done if already created)
	locEventLoop->GetSingle(dEventWriterROOT);
	dEventWriterROOT->Create_DataTrees(locEventLoop);
	*/

	return NOERROR;
}

//------------------
// evnt
//------------------




jerror_t DEventProcessor_BCAL_gainmatrix::evnt(jana::JEventLoop* locEventLoop, int locEventNumber)
{

	eventnum = locEventNumber;
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
	vector<const DBCALCluster*> locBCALClusters;
	vector<const DBCALPoint*> locBCALPoints;
	vector<const DVertex*> locVertex;
	locEventLoop->Get(locBCALShowers);
	locEventLoop->Get(locBCALClusters);
	locEventLoop->Get(locBCALPoints);
	locEventLoop->Get(locVertex);

	vector<const DTrackTimeBased*> locTrackTimeBased;
	locEventLoop->Get(locTrackTimeBased);

	vector <const DBCALShower *> matchedShowers;
	DVector3 mypos(0.0,0.0,0.0);
	for (unsigned int i=0; i < locTrackTimeBased.size() ; ++i){
		for (unsigned int j=0; j< locBCALShowers.size(); ++j){
	  	double x = locBCALShowers[j]->x;
		double y = locBCALShowers[j]->y;
		double z = locBCALShowers[j]->z;
		DVector3 pos_bcal(x,y,z);
		double R = pos_bcal.Perp();
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(R, mypos);
		
		double dPhi = TMath::Abs(mypos.Phi()-pos_bcal.Phi());
		double dZ = TMath::Abs(mypos.Z() - z);
		
		if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R) matchedShowers.push_back(locBCALShowers[j]);
	 	}
	}

	japp->RootWriteLock();
	double vertexX = 0.0;
	double vertexY = 0.0;
	double vertexZ = 0.0;
	for (unsigned int i = 0 ; i < locVertex.size(); i++)
	{
		vertexX = locVertex[i]->dSpacetimeVertex.X();
		vertexY = locVertex[i]->dSpacetimeVertex.Y();
		vertexZ = locVertex[i]->dSpacetimeVertex.Z();
	}
	
//OoOoOoOoOoOoOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO CALCULATING THE PI0 INVARIANT MASS PORTION O0O0O0O0O0O0O0O0OOooooooooooOooOoOoooooOOooooOOo

	for(unsigned int i=0; i<locBCALShowers.size(); i++){
	  	double pi0_mass = 0.1349766;
		if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[i]) != matchedShowers.end()) continue;
		const DBCALShower *s1 = locBCALShowers[i];
		vector<const DBCALCluster*> associated_clusters1;
		s1->Get(associated_clusters1);
		double dx1 = s1->x - vertexX;
		double dy1 = s1->y - vertexY;
		double dz1 = s1->z - vertexZ;
		double dt1 = s1->t;
		t1 = dt1;
		z1 = s1->z;
		x1 = s1->x;
		y1 = s1->y;
		double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
		E1 = s1->E;
		E1_raw = s1->E_raw;
		TLorentzVector sh1_p(E1*dx1/R1, E1*dy1/R1, E1*dz1/R1, E1);
		TLorentzVector sh1_p_raw(E1_raw*dx1/R1, E1_raw*dy1/R1, E1_raw*dz1/R1, E1_raw);
		vertexz = vertexZ;
			for(unsigned int j=i+1; j<locBCALShowers.size(); j++){
			const DBCALShower *s2 = locBCALShowers[j];
		     	if (find(matchedShowers.begin(), matchedShowers.end(),s2) != matchedShowers.end()) continue;
			vector<const DBCALCluster*> associated_clusters2;
			s2->Get(associated_clusters2);
			double dx2 = s2->x - vertexX;
			double dy2 = s2->y - vertexY;
			double dz2 = s2->z - vertexZ; // shift to coordinate relative to center of target
			double dt2 = s2->t;
			t2=dt2;
			z2 = s2->z;
			x2 = s2->x;
			y2 = s2->y;
			double R2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
			E2 = s2->E;
			E2_raw = s2->E_raw;
			TLorentzVector sh2_p(E2*dx2/R2, E2*dy2/R2, E2*dz2/R2, E2);	
			TLorentzVector sh2_p_raw(E2_raw*dx2/R2, E2_raw*dy2/R2, E2_raw*dz2/R2, E2_raw);	
			double cospsi = (dx1*dx2+dy1*dy2+dz1*dz2)/R1/R2;
			psi=acos(cospsi);
			TLorentzVector ptot = sh1_p+sh2_p;
			TLorentzVector ptot_raw = sh1_p_raw + sh2_p_raw ;
			inv_mass = ptot.M();
			inv_mass_raw = ptot_raw.M();
			point1_energy_calib.clear();
			point2_energy_calib.clear();
			frac_en.clear();
			point1_channel.clear();
			point2_channel.clear();
			point_energy.clear();
			logical1=0;
			logical2=0;
			for(unsigned int loc_j = 0; loc_j < associated_clusters1.size(); loc_j++){
				const DBCALCluster* a_cluster1 = associated_clusters1[loc_j];
				vector<const DBCALPoint*> associated_points1;
				a_cluster1->Get(associated_points1); 
				for(unsigned int loc_jj = 0; loc_jj < associated_clusters2.size(); loc_jj++){
					const DBCALCluster* a_cluster2 = associated_clusters2[loc_jj];
					vector<const DBCALPoint*> associated_points2;
					a_cluster2->Get(associated_points2); 
					for(unsigned int loc_k = 0; loc_k < associated_points1.size(); loc_k++){
						int module1 = associated_points1[loc_k]->module(); 
						int layer1 = associated_points1[loc_k]->layer();
						int sector1 = associated_points1[loc_k]->sector();
						double energy1_US = associated_points1[loc_k]->E_US();
						double energy1_DS = associated_points1[loc_k]->E_DS();
						int channel1 = 16*(module1-1)+4*(layer1-1)+(sector1-1);
						if(layer1==4) logical1=1;
						point1_energy_calib.push_back(sqrt(energy1_US*energy1_DS/E1_raw/E1_raw));
						point1_channel.push_back(channel1);
						for(unsigned int loc_p = 0 ; loc_p < associated_points2.size(); loc_p++){
							int module2 = associated_points2[loc_p]->module(); 
							int layer2 = associated_points2[loc_p]->layer();
							int sector2 = associated_points2[loc_p]->sector();
							double energy2_US = associated_points2[loc_p]->E_US();
							double energy2_DS = associated_points2[loc_p]->E_DS();
							int channel2 = 16*(module2-1)+4*(layer2-1)+(sector2-1);
							if(layer2==4) logical2=1;
							point2_energy_calib.push_back(sqrt(energy2_US*energy2_DS/E2_raw/E2_raw));
							point2_channel.push_back(channel2);
						}
					} 
				}
			}   
			for ( unsigned int f = 0 ; f < point1_energy_calib.size() ; ++f){
				frac_en.push_back(make_pair(point1_energy_calib[f],point1_channel[f]));
			}

			for ( unsigned int h = 0 ; h < point2_energy_calib.size() ; ++h){
				frac_en.push_back(make_pair(point2_energy_calib[h],point2_channel[h]));
			}

			for (unsigned int a = 0 ; a < frac_en.size() ; ++a){
				point_energy.push_back(frac_en[a].first);
				if(inv_mass_raw>0.12&&inv_mass_raw<.145&&E1_raw>.9&&E2_raw>.9&&vertexZ!=65.0&&psi>.07&&psi<.14 &&frac_en[a].first > 0.0 && frac_en[a].first < 30.0)  m_mD[frac_en[a].second][0] += -(inv_mass_raw*inv_mass_raw - pi0_mass*pi0_mass)*(inv_mass_raw*inv_mass_raw)*frac_en[a].first;
				if(inv_mass_raw>0.12&&inv_mass_raw<.145&&E1_raw>.9&&E2_raw>.9&&vertexZ!=65.0&&psi>.07&&psi<.14 &&frac_en[a].first > 0.0 && frac_en[a].first < 30.0 ) m_mL[frac_en[a].second][0] += (inv_mass_raw*inv_mass_raw)*frac_en[a].first;
				if(inv_mass_raw > 0.12 && inv_mass_raw <.145 && E1_raw>0.9&&E2_raw>0.9&&vertexZ!=65.0&&psi>.07&&psi<.14 && frac_en[a].first > 0.0 && frac_en[a].first < 30.0 ) m_nhits[frac_en[a].second] += 1;

				for(unsigned int b = 0 ; b < frac_en.size(); ++b)
				{
			     	if(inv_mass_raw>.12&&inv_mass_raw<.145&&E1_raw>.9&&E2_raw>.9&&vertexZ!=65.0&&psi>.07&&psi<.14 && frac_en[a].first > 0.0 && frac_en[a].first < 30.0&& frac_en[b].first > 0.0 && frac_en[b].first < 30.0 ) m_mC[frac_en[a].second][frac_en[b].second] += (inv_mass_raw*inv_mass_raw*inv_mass_raw*inv_mass_raw)*(frac_en[a].first)*(frac_en[b].first);
				}
		        }
			for(unsigned int h=0;h<point_energy.size();h++){
			if(inv_mass_raw>.12&&inv_mass_raw<.145&&E1_raw>.9&&E2_raw>.9&&vertexZ!=65.0&&psi>.07&&psi<.14&& frac_en[h].first > 0.0 && frac_en[h].first < 30.0 ) m_massbias += (inv_mass_raw*inv_mass_raw - pi0_mass*pi0_mass);
			}	

			BCAL_Neutrals->Fill();
			}
		
	}   

	japp->RootUnLock();

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_BCAL_gainmatrix::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

  japp->RootWriteLock();
  int n_channels = 768;

  h1D_nhits = new TH1F("h1D_nhits", "Number of hits as function of channel number", n_channels,0,767);
  h1D_nhits->SetYTitle("Number of Hits");
  h1D_massbias->SetBinContent(2,m_massbias);

    for(int i=0; i<n_channels; i++){
    h1D_mD->SetBinContent(i+1,m_mD[i][0]);
    h1D_mL->SetBinContent(i+1,m_mL[i][0]);
    h1D_nhits->SetBinContent(i+1,m_nhits[i][0]);


   	 for(int j=0; j<n_channels; j++){
     	 h2D_mC->SetBinContent(i+1,j+1,m_mC[i][j]);
   	 }
    }

    japp->RootUnLock();
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_BCAL_gainmatrix::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}


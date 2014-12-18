// $Id$
//
//    File: DEventProcessor_BCAL_Shower.cc
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "DEventProcessor_BCAL_Shower.h"

#include <TLorentzVector.h>
#include "TMath.h"

#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALHit.h"
#include "FCAL/DFCALShower.h"
#include "TRACKING/DMCThrown.h"
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"

#include <vector>
#include <deque>
#include <string>
#include <iostream>

// Routine used to create our DEventProcessor

extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new DEventProcessor_BCAL_Shower()); //register this plugin
	}
} // "C"

#define FCAL_Z_OFFSET 640.0-65.0
//------------------
// init
//------------------
jerror_t DEventProcessor_BCAL_Shower::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	 japp->RootWriteLock();
	//  ... create historgrams or trees ...

	TDirectory *dir = new TDirectoryFile("BCAL","BCAL");
	dir->cd();
	
	two_gamma_mass = new TH1F("two_gamma_mass","two_gamma_mass",2000, 0.0, 3.00);
	bcal_diphoton_mass = new TH1F("bcal_diphoton_mass","bcal_diphoton_mass",2000,0.0,3.0);
	bcal_diphoton_mass_highE = new TH1F("bcal_diphoton_mass_highE","bcal_diphoton_mass_highE",2000,0.0,3.0);
	//two_gamma_mass_corr = new TH1F("two_gamma_mass_corr","two_gamma_mass_corr",100, 0.0, 1.00);
	//two_gamma_mass_cut = new TH1F("two_gamma_mass_cut","two_gamma_mass_cut",100, 0.0, 0.300);
	bcal_fcal_two_gamma_mass = new TH1F("one_fcal_one_bcal_gamma_mass","bcal_fcal_two_gamma_mass",1000, 0.0, 1.00);
	//bcal_fcal_two_gamma_mass_cut = new TH1F("fcal_two_gamma_mass_cut","bcal_fcal_two_gamma_mass_cut",100, 0.0, 0.300);
	two_fcal_gamma_mass = new TH1F("two_fcal_gamma_mass","two_fcal_gamma_mass",1000,0.0,1.0);
	two_fcal_gamma_mass_test = new TH1F("two_fcal_gamma_mass_test","two_fcal_gamma_mass_test",1000,0.0,1.0);
	mass_summed = new TH1F("inv_mass","inv_mass",100,0.0,1.0);
	//xy_shower = new TH2F("xy_shower","xy_shower",100, -100.0, 100., 100 , -100.0, 100.0);
	z_shower = new TH1F("z_shower","z_shower",1600, 0.0, 6.0);
	E_shower = new TH1F("E_shower","E_shower", 1600, 0.0, 6.0);
	Neutral_E = new TH1F("Neutral_E","Neutral_E",1000,0.0, 400.0);
	VertexZ = new TH1F("VertexZ","VertexZ",1000,-200,300);
	goodVertexZ = new TH1F("goodVertexZ","goodVertexZ",1000,-200,300);
	Theta_Hist_Both_BCAL = new TH1F("Theta_Both_BCAL","Theta_Both_BCAL",300,0.0,150.0);
	Phi_Hist_Both_BCAL = new TH1F("Phi Both BCAL","Phi Both BCAL",720,-180.0,180.0);
	Psi_Hist_Both_BCAL = new TH1F("Psi Both BCAL","Psi Both BCAL",1000,0.0,150.0);
	Theta_Hist_Split_Gammas = new TH1F("Theta Split Gammas","Theta Split Gammas",300,0.0,150.0);
	Phi_Hist_Split_Gammas = new TH1F("Phi Split Gammas","Phi Split Gammas",720,-180.0,180.0);
	Psi_Hist_Split_Gammas = new TH1F("Psi Split Gammas","Psi Split Gammas",1000,0.0,150.0);
	Theta_Hist_Both_FCAL = new TH1F("Theta Both FCAL","Theta Both FCAL",300,0.0,150.0);
	Phi_Hist_Both_FCAL = new TH1F("Phi Both FCAL","Phi Both FCAL",720,-180.0,180.0);
	Psi_Hist_Both_FCAL = new TH1F("Psi Both FCAL","Psi Both FCAL",1000,0.0,150.0);
	E1_Vs_E2CosTheta = new TH2F("E1 Vs E2(1-cosTheta)","E1 Vs E2(1-cosTheta)",2000,0.0,1.5,2000.0,0.0,5.0);
	Thrown_Gamma_Theta = new TH1F("Thrown Gamma Theta","Thrown Gamma Theta",300,0.0,150.0);	
	Thrown_Inv_Mass = new TH1F("Thrown Gamma Inv Mass","Thrown Gamma Inv Mass",1000,0.0,1.50);
	Thrown_Vertex = new TH1F("Thrown Gamma Vertex","Thrown Gamma Vertex",2000,-300.0,400.0);
	Thrown_Vertex_Frequency = new TH1F("Vertex Counting","Vertex Counting",100,0.0,10.0);
	Pi0_velocity = new TH1F("v","v",1000,0.0,2.0);
	Cluster_Energy = new TH1F("Cluster E","Cluster E",1000,0.0,2.0);
	Point_Energy = new TH1F("Point Energy","Point Energy",1000,0.0,2.0);
	Point_Module = new TH1F("Module","Module",48,0,48.0);
	Point_z = new TH1F("Point z","Point z",1000,-100.0,600.0);
	Assoc_E = new TH1F("Assoc E","Assoc E",1000,0,2.0);
	Layer1_Energy_vs_Channel = new TH2F("Layer1 Energy vs channel","Layer1 Energy vs channel",768,0.0,768.0,1000,0.0,0.5);
	Layer2_Energy_vs_Channel = new TH2F("Layer2 Energy vs channel","Layer2 Energy vs channel",768,0.0,768.0,1000,0.0,0.5);
	Layer3_Energy_vs_Channel = new TH2F("Layer3 Energy vs channel","Layer3 Energy vs channel",768,0.0,768.0,1000,0.0,0.5);
	Layer4_Energy_vs_Channel = new TH2F("Layer4 Energy vs channel","Layer4 Energy vs channel",768,0.0,768.0,1000,0.0,0.5);
	Layer1_Energy = new TH1F("Layer1 Energy","Layer1 Energy",1500,0.0,1.5);
	Layer2_Energy = new TH1F("Layer2 Energy","Layer2 Energy",1500,0.0,1.5);
	Layer3_Energy = new TH1F("Layer3 Energy","Layer3 Energy",1500,0.0,1.5);
	Layer4_Energy = new TH1F("Layer4 Energy","Layer4 Energy",1500,0.0,1.5);	
	Layer1_Energy_v2 = new TH1F("Layer1 Energy v2","Layer1 Energy v2",1500,0.0,1.5);
	Layer2_Energy_v2 = new TH1F("Layer2 Energy v2","Layer2 Energy v2",1500,0.0,1.5);
	Layer3_Energy_v2 = new TH1F("Layer3 Energy v2","Layer3 Energy v2",1500,0.0,1.5);
	Layer4_Energy_v2 = new TH1F("Layer4 Energy v2","Layer4 Energy v2",1500,0.0,1.5);	
	Layer1_Energy_pos = new TH1F("Layer1 Energy pos","Layer1 Energy pos",1500,0.0,0.1);
	Layer2_Energy_pos = new TH1F("Layer2 Energy pos","Layer2 Energy pos",1500,0.0,0.1);
	Layer3_Energy_pos = new TH1F("Layer3 Energy pos","Layer3 Energy pos",1500,0.0,0.15);
	Layer4_Energy_pos = new TH1F("Layer4 Energy pos","Layer4 Energy pos",1500,0.0,0.15);	
	Point_E_M1S1L1 = new TH1F("M1S1L1 E","M1S1L1 E",1500,0.0,1.0);
	Point_E_M12S2L2 = new TH1F("M12S2L2 E","M12S2L2 E",1500,0.0,1.0);
	Point_E_M25S3L3 = new TH1F("M25S3L3 E","M25S3L3 E",1500,0.0,1.0);
	Point_E_M37S4L4 = new TH1F("M37S4L4 E","M37S4L4 E",1500,0.0,1.0);
	Time_Diff = new TH1F("Time Diff", "Time Diff",3000,0.0,500);
	Eoverp_plus_nocuts = new TH1F(" E over p plus no cuts "," E over p plus no cuts ", 1000,0.0,5.0);
	Eoverp_minus_nocuts = new TH1F(" E over p minus no cuts ", " E over p minus no cuts ",1000,0.0,5.0);
	Eoverp_plus_cuts = new TH1F(" E over p plus cuts "," E over p plus cuts ", 1000,0.0,5.0);
	Eoverp_minus_cuts = new TH1F(" E over p minus cuts ", " E over p minus cuts ",1000,0.0,5.0);
	Evsp_plus = new TH2F(" E vs p plus ", " E vs p plus", 1000,0.0,10.0, 1000, 0.0, 10.0);
	Evsp_minus = new TH2F(" E vs p minus ", " E vs p minus", 1000,0.0,10.0, 1000, 0.0, 10.0);
	Eoverpvsp_plus = new TH2F(" E over p vs p plus " , " E over p vs p plus " , 1000, 0.0, 10.0, 1000, 0.0, 10.0);
	Eoverpvsp_minus = new TH2F(" E over p vs p minus " , " E over p vs p minus " , 1000,0.0,10.0,1000,0.0,10.0);
	FCAL_Evsp = new TH2F(" FCAL E vs p " , " FCAL E vs p " ,500,0.0,5.0,250,0.0,2.5);
	FCAL_Eoverpvsp = new TH2F( " FCAL Eoverp vs p " , " FCAL Eoverp vs p " , 500, 0.0, 5.0, 120, 0.0, 1.2);
	FCAL_Eoverp_cuts = new TH1F ( " FCAL E over p cuts " , " FCAL E over p cuts " , 1000, 0.0, 1.0);
	FCAL_Eoverp_nocuts = new TH1F( " FCAL E over p no cuts " , " FCAL E over p no cuts " , 1000,0.0,1.0);

	BCALShowerTrack_Energy = new TH1F("Charged energy","Charged energy",1000,0.0,2.0);
	All_Layer_Energy = new TH1F("All Layer Point Energy","All Layer Point Energy",1500,0.0,1.0);
	
	//Erec_over_Ethrown_vs_z = new TH2F("Erec_over_Ethrown_vs_z","Erec_over_Ethrown_vs_z", 200, -50.0, 600.0, 200, 0.0, 2.0);
	//Ecorr_over_Erec_vs_z = new TH2F("Ecorr_over_Erec_vs_z","Ecorr_over_Erec_vs_z", 200, -50.0, 600.0, 200, 0.0, 4.0);
	//Ereconstructed_vs_Ethrown = new TH2F("Ereconstructed_vs_Ethrown","BCAL total reconstructed E to total thrown E", 200, 0.0, 6.0, 200, 0.0, 6.0);

	//Etot_truth = new TH1F("Etot_truth", "Sum of all truth showers (GeV)", 200, 0.0, 6.0);
	Etot_hits = new TH1F("Etot_hits", "Sum of all hits (GeV)", 200, 0.0, 6.0);
	//Etruth_over_Ethrown_vs_z = new TH2F("Etruth_over_Ethrown_vs_z","Etruth_over_Ethrown_vs_z", 200, -50.0, 600.0, 200, 0.0, 2.0);

	//Edeposited_over_Ethrown_vs_z = new TH2F("Edeposited_over_Ethrown_vs_z","Edeposited_over_Ethrown_vs_z", 80, -50.0, 450.0, 200, 0.0, 2.0);

	// Go back up to the parent directory
	dir->cd("../");

	
	BCALPoint_Charged_neg = new TTree("BCALPointChargedNeg"," from DBCALPoint object");
	//BCAL_Energy->Branch("shower_energy",&shower_energy);
	BCALPoint_Charged_neg->Branch("energy_point",&energy_point);
	BCALPoint_Charged_neg->Branch("energy_shower",&energy_shower);
	BCALPoint_Charged_neg->Branch("energy_raw_shower",&energy_raw_shower);
	BCALPoint_Charged_neg->Branch("track_momentum",&track_momentum);
	BCALPoint_Charged_neg->Branch("layer1_energysum",&layer1_energysum);
	BCALPoint_Charged_neg->Branch("layer2_energysum",&layer2_energysum);
	BCALPoint_Charged_neg->Branch("layer3_energysum",&layer3_energysum);
	BCALPoint_Charged_neg->Branch("layer4_energysum",&layer4_energysum);
	BCALPoint_Charged_neg->Branch("module",&module,"module/i");
	BCALPoint_Charged_neg->Branch("sector",&sector,"sector/i");
	BCALPoint_Charged_neg->Branch("layer",&layer,"layer/i");
	BCALPoint_Charged_neg->Branch("eventnum",&eventnum,"eventnum/i");
	BCALPoint_Charged_neg->Branch("theta",&theta);
	BCALPoint_Charged_neg->Branch("channel",&channel,"channel/i");
	//BCALPoint_Charged_neg->Branch("charge",&charge);
	BCALPoint_Charged_neg->Branch("z", &z);
	BCALPoint_Charged_neg->Branch("r", &r);
	BCALPoint_Charged_neg->Branch("phi", &phi);
	BCALPoint_Charged_neg->Branch("L1_pathlength", &L1_pathlength);
	BCALPoint_Charged_neg->Branch("L2_pathlength", &L2_pathlength);
	BCALPoint_Charged_neg->Branch("L3_pathlength", &L3_pathlength);
	BCALPoint_Charged_neg->Branch("L4_pathlength", &L4_pathlength);

	
	BCALPoint_Charged_pos = new TTree("BCALPointChargedPos"," from DBCALPoint object pos");
	//BCAL_Energy->Branch("shower_energy",&shower_energy);
	BCALPoint_Charged_pos->Branch("energy_point",&energy_point);
	BCALPoint_Charged_pos->Branch("energy_shower",&energy_shower);
	BCALPoint_Charged_pos->Branch("energy_raw_shower",&energy_raw_shower);
	BCALPoint_Charged_pos->Branch("track_momentum",&track_momentum);
	BCALPoint_Charged_pos->Branch("layer1_energysum",&layer1_energysum);
	BCALPoint_Charged_pos->Branch("layer2_energysum",&layer2_energysum);
	BCALPoint_Charged_pos->Branch("layer3_energysum",&layer3_energysum);
	BCALPoint_Charged_pos->Branch("layer4_energysum",&layer4_energysum);
	BCALPoint_Charged_pos->Branch("module",&module,"module/i");
	BCALPoint_Charged_pos->Branch("sector",&sector,"sector/i");
	BCALPoint_Charged_pos->Branch("layer",&layer,"layer/i");
	BCALPoint_Charged_pos->Branch("eventnum",&eventnum,"eventnum/i");
	BCALPoint_Charged_pos->Branch("theta",&theta);
	BCALPoint_Charged_pos->Branch("channel",&channel,"channel/i");
	//BCALPoint_Charged_neg->Branch("charge",&charge);
	BCALPoint_Charged_pos->Branch("z", &z);
	BCALPoint_Charged_pos->Branch("r", &r);
	BCALPoint_Charged_pos->Branch("phi", &phi);
	BCALPoint_Charged_pos->Branch("L1_pathlength", &L1_pathlength);
	BCALPoint_Charged_pos->Branch("L2_pathlength", &L2_pathlength);
	BCALPoint_Charged_pos->Branch("L3_pathlength", &L3_pathlength);
	BCALPoint_Charged_pos->Branch("L4_pathlength", &L4_pathlength);

	BCALPoint_Neutral = new TTree("BCAL Point Neutral", " neutral bcal points");
     //  	BCALPoint_Neutral->Branch("energy_point",&energy_point);
	BCALPoint_Neutral->Branch("module",&module,"module/i");
	BCALPoint_Neutral->Branch("sector",&sector,"sector/i");
	BCALPoint_Neutral->Branch("layer",&layer,"layer/i");
	BCALPoint_Neutral->Branch("eventnum",&eventnum,"eventnum/i");
	BCALPoint_Neutral->Branch("theta",&theta);
	BCALPoint_Neutral->Branch("channel",&channel,"channel/i");

	BCAL_Neutrals = new TTree("BCALNeutrals","BCALNeutrals");
	BCAL_Neutrals->Branch("eventnum",&eventnum,"eventnum/i");
	BCAL_Neutrals->Branch("E1",&E1);
	BCAL_Neutrals->Branch("E1_raw",&E1_raw);
	BCAL_Neutrals->Branch("E2",&E2);
	BCAL_Neutrals->Branch("E2_raw",&E2_raw);
	BCAL_Neutrals->Branch("p1_mag",&p1_mag);
	BCAL_Neutrals->Branch("p1_raw_mag",&p1_raw_mag);
	BCAL_Neutrals->Branch("p2_mag",&p2_mag);
	BCAL_Neutrals->Branch("p2_raw_mag",&p2_raw_mag);
	BCAL_Neutrals->Branch("inv_mass",&inv_mass);
	BCAL_Neutrals->Branch("inv_mass_raw",&inv_mass_raw);
	BCAL_Neutrals->Branch("x1",&x1);
	BCAL_Neutrals->Branch("y1",&y1);
	BCAL_Neutrals->Branch("x2",&x2);
	BCAL_Neutrals->Branch("y2",&y2);
	BCAL_Neutrals->Branch("z1",&z1);
	BCAL_Neutrals->Branch("z2",&z2);
	BCAL_Neutrals->Branch("phi1",&phi1);
	BCAL_Neutrals->Branch("phi2",&phi2);
	BCAL_Neutrals->Branch("t1",&t1);
	BCAL_Neutrals->Branch("t2",&t2);
	BCAL_Neutrals->Branch("vertexZ",&vertexZ);
	BCAL_Neutrals->Branch("vertexX",&vertexX);
	BCAL_Neutrals->Branch("vertexY",&vertexY);
	BCAL_Neutrals->Branch("BCALShowers_per_event",&BCALShowers_per_event,"BCALShowers_per_event/i");
	BCAL_Neutrals->Branch("channel",&channel);

	FCAL_Neutrals = new TTree("FCALNeutrals","FCALNeutrals");
	FCAL_Neutrals->Branch("eventnum",&eventnum,"eventnum/i");
	FCAL_Neutrals->Branch("E1",&E1);
	FCAL_Neutrals->Branch("E2",&E2);
	FCAL_Neutrals->Branch("p1_mag",&p1_mag);
	FCAL_Neutrals->Branch("p2_mag",&p2_mag);
	FCAL_Neutrals->Branch("inv_mass",&inv_mass);
	FCAL_Neutrals->Branch("x1",&x1);
	FCAL_Neutrals->Branch("x2",&x2);
	FCAL_Neutrals->Branch("y1",&y1);
	FCAL_Neutrals->Branch("y2",&y2);
	FCAL_Neutrals->Branch("z1",&z1);
	FCAL_Neutrals->Branch("z2",&z2);
	FCAL_Neutrals->Branch("t1",&t1);
	FCAL_Neutrals->Branch("t2",&t2);
	FCAL_Neutrals->Branch("vertexZ",&vertexZ);
	FCAL_Neutrals->Branch("vertexX",&vertexX);
	FCAL_Neutrals->Branch("vertexY",&vertexY);
	FCAL_Neutrals->Branch("FCALShowers_per_event",&FCALShowers_per_event,"FCALShowers_per_event/i");

	FCALClusterNeutrals = new TTree("FCALClusterNeutrals","FCALClusterNeutrals");
	FCALClusterNeutrals->Branch("eventnum",&eventnum,"eventnum/i");
	FCALClusterNeutrals->Branch("E1",&E1);
	FCALClusterNeutrals->Branch("E2",&E2);
	FCALClusterNeutrals->Branch("p1_mag",&p1_mag);
	FCALClusterNeutrals->Branch("p2_mag",&p2_mag);
	FCALClusterNeutrals->Branch("inv_mass",&inv_mass);
	FCALClusterNeutrals->Branch("x1",&x1);
	FCALClusterNeutrals->Branch("x2",&x2);
	FCALClusterNeutrals->Branch("y1",&y1);
	FCALClusterNeutrals->Branch("y2",&y2);
	FCALClusterNeutrals->Branch("z1",&z1);
	FCALClusterNeutrals->Branch("z2",&z2);
	FCALClusterNeutrals->Branch("t1",&t1);
	FCALClusterNeutrals->Branch("t2",&t2);
	FCALClusterNeutrals->Branch("vertexZ",&vertexZ);
	FCALClusterNeutrals->Branch("vertexX",&vertexX);
	FCALClusterNeutrals->Branch("vertexY",&vertexY);
	FCALClusterNeutrals->Branch("FCALClusters_per_event",&FCALShowers_per_event,"FCALClusters_per_event/i");


	Triple_FCAL_Neutrals = new TTree("TripleFCALNeutrals","TripleFCALNeutrals");
	Triple_FCAL_Neutrals->Branch("eventnum",&eventnum,"eventnum/i");
	Triple_FCAL_Neutrals->Branch("E1",&E1);
	Triple_FCAL_Neutrals->Branch("E2",&E2);
	Triple_FCAL_Neutrals->Branch("E3",&E3);
	Triple_FCAL_Neutrals->Branch("p1_mag",&p1_mag);
	Triple_FCAL_Neutrals->Branch("p2_mag",&p2_mag);
	Triple_FCAL_Neutrals->Branch("p3_mag",&p3_mag);
	Triple_FCAL_Neutrals->Branch("inv_mass",&inv_mass);
	Triple_FCAL_Neutrals->Branch("x1",&x1);
	Triple_FCAL_Neutrals->Branch("x2",&x2);
	Triple_FCAL_Neutrals->Branch("x3",&x3);
	Triple_FCAL_Neutrals->Branch("y1",&y1);
	Triple_FCAL_Neutrals->Branch("y2",&y2);
	Triple_FCAL_Neutrals->Branch("y3",&y3);
	Triple_FCAL_Neutrals->Branch("z1",&z1);
	Triple_FCAL_Neutrals->Branch("z2",&z2);
	Triple_FCAL_Neutrals->Branch("z3",&z3);
	Triple_FCAL_Neutrals->Branch("t1",&t1);
	Triple_FCAL_Neutrals->Branch("t2",&t2);
	Triple_FCAL_Neutrals->Branch("t3",&t3);
	Triple_FCAL_Neutrals->Branch("vertexZ",&vertexZ);
	Triple_FCAL_Neutrals->Branch("vertexX",&vertexX);
	Triple_FCAL_Neutrals->Branch("vertexY",&vertexY);
	Triple_FCAL_Neutrals->Branch("FCALShowers_per_event",&FCALShowers_per_event,"FCALShowers_per_event/i");


	Split_Gamma_Neutrals = new TTree ("SplitGammasNeutrals","SplitGammasNeutrals");
	Split_Gamma_Neutrals->Branch("eventnum",&eventnum,"eventnum/i");
	Split_Gamma_Neutrals->Branch("bcal_E",&bcal_E);
	Split_Gamma_Neutrals->Branch("bcal_x",&bcal_x);
	Split_Gamma_Neutrals->Branch("bcal_y",&bcal_y);
	Split_Gamma_Neutrals->Branch("bcal_z",&bcal_z);
	Split_Gamma_Neutrals->Branch("bcal_phi",&bcal_phi);
	Split_Gamma_Neutrals->Branch("bcal_t",&bcal_t);
	Split_Gamma_Neutrals->Branch("fcal_E",&fcal_E);
	Split_Gamma_Neutrals->Branch("fcal_x",&fcal_x);
	Split_Gamma_Neutrals->Branch("fcal_y",&fcal_y);
	Split_Gamma_Neutrals->Branch("fcal_z",&fcal_z);
	Split_Gamma_Neutrals->Branch("fcal_t",&fcal_t);
	Split_Gamma_Neutrals->Branch("bcal_p",&bcal_p);
	Split_Gamma_Neutrals->Branch("fcal_p",&fcal_p);
	Split_Gamma_Neutrals->Branch("inv_mass",&inv_mass);
	Split_Gamma_Neutrals->Branch("BCALShowers_per_event",&BCALShowers_per_event);
	Split_Gamma_Neutrals->Branch("FCALShowers_per_event",&FCALShowers_per_event);
	Split_Gamma_Neutrals->Branch("vertexZ",&vertexZ);
	Split_Gamma_Neutrals->Branch("vertexX",&vertexX);
	Split_Gamma_Neutrals->Branch("vertexY",&vertexY);


	Split_Gamma_Neutrals_raw = new TTree ("SplitGammasNeutralsRaw","SplitGammasNeutralsRaw");
	Split_Gamma_Neutrals_raw->Branch("eventnum",&eventnum,"eventnum/i");
	Split_Gamma_Neutrals_raw->Branch("bcal_E",&bcal_E);
	Split_Gamma_Neutrals_raw->Branch("bcal_x",&bcal_x);
	Split_Gamma_Neutrals_raw->Branch("bcal_y",&bcal_y);
	Split_Gamma_Neutrals_raw->Branch("bcal_z",&bcal_z);
	Split_Gamma_Neutrals_raw->Branch("bcal_phi",&bcal_phi);
	Split_Gamma_Neutrals_raw->Branch("bcal_t",&bcal_t);
	Split_Gamma_Neutrals_raw->Branch("fcal_E",&fcal_E);
	Split_Gamma_Neutrals_raw->Branch("fcal_x",&fcal_x);
	Split_Gamma_Neutrals_raw->Branch("fcal_y",&fcal_y);
	Split_Gamma_Neutrals_raw->Branch("fcal_z",&fcal_z);
	Split_Gamma_Neutrals_raw->Branch("fcal_t",&fcal_t);
	Split_Gamma_Neutrals_raw->Branch("bcal_p",&bcal_p);
	Split_Gamma_Neutrals_raw->Branch("fcal_p",&fcal_p);
	Split_Gamma_Neutrals_raw->Branch("inv_mass",&inv_mass);
	Split_Gamma_Neutrals_raw->Branch("BCALShowers_per_event",&BCALShowers_per_event);
	Split_Gamma_Neutrals_raw->Branch("FCALClusters_per_event",&FCALShowers_per_event);
	Split_Gamma_Neutrals_raw->Branch("vertexZ",&vertexZ);
	Split_Gamma_Neutrals_raw->Branch("vertexX",&vertexX);
	Split_Gamma_Neutrals_raw->Branch("vertexY",&vertexY);
	

	 japp->RootUnLock();
	

	dEventWriterROOT = NULL;
	dEventWriterREST = NULL;

	return NOERROR;
}


//------------------
// brun
//------------------
jerror_t DEventProcessor_BCAL_Shower::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes

	/*
	//Optional: Retrieve REST writer for writing out skims
	locEventLoop->GetSingle(dEventWriterREST);
	*/

	//vector<const DTrackFinder *> finders;
	//locEventLoop->Get(finders);
	//finder = const_cast<DTrackFinder*>(finders[0]);

	/*
	//Recommeded: Create output ROOT TTrees (nothing is done if already created)
	locEventLoop->GetSingle(dEventWriterROOT);
	dEventWriterROOT->Create_DataTrees(locEventLoop);
	*/

	return NOERROR;
}

//------------------
// evnt
//------------------




jerror_t DEventProcessor_BCAL_Shower::evnt(jana::JEventLoop* locEventLoop, int locEventNumber)
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
	vector<const DFCALShower*> locFCALShowers;
	vector<const DBCALTruthShower*> truthshowers;	
	vector<const DMCThrown*> mcthrowns;
	vector<const DBCALHit*> bcalhits;
	vector<const DBCALCluster*> locBCALClusters;
	vector<const DFCALCluster*> locFCALClusters;
	vector<const DBCALPoint*> locBCALPoints;
	vector<const DVertex*> kinfitVertex;
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);
	locEventLoop->Get(locBCALShowers);
	locEventLoop->Get(locFCALShowers);
	locEventLoop->Get(bcalhits);
	locEventLoop->Get(mcthrowns);
	locEventLoop->Get(truthshowers);
	locEventLoop->Get(locBCALClusters);
	locEventLoop->Get(locFCALClusters);
	locEventLoop->Get(locBCALPoints);
	locEventLoop->Get(kinfitVertex);

	vector<const DTrackTimeBased*> locTrackTimeBased;
	locEventLoop->Get(locTrackTimeBased);

	vector <const DBCALShower *> matchedShowers;
	vector < const DBCALShower *> matchedShowersneg;
	vector < const DBCALShower *> matchedShowerspos;
	vector <const DTrackTimeBased* > matchedTracks;
	vector <const DFCALShower *> matchedFCALShowers;
	vector < const DFCALCluster *> matchedFCALClusters;
	DVector3 mypos(0.0,0.0,0.0);
	DVector3 myposL1(0.0,0.0,0.0);
	DVector3 myposL2(0.0,0.0,0.0);
	DVector3 myposL3(0.0,0.0,0.0);
	DVector3 myposL4(0.0,0.0,0.0);
	DVector3 myposL5(0.0,0.0,0.0);
	DVector3 myposL41(0.0,0.0,0.0);
	DVector3 myposL42(0.0,0.0,0.0);
	DVector3 myposL43(0.0,0.0,0.0);
	double p;
	for (unsigned int i=0; i < locTrackTimeBased.size() ; ++i){
	  for (unsigned int j=0; j< locBCALShowers.size(); ++j){
	
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
		double L41 = L4 + 0.97*2.54*4*1/4;
		double L42 = L4 + 0.97*2.54*4*2/4;
		double L43 = L4 + 0.97*2.54*4*3/4;
		double FOM = TMath::Prob(locTrackTimeBased[i]->chisq, locTrackTimeBased[i]->Ndof);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(R, mypos);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(65.0,myposL1);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(L2,myposL2);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(L3,myposL3);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(L4,myposL4);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(L5,myposL5);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(L41,myposL41);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(L42,myposL42);
		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(L43,myposL43);

		double q = locTrackTimeBased[i]->rt->q;
		 p = locTrackTimeBased[i]->momentum().Mag();
		double dPhi = TMath::Abs(mypos.Phi()-pos_bcal.Phi());
		double dZ = TMath::Abs(mypos.Z() - z);
	//	cout << " energy before matching = " << E << endl;
	 	if( mypos.Perp() == R && p > 0.8 && dZ < 30.0 && dPhi < 0.18 && q == 1.0) Eoverp_plus_cuts->Fill(E/p);
		if( mypos.Perp() == R && p > 0.8 && dZ < 30.0 && dPhi < 0.18 && q == -1.0) Eoverp_minus_cuts->Fill(E/p);
		if( mypos.Perp() == R && dZ < 30.0 && dPhi < 0.18 && q ==1.0) Eoverp_plus_nocuts->Fill(E/p);
		if( mypos.Perp() == R && dZ < 30.0 && dPhi < 0.18 && q ==-1.0) Eoverp_minus_nocuts->Fill(E/p);
		if( mypos.Perp() == R && dZ < 30.0 && dPhi < 0.18 && q ==1.0) Evsp_plus->Fill(p,E);
		if( mypos.Perp() == R && dZ < 30.0 && dPhi < 0.18 && q ==-1.0) Evsp_minus->Fill(p,E);
		if( mypos.Perp() == R && dZ < 30.0 && dPhi < 0.18 && q ==1.0) Eoverpvsp_plus->Fill(p,E/p);
		if( mypos.Perp() == R && dZ < 30.0 && dPhi < 0.18 && q ==-1.0) Eoverpvsp_minus->Fill(p,E/p);
		
		if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R) {
		  matchedShowers.push_back(locBCALShowers[j]);
	          matchedTracks.push_back(locTrackTimeBased[i]);
			z_shower->Fill(E);
			//cout << " energy after matching = " << E << " q = " << q << endl;
		}
		if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R && q == -1.0){
		matchedShowersneg.push_back(locBCALShowers[j]);
		}
		if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R && q == 1.0){
		matchedShowerspos.push_back(locBCALShowers[j]);
		}
	  }
	}

	for(unsigned int i = 0 ; i < locTrackTimeBased.size() ; ++i)
	{
		for(unsigned int j = 0 ; j < locFCALShowers.size(); ++j)
		{
			const DFCALShower *s2 = locFCALShowers[j];
			//double dx, dy, dz, R, thetarad1, phirad1,
			double x = s2->getPosition().X();
			double y = s2->getPosition().Y();
			double z = s2->getPosition().Z();
			DVector3 fcalpos(x,y,z);
			//cout << " x = " << x << " y = " << y << endl;
			DVector3 norm(0.0,0.0,-1);
			DVector3 pos;
			locTrackTimeBased[i]->rt->GetIntersectionWithPlane(fcalpos,norm,pos);
			
				double diffX = TMath::Abs(x - pos.X());
				double diffY = TMath::Abs(y - pos.Y());
				if(diffX < 3.0 && diffY < 3.0) 
				{
					matchedFCALShowers.push_back(locFCALShowers[j]);
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


	//What to do later
	// if (matchedShowes.find( <<< cost BCALShower * >>> ) != matchedShower.end()) continue;

	const DKinematicData* locKinematicData;
	//locEventLoop->Get(locKinematicData);


	DVector3 crude_vertex;
	vector <const DChargedTrackHypothesis*> locChargedTrackHypothesis;
	locEventLoop->Get(locChargedTrackHypothesis);

//	vector <const DTrackCandidate *> locTrackCandidate;
//	locEventLoop->Get(locTrackCandidate,"StraightLine");


	       japp->RootWriteLock();
// = static_cast<const DChargedTrackHypothesis* (locKinematicData);
	//DChargedTrackHypothesis * locChargedTrack = const_cast<const DChargedTrackHypothesis *> (locChargedTrackHypothesis);
	//const DChargedTrackHypothesis* locChargedTrack = static_cast<const DChargedTrackHypothesis* (locKinematicData);
	//deque <const DChargedTrackHypothesis *> locChargedTrack = static_cast <const DChargedTrackHypothesis*>(locKinematicData);

/*	vector <const DNeutralParticleHypothesis *> locNeutralParticleHypothesis;
	locEventLoop->Get(locNeutralParticleHypothesis);
	double locTheta;
	locTheta = locNeutralParticleHypothesis->momentum().Theta*180.0/TMath::Pi();
	Theta->Fill(locTheta);*/


 	vector <const DChargedTrackHypothesis*> locParticles;
	locEventLoop->Get(locParticles);
	//vector <const DKinematicData*> locParticles;
	//locEventLoop->Get(locParticles);
  // DVector3 Calc_CrudeVertex(const deque<const DKinematicData*>& locParticles) const
   //{
	DVector3 locVertex(0.0,0.0,dTargetZCenter);
	double locVertexZ, locVertexY, locVertexX;
	//cout << " vertex before = " << locVertexZ << " dTargetZcenter = " << dTargetZCenter << endl;
	//locVertexZ = locVertex.Z();
	if(locParticles.size() == 0) 
		return NOERROR;
	//if(locParticles.size() == 1)
	//	return locParticles[0]->position();

	double locDOCA, locDOCA2, locSmallestDOCA, locTime0;
	double locAverageTime = 0.0;
	DVector3 locTempVertex;
	DVector3 locPOCA;
	DVector3 locDeltaVertex;
	DVector3 locMomentum;

	locSmallestDOCA = 9.9E9;
	if(locParticles.size()>1){
	//cout << " locParticles.size() = " << locParticles.size() << endl;
	for(int j = 0; j < (int(locParticles.size())-1); ++j) 
	{
			//cout << " please be this far lol " << endl;
		for(size_t k= j+1; k < locParticles.size(); ++k)
		{
			locDOCA = dAnalysisUtilities->Calc_DOCAVertex(locParticles[j],locParticles[k], locTempVertex);

				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
				locVertexZ = locVertex.Z();
				locVertexY = locVertex.Y();
				locVertexX = locVertex.X();
			
		}
	}
	


        }

	double kinfitVertexX, kinfitVertexY, kinfitVertexZ, kinfitVertexT;
	for (int i = 0 ; i < kinfitVertex.size(); i++)
	{
		kinfitVertexX = kinfitVertex[i]->dSpacetimeVertex.X();
		kinfitVertexY = kinfitVertex[i]->dSpacetimeVertex.Y();
		kinfitVertexZ = kinfitVertex[i]->dSpacetimeVertex.Z();
		kinfitVertexT = kinfitVertex[i]->dSpacetimeVertex.T();
		goodVertexZ->Fill(kinfitVertexZ);
	}
	

	//DVector3 poca;
	//double d = finder->FindDoca(position1,momentum1,position2,momentum2, &poca);
	

/*	for(int j = 0; j < locParticles.size(); ++j) 
	{
		locDOCA2 = dAnalysisUtilities->Calc_DOCAToVertex(locParticles[j], locVertex, locPOCA);
		locDeltaVertex = locPOCA - locParticles[j]->position();
		locMomentum = locParticles[j]->momentum();
		double locTime = locParticles[j]->time() + locDeltaVertex.Dot(locMomentum)*locParticles[j]->energy()/(29.9792458*locMomentum.Mag2());
		locAverageTime += locTime;
	}
	locTime0 = locAverageTime/(double(locParticles.size())); */
		//cout << " vertex after = " << locVertexZ << endl;
		VertexZ->Fill(locVertexZ);
	//if (FUCKTHISDEQUE.size() > 0){
		//crude_vertex = dAnalysisUtilities->Calc_CrudeVertex(FUCKTHISDEQUE);
		double VertexZ;
		//VertexZ = crude_vertex.Z();
		//cout << " vertex = " << locVertex << endl;
//	}
//	cout << "check " << endl;
	//single photon details
	double Etot_reconstructed = 0.0;
	for(unsigned int i=0; i<locBCALShowers.size(); i++)
	{
               // if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[i]))
              //  continue;
		const DBCALShower *s = locBCALShowers[i];
		//xy_shower->Fill(s->x, s->y);
		//if (locVertexZ > 50 && locVertexZ < 80) 
		//z_shower->Fill(s->E);
		 E_shower->Fill(s->E);
		Etot_reconstructed += s->E;
	}

	// test to reproduce neutral shower energy histos

   

  for(unsigned int loc_i = 0; loc_i < locBCALShowers.size(); loc_i++)
    {
      // if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[loc_i]))
       //   continue;
	if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[loc_i]) != matchedShowers.end()) continue;
	const DBCALShower *s1 = locBCALShowers[loc_i];
	const DBCALShower* a_shower = locBCALShowers[loc_i];
	vector<const DBCALCluster*> associated_clusters;
	a_shower->Get(associated_clusters);
	double E1 = s1->E;
	double z_one = s1->z;
	Neutral_E->Fill(z_one);
	
	double E_shower = a_shower->E;	
	double E_shower_raw = a_shower->E_raw;
//	if(associated_clusters.size() == 0) cout << " assoc cluster is empty " << endl;
	//if(associated_clusters.size() != 0) cout << " assoc cluster is not empty " << endl;

	for(unsigned int loc_j = 0; loc_j < associated_clusters.size(); loc_j++)
	{

		const DBCALCluster* a_cluster = associated_clusters[loc_j];
		vector<const DBCALPoint*> associated_points;
		a_cluster->Get(associated_points); 

		if(associated_points.size() == 0) cout << " assoc points is empty " << endl;

		for(unsigned int loc_k = 0; loc_k < associated_points.size(); loc_k++)
		{
			int module2 = associated_points[loc_k]->module(); 
			int layer = associated_points[loc_k]->layer();
			double point_energy = associated_points[loc_k]->E();
		
			//cout << " associated point vector size = " << associated_points.size() << " module = " << module2 << " layer = " << associated_points[loc_k]->layer() << " sector = " << associated_points[loc_k]->sector() << " shower energy = " << E_shower << " raw shower energy = " << E_shower_raw << " point energy = " << point_energy << " associated cluster size = " << associated_clusters.size() << " shower sizse = " << locBCALShowers.size() << endl;

			//double theta = associated_clusters[loc_j]->theta();

			//double theta = associated_clusters.rho();	
 			//if(layer==1) Layer1_Energy->Fill(point_energy);
			//if(layer==2) Layer2_Energy->Fill(point_energy);
			//if(layer==3) Layer3_Energy->Fill(point_energy);
			//if(layer==4) Layer4_Energy->Fill(point_energy);

			
			//Assoc_E->Fill(point_energy);
		}
 	}
    } 

// investigating DBCALPoint members

	for(unsigned int i=0; i< locBCALPoints.size(); i++)
	{
		const DBCALPoint *p1 = locBCALPoints[i];
		double E = p1->E();
		double r = p1->r();
		double phi = p1->phi();
		double theta = p1->theta();
		double z = p1->z();
		int module = p1->module();
		int layer = p1->layer();
		int sector = p1->sector();
		Point_Energy->Fill(E);
		Point_Module->Fill(module);
		Point_z->Fill(z);
		//cout << " Point vector size = " << locBCALPoints.size() << " module = " << module << endl;
	}
	for(unsigned int i=0; i< locBCALClusters.size(); i++)
	{
		const DBCALCluster *c1 = locBCALClusters[i];
		double E = c1->E();
		Cluster_Energy->Fill(E);
	}


/*


	for(unsigned int i=0; i<locBCALShowers.size(); i++)
	{
	        if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[i]))
                continue;
		const DBCALShower *s1 = locBCALShowers[i];
		vector<const DBCALCluster*> associated_clusters1;
		s1->Get(associated_clusters1);
		double dx1 = s1->x - locVertexX;
		double dy1 = s1->y - locVertexY;
		double dz1 = s1->z - locVertexZ;
		double dt1 = s1->t;
		double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
		//double path1 = sqrt((dx1-locVertexX)*(dx1-locVertexX)+(dy1-locVertexY)*(dy1-locVertexY)+(dz1)*(dz1));
		double E1 = s1->E;
		double ECalc = s1->E*(1.106+(dz1+65.0-208.4)*(dz1+65.0-208.4)*6.851E-6);
		TLorentzVector p1(E1*dx1/R1, E1*dy1/R1, E1*dz1/R1, E1);
		double thetadeg1, thetarad1, phideg1, phirad1;
		thetadeg1 = p1.Theta()*180.0/TMath::Pi();
		thetarad1 = p1.Theta();
		phideg1 = p1.Phi()*180.0/TMath::Pi();
		phirad1 = p1.Phi();	
		TLorentzVector p1Calc(ECalc*dx1/R1, ECalc*dy1/R1, ECalc*dz1/R1, ECalc);
		for(unsigned int loc_j = 0; loc_j < associated_clusters1.size(); loc_j++)
		{
			const DBCALCluster* a_cluster1 = associated_clusters1[loc_j];
			vector<const DBCALPoint*> associated_points1;
			a_cluster1->Get(associated_points1); 
			for(unsigned int loc_k = 0; loc_k < associated_points1.size(); loc_k++)
			{
				int module1 = associated_points1[loc_k]->module(); 
			}
		}
			for(unsigned int j=i+1; j<locBCALShowers.size(); j++){
	                if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[j]))
                        continue;
			const DBCALShower *s2 = locBCALShowers[j];
			vector<const DBCALCluster*> associated_clusters2;
			s2->Get(associated_clusters2);
			double dx2 = s2->x - locVertexX;
			double dy2 = s2->y - locVertexY;
			double dz2 = s2->z - locVertexZ; // shift to coordinate relative to center of target
			double dt2 = s2->t;
			double R2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
			//double path2 = sqrt((dx2-locVertexX)*(dx2-locVertexX)+(dy2-locVertexY)*(dy2-locVertexY)+(dz2)*(dz2));
			double E2 = s2->E;
			double ECalc = s2->E*(1.106+(dz2+65.0-208.4)*(dz2+65.0-208.4)*6.851E-6);
			TLorentzVector p2(E2*dx2/R2, E2*dy2/R2, E2*dz2/R2, E2);		
			TLorentzVector p2Calc(ECalc*dx2/R2, ECalc*dy2/R2, ECalc*dz2/R2, ECalc);	
			double thetadeg2, thetarad2, phideg2, phirad2, cospsi, psi;
			thetarad2 = p2.Theta();
			phirad2 = p2.Phi();
			thetadeg2 = p2.Theta()*180.0/TMath::Pi();
			phideg2 = p2.Phi()*180.0/TMath::Pi();					
			cospsi = sin(thetarad1)*sin(thetarad2)*cos(phirad1-phirad2)+cos(thetarad1)*cos(thetarad2);
			psi=acos(cospsi)*180/TMath::Pi();
			TLorentzVector ptot = p1+p2;

			Theta_Hist_Both_BCAL->Fill(thetadeg1);
			Theta_Hist_Both_BCAL->Fill(thetadeg2);
			Phi_Hist_Both_BCAL->Fill(phideg2);
			Phi_Hist_Both_BCAL->Fill(phideg2);
			Psi_Hist_Both_BCAL->Fill(psi);
			double makes_sense = 0;
			if (R1 > R2 && dt1 > dt2) makes_sense = 1;
			if (R1 < R2 && dt1 < dt2) makes_sense = 1;
			if (R1 < R2 && dt1 > dt2) makes_sense = 0;
			if (R2 < R1 && dt1 > dt1) makes_sense = 0;

				for(unsigned int loc_jj = 0; loc_jj < associated_clusters1.size(); loc_jj++)
				{
					const DBCALCluster* a_cluster2 = associated_clusters2[loc_jj];
					vector<const DBCALPoint*> associated_points2;
					a_cluster2->Get(associated_points2); 

						for(unsigned int loc_kk = 0; loc_kk < associated_points2.size(); loc_kk++)
						{
							//int module2 = associated_points2[loc_kk]->module(); 
		
							if ( locVertexZ > 55.0 && locVertexZ < 75.0 && E1 > 0.2 && E2 > 0.2 && makes_sense==1 && locVertexX < 15.0 && locVertexX > -15.0 && locVertexY > -15.0 && locVertexY < 15.0) two_gamma_mass->Fill(ptot.M());
						}
					}
				
			

		}
	} 

*/

			





// MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM MIN IONIZE REAL DATA MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM



			for(unsigned int loc_i = 0; loc_i < locBCALShowers.size(); loc_i++)
   			 {
				const DBCALShower* locBCALShower = locBCALShowers[loc_i];
				double E_before = locBCALShower->E;
				//cout << " E before = " << E_before << endl;
				if(find(matchedShowers.begin(), matchedShowers.end(), locBCALShower) == matchedShowers.end()) continue;
  				
				double x = locBCALShower->x;
				double y = locBCALShower->y;
				double z = locBCALShower->z;
				double E1 = locBCALShower->E;
				double first_pos = TMath::Sqrt( TMath::Power(myposL1.X(),2)+TMath::Power(myposL1.Y(),2)+TMath::Power(myposL1.Z(),2));
				double second_pos = TMath::Sqrt( TMath::Power(myposL2.X(),2)+TMath::Power(myposL2.Y(),2)+TMath::Power(myposL2.Z(),2));
				double third_pos = TMath::Sqrt( TMath::Power(myposL3.X(),2)+TMath::Power(myposL3.Y(),2)+TMath::Power(myposL3.Z(),2));
				double fourth_pos = TMath::Sqrt( TMath::Power(myposL4.X(),2)+TMath::Power(myposL4.Y(),2)+TMath::Power(myposL4.Z(),2));
				double fifth_pos = TMath::Sqrt( TMath::Power(myposL5.X(),2)+TMath::Power(myposL5.Y(),2)+TMath::Power(myposL5.Z(),2));
				double fourth_pos2_test = TMath::Sqrt( TMath::Power(myposL4.X()-myposL41.X(),2)+TMath::Power(myposL4.Y()-myposL41.Y(),2)+TMath::Power(myposL4.Z()-myposL41.Z(),2));
				double fourth_pos3_test = TMath::Sqrt( TMath::Power(myposL41.X()-myposL42.X(),2)+TMath::Power(myposL41.Y()-myposL42.Y(),2)+TMath::Power(myposL41.Z()-myposL42.Z(),2));
				double fourth_pos4_test = TMath::Sqrt( TMath::Power(myposL42.X()-myposL43.X(),2)+TMath::Power(myposL42.Y()-myposL43.Y(),2)+TMath::Power(myposL42.Z()-myposL43.Z(),2));
				double fourth_pos5_test = TMath::Sqrt( TMath::Power(myposL43.X()-myposL5.X(),2)+TMath::Power(myposL43.Y()-myposL5.Y(),2)+TMath::Power(myposL43.Z()-myposL5.Z(),2));
	
				double less_approx_dist = fourth_pos2_test + fourth_pos3_test + fourth_pos4_test + fourth_pos5_test;
	
				double L1_pathlength = second_pos - first_pos;
				double L2_pathlength = third_pos - second_pos;
				double L3_pathlength = fourth_pos - third_pos;
				double L4_pathlength = fifth_pos - fourth_pos;
				 
				
				//cout << " first pos = " << first_pos << " second_pos = " << second_pos << " third pos = " << third_pos << " fourth pos = " << fourth_pos << " fifth pos = " << fifth_pos << " event num = " << eventnum << endl;
				//cout << " E after = " << E1 << endl;
				double E_shower = locBCALShower->E;	
				double E_shower_raw = locBCALShower->E_raw;
				double x1 = locBCALShower->x - kinfitVertexX;
				double y1 = locBCALShower->y - kinfitVertexY;
				double z1 = locBCALShower->z - kinfitVertexZ;
				double t1 = locBCALShower->t;
				double R1 = sqrt(x1*x1 + y1*y1 + z1*z1);
				//double path1 = sqrt((dx1-kinfitVertexX)*(dx1-kinfitVertexX)+(dy1-kinfitVertexY)*(dy1-kinfitVertexY)+(dz1)*(dz1));
				TLorentzVector p1(E1*x1/R1, E1*y1/R1, E1*z1/R1, E1);
				vector<const DBCALCluster*> associated_clusters;
				locBCALShower->Get(associated_clusters);

				if(associated_clusters.size() == 0) cout << " assoc cluster is empty " << endl;

				for(unsigned int loc_j = 0; loc_j < associated_clusters.size(); loc_j++)
				{

					const DBCALCluster* a_cluster = associated_clusters[loc_j];
					vector<const DBCALPoint*> associated_points;
					a_cluster->Get(associated_points); 

					if(associated_points.size() == 0) cout << " assoc points is empty " << endl;
						double Layer1_Energy_Sum=0.0;
						double Layer2_Energy_Sum=0.0;
						double Layer3_Energy_Sum=0.0;
						double Layer4_Energy_Sum=0.0;

						double Layer1_Energy_Sum2=0.0;
						double Layer2_Energy_Sum2=0.0;
						double Layer3_Energy_Sum2=0.0;
						double Layer4_Energy_Sum2=0.0;

						double Layer1_Energy_Sumpos = 0.0;
						double Layer2_Energy_Sumpos = 0.0;
						double Layer3_Energy_Sumpos = 0.0;
						double Layer4_Energy_Sumpos = 0.0;
					for(unsigned int loc_k = 0; loc_k < associated_points.size(); loc_k++)
					{
						const DBCALPoint* a_point = associated_points[loc_k];
						vector<const DBCALUnifiedHit*> associated_unifiedhits;
						a_point->Get(associated_unifiedhits);
						int point_size = associated_points.size();
						int module = associated_points[loc_k]->module(); 
						int layer = associated_points[loc_k]->layer();
						int sector = associated_points[loc_k]->sector();
						double theta = associated_points[loc_k]->theta();
						double sintheta = sin(theta);
						double point_energy = associated_points[loc_k]->E();
						int channel_per_module;
					//	int charge = locChargedTrackHypothesis->charge();
						double z = associated_points[loc_k]->z();
						double r = associated_points[loc_k]->r();
						double theta_wrt_vertex = atan(r/(z-kinfitVertexZ));
						double point_energy_sum;
						double point_energy_over_angle = TMath::Abs(point_energy/sin(theta_wrt_vertex));
						point_energy_sum += point_energy_over_angle;
						int logical = 1;
					//	if( layer == 1 && point_energy/sintheta/point_energy_sum > 0.13) logical = 0;
				//		if( layer == 2 && point_energy/sintheta/point_energy_sum > 2*1.5/associated_points.size()) logical = 0;
				//		if( layer == 3 && point_energy/sintheta/point_energy_sum > 3*1.5/associated_points.size()) logical = 0;
				//		if( layer == 4 && point_energy/sintheta/point_energy_sum > 4*1.5/associated_points.size()) logical = 0;
						//if(E_shower>0.3) logical = 0;
						if(associated_points.size()<4) logical = 0;
						//if(layer==1 && point_energy/sintheta >0.4) logical=0;
						//if(layer==2 && (point_energy/sintheta < 0.2 || point_energy/sintheta > 0.6)) logical=0;
						//if(layer==3 && (point_energy/sintheta < 0.4 || point_energy/sintheta > 0.8)) logical=0;
						//if(layer==4 && (point_energy/sintheta < 0.6 || point_energy/sintheta > 1.0)) logical=0;
						if (layer == 1 && sector == 1) channel_per_module = 0;
						if (layer == 1 && sector == 2) channel_per_module = 1;
						if (layer == 1 && sector == 3) channel_per_module = 2;	
						if (layer == 1 && sector == 4) channel_per_module = 3;	
						if (layer == 2 && sector == 1) channel_per_module = 4;
						if (layer == 2 && sector == 2) channel_per_module = 5;	
						if (layer == 2 && sector == 3) channel_per_module = 6;
						if (layer == 2 && sector == 4) channel_per_module = 7;
						if (layer == 3 && sector == 1) channel_per_module = 8;
						if (layer == 3 && sector == 2) channel_per_module = 9;	
						if (layer == 3 && sector == 3) channel_per_module = 10;	
						if (layer == 3 && sector == 4) channel_per_module = 11;	
						if (layer == 4 && sector == 1) channel_per_module = 12;
						if (layer == 4 && sector == 2) channel_per_module = 13;
						if (layer == 4 && sector == 3) channel_per_module = 14;
						if (layer == 4 && sector == 4) channel_per_module = 15;
						int channel = channel_per_module + (module-1)*16;

						//cout << " point energy = " << point_energy << " E shower = " << E1 << " logical = " << logical << " sine = " << sin(theta_wrt_vertex) << " points size = " << associated_points.size() << endl;
						if(associated_points.size()>3)
						{

						if(layer==1)Layer1_Energy_Sum += point_energy;
						if(layer==2)Layer2_Energy_Sum += point_energy;
						if(layer==3)Layer3_Energy_Sum += point_energy;
						if(layer==4)Layer4_Energy_Sum += point_energy;
						//cout << " Layer 1 Energy Sum = " << Layer1_Energy_Sum << " layer 2 e sum = " << Layer2_Energy_Sum << " Layer3_Energy_Sum = " << Layer3_Energy_Sum << " Layer 4 e sum = " << Layer4_Energy_Sum << " point size = " << associated_points.size() << " event num = " << eventnum << " iteration = " << loc_k << endl;

						if(layer==1)Layer1_Energy_Sum2 += point_energy/sintheta;
						if(layer==2)Layer2_Energy_Sum2 += point_energy/sintheta;
						if(layer==3)Layer3_Energy_Sum2 += point_energy/sintheta;
						if(layer==4)Layer4_Energy_Sum2 += point_energy/sintheta;

						if(layer==1)Layer1_Energy_Sumpos += point_energy/L1_pathlength;
						if(layer==2)Layer2_Energy_Sumpos += point_energy/L2_pathlength;
						if(layer==3)Layer3_Energy_Sumpos += point_energy/L3_pathlength;
						if(layer==4)Layer4_Energy_Sumpos += point_energy/L4_pathlength;
						All_Layer_Energy->Fill(point_energy/sin(theta_wrt_vertex));
						
								       if(layer==1) Layer1_Energy_vs_Channel->Fill(channel, point_energy/sin(theta_wrt_vertex));
								       if(layer==2) Layer2_Energy_vs_Channel->Fill(channel, point_energy/sin(theta_wrt_vertex));
								       if(layer==3) Layer3_Energy_vs_Channel->Fill(channel, point_energy/sin(theta_wrt_vertex));
								       if(layer==4) Layer4_Energy_vs_Channel->Fill(channel, point_energy/sin(theta_wrt_vertex));
						if(loc_k+1 == associated_points.size() ){
						  if(Layer1_Energy_Sum < Layer4_Energy_Sum){
							  if(Layer1_Energy_Sum > 0.001 && Layer2_Energy_Sum > 0.001 && Layer3_Energy_Sum > 0.001 && Layer4_Energy_Sum > 0.001) Layer1_Energy_v2->Fill(Layer1_Energy_Sum);
							  if(Layer1_Energy_Sum > 0.001 && Layer2_Energy_Sum > 0.001 && Layer3_Energy_Sum > 0.001 && Layer4_Energy_Sum > 0.001) Layer2_Energy_v2->Fill(Layer2_Energy_Sum);
							  if(Layer1_Energy_Sum > 0.001 && Layer2_Energy_Sum > 0.001 && Layer3_Energy_Sum > 0.001 && Layer4_Energy_Sum > 0.001) Layer3_Energy_v2->Fill(Layer3_Energy_Sum);
							  if(Layer1_Energy_Sum > 0.001 && Layer2_Energy_Sum > 0.001 && Layer3_Energy_Sum > 0.001 && Layer4_Energy_Sum > 0.001) Layer4_Energy_v2->Fill(Layer4_Energy_Sum);
							}
						  if(Layer1_Energy_Sum2 < Layer4_Energy_Sum2 && Layer1_Energy_Sum2 > 0.005){
						    		        Layer1_Energy->Fill(Layer1_Energy_Sum2);
								        Layer2_Energy->Fill(Layer2_Energy_Sum2);
								        Layer3_Energy->Fill(Layer3_Energy_Sum2);
								        Layer4_Energy->Fill(Layer4_Energy_Sum2);
						}
						if(Layer1_Energy_Sumpos < Layer4_Energy_Sumpos && Layer1_Energy_Sumpos > 0.005 && Layer4_Energy_Sumpos > 0.005){
									Layer1_Energy_pos->Fill(Layer1_Energy_Sumpos);
									Layer2_Energy_pos->Fill(Layer2_Energy_Sumpos);
									Layer3_Energy_pos->Fill(Layer3_Energy_Sumpos);
									Layer4_Energy_pos->Fill(Layer4_Energy_Sumpos);
								//cout << " approx dist = " << L4_pathlength << " less approx dist = " << less_approx_dist << " ratio = " << less_approx_dist/L4_pathlength << endl;
						}
						//cout << " layer 1 e sum = " << Layer1_Energy_Sum2 << " layer 2 e sum = " << Layer2_Energy_Sum2 << " layer 3 e sum = " << Layer3_Energy_Sum2 << " layer 4 e sum = " << Layer4_Energy_Sum2 << " layer 1 path = " << L1_pathlength << " layer 2 path = " << L2_pathlength << " layer 3 path = " << L3_pathlength << " layer 4 path = " << L4_pathlength << " event num = " << eventnum << endl;
						}
						}
						//cout << " charge = " << locChargedTrackHypothesis->charge() << " point energy/sintheta= " << point_energy/sintheta1 << " 'shower energy' = " << E1 << " shower momentum = " << p1.M() << " channel = " << channel << " module = " << module1 << " layer = " << layer1 << " point size = " << associated_points.size() << endl;
						for(unsigned int loc_m = 0; loc_m < associated_unifiedhits.size(); loc_m++)
						{
							int modulehit = associated_unifiedhits[loc_m]->module;
							int layerhit = associated_unifiedhits[loc_m]->layer;
							int sectorhit = associated_unifiedhits[loc_m]->sector;
							int end = associated_unifiedhits[loc_m]->end;
							double unifiedhit_energy = associated_unifiedhits[loc_m]->E;
							//cout << " point E = " << point_energy << " hit E = " << unifiedhit_energy << " module hit = " << modulehit << " module point = " << module1 << " layer hit = " << layerhit << " layer point = " << layer1 << " end = " << end << " point size = " << associated_points.size() << " hit size = " << associated_unifiedhits.size() << endl;


						}

	
						
					}
				}
			//}

		}
		//}
	//}


/*
	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	for(size_t i = 0 ; i < locChargedTracks.size(); ++i)
	{
	
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[i]->Get_BestFOM();
		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		const DShowerMatchParams& locBCALShowerMatchParams = locChargedTrackHypothesis->dBCALShowerMatchParams;
		 
		if(locBCALShowerMatchParams.dTrackTimeBased != NULL)
		{
			const DBCALShower* locBCALShower = dynamic_cast <const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
			if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShower)){ 
			BCALShowerTrack_Energy->Fill(locBCALShower->E);
			
			//for(unsigned int loc_i = 0; loc_i < locBCALShowers.size(); loc_i++)
   			 //{
  				//const DBCALShower* s1 = locBCALShowers;
				//const DBCALShower* a_shower = locBCALShowers;
	
				double E1 = locBCALShower->E;
				double E_shower = locBCALShower->E;	
				double E_shower_raw = locBCALShower->E_raw;
				double x1 = locBCALShower->x - kinfitVertexX;
				double y1 = locBCALShower->y - kinfitVertexY;
				double z1 = locBCALShower->z - kinfitVertexZ;
				double t1 = locBCALShower->t;
				double R1 = sqrt(x1*x1 + y1*y1 + z1*z1);
				//double path1 = sqrt((dx1-kinfitVertexX)*(dx1-kinfitVertexX)+(dy1-kinfitVertexY)*(dy1-kinfitVertexY)+(dz1)*(dz1));
				TLorentzVector p1(E1*x1/R1, E1*y1/R1, E1*z1/R1, E1);
				vector<const DBCALCluster*> associated_clusters;
				locBCALShower->Get(associated_clusters);

				if(associated_clusters.size() == 0) cout << " assoc cluster is empty " << endl;

				for(unsigned int loc_j = 0; loc_j < associated_clusters.size(); loc_j++)
				{

					const DBCALCluster* a_cluster = associated_clusters[loc_j];
					vector<const DBCALPoint*> associated_points;
					a_cluster->Get(associated_points); 

					if(associated_points.size() == 0) cout << " assoc points is empty " << endl;

					for(unsigned int loc_k = 0; loc_k < associated_points.size(); loc_k++)
					{
						const DBCALPoint* a_point = associated_points[loc_k];
						vector<const DBCALUnifiedHit*> associated_unifiedhits;
						a_point->Get(associated_unifiedhits);
						int point_size = associated_points.size();
						int module = associated_points[loc_k]->module(); 
						int layer = associated_points[loc_k]->layer();
						int sector = associated_points[loc_k]->sector();
						double theta = associated_points[loc_k]->theta();
						double sintheta = sin(theta);
						double point_energy = associated_points[loc_k]->E();
						int channel_per_module;
						int charge = locChargedTrackHypothesis->charge();
						double z = associated_points[loc_k]->z();
						double r = associated_points[loc_k]->r();
						double theta_wrt_vertex = atan(r/(z-kinfitVertexZ));
						double point_energy_sum;
						point_energy_sum += point_energy/sin(theta_wrt_vertex);
						int logical = 1;
					//	if( layer == 1 && point_energy/sintheta/point_energy_sum > 0.13) logical = 0;
				//		if( layer == 2 && point_energy/sintheta/point_energy_sum > 2*1.5/associated_points.size()) logical = 0;
				//		if( layer == 3 && point_energy/sintheta/point_energy_sum > 3*1.5/associated_points.size()) logical = 0;
				//		if( layer == 4 && point_energy/sintheta/point_energy_sum > 4*1.5/associated_points.size()) logical = 0;
						if(E_shower>0.3) logical = 0;
						if(associated_points.size()<4) logical = 0;
						//if(layer==1 && point_energy/sintheta >0.4) logical=0;
						//if(layer==2 && (point_energy/sintheta < 0.2 || point_energy/sintheta > 0.6)) logical=0;
						//if(layer==3 && (point_energy/sintheta < 0.4 || point_energy/sintheta > 0.8)) logical=0;
						//if(layer==4 && (point_energy/sintheta < 0.6 || point_energy/sintheta > 1.0)) logical=0;
						if (layer == 1 && sector == 1) channel_per_module = 0;
						if (layer == 1 && sector == 2) channel_per_module = 1;
						if (layer == 1 && sector == 3) channel_per_module = 2;	
						if (layer == 1 && sector == 4) channel_per_module = 3;	
						if (layer == 2 && sector == 1) channel_per_module = 4;
						if (layer == 2 && sector == 2) channel_per_module = 5;	
						if (layer == 2 && sector == 3) channel_per_module = 6;
						if (layer == 2 && sector == 4) channel_per_module = 7;
						if (layer == 3 && sector == 1) channel_per_module = 8;
						if (layer == 3 && sector == 2) channel_per_module = 9;	
						if (layer == 3 && sector == 3) channel_per_module = 10;	
						if (layer == 3 && sector == 4) channel_per_module = 11;	
						if (layer == 4 && sector == 1) channel_per_module = 12;
						if (layer == 4 && sector == 2) channel_per_module = 13;
						if (layer == 4 && sector == 3) channel_per_module = 14;
						if (layer == 4 && sector == 4) channel_per_module = 15;
						int channel = channel_per_module + (module-1)*16;
						if(charge < 0 && logical == 1)
						{
						double Layer1_Energy_Sum;
						double Layer2_Energy_Sum;
						double Layer3_Energy_Sum;
						double Layer4_Energy_Sum;
						if(layer==1)Layer1_Energy_Sum += point_energy/sin(theta_wrt_vertex);
						if(layer==2)Layer2_Energy_Sum += point_energy/sin(theta_wrt_vertex);
						if(layer==3)Layer3_Energy_Sum += point_energy/sin(theta_wrt_vertex);
						if(layer==4)Layer4_Energy_Sum += point_energy/sin(theta_wrt_vertex);
						double Layer1_Energy_Sum2;
						double Layer2_Energy_Sum2;
						double Layer3_Energy_Sum2;
						double Layer4_Energy_Sum2;
						if(layer==1)Layer1_Energy_Sum2 += point_energy/sintheta;
						if(layer==2)Layer2_Energy_Sum2 += point_energy/sintheta;
						if(layer==3)Layer3_Energy_Sum2 += point_energy/sintheta;
						if(layer==4)Layer4_Energy_Sum2 += point_energy/sintheta;
						All_Layer_Energy->Fill(point_energy/sin(theta_wrt_vertex));
						
								       if(layer==1) Layer1_Energy_vs_Channel->Fill(channel, point_energy/sin(theta_wrt_vertex));
								       if(layer==2) Layer2_Energy_vs_Channel->Fill(channel, point_energy/sin(theta_wrt_vertex));
								       if(layer==3) Layer3_Energy_vs_Channel->Fill(channel, point_energy/sin(theta_wrt_vertex));
								       if(layer==4) Layer4_Energy_vs_Channel->Fill(channel, point_energy/sin(theta_wrt_vertex));
						if(loc_k+1 == associated_points.size() ){
						  if(Layer1_Energy_Sum < Layer4_Energy_Sum){
							  if(Layer1_Energy_Sum > 0.001) Layer1_Energy_v2->Fill(Layer1_Energy_Sum);
							  if(Layer2_Energy_Sum > 0.001) Layer2_Energy_v2->Fill(Layer2_Energy_Sum);
							  if(Layer3_Energy_Sum > 0.001) Layer3_Energy_v2->Fill(Layer3_Energy_Sum);
							  if(Layer4_Energy_Sum > 0.001) Layer4_Energy_v2->Fill(Layer4_Energy_Sum);
							}
						  if(Layer1_Energy_Sum2 < Layer4_Energy_Sum2){
						    		        Layer1_Energy->Fill(Layer1_Energy_Sum2);
								        Layer2_Energy->Fill(Layer2_Energy_Sum2);
								        Layer3_Energy->Fill(Layer3_Energy_Sum2);
								        Layer4_Energy->Fill(Layer4_Energy_Sum2);
						}
						}
						cout << " point energy/sintheta= " << point_energy/sintheta << " layer 1 e sum = " << Layer1_Energy_Sum2 << " layer 2 e sum = " << Layer2_Energy_Sum2 << " layer 3 e sum = " << Layer3_Energy_Sum2 << " layer 4 e sum = " << Layer4_Energy_Sum2 << "evenr num = " << eventnum << endl;
						}
						//cout << " charge = " << locChargedTrackHypothesis->charge() << " point energy/sintheta= " << point_energy/sintheta1 << " 'shower energy' = " << E1 << " shower momentum = " << p1.M() << " channel = " << channel << " module = " << module1 << " layer = " << layer1 << " point size = " << associated_points.size() << endl;
						for(unsigned int loc_m = 0; loc_m < associated_unifiedhits.size(); loc_m++)
						{
							int modulehit = associated_unifiedhits[loc_m]->module;
							int layerhit = associated_unifiedhits[loc_m]->layer;
							int sectorhit = associated_unifiedhits[loc_m]->sector;
							int end = associated_unifiedhits[loc_m]->end;
							double unifiedhit_energy = associated_unifiedhits[loc_m]->E;
							//cout << " point E = " << point_energy << " hit E = " << unifiedhit_energy << " module hit = " << modulehit << " module point = " << module1 << " layer hit = " << layerhit << " layer point = " << layer1 << " end = " << end << " point size = " << associated_points.size() << " hit size = " << associated_unifiedhits.size() << endl;
							if(locChargedTrackHypothesis->charge() < 0.0 && module == 1 && sector == 1 && layer == 1  && logical ==1) Point_E_M1S1L1->Fill(point_energy/sintheta);
							if(locChargedTrackHypothesis->charge() < 0.0 && module == 12 && sector == 2 && layer == 2 && logical ==1) Point_E_M12S2L2->Fill(point_energy/sintheta);
							if(locChargedTrackHypothesis->charge() < 0.0 && module == 25 && sector == 3 && layer == 3 && logical == 1) Point_E_M25S3L3->Fill(point_energy/sintheta);
							if(locChargedTrackHypothesis->charge() < 0.0 && module == 37 && sector == 4 && layer == 4 && logical == 1) Point_E_M37S4L4->Fill(point_energy/sintheta);

						}

	
						
					}
				}
			//}

		}
		}
	}


*/






/*	for(size_t i = 0 ; i < locChargedTracks.size(); ++i)
	{
	
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[i]->Get_BestFOM();
		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		const DShowerMatchParams& locBCALShowerMatchParams = locChargedTrackHypothesis->dBCALShowerMatchParams;
		 
		if(locBCALShowerMatchParams.dTrackTimeBased != NULL)
		{
			const DBCALShower* locBCALShower = dynamic_cast <const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
			if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShower)){
			BCALShowerTrack_Energy->Fill(locBCALShower->E);
			
			//for(unsigned int loc_i = 0; loc_i < locBCALShowers.size(); loc_i++)
   			 //{
  				//const DBCALShower* s1 = locBCALShowers;
				//const DBCALShower* a_shower = locBCALShowers;
	
				double E1 = locBCALShower->E;
				double E_shower = locBCALShower->E;	
				double E_shower_raw = locBCALShower->E_raw;
				double x1 = locBCALShower->x - kinfitVertexX;
				double y1 = locBCALShower->y - kinfitVertexY;
				double z1 = locBCALShower->z - kinfitVertexZ;
				double t1 = locBCALShower->t;
				double R1 = sqrt(x1*x1 + y1*y1 + z1*z1);
				//double path1 = sqrt((dx1-kinfitVertexX)*(dx1-kinfitVertexX)+(dy1-kinfitVertexY)*(dy1-kinfitVertexY)+(dz1)*(dz1));
				TLorentzVector p1(E1*x1/R1, E1*y1/R1, E1*z1/R1, E1);
				vector<const DBCALCluster*> associated_clusters;
				locBCALShower->Get(associated_clusters);

				if(associated_clusters.size() == 0) cout << " assoc cluster is empty " << endl;

				for(unsigned int loc_j = 0; loc_j < associated_clusters.size(); loc_j++)
				{

					const DBCALCluster* a_cluster = associated_clusters[loc_j];
					vector<const DBCALPoint*> associated_points;
					a_cluster->Get(associated_points); 

					if(associated_points.size() == 0) cout << " assoc points is empty " << endl;

					for(unsigned int loc_k = 0; loc_k < associated_points.size(); loc_k++)
					{
						const DBCALPoint* a_point = associated_points[loc_k];
						vector<const DBCALUnifiedHit*> associated_unifiedhits;
						a_point->Get(associated_unifiedhits);
						int module = associated_points[loc_k]->module(); 
						int layer = associated_points[loc_k]->layer();
						int sector = associated_points[loc_k]->sector();
						double theta = associated_points[loc_k]->theta();
						double sintheta = sin(theta);
						double point_energy = associated_points[loc_k]->E();
						int channel_per_module;
						int charge = locChargedTrackHypothesis->charge();
						float point_energy_sum;
						point_energy_sum += point_energy;
						int logical = 1;
						if( layer == 1 && point_energy/sintheta/point_energy_sum > 1.5/associated_points.size()) logical = 0;
						if( layer == 2 && point_energy/sintheta/point_energy_sum > 2*1.5/associated_points.size()) logical = 0;
						if( layer == 3 && point_energy/sintheta/point_energy_sum > 3*1.5/associated_points.size()) logical = 0;
						if( layer == 4 && point_energy/sintheta/point_energy_sum > 4*1.5/associated_points.size()) logical = 0;
						if(E_shower>0.3) logical = 0;
						if(associated_points.size()<4) logical = 0;
						if(layer==1 && point_energy/sintheta >0.4) logical=0;
						if(layer==2 && point_energy/sintheta < 0.2) logical=0;
						if(layer==2 && point_energy/sintheta > 0.6) logical=0;
						if(layer==3 && point_energy/sintheta < 0.4) logical=0;
						if(layer==3 && point_energy/sintheta > 0.8) logical=0;
						if(layer==4 && point_energy/sintheta < 0.6) logical=0;
						if(layer==4 && point_energy/sintheta > 1.0) logical=0; 
						if (layer == 1 && sector == 1) channel_per_module = 0;
						if (layer == 1 && sector == 2) channel_per_module = 1;
						if (layer == 1 && sector == 3) channel_per_module = 2;	
						if (layer == 1 && sector == 4) channel_per_module = 3;	
						if (layer == 2 && sector == 1) channel_per_module = 4;
						if (layer == 2 && sector == 2) channel_per_module = 5;	
						if (layer == 2 && sector == 3) channel_per_module = 6;
						if (layer == 2 && sector == 4) channel_per_module = 7;
						if (layer == 3 && sector == 1) channel_per_module = 8;
						if (layer == 3 && sector == 2) channel_per_module = 9;	
						if (layer == 3 && sector == 3) channel_per_module = 10;	
						if (layer == 3 && sector == 4) channel_per_module = 11;	
						if (layer == 4 && sector == 1) channel_per_module = 12;
						if (layer == 4 && sector == 2) channel_per_module = 13;
						if (layer == 4 && sector == 3) channel_per_module = 14;
						if (layer == 4 && sector == 4) channel_per_module = 15;
						int channel = channel_per_module + (module-1)*16;
						if(charge < 0 && logical == 1)
						{
				
						//if(layer==1) Layer1_Energy_v2->Fill(point_energy/sintheta);
						//if(layer==2) Layer2_Energy_v2->Fill(point_energy/sintheta);
						//if(layer==3) Layer3_Energy_v2->Fill(point_energy/sintheta);
						//if(layer==4) Layer4_Energy_v2->Fill(point_energy/sintheta);
						//cout << " charge = " << locChargedTrackHypothesis->charge() << " point energy/sintheta= " << point_energy/sintheta << " 'shower energy' = " << E1 << " energy sum = " << point_energy_sum << " point energy = " << point_energy << " theta = " << theta << " channel = " << channel << " module = " << module << " layer = " << layer  << " point size = " << associated_points.size() << endl;
						}
						//cout << " charge = " << locChargedTrackHypothesis->charge() << " point energy/sintheta= " << point_energy/sintheta1 << " 'shower energy' = " << E1 << " shower momentum = " << p1.M() << " channel = " << channel << " module = " << module1 << " layer = " << layer1 << " point size = " << associated_points.size() << endl;
				

						

	
						
					}
				}
			//}

		}
		}
	}


*/


			for(unsigned int loc_i = 0; loc_i < locBCALShowers.size(); loc_i++)
   			 {
				const DBCALShower* locBCALShower = locBCALShowers[loc_i];
				if(find(matchedShowersneg.begin(), matchedShowersneg.end(), locBCALShower) == matchedShowersneg.end()) continue;
				 energy_shower = locBCALShower->E;	
				energy_raw_shower = locBCALShower->E_raw;
				 energy_raw_shower = locBCALShower->E_raw;
				track_momentum = p;
				double shower_x = locBCALShower->x;
				double shower_y = locBCALShower->y;

				 r = TMath::Sqrt(TMath::Power(shower_x,2)+TMath::Power(shower_y,2));
				z = locBCALShower->z;
				//phi = locBCALShower->phi;
				//theta = locBCALShower->theta;
				vector<const DBCALCluster*> associated_clusters;
				locBCALShower->Get(associated_clusters);

				if(associated_clusters.size() == 0) cout << " assoc cluster is empty " << endl;

				for(unsigned int loc_j = 0; loc_j < associated_clusters.size(); loc_j++)
				{

					const DBCALCluster* a_cluster = associated_clusters[loc_j];
					vector<const DBCALPoint*> associated_points;
					a_cluster->Get(associated_points); 
					layer1_energysum = 0.0;
					layer2_energysum = 0.0;
					layer3_energysum = 0.0;
					layer4_energysum = 0.0;
					double layer1_energysum_inter = 0.0;
					double layer2_energysum_inter = 0.0;
					double layer3_energysum_inter = 0.0;
					double layer4_energysum_inter = 0.0;
					phi = associated_clusters[loc_j]->phi();
					theta = associated_clusters[loc_j]->theta();

					for(unsigned int loc_k = 0 ; loc_k < associated_points.size(); loc_k++)
					{
						energy_point = associated_points[loc_k]->E();
						module = associated_points[loc_k]->module(); 
						layer = associated_points[loc_k]->layer();
						sector = associated_points[loc_k]->sector();
						//theta = associated_points[loc_k]->theta();
					       //	 z = associated_points[loc_k]->z();
			       			// phi = associated_points[loc_k]->phi();
		       				// r = associated_points[loc_k]->r();
						double first_pos = TMath::Sqrt( TMath::Power(myposL1.X(),2)+TMath::Power(myposL1.Y(),2)+TMath::Power(myposL1.Z(),2));
						double second_pos = TMath::Sqrt( TMath::Power(myposL2.X(),2)+TMath::Power(myposL2.Y(),2)+TMath::Power(myposL2.Z(),2));
						double third_pos = TMath::Sqrt( TMath::Power(myposL3.X(),2)+TMath::Power(myposL3.Y(),2)+TMath::Power(myposL3.Z(),2));
						double fourth_pos = TMath::Sqrt( TMath::Power(myposL4.X(),2)+TMath::Power(myposL4.Y(),2)+TMath::Power(myposL4.Z(),2));
						double fifth_pos = TMath::Sqrt( TMath::Power(myposL5.X(),2)+TMath::Power(myposL5.Y(),2)+TMath::Power(myposL5.Z(),2));
						 L1_pathlength = second_pos - first_pos;
						 L2_pathlength = third_pos - second_pos;
						 L3_pathlength = fourth_pos - third_pos;
						 L4_pathlength = fifth_pos - fourth_pos;

						if(layer==1) layer1_energysum_inter += energy_point;
						if(layer==2) layer2_energysum_inter += energy_point;
						if(layer==3) layer3_energysum_inter += energy_point;
						if(layer==4) layer4_energysum_inter += energy_point;
						if(loc_k+1 == associated_points.size()){
							layer1_energysum = layer1_energysum_inter;
							layer2_energysum = layer2_energysum_inter;
							layer3_energysum = layer3_energysum_inter;
							layer4_energysum = layer4_energysum_inter;
						}
						//cout << " shower E = " << energy_shower << " raw shower E = " << energy_raw_shower << " layer1_energysum = " << layer1_energysum << " layer 2 energy sum = " << layer2_energysum << " layer 3 energy sum = " << layer3_energysum << " layer4 energy sum = " << layer4_energysum << " p = " << track_momentum << endl;
						//if(layer1_energysum > 0.0 || layer2_energysum > 0.0 || layer3_energysum > 0.0 || layer4_energysum > 0.0 ) cout << " shower E = " << energy_shower << " raw shower E = " << energy_raw_shower << " layer1_energysum = " << layer1_energysum << " layer 2 energy sum = " << layer2_energysum << " layer 3 energy sum = " << layer3_energysum << " layer4 energy sum = " << layer4_energysum << " p = " << track_momentum << " REPEATER " << endl;
						int channel_per_module;
						if (layer == 1 && sector == 1) channel_per_module = 0;
						if (layer == 1 && sector == 2) channel_per_module = 1;
						if (layer == 1 && sector == 3) channel_per_module = 2;	
						if (layer == 1 && sector == 4) channel_per_module = 3;	
						if (layer == 2 && sector == 1) channel_per_module = 4;
						if (layer == 2 && sector == 2) channel_per_module = 5;	
						if (layer == 2 && sector == 3) channel_per_module = 6;
						if (layer == 2 && sector == 4) channel_per_module = 7;
						if (layer == 3 && sector == 1) channel_per_module = 8;
						if (layer == 3 && sector == 2) channel_per_module = 9;	
						if (layer == 3 && sector == 3) channel_per_module = 10;	
						if (layer == 3 && sector == 4) channel_per_module = 11;	
						if (layer == 4 && sector == 1) channel_per_module = 12;
						if (layer == 4 && sector == 2) channel_per_module = 13;
						if (layer == 4 && sector == 3) channel_per_module = 14;
						if (layer == 4 && sector == 4) channel_per_module = 15;
						 channel = channel_per_module + (module-1)*16;
						//double sintheta = sin(theta);
						
					//	 cout << " energy = " << energy << " charge = " << charge << " module = " << module << " layer = " << layer << " point size = " << associated_points.size() << " event num = " << eventnum << " theta = " << theta << " z = " << z << " r = " << r << " vertexZ = " << kinfitVertexZ << " theta wrt vertex = " << theta_wrt_vertex << endl;


						if(layer1_energysum > 0.0 || layer2_energysum > 0.0 || layer3_energysum > 0.0 || layer4_energysum > 0.0 ) BCALPoint_Charged_neg->Fill();
					}
				}
			}						



			for(unsigned int loc_i = 0; loc_i < locBCALShowers.size(); loc_i++)
   			 {
				const DBCALShower* locBCALShower = locBCALShowers[loc_i];
				if(find(matchedShowerspos.begin(), matchedShowerspos.end(), locBCALShower) == matchedShowerspos.end()) continue;
				 energy_shower = locBCALShower->E;	
				 energy_raw_shower = locBCALShower->E_raw;
				double shower_x = locBCALShower->x;
				double shower_y = locBCALShower->y;
				track_momentum = p;
				z = locBCALShower->z;
				//phi = locBCALShower->phi;
				//theta = locBCALShower->theta;
				r = TMath::Sqrt(TMath::Power(shower_x,2)+TMath::Power(shower_y,2));
				vector<const DBCALCluster*> associated_clusters;
				locBCALShower->Get(associated_clusters);

				if(associated_clusters.size() == 0) cout << " assoc cluster is empty " << endl;

				for(unsigned int loc_j = 0; loc_j < associated_clusters.size(); loc_j++)
				{

					const DBCALCluster* a_cluster = associated_clusters[loc_j];
					vector<const DBCALPoint*> associated_points;
					a_cluster->Get(associated_points); 

					phi = associated_clusters[loc_j]->phi();
					theta = associated_clusters[loc_j]->theta();

					 layer1_energysum = 0.0;
                                        layer2_energysum = 0.0;
                                        layer3_energysum = 0.0;
                                        layer4_energysum = 0.0;
                                        double layer1_energysum_inter = 0.0;
                                        double layer2_energysum_inter = 0.0;
                                        double layer3_energysum_inter = 0.0;
                                        double layer4_energysum_inter = 0.0;

					for(unsigned int loc_k = 0 ; loc_k < associated_points.size(); loc_k++)
					{
						energy_point = associated_points[loc_k]->E();
						module = associated_points[loc_k]->module(); 
						layer = associated_points[loc_k]->layer();
						sector = associated_points[loc_k]->sector();
						//theta = associated_points[loc_k]->theta();
					    //   	 z = associated_points[loc_k]->z();
			       			// phi = associated_points[loc_k]->phi();
		       				// r = associated_points[loc_k]->r();
						double first_pos = TMath::Sqrt( TMath::Power(myposL1.X(),2)+TMath::Power(myposL1.Y(),2)+TMath::Power(myposL1.Z(),2));
						double second_pos = TMath::Sqrt( TMath::Power(myposL2.X(),2)+TMath::Power(myposL2.Y(),2)+TMath::Power(myposL2.Z(),2));
						double third_pos = TMath::Sqrt( TMath::Power(myposL3.X(),2)+TMath::Power(myposL3.Y(),2)+TMath::Power(myposL3.Z(),2));
						double fourth_pos = TMath::Sqrt( TMath::Power(myposL4.X(),2)+TMath::Power(myposL4.Y(),2)+TMath::Power(myposL4.Z(),2));
						double fifth_pos = TMath::Sqrt( TMath::Power(myposL5.X(),2)+TMath::Power(myposL5.Y(),2)+TMath::Power(myposL5.Z(),2));
						 L1_pathlength = second_pos - first_pos;
						 L2_pathlength = third_pos - second_pos;
						 L3_pathlength = fourth_pos - third_pos;
						 L4_pathlength = fifth_pos - fourth_pos;

						int channel_per_module;
						if (layer == 1 && sector == 1) channel_per_module = 0;
						if (layer == 1 && sector == 2) channel_per_module = 1;
						if (layer == 1 && sector == 3) channel_per_module = 2;	
						if (layer == 1 && sector == 4) channel_per_module = 3;	
						if (layer == 2 && sector == 1) channel_per_module = 4;
						if (layer == 2 && sector == 2) channel_per_module = 5;	
						if (layer == 2 && sector == 3) channel_per_module = 6;
						if (layer == 2 && sector == 4) channel_per_module = 7;
						if (layer == 3 && sector == 1) channel_per_module = 8;
						if (layer == 3 && sector == 2) channel_per_module = 9;	
						if (layer == 3 && sector == 3) channel_per_module = 10;	
						if (layer == 3 && sector == 4) channel_per_module = 11;	
						if (layer == 4 && sector == 1) channel_per_module = 12;
						if (layer == 4 && sector == 2) channel_per_module = 13;
						if (layer == 4 && sector == 3) channel_per_module = 14;
						if (layer == 4 && sector == 4) channel_per_module = 15;
						 channel = channel_per_module + (module-1)*16;
						//double sintheta = sin(theta);

                                                if(layer==1) layer1_energysum_inter += energy_point;
                                                if(layer==2) layer2_energysum_inter += energy_point;
                                                if(layer==3) layer3_energysum_inter += energy_point;
                                                if(layer==4) layer4_energysum_inter += energy_point;
                                                if(loc_k+1 == associated_points.size()){
                                                        layer1_energysum = layer1_energysum_inter;
                                                        layer2_energysum = layer2_energysum_inter;
                                                        layer3_energysum = layer3_energysum_inter;
                                                        layer4_energysum = layer4_energysum_inter;
						}						
					//	 cout << " energy = " << energy << " charge = " << charge << " module = " << module << " layer = " << layer << " point size = " << associated_points.size() << " event num = " << eventnum << " theta = " << theta << " z = " << z << " r = " << r << " vertexZ = " << kinfitVertexZ << " theta wrt vertex = " << theta_wrt_vertex << endl;


						if(layer1_energysum > 0.0 || layer2_energysum > 0.0 || layer3_energysum > 0.0 || layer4_energysum > 0.0 ) BCALPoint_Charged_pos->Fill();
					}
				}
			}						




//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL LOOKING FOR E OVER P LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

	


	for (unsigned int i=0; i < locTrackTimeBased.size() ; i++)
	{
	if (find(matchedTracks.begin(), matchedTracks.end(), locTrackTimeBased[i]) == matchedTracks.end()) continue;
	double charge = locTrackTimeBased[i]->rt->q;
	double p = locTrackTimeBased[i]->momentum().Mag();
	//cout << " p = " << p << endl;

		for(unsigned int i=0; i<locBCALShowers.size(); i++)
		{
			const DBCALShower *s1 = locBCALShowers[i];
			if (find(matchedShowers.begin(), matchedShowers.end(),s1) == matchedShowers.end()) continue;
			double E_MatchedShower = s1->E;
			//cout << " shower E = " << E_MatchedShower << " track p = " << p << endl;
			//if(charge==1.0)Eoverp_plus->Fill(E_MatchedShower/p);
			//if(charge==-1.0)Eoverp_minus->Fill(E_MatchedShower/p);
			//if(charge==1.0)Evsp_plus->Fill(E_MatchedShower,p);
			//if(charge==-1.0)Evsp_minus->Fill(E_MatchedShower,p);
		
		}	
	
	}
	











// MMMMMMMMMMMMMMMMMMMMMMMMMMMM MATHCING TEST MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
/*
double locDistance=99999;
const DBCALShower* locBCALShowerz;
for(size_t i = 0 ; i < locBCALShowers.size() ; ++i)
{
locDetectorMatches->Get_DistanceToNearestTrack(locBCALShowerz, locDistance);
cout << " shower track doca = " << locDistance << endl;
}
*/
  


	//cons locTrackTimeBased[j]->DReferenceTrajectory;
	//if(rt == NULL) continue;

  //	DVector3 mypos(0.0,0.0,0.0);
	
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
	{	
		double x = locBCALShowers[loc_i]->x;
		double y = locBCALShowers[loc_i]->y;
		double z = locBCALShowers[loc_i]->z;
		double E = locBCALShowers[loc_i]->E;
		DVector3 pos_bcal(x,y,z);
		double R = pos_bcal.Perp();
		double R1 = sqrt((x-kinfitVertexX)*(x-kinfitVertexX)+(y-kinfitVertexY)*(y-kinfitVertexY)+(z-kinfitVertexZ)*(z-kinfitVertexZ));
		double phi = pos_bcal.Phi();
	  for(size_t j = 0 ; j < locTrackTimeBased.size(); ++j)
	    {
	
		locTrackTimeBased[j]->rt->GetIntersectionWithRadius(R, mypos);
		//double dPhi = sqrt((mypos.Phi() - pos_bcal.Phi())*(mypos.Phi() - pos_bcal().Phi()));
		double dPhi = TMath::Abs(mypos.Phi()-pos_bcal.Phi());
		double dZ = TMath::Abs(mypos.Z() - z);
		TLorentzVector p1(E*(x-kinfitVertexX)/R, E*(y-kinfitVertexY)/R1, E*(z-kinfitVertexZ)/R1, E);
	//	if(dZ < 10.0 && dPhi < 1.0) EoverP->Fill(E/p1.M());
		//	if(R==mypos.Perp())cout << " radius = " << mypos.Perp() << " R = " << R << " phi track = " << mypos.Phi() << " bcal phi = " << pos_bcal.Phi() << " dPhi = " << dPhi << " dZ = " << dZ << endl;
		//cout << " pos_bcal.perp = " << pos_bcal.Perp() << " pos_bcal.phi = " << pos_bcal.Phi() << endl;
	}
}

//double x = locBCALShowers->x;
//double y = locBCALShowers->y;
//double z = locBCALShowers->z;
//DVector3 bcal_pos(x, y, z);
//DVector3 proj_pos = rt->GetLastDOCAPoint();
//cout << " z pos = " << proj_pos.Perp() << endl;



for(unsigned int j = 0 ; j < locTrackTimeBased.size(); ++j)
	{		   
    for(unsigned int i = 0 ; i < locFCALShowers.size(); ++i)
      {
	double E = locFCALShowers[i]->getEnergy();
	double x = locFCALShowers[i]->getPosition().X();
	double y = locFCALShowers[i]->getPosition().Y();
	double z = locFCALShowers[i]->getPosition().Z();
	DVector3 origin(0.0,0.0,z);
	DVector3 norm(0.0,0.0,-1.0);
	DVector3 mypos2(0.0,0.0,0.0);
	DVector3 mymom(0.0,0.0,0.0);
	DVector3 pos_fcal(x,y,z);
	double myorigin = pos_fcal.Perp();
	

	locTrackTimeBased[j]->rt->GetIntersectionWithPlane(pos_fcal,norm,mypos2,mymom);
	double pmag = mymom.Mag();
	//cout << " p mag = " << pmag << " my pos perp = " << mypos2.Perp() << " my pos x = " << mypos2.x() << " my pos y = " << mypos2.y() << "my pos z = " << mypos2.z() << " pos fcal x = " << pos_fcal.x() << " pos fcal y = " << pos_fcal.y() << " pos fcal z = " << pos_fcal.z() << " pos fcal perp = " << pos_fcal.Perp() << " energy = " << E << " shower size = " << locFCALShowers.size() << " track size = " << locTrackTimeBased.size () << endl;

  
	double dX = TMath::Abs(mypos2.X() - x);
	double dY = TMath::Abs(mypos2.Y() - y);
	double dZ = TMath::Abs(mypos2.Z() - z);
	double dR = TMath::Abs(mypos2.Perp() - pos_fcal.Perp());
	//if(pos_fcal.z() == mypos2.z() && pos_fcal.Perp() > 20 && dR < 3.0) FCAL_Evsp->Fill(pmag,E);
	//if(pos_fcal.z() == mypos2.z() && pos_fcal.Perp() > 20 && dR < 3.0) FCAL_Eoverpvsp->Fill(pmag,E/pmag);
	//if(pos_fcal.z() == mypos2.z() && pos_fcal.Perp() > 20 && dR < 3.0 && pmag > 1.3 && pos_fcal.Perp() > 60.0) FCAL_Eoverp_nocuts->Fill(E/pmag);
	//if(pos_fcal.z() == mypos2.z() && dR < 3.0 && pmag > 1.3 && pos_fcal.Perp() > 20.0 && pos_fcal.Perp() < 60.0) FCAL_Eoverp_cuts->Fill(E/pmag);
	 //cout << " E = " << E << " p mag = " << pmag << " my pos 2 perp = " << mypos2.Perp() << " dx = " << dX << " dy = " << dY << " eventnum = " << eventnum << endl;
      }
  
}





//vector<double>  evsp_vec, p_vec;
for(size_t j = 0 ; j < locTrackTimeBased.size(); ++j)
{
 for(size_t loc_i = 0; loc_i < locFCALClusters.size(); ++loc_i)
{
//double x = locFCALShowers[loc_i]->getPosition().X();
//double y = locFCALShowers[loc_i]->getPosition().Y();
//double z = locFCALShowers[loc_i]->getPosition().Z();
DVector3 fcalpos=locFCALClusters[loc_i]->getCentroid();
//cout << " cluster x = " << fcalpos.X() << " cluster y = " << fcalpos.Y() << endl;
      //          fcalpos.SetZ( fcalpos.Z() + m_targetZ );
//DVector3 fcalpos(x,y,z);
DVector3 norm(0.0,0.0,-1);
DVector3 pos,mom;

//locTrackTimeBased[j]->rt->GetIntersectionWithPlane(origin,norm,pos);
if(locTrackTimeBased[j]->rt->GetIntersectionWithPlane(fcalpos,norm,pos,mom)==NOERROR)
{
double diffX = TMath::Abs(fcalpos.X() - pos.X());
double diffY = TMath::Abs(fcalpos.Y() - pos.Y());
double track_mom = TMath::Sqrt( TMath::Power(mom.X(),2) +
TMath::Power(mom.Y(),2) + TMath::Power(mom.Z(),2));
double cluster_energy = locFCALClusters[loc_i]->getEnergy();
double EvsP = cluster_energy/track_mom;
double cluster_radius = TMath::Sqrt( TMath::Power(fcalpos.X(),2) +
TMath::Power(fcalpos.Y(),2));
//cout << " cluster E = " << cluster_energy << " cluster radius = " << cluster_radius << " track momentum = " << track_mom << " diff x = " << diffX << " diff y = " << diffY << " event num = " << eventnum << endl;
//m_evsp->Fill( EvsP );
//if(diffX < 10 && diffY < 10) cout << " clusrer E = " << cluster_energy << " track p = " << track_mom << " cluster radius = " << cluster_radius << endl;
if(diffX < 10 && diffY < 10){
 FCAL_Eoverp_nocuts->Fill( EvsP );
 FCAL_Evsp->Fill(track_mom,cluster_energy);
 FCAL_Eoverpvsp->Fill(track_mom,cluster_energy/track_mom);
}
//std::cout<<fcalpos.X()<<"\t"<<fcalpos.Y()<<"\t"<<fcalpos.Z()<<endl;
//std::cout<<pos.X()<<"\t"<<pos.Y()<<"\t"<<pos.Z()<<endl;
//std::cout<<endl;
if (diffX < 5 && diffY < 5 && track_mom >1.3 && cluster_radius > 20.0 && cluster_radius < 60.0)
{
//m_posX->Fill( fcalpos.X() - pos.X() );
//m_posY->Fill( fcalpos.Y() - pos.Y() );
//m_energy->Fill(cluster_energy);
//m_mom->Fill(track_mom);
//m_evspcut->Fill( EvsP);
FCAL_Eoverp_cuts->Fill( EvsP );
//m_radius->Fill(cluster_radius);
//evsp_vec.push_back(EvsP);
//p_vec.push_back(track_mom);
//m_ebypvsp->Fill( EvsP, track_mom );
//std::cout<<"Real Tracks:"<<endl;
//std::cout<<fcalpos.X()<<"\t"<<fcalpos.Y()<<"\t"<<fcalpos.Z()<<endl;
//std::cout<<pos.X()<<"\t"<<pos.Y()<<"\t"<<pos.Z()<<endl;
//std::cout<<endl;
}
  }
}
}


  





//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw Fill ROOT Tree Test WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

/*
for(size_t i = 0 ; i < locChargedTracks.size(); ++i)
	{
	
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[i]->Get_BestFOM();
	//	cout << " Best FOM = " << locChargedTracks[i]->Get_BestFOM() << endl;
		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		const DShowerMatchParams& locBCALShowerMatchParams = locChargedTrackHypothesis->dBCALShowerMatchParams;
		 
	//	if(locBCALShowerMatchParams.dTrackTimeBased != NULL)
	//	{
			double d_min = locBCALShowerMatchParams.dDOCAToShower;
			//if(d_min != 0.0) cout <<  " d min = " << d_min << endl;
			const DBCALShower* locBCALShower = dynamic_cast <const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
			if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShower)){
			BCALShowerTrack_Energy->Fill(locBCALShower->E);
			charge = locChargedTrackHypothesis->charge();
			//cout << " shower track E = " << BCALShowerTrack_Energy << endl;		
			//for(unsigned int loc_i = 0; loc_i < locBCALShowers.size(); loc_i++)
   			 //{
  				//const DBCALShower* s1 = locBCALShowers;
				//const DBCALShower* a_shower = locBCALShowers;

				//double shower_energy = locBCALShower->E;	
				//double path1 = sqrt((dx1-kinfitVertexX)*(dx1-kinfitVertexX)+(dy1-kinfitVertexY)*(dy1-kinfitVertexY)+(dz1)*(dz1));
				vector<const DBCALCluster*> associated_clusters;
				locBCALShower->Get(associated_clusters);

				if(associated_clusters.size() == 0) cout << " assoc cluster is empty " << endl;

				for(unsigned int loc_j = 0; loc_j < associated_clusters.size(); loc_j++)
				{

					const DBCALCluster* a_cluster = associated_clusters[loc_j];
					vector<const DBCALPoint*> associated_points;
					a_cluster->Get(associated_points); 

					if(associated_points.size() == 0) cout << " assoc points is empty " << endl;

					for(unsigned int loc_k = 0; loc_k < associated_points.size(); loc_k++)
					{
						//cout << " charged point size = " << associated_points.size() << endl;
						energy = associated_points[loc_k]->E();
						module = associated_points[loc_k]->module(); 
						layer = associated_points[loc_k]->layer();
						sector = associated_points[loc_k]->sector();
						theta = associated_points[loc_k]->theta();
					       	 z = associated_points[loc_k]->z();
			       			 phi = associated_points[loc_k]->phi();
		       				 r = associated_points[loc_k]->r();
					       	theta_wrt_vertex =atan(r/(z-kinfitVertexZ));

						int channel_per_module;
						if (layer == 1 && sector == 1) channel_per_module = 0;
						if (layer == 1 && sector == 2) channel_per_module = 1;
						if (layer == 1 && sector == 3) channel_per_module = 2;	
						if (layer == 1 && sector == 4) channel_per_module = 3;	
						if (layer == 2 && sector == 1) channel_per_module = 4;
						if (layer == 2 && sector == 2) channel_per_module = 5;	
						if (layer == 2 && sector == 3) channel_per_module = 6;
						if (layer == 2 && sector == 4) channel_per_module = 7;
						if (layer == 3 && sector == 1) channel_per_module = 8;
						if (layer == 3 && sector == 2) channel_per_module = 9;	
						if (layer == 3 && sector == 3) channel_per_module = 10;	
						if (layer == 3 && sector == 4) channel_per_module = 11;	
						if (layer == 4 && sector == 1) channel_per_module = 12;
						if (layer == 4 && sector == 2) channel_per_module = 13;
						if (layer == 4 && sector == 3) channel_per_module = 14;
						if (layer == 4 && sector == 4) channel_per_module = 15;
						 channel = channel_per_module + (module-1)*16;
						//double sintheta = sin(theta);
						
					//	 cout << " energy = " << energy << " charge = " << charge << " module = " << module << " layer = " << layer << " point size = " << associated_points.size() << " event num = " << eventnum << " theta = " << theta << " z = " << z << " r = " << r << " vertexZ = " << kinfitVertexZ << " theta wrt vertex = " << theta_wrt_vertex << endl;


						BCALPoint_Charged->Fill();
						
					}
				}
			//}

		}
	//	}
	} 


*/

                                      

	for(unsigned int i=0; i<locBCALShowers.size(); i++)
	{
	     //   if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[i]))
               // continue;
		if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[i]) != matchedShowers.end()) continue;
		const DBCALShower *s1 = locBCALShowers[i];
		vector<const DBCALCluster*> associated_clusters1;
		
		s1->Get(associated_clusters1);
			double shower_energy = s1->E;
			double raw_shower_energy = s1->E_raw;
		//cout << " shower E = " << shower_energy << endl;
		//cout << " shower raw E = " << raw_shower_energy << endl;
		for(unsigned int loc_j = 0; loc_j < associated_clusters1.size(); loc_j++)
		{
			double cluster_energy = associated_clusters1[loc_j]->E();
			//cout << " cluster E = " << cluster_energy << endl;
			const DBCALCluster* a_cluster1 = associated_clusters1[loc_j];
			vector<const DBCALPoint*> associated_points1;
			a_cluster1->Get(associated_points1); 
			for(unsigned int loc_k = 0; loc_k < associated_points1.size(); loc_k++)
			{
				 module = associated_points1[loc_k]->module(); 
				 layer = associated_points1[loc_k]->layer();
				 sector = associated_points1[loc_k]->sector();
			      	int channel_per_module;
	       			if (layer == 1 && sector == 1) channel_per_module = 0;
       				if (layer == 1 && sector == 2) channel_per_module = 1;
	       			if (layer == 1 && sector == 3) channel_per_module = 2;	
	       			if (layer == 1 && sector == 4) channel_per_module = 3;	
       				if (layer == 2 && sector == 1) channel_per_module = 4;
			      	if (layer == 2 && sector == 2) channel_per_module = 5;	
	      			if (layer == 2 && sector == 3) channel_per_module = 6;
      			        if (layer == 2 && sector == 4) channel_per_module = 7;
	       			if (layer == 3 && sector == 1) channel_per_module = 8;
  		      		if (layer == 3 && sector == 2) channel_per_module = 9;	
	      			if (layer == 3 && sector == 3) channel_per_module = 10;	
				if (layer == 3 && sector == 4) channel_per_module = 11;	
			       	if (layer == 4 && sector == 1) channel_per_module = 12;
		       		if (layer == 4 && sector == 2) channel_per_module = 13;
	       			if (layer == 4 && sector == 3) channel_per_module = 14;
       				if (layer == 4 && sector == 4) channel_per_module = 15;
		       		 channel = channel_per_module + (module-1)*16;
				 energy_point = associated_points1[loc_k]->E();
				 theta = associated_points1[loc_k]->theta();
				 	// cout  << " point E = " << energy << " point size = " << associated_points1.size() << " shower size = " << locBCALShowers.size() << " channel = " << channel << " event num = " << eventnum << endl;
				 BCALPoint_Neutral->Fill();

			}
		}
	}



 


/*			
			vector<const DBCALPoint* > bcalpoints;
			locEventLoop->Get(bcalpoints);
				for(unsigned int i = 0; i < bcalpoints.size(); i++)
				  {
			           if(locDetectorMatches->Get_IsMatchedToTrack(bcalpoints[i]))
				     continue;
					
					 theta = bcalpoints[i]->theta();
					 energy = bcalpoints[i]->E();
					 module = bcalpoints[i]->module(); 
					 layer = bcalpoints[i]->layer();
					 sector = bcalpoints[i]->sector();	
						int channel_per_module;
						if (layer == 1 && sector == 1) channel_per_module = 0;
						if (layer == 1 && sector == 2) channel_per_module = 1;
						if (layer == 1 && sector == 3) channel_per_module = 2;	
						if (layer == 1 && sector == 4) channel_per_module = 3;	
						if (layer == 2 && sector == 1) channel_per_module = 4;
						if (layer == 2 && sector == 2) channel_per_module = 5;	
						if (layer == 2 && sector == 3) channel_per_module = 6;
						if (layer == 2 && sector == 4) channel_per_module = 7;
						if (layer == 3 && sector == 1) channel_per_module = 8;
						if (layer == 3 && sector == 2) channel_per_module = 9;	
						if (layer == 3 && sector == 3) channel_per_module = 10;	
						if (layer == 3 && sector == 4) channel_per_module = 11;	
						if (layer == 4 && sector == 1) channel_per_module = 12;
						if (layer == 4 && sector == 2) channel_per_module = 13;
						if (layer == 4 && sector == 3) channel_per_module = 14;
						if (layer == 4 && sector == 4) channel_per_module = 15;
						 channel = channel_per_module + (module-1)*16;
					 BCALPoint_Neutral->Fill();
					
				}






*/










//OoOoOoOoOoOoOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO CALCULATING THE PI0 INVARIANT MASS PORTION O0O0O0O0O0O0O0O0OOooooooooooOooOoOoooooOOooooOOo





	//2 gamma inv mass 




	for(unsigned int i = 0; i < matchedShowers.size(); i++)
	{
	double E = matchedShowers[i]->E;
	//cout << " E = " << E << endl;
	const DBCALShower *ms = matchedShowers[i];
	vector<const DBCALCluster*> associated_clusters;
	ms->Get(associated_clusters);
	for(unsigned int j = 0 ; j < associated_clusters.size(); j++)
		{
		double cluster_E = associated_clusters[j]->E();
		//cout << " cluster E = " << cluster_E << endl;
		}
	
	}



	for(unsigned int i=0; i<locBCALShowers.size(); i++)
	{
	     //   if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[i]))
               // continue;
		double E_before = locBCALShowers[i]->E;
		BCALShowers_per_event = locBCALShowers.size();
		//cout << " E before = " << E_before << endl;
		if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[i]) != matchedShowers.end()) continue;
		const DBCALShower *s1 = locBCALShowers[i];
		vector<const DBCALCluster*> associated_clusters1;
		s1->Get(associated_clusters1);
		double dx1 = s1->x - kinfitVertexX;
		double dy1 = s1->y - kinfitVertexY;
		double dz1 = s1->z - kinfitVertexZ;
		double dt1 = s1->t;
		t1 = dt1;
		z1 = s1->z;
		x1 = s1->x;
		y1 = s1->y;
	
		double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
		//double path1 = sqrt((dx1-kinfitVertexX)*(dx1-kinfitVertexX)+(dy1-kinfitVertexY)*(dy1-kinfitVertexY)+(dz1)*(dz1));
		 E1 = s1->E;
		 E1_raw = s1->E_raw;
		
		//cout << " E after = " << E1 << endl;
		double ECalc = s1->E*(1.106+(dz1+65.0-208.4)*(dz1+65.0-208.4)*6.851E-6);
		TLorentzVector p1(E1*dx1/R1, E1*dy1/R1, E1*dz1/R1, E1);
		TLorentzVector p1_raw(E1_raw*dx1/R1, E1_raw*dy1/R1, E1_raw*dz1/R1, E1_raw);
		double thetadeg1, thetarad1, phideg1, phirad1;
		thetadeg1 = p1.Theta()*180.0/TMath::Pi();
		thetarad1 = p1.Theta();
		phideg1 = p1.Phi()*180.0/TMath::Pi();
		phirad1 = p1.Phi();	
		phi1 = phirad1;
		p1_mag = p1.M();
		p1_raw_mag = p1_raw.M();
		vertexZ = kinfitVertexZ;
		vertexX = kinfitVertexX;
		vertexY = kinfitVertexY;
		TLorentzVector p1Calc(ECalc*dx1/R1, ECalc*dy1/R1, ECalc*dz1/R1, ECalc);


			for(unsigned int j=i+1; j<locBCALShowers.size(); j++){
	         //       if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[j]))
                   //     continue;
			const DBCALShower *s2 = locBCALShowers[j];
			//cout << " energy before matching = " << s2->E << endl;
		     	if (find(matchedShowers.begin(), matchedShowers.end(),s2) != matchedShowers.end()) continue;
			//cout << " energy after matching = " << s2->E << endl;
			vector<const DBCALCluster*> associated_clusters2;
			s2->Get(associated_clusters2);
			double dx2 = s2->x - kinfitVertexX;
			double dy2 = s2->y - kinfitVertexY;
			double dz2 = s2->z - kinfitVertexZ; // shift to coordinate relative to center of target
			double dt2 = s2->t;
			t2=dt2;
			z2 = s2->z;
			x2 = s2->x;
			y2 = s2->y;
			double deltaT = TMath::Abs(dt1-dt2);
			Time_Diff->Fill(deltaT);
			double R2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
			//double path2 = sqrt((dx2-kinfitVertexX)*(dx2-kinfitVertexX)+(dy2-kinfitVertexY)*(dy2-kinfitVertexY)+(dz2)*(dz2));
			 E2 = s2->E;
			 E2_raw = s2->E_raw;
			double ECalc = s2->E*(1.106+(dz2+65.0-208.4)*(dz2+65.0-208.4)*6.851E-6);
			TLorentzVector p2(E2*dx2/R2, E2*dy2/R2, E2*dz2/R2, E2);	
			TLorentzVector p2_raw(E2_raw*dx2/R2, E2_raw*dy2/R2, E2_raw*dz2/R2, E2_raw);	
			p2_mag = p2.M();
			p2_raw_mag = p2_raw.M();
			TLorentzVector p2Calc(ECalc*dx2/R2, ECalc*dy2/R2, ECalc*dz2/R2, ECalc);	
			double thetadeg2, thetarad2, phideg2, phirad2, cospsi, psi;
			thetarad2 = p2.Theta();
			phirad2 = p2.Phi();
			phi2 = phirad2;
			thetadeg2 = p2.Theta()*180.0/TMath::Pi();
			phideg2 = p2.Phi()*180.0/TMath::Pi();					
			cospsi = sin(thetarad1)*sin(thetarad2)*cos(phirad1-phirad2)+cos(thetarad1)*cos(thetarad2);
			psi=acos(cospsi)*180/TMath::Pi();
			TLorentzVector ptot = p1+p2;
			TLorentzVector ptot_raw = p1_raw + p2_raw ;
			inv_mass = ptot.M();
			inv_mass_raw = ptot_raw.M();

			Theta_Hist_Both_BCAL->Fill(thetadeg1);
			Theta_Hist_Both_BCAL->Fill(thetadeg2);
			Phi_Hist_Both_BCAL->Fill(phideg2);
			Phi_Hist_Both_BCAL->Fill(phideg2);
			Psi_Hist_Both_BCAL->Fill(psi);
			
			double makes_sense = 0;
			if (R1 > R2 && dt1 > dt2) makes_sense = 1;
			if (R1 < R2 && dt1 < dt2) makes_sense = 1;
			if (R1 < R2 && dt1 > dt2) makes_sense = 0;
			if (R2 < R1 && dt1 > dt1) makes_sense = 0;
		/*	for(unsigned int loc_j = 0; loc_j < associated_clusters1.size(); loc_j++)
			{
				const DBCALCluster* a_cluster1 = associated_clusters1[loc_j];
				vector<const DBCALPoint*> associated_points1;
				a_cluster1->Get(associated_points1); 
					for(unsigned int loc_jj = 0; loc_jj < associated_clusters2.size(); loc_jj++)
				{
					const DBCALCluster* a_cluster2 = associated_clusters2[loc_jj];
					vector<const DBCALPoint*> associated_points2;
					a_cluster2->Get(associated_points2); 
					for(unsigned int loc_k = 0; loc_k < associated_points1.size(); loc_k++)
					{
						int module1 = associated_points1[loc_k]->module(); 
						int layer1 = associated_points1[loc_k]->layer();
						int sector1 = associated_points1[loc_k]->sector();
						double energy1 = associated_points1[loc_k]->E();
						int channel_per_module1;
						if (layer1 == 1 && sector1 == 1) channel_per_module1 = 0;
						if (layer1 == 1 && sector1 == 2) channel_per_module1 = 1;
						if (layer1 == 1 && sector1 == 3) channel_per_module1 = 2;	
						if (layer1 == 1 && sector1 == 4) channel_per_module1 = 3;	
						if (layer1 == 2 && sector1 == 1) channel_per_module1 = 4;
						if (layer1 == 2 && sector1 == 2) channel_per_module1 = 5;	
						if (layer1 == 2 && sector1 == 3) channel_per_module1 = 6;
						if (layer1 == 2 && sector1 == 4) channel_per_module1 = 7;
						if (layer1 == 3 && sector1 == 1) channel_per_module1 = 8;
						if (layer1 == 3 && sector1 == 2) channel_per_module1 = 9;	
						if (layer1 == 3 && sector1 == 3) channel_per_module1 = 10;	
						if (layer1 == 3 && sector1 == 4) channel_per_module1 = 11;	
						if (layer1 == 4 && sector1 == 1) channel_per_module1 = 12;
						if (layer1 == 4 && sector1 == 2) channel_per_module1 = 13;
						if (layer1 == 4 && sector1 == 3) channel_per_module1 = 14;
						if (layer1 == 4 && sector1 == 4) channel_per_module1 = 15;
						 channel = channel_per_module1 + (module1-1)*16;
						
						//cout << " moodule = " << module1 << " layer = " << layer1 << " sector = " << sector1 << " channel = " << channel << endl;
					
						//if( ptot.M() > 0.12 && ptot.M() < 0.18) cout << "shower1 energy = " << E1  << " shower2 energy =  " << E2 << " shower1 time = " << dt1 << " shower 2 time = " << dt2 << " shower 1 z = " << dz1 << " shower 2 z = " << dz2 << " ptot = " << ptot.M() <<  " event num = " << eventnum << endl;
					
							
							for(unsigned int loc_kk = 0; loc_kk < associated_points1.size(); loc_kk++)
							{
							  //	int module2 = associated_points2[loc_kk]->module(); 
							  //	int layer2 = associated_points2[loc_kk]->layer();
							  //	int sector2 = associated_points2[loc_kk]->sector();
							  //	int energy2 = associated_points2[loc_kk]->E();
							//	cout << "shower1 energy = " << E1  << " shower2 energy =  " << E2 <<"  point1 energy = " << energy1 << " point2 energy = " << energy2 << " channel = " << channel << " point size = " << associated_points1.size() << " shower size = " << locBCALShowers.size()  << " ptot = " << ptot.M() <<  " event num = " << eventnum << endl;
							
							}
					}
					}
			}   */
			//cout << " z1 = " << z1 << " z2 = " << z2 << " shower size = " << locBCALShowers.size() << " eventnum = " << eventnum << endl;
			BCAL_Neutrals->Fill();
			if (E1 > 0.5 && E2 > 0.5 && kinfitVertexZ > 63.0 && kinfitVertexZ < 67.0 && kinfitVertexX > -15.0 && kinfitVertexX < 15.0 && kinfitVertexY > -15.0 && kinfitVertexY < 15.0) two_gamma_mass->Fill(ptot.M());
			if (E1 > 0.2 && E2 > 0.2 && kinfitVertexZ > 63.0 && kinfitVertexZ < 67.0) bcal_diphoton_mass->Fill(ptot.M());
			if (E1 > 0.7 && E2 > 0.7 && kinfitVertexZ > 63.0 && kinfitVertexZ < 67.0 && kinfitVertexX > -15.0 && kinfitVertexX < 15.0 && kinfitVertexY > -15.0 && kinfitVertexY < 15.0) bcal_diphoton_mass_highE->Fill(ptot.M());
			}
	}   


	/*	for(unsigned int j=0; j<locFCALShowers.size(); j++)
		{
			const DFCALShower *s3 = locFCALShowers[j];
 
			double dx3 = s3->getPosition().X() - kinfitVertexX;
			double dy3 = s3->getPosition().Y() - kinfitVertexY;
			double dz3 = 620.0 - kinfitVertexZ; // shift to coordinate relative to center of target
			double dt3 = s3->getTime();
			double R3 = sqrt(dx3*dx3 + dy3*dy3 + dz3*dz3);
			double E3 = s3->getEnergy();
			//double path = sqrt((dx-kinfitVertexX)*(dx-kinfitVertexX)+(dy-kinfitVertexY)*(dy-kinfitVertexY)+(dz)*(dz));
			TLorentzVector p3(E3*dx3/R3, E3*dy3/R3, E3*dz3/R3, E3);
			int makes_sense2 = 0;
			if(R3 > R1 && dt3 > dt1) makes_sense2=1;
			if(R3 < R1 && dt3 < dt1) makes_sense2=1;
			if(R3 > R1 && dt3 < dt1) makes_sense2=0;
			if(R3 < R1 && dt3 > dt1) makes_sense2=0;
			double thetadeg3, thetarad3, phideg3, phirad3, cospsi3, psi3;
			thetarad3 = p3.Theta();
			phirad3 = p3.Phi();
			thetadeg3 = p3.Theta()*180.0/TMath::Pi();
			phideg3 = p3.Phi()*180.0/TMath::Pi();					
			cospsi3 = sin(thetarad1)*sin(thetarad3)*cos(phirad1-phirad3)+cos(thetarad1)*cos(thetarad3);
			psi3=acos(cospsi3)*180/TMath::Pi();		
			
			TLorentzVector ptot = p1+p3;
			//if(E1>0.3 && E2>0.3 && kinfitVertexZ > 50 && kinfitVertexZ < 80) 
			if( kinfitVertexZ > 62.0 && kinfitVertexZ < 68.0 && E1 > 0.4 && E3 > 0.4 && makes_sense2==1 && kinfitVertexX < 15.0 && kinfitVertexX > -15.0 && kinfitVertexY > -15.0 && kinfitVertexY < 15.0)  bcal_fcal_two_gamma_mass->Fill(ptot.M());
			//if(E1>0.2 && E>0.2) bcal_fcal_two_gamma_mass->Fill(ptot.M());

			Theta_Hist_Split_Gammas->Fill(thetadeg1);
			Theta_Hist_Split_Gammas->Fill(thetadeg3);
			Phi_Hist_Split_Gammas->Fill(phideg1);
			Phi_Hist_Split_Gammas->Fill(phideg3);
			Psi_Hist_Split_Gammas->Fill(psi3);

			if(locBCALShowers.size()==1 && locFCALShowers.size()==1){
				//bcal_fcal_two_gamma_mass_cut->Fill(ptot.M());
			}
			
		
		}  */
		


		for(unsigned int j=0; j<locFCALShowers.size(); j++)
		{
		  // if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[j]))
		  //    continue;

			if (find(matchedFCALShowers.begin(), matchedFCALShowers.end(),locFCALShowers[j]) != matchedFCALShowers.end()) continue;
			const DFCALShower *s3 = locFCALShowers[j];
 
			double dx3 = s3->getPosition().X() - kinfitVertexX;
			double dy3 = s3->getPosition().Y() - kinfitVertexY;
			double dz3 = s3->getPosition().Z() - kinfitVertexZ; // shift to coordinate relative to center of target
			double dt3 = s3->getTime();
			double R3 = sqrt(dx3*dx3 + dy3*dy3 + dz3*dz3);
			double E3 = s3->getEnergy();
			fcal_x = s3->getPosition().X();
			fcal_y = s3->getPosition().Y();
			fcal_z = s3->getPosition().Z();
			fcal_E = s3->getEnergy();
			fcal_t = s3->getTime();
			FCALShowers_per_event = locFCALShowers.size();
			//double path = sqrt((dx-kinfitVertexX)*(dx-kinfitVertexX)+(dy-kinfitVertexY)*(dy-kinfitVertexY)+(dz)*(dz));
			TLorentzVector p3(E3*dx3/R3, E3*dy3/R3, E3*dz3/R3, E3);
			fcal_p=p3.M();
			double thetadeg3, thetarad3, phideg3, phirad3, cospsi3, psi3;
			thetarad3 = p3.Theta();
			phirad3 = p3.Phi();
			thetadeg3 = p3.Theta()*180.0/TMath::Pi();
			phideg3 = p3.Phi()*180.0/TMath::Pi();					
			


			for( unsigned int i=0; i<locBCALShowers.size(); i++)
			{
			  //   if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[i]))
			  //    continue;
		
				if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[i]) != matchedShowers.end()) continue;
				const DBCALShower *s1 = locBCALShowers[i];
				double dx1 = s1->x - kinfitVertexX;
				double dy1 = s1->y - kinfitVertexY;
				double dz1 = s1->z - kinfitVertexZ;
				double dt1 = s1->t;
				double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
				bcal_x = s1->x;
				bcal_y = s1->y;
				bcal_z = s1->z;
				bcal_E = s1->E;
				bcal_t = s1->t;
				vertexZ = kinfitVertexZ;
				vertexX = kinfitVertexX;
				vertexY = kinfitVertexY;
				BCALShowers_per_event = locBCALShowers.size();
				//double path1 = sqrt((dx1-kinfitVertexX)*(dx1-kinfitVertexX)+(dy1-kinfitVertexY)*(dy1-kinfitVertexY)+(dz1)*(dz1));
				double E1 = s1->E;
				double ECalc = s1->E*(1.106+(dz1+65.0-208.4)*(dz1+65.0-208.4)*6.851E-6);
				TLorentzVector p1(E1*dx1/R1, E1*dy1/R1, E1*dz1/R1, E1);
				bcal_p = p1.M();
				double thetadeg1, thetarad1, phideg1, phirad1;
				thetadeg1 = p1.Theta()*180.0/TMath::Pi();
				thetarad1 = p1.Theta();
				phideg1 = p1.Phi()*180.0/TMath::Pi();
				phirad1 = p1.Phi();	
				bcal_phi = phirad1;
				TLorentzVector p1Calc(ECalc*dx1/R1, ECalc*dy1/R1, ECalc*dz1/R1, ECalc);
				TLorentzVector ptot = p1+p3;
				inv_mass = ptot.M();
				
				cospsi3 = sin(thetarad1)*sin(thetarad3)*cos(phirad1-phirad3)+cos(thetarad1)*cos(thetarad3);
				psi3=acos(cospsi3)*180/TMath::Pi();		
				//if(E1>0.3 && E2>0.3 && kinfitVertexZ > 50 && kinfitVertexZ < 80) 
				if( kinfitVertexZ > 63.0 && kinfitVertexZ < 67.0 && E1 > 0.5 && E3 > 0.5)  bcal_fcal_two_gamma_mass->Fill(ptot.M());
				//if(E1>0.2 && E>0.2) bcal_fcal_two_gamma_mass->Fill(ptot.M());

				Theta_Hist_Split_Gammas->Fill(thetadeg1);
				Theta_Hist_Split_Gammas->Fill(thetadeg3);
				Phi_Hist_Split_Gammas->Fill(phideg1);
				Phi_Hist_Split_Gammas->Fill(phideg3);
				Psi_Hist_Split_Gammas->Fill(psi3);
				
		
				Split_Gamma_Neutrals->Fill();
			}
		
		}

		

		for(unsigned int j=0; j<locFCALClusters.size(); j++)
		{
		  // if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[j]))
		  //    continue;

			if (find(matchedFCALClusters.begin(), matchedFCALClusters.end(),locFCALClusters[j]) != matchedFCALClusters.end()) continue;
			const DFCALCluster *s3 = locFCALClusters[j];
 
			double dx3 = s3->getCentroid().X() - kinfitVertexX;
			double dy3 = s3->getCentroid().Y() - kinfitVertexY;
			double dz3 = s3->getCentroid().Z() - kinfitVertexZ; // shift to coordinate relative to center of target
			double dt3 = s3->getTime();
			double R3 = sqrt(dx3*dx3 + dy3*dy3 + dz3*dz3);
			double E3 = s3->getEnergy();
			fcal_x = s3->getCentroid().X();
			fcal_y = s3->getCentroid().Y();
			fcal_z = s3->getCentroid().Z();
			fcal_E = s3->getEnergy();
			fcal_t = s3->getTime();
			FCALClusters_per_event = locFCALClusters.size();
			//double path = sqrt((dx-kinfitVertexX)*(dx-kinfitVertexX)+(dy-kinfitVertexY)*(dy-kinfitVertexY)+(dz)*(dz));
			TLorentzVector p3(E3*dx3/R3, E3*dy3/R3, E3*dz3/R3, E3);
			fcal_p=p3.M();
			double thetadeg3, thetarad3, phideg3, phirad3, cospsi3, psi3;
			thetarad3 = p3.Theta();
			phirad3 = p3.Phi();
			thetadeg3 = p3.Theta()*180.0/TMath::Pi();
			phideg3 = p3.Phi()*180.0/TMath::Pi();					
			


			for( unsigned int i=0; i<locBCALShowers.size(); i++)
			{
			  //   if(locDetectorMatches->Get_IsMatchedToTrack(locBCALShowers[i]))
			  //    continue;
		
				if (find(matchedShowers.begin(), matchedShowers.end(),locBCALShowers[i]) != matchedShowers.end()) continue;
				const DBCALShower *s1 = locBCALShowers[i];
				double dx1 = s1->x - kinfitVertexX;
				double dy1 = s1->y - kinfitVertexY;
				double dz1 = s1->z - kinfitVertexZ;
				double dt1 = s1->t;
				double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
				bcal_x = s1->x;
				bcal_y = s1->y;
				bcal_z = s1->z;
				bcal_E = s1->E_raw;
				bcal_t = s1->t;
				vertexZ = kinfitVertexZ;
				vertexX = kinfitVertexX;
				vertexY = kinfitVertexY;
				BCALShowers_per_event = locBCALShowers.size();
				//double path1 = sqrt((dx1-kinfitVertexX)*(dx1-kinfitVertexX)+(dy1-kinfitVertexY)*(dy1-kinfitVertexY)+(dz1)*(dz1));
				double E1_raw = s1->E_raw;
				double ECalc = s1->E*(1.106+(dz1+65.0-208.4)*(dz1+65.0-208.4)*6.851E-6);
				TLorentzVector p1(E1_raw*dx1/R1, E1_raw*dy1/R1, E1_raw*dz1/R1, E1_raw);
				bcal_p = p1.M();
				double thetadeg1, thetarad1, phideg1, phirad1;
				thetadeg1 = p1.Theta()*180.0/TMath::Pi();
				thetarad1 = p1.Theta();
				phideg1 = p1.Phi()*180.0/TMath::Pi();
				phirad1 = p1.Phi();	
				bcal_phi = phirad1;
				TLorentzVector p1Calc(ECalc*dx1/R1, ECalc*dy1/R1, ECalc*dz1/R1, ECalc);
				TLorentzVector ptot = p1+p3;
				inv_mass = ptot.M();
				
				cospsi3 = sin(thetarad1)*sin(thetarad3)*cos(phirad1-phirad3)+cos(thetarad1)*cos(thetarad3);
				psi3=acos(cospsi3)*180/TMath::Pi();		
				//if(E1>0.3 && E2>0.3 && kinfitVertexZ > 50 && kinfitVertexZ < 80) 
				//if( kinfitVertexZ > 63.0 && kinfitVertexZ < 67.0 && E1 > 0.5 && E3 > 0.5)  bcal_fcal_two_gamma_mass->Fill(ptot.M());
				//if(E1>0.2 && E>0.2) bcal_fcal_two_gamma_mass->Fill(ptot.M());

				Theta_Hist_Split_Gammas->Fill(thetadeg1);
				Theta_Hist_Split_Gammas->Fill(thetadeg3);
				Phi_Hist_Split_Gammas->Fill(phideg1);
				Phi_Hist_Split_Gammas->Fill(phideg3);
				Psi_Hist_Split_Gammas->Fill(psi3);
				
		
				Split_Gamma_Neutrals_raw->Fill();
			}
		
		}

		




		for(unsigned int ii=0; ii<locFCALShowers.size(); ii++)
		{

		  //    if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[ii]))
		  //     continue;
			FCALShowers_per_event = locFCALShowers.size();
			if (find(matchedFCALShowers.begin(), matchedFCALShowers.end(),locFCALShowers[ii]) != matchedFCALShowers.end()) continue;
			const DFCALShower *s2 = locFCALShowers[ii];
			
			//double dx, dy, dz, R, thetarad1, phirad1,
			double dx1 = s2->getPosition().X() - kinfitVertexX;
			double dy1 = s2->getPosition().Y() - kinfitVertexY;
			double dz1 = s2->getPosition().Z() - kinfitVertexZ; // shift to coordinate relative to center of target
			double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
			 E1 = s2->getEnergy();
			double dt1 = s2->getTime();
			//double path1 = sqrt((dx-kinfitVertexX)*(dx-kinfitVertexX)+(dy-kinfitVertexY)*(dy-kinfitVertexY)+(dz)*(dz));g
			TLorentzVector p1(E1*dx1/R1, E1*dy1/R1, E1*dz1/R1, E1);
			x1 = s2->getPosition().X();
			y1 = s2->getPosition().Y();
			z1 = s2->getPosition().Z();
			t1 = s2->getTime();
			p1_mag = p1.M();

			vertexZ = kinfitVertexZ;
			vertexX = kinfitVertexX;
			vertexY = kinfitVertexY;
			//double thetadeg2, thetarad2, phideg2, phirad2, cospsi, psi;
			
			double thetarad2 = p1.Theta();
			double phirad2 = p1.Phi();
			double thetadeg2 = p1.Theta()*180.0/TMath::Pi();
			double phideg2 = p1.Phi()*180.0/TMath::Pi();					
			//double cospsi = sin(thetarad1)*sin(thetarad2)*cos(phirad1-phirad2)+cos(thetarad1)*cos(thetarad2);
			//double psi=acos(cospsi)*180/TMath::Pi();		
			

			
			for(unsigned int k=ii+1; k<locFCALShowers.size(); k++)
			{
			  //        if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[k]))
			  //    continue;
				if (find(matchedFCALShowers.begin(), matchedFCALShowers.end(),locFCALShowers[k]) != matchedFCALShowers.end()) continue;
				const DFCALShower *s3 = locFCALShowers[k];
				double dx2 = s3->getPosition().X() - kinfitVertexX;
				double dy2 = s3->getPosition().Y() - kinfitVertexY;
				double dz2 = s3->getPosition().Z()-kinfitVertexZ; // shift to coordinate relative to center of target
				double R2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
				 E2 = s3->getEnergy();
				double dt2 = s3->getTime();
				x2 = s3->getPosition().X();
				y2 = s3->getPosition().Y();
				z2 = s3->getPosition().Z();
				t2 = s3->getTime();
				
				
				TLorentzVector p2(E2*dx2/R2, E2*dy2/R2, E2*dz2/R2, E2);
				p2_mag = p2.M();
				double thetadeg3,thetarad3,phideg3,phirad3,cospsi3,psi3;
				int makes_sense3 = 0;
				if(R2 > R1 && dt2 > dt1) makes_sense3=1;
				if(R2 < R1 && dt2 < dt1) makes_sense3=1; 
				if(R2 > R1 && dt2 < dt1) makes_sense3=0;
				if(R2 < R1 && dt2 > dt1) makes_sense3=0;
				TLorentzVector ptot = p2+p1;
				inv_mass = ptot.M();
				thetadeg3 = p2.Theta()*180.0/TMath::Pi();
				thetarad3 = p2.Theta();
				phideg3 = p2.Phi()*180.0/TMath::Pi();
				phirad3 = p2.Phi();
				cospsi3 = sin(thetarad2)*sin(thetarad3)*cos(phirad2-phirad3)+cos(thetarad2)*cos(thetarad3);
				psi3=acos(cospsi3)*180/TMath::Pi();
				Theta_Hist_Both_FCAL->Fill(thetadeg3);
				Theta_Hist_Both_FCAL->Fill(thetadeg2);
				Phi_Hist_Both_FCAL->Fill(phideg3);
				Phi_Hist_Both_FCAL->Fill(phideg2);
				Psi_Hist_Both_FCAL->Fill(psi3);
				//if(E2>0.3 && E3>0.3 && kinfitVertexZ > 50 && kinfitVertexZ < 80)
				//double deltaT = TMath::Abs(dt1-dt2);
				//if(ptot.M() > 0.05 && ptot.M() < 1.0)Time_Diff->Fill(deltaT);
				//	cout << " E1 = " << E1 << " E2 = " << E2 << " d1 = " << dt1 << " t2 = " << dt2 << " ptot = " << ptot.M() << endl;
				FCAL_Neutrals->Fill();
			       if( kinfitVertexZ > 62.0 && kinfitVertexZ < 68.0 && E2 > 0.5 && E1 > 0.5 && kinfitVertexX < 15.0 && kinfitVertexX > -15.0 && kinfitVertexY > -15.0 && kinfitVertexY < 15.0) two_fcal_gamma_mass->Fill(ptot.M());
				//cout << " shower1 energy =  " << E1 << " shower2 energy = " << E2 << " ptot = " << ptot.M() << " fcal shower size = " << locFCALShowers.size() << endl;
				//if(E2>0.2 && E1>0.2) two_fcal_gamma_mass->Fill(ptot.M());
				
			}
		}


		for(unsigned int ii=0; ii<locFCALShowers.size(); ii++)
		{

		  //    if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[ii]))
		  //     continue;
			FCALShowers_per_event = locFCALShowers.size();
			if (find(matchedFCALShowers.begin(), matchedFCALShowers.end(),locFCALShowers[ii]) != matchedFCALShowers.end()) continue;
			const DFCALShower *s2 = locFCALShowers[ii];
			
			//double dx, dy, dz, R, thetarad1, phirad1,
			double dx1 = s2->getPosition().X() - kinfitVertexX;
			double dy1 = s2->getPosition().Y() - kinfitVertexY;
			double dz1 = s2->getPosition().Z() - kinfitVertexZ; // shift to coordinate relative to center of target
			double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
			 E1 = s2->getEnergy();
			double dt1 = s2->getTime();
			//double path1 = sqrt((dx-kinfitVertexX)*(dx-kinfitVertexX)+(dy-kinfitVertexY)*(dy-kinfitVertexY)+(dz)*(dz));g
			TLorentzVector p1(E1*dx1/R1, E1*dy1/R1, E1*dz1/R1, E1);
			x1 = s2->getPosition().X();
			y1 = s2->getPosition().Y();
			z1 = s2->getPosition().Z();
			t1 = s2->getTime();
			p1_mag = p1.M();

			vertexZ = kinfitVertexZ;
			vertexX = kinfitVertexX;
			vertexY = kinfitVertexY;
			//double thetadeg2, thetarad2, phideg2, phirad2, cospsi, psi;
			
			double thetarad2 = p1.Theta();
			double phirad2 = p1.Phi();
			double thetadeg2 = p1.Theta()*180.0/TMath::Pi();
			double phideg2 = p1.Phi()*180.0/TMath::Pi();					
			//double cospsi = sin(thetarad1)*sin(thetarad2)*cos(phirad1-phirad2)+cos(thetarad1)*cos(thetarad2);
			//double psi=acos(cospsi)*180/TMath::Pi();		
			

			
			for(unsigned int k=ii+1; k<locFCALShowers.size(); k++)
			{
			  //        if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[k]))
			  //    continue;
				if (find(matchedFCALShowers.begin(), matchedFCALShowers.end(),locFCALShowers[k]) != matchedFCALShowers.end()) continue;
				const DFCALShower *s3 = locFCALShowers[k];
				double dx2 = s3->getPosition().X() - kinfitVertexX;
				double dy2 = s3->getPosition().Y() - kinfitVertexY;
				double dz2 = s3->getPosition().Z() - kinfitVertexZ; // shift to coordinate relative to center of target
				double R2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
				 E2 = s3->getEnergy();
				double dt2 = s3->getTime();
				x2 = s3->getPosition().X();
				y2 = s3->getPosition().Y();
				z2 = s3->getPosition().Z();
				t2 = s3->getTime();
				
				
				TLorentzVector p2(E2*dx2/R2, E2*dy2/R2, E2*dz2/R2, E2);
				p2_mag = p2.M();
				double thetadeg3,thetarad3,phideg3,phirad3,cospsi3,psi3;
				int makes_sense3 = 0;
				if(R2 > R1 && dt2 > dt1) makes_sense3=1;
				if(R2 < R1 && dt2 < dt1) makes_sense3=1; 
				if(R2 > R1 && dt2 < dt1) makes_sense3=0;
				if(R2 < R1 && dt2 > dt1) makes_sense3=0;
				thetadeg3 = p2.Theta()*180.0/TMath::Pi();
				thetarad3 = p2.Theta();
				phideg3 = p2.Phi()*180.0/TMath::Pi();
				phirad3 = p2.Phi();
				cospsi3 = sin(thetarad2)*sin(thetarad3)*cos(phirad2-phirad3)+cos(thetarad2)*cos(thetarad3);
				psi3=acos(cospsi3)*180/TMath::Pi();
				Theta_Hist_Both_FCAL->Fill(thetadeg3);
				Theta_Hist_Both_FCAL->Fill(thetadeg2);
				Phi_Hist_Both_FCAL->Fill(phideg3);
				Phi_Hist_Both_FCAL->Fill(phideg2);
				Psi_Hist_Both_FCAL->Fill(psi3);
				for(unsigned int m = ii+2; m < locFCALShowers.size(); m ++)
				{
					if (find(matchedFCALShowers.begin(), matchedFCALShowers.end(),locFCALShowers[m]) != matchedFCALShowers.end()) continue;
					double dx3 = locFCALShowers[m]->getPosition().X() - kinfitVertexX;
					double dy3 = locFCALShowers[m]->getPosition().Y() - kinfitVertexY;
					double dz3 = locFCALShowers[m]->getPosition().Z() - kinfitVertexZ;
					double R3 = sqrt(dx3*dx3+dy3*dy3+dz3*dz3);
					E3 = locFCALShowers[m]->getEnergy();
					double dt3 = locFCALShowers[m]->getTime();
					x3 = locFCALShowers[m]->getPosition().X();
					y3 = locFCALShowers[m]->getPosition().Y();
					z3 = locFCALShowers[m]->getPosition().Z();
					t3 = locFCALShowers[m]->getTime();
					TLorentzVector p3(E3*dx3/R3, E3*dy3/R3, E3*dz3/R3, E3);
					p3_mag = p3.M();
					TLorentzVector ptot = p1 + p2 + p3;
					inv_mass = ptot.M();
					Triple_FCAL_Neutrals->Fill();
				}
				
			}
		}






	// time to look at clusters lmao

		for(unsigned int ii=0; ii<locFCALClusters.size(); ii++)
		{

		  //    if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[ii]))
		  //     continue;
			FCALClusters_per_event = locFCALClusters.size();
			if (find(matchedFCALClusters.begin(), matchedFCALClusters.end(),locFCALClusters[ii]) != matchedFCALClusters.end()) continue;
			const DFCALCluster *s2 = locFCALClusters[ii];
			
			//double dx, dy, dz, R, thetarad1, phirad1,
			double dx1 = s2->getCentroid().X() - kinfitVertexX;
			double dy1 = s2->getCentroid().Y() - kinfitVertexY;
			double dz1 = s2->getCentroid().Z() - kinfitVertexZ; // shift to coordinate relative to center of target
			double R1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
			 E1 = s2->getEnergy();
			double dt1 = s2->getTime();
			//double path1 = sqrt((dx-kinfitVertexX)*(dx-kinfitVertexX)+(dy-kinfitVertexY)*(dy-kinfitVertexY)+(dz)*(dz));g
			TLorentzVector p1(E1*dx1/R1, E1*dy1/R1, E1*dz1/R1, E1);
			x1 = s2->getCentroid().X();
			y1 = s2->getCentroid().Y();
			z1 = s2->getCentroid().Z();
			t1 = s2->getTime();
			p1_mag = p1.M();

			vertexZ = kinfitVertexZ;
			vertexX = kinfitVertexX;
			vertexY = kinfitVertexY;
			//double thetadeg2, thetarad2, phideg2, phirad2, cospsi, psi;
			
			double thetarad2 = p1.Theta();
			double phirad2 = p1.Phi();
			double thetadeg2 = p1.Theta()*180.0/TMath::Pi();
			double phideg2 = p1.Phi()*180.0/TMath::Pi();					
			//double cospsi = sin(thetarad1)*sin(thetarad2)*cos(phirad1-phirad2)+cos(thetarad1)*cos(thetarad2);
			//double psi=acos(cospsi)*180/TMath::Pi();		
			

			
			for(unsigned int k=ii+1; k<locFCALClusters.size(); k++)
			{
			  //        if(locDetectorMatches->Get_IsMatchedToTrack(locFCALShowers[k]))
			  //    continue;
				if (find(matchedFCALClusters.begin(), matchedFCALClusters.end(),locFCALClusters[k]) != matchedFCALClusters.end()) continue;
				const DFCALCluster *s3 = locFCALClusters[k];
				double dx2 = s3->getCentroid().X() - kinfitVertexX;
				double dy2 = s3->getCentroid().Y() - kinfitVertexY;
				double dz2 = s3->getCentroid().Z() - kinfitVertexZ; // shift to coordinate relative to center of target
				double R2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
				 E2 = s3->getEnergy();
				double dt2 = s3->getTime();
				x2 = s3->getCentroid().X();
				y2 = s3->getCentroid().Y();
				z2 = s3->getCentroid().Z();
				t2 = s3->getTime();
				
				
				TLorentzVector p2(E2*dx2/R2, E2*dy2/R2, E2*dz2/R2, E2);
				p2_mag = p2.M();
				double thetadeg3,thetarad3,phideg3,phirad3,cospsi3,psi3;
				int makes_sense3 = 0;
				if(R2 > R1 && dt2 > dt1) makes_sense3=1;
				if(R2 < R1 && dt2 < dt1) makes_sense3=1; 
				if(R2 > R1 && dt2 < dt1) makes_sense3=0;
				if(R2 < R1 && dt2 > dt1) makes_sense3=0;
				thetadeg3 = p2.Theta()*180.0/TMath::Pi();
				thetarad3 = p2.Theta();
				phideg3 = p2.Phi()*180.0/TMath::Pi();
				phirad3 = p2.Phi();
				cospsi3 = sin(thetarad2)*sin(thetarad3)*cos(phirad2-phirad3)+cos(thetarad2)*cos(thetarad3);
				psi3=acos(cospsi3)*180/TMath::Pi();
				Theta_Hist_Both_FCAL->Fill(thetadeg3);
				Theta_Hist_Both_FCAL->Fill(thetadeg2);
				Phi_Hist_Both_FCAL->Fill(phideg3);
				Phi_Hist_Both_FCAL->Fill(phideg2);
				Psi_Hist_Both_FCAL->Fill(psi3);
				TLorentzVector ptot = p1 + p2;
				inv_mass = ptot.M();
				FCALClusterNeutrals->Fill();
			}
		}












	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);
	//const DMCThrown* locMCThrown;
	//const DMCThrown* locMCThrown2




/*	for(unsigned int i = 0; i < locMCThrowns.size(); i++)
	{
		const DMCThrown *locMCThrown = locMCThrowns[i];
		Particle_t locPID; 
		locPID = (Particle_t) locMCThrown->type;
		if (locPID != Gamma) continue;
		int numZ = 0;
		double dz1;
		//if (locPID == Gamma)
	//	{
			 dz1 = locMCThrown->position().Z();
			//Thrown_Vertex->Fill(dz1);
	//	}
			for(unsigned int j=i+1; j < locMCThrowns.size(); j++)
			{
				Particle_t locPID2; 
				const DMCThrown *locMCThrown2 = locMCThrowns[j];
				locPID2 = (Particle_t) locMCThrown2->type;
				if (locPID2 != Gamma) continue;
				//if(locPID2==Gamma)
				//{
					double dz2 = locMCThrown2->position().Z();
					if(dz1 < (dz2 + 1) && dz1 > (dz2 - 1)) numZ = 1;
					else numZ = 2;
					//if(dz1 != dz2) numZ = 2;
					Thrown_Vertex_Frequency->Fill(numZ);
				//	cout << " thrown size = " << locMCThrowns.size() << " vertex of gamma 1 = " << dz1 << " vertex of gamma 2 = " << dz2 << " numz = " << numZ << endl;
				//}
			}
		
				
		
	} 



	Particle_t locPID2;
//=(Particle_t) locMCThrown->type;
	for(unsigned int i = 0; i < locMCThrowns.size(); i++)
	{
		const DMCThrown *locMCThrown = locMCThrowns[i];
		Particle_t locPID; 
		locPID = (Particle_t) locMCThrown->type;
		DVector3 locMomentum1;
		DLorentzVector lorentzmomentum1;
		TLorentzVector p1;
		if (locPID == Gamma)
		{
			DVector3 locMomentum1 = locMCThrown->momentum();
			DLorentzVector lorentzmomentum1 = locMCThrown->lorentzMomentum();
			double dx = locMCThrown->position().X();
			double dy = locMCThrown->position().Y();
			double dz = locMCThrown->position().Z();
			double R = sqrt(dx*dx + dy*dy + dz*dz);
			double E1 = locMCThrown->energy();
			TLorentzVector p1(E1*dx/R, E1*dy/R, E1*dz/R,E1);
			double locTheta = locMomentum1.Theta()*180.0/TMath::Pi();
			//Thrown_Vertex->Fill(dz);
			//cout << " vertex in z = " << dz << " locMCThrownSize = " << locMCThrowns.size() << endl;


		//}	
		for(unsigned int j =i+1; j < locMCThrowns.size(); j++)
		{
			const DMCThrown *locMCThrown2 = locMCThrowns[j];
			locPID2 = (Particle_t) locMCThrown2->type;
			//if (locPID2 == 1)
			//{
				double dx = locMCThrown2->position().X();
				double dy = locMCThrown2->position().Y();
				double dz = locMCThrown2->position().Z();
				double R = sqrt(dx*dx + dy*dy + dz*dz);
				double E2 = locMCThrown2->energy();
				TLorentzVector p2(E2*dx/R, E2*dy/R, E2*dz/R,E2);
				DVector3 locMomentum2 = locMCThrown2->momentum();
				double momentum2Mag = locMomentum2.Mag();
				DLorentzVector lorentzmomentum2 = locMCThrown2->lorentzMomentum();
				DVector3 locPSum = locMomentum1 + locMomentum2;
				DLorentzVector LorentzPSum = lorentzmomentum1 + lorentzmomentum2;
				TLorentzVector ptot = p1+p2;
				
				Thrown_Gamma_Theta->Fill(lorentzmomentum1.Theta()*180.0/TMath::Pi());
				Thrown_Gamma_Theta->Fill(lorentzmomentum2.Theta()*180.0/TMath::Pi());
				 Thrown_Inv_Mass->Fill(LorentzPSum.M());

			}
		}
	}
*/
	// Total energy in hits
	double Ehit_tot = 0.0;
	for(unsigned int i=0; i<bcalhits.size(); i++){
		Ehit_tot += bcalhits[i]->E;
	}
	Etot_hits->Fill(Ehit_tot);
	
	// Truth values
	double Etruth_tot = 0.0;
	double z_truth = 0.0;
	for(unsigned int i=0; i<truthshowers.size(); i++){
		Etruth_tot += truthshowers[i]->E;
		z_truth += truthshowers[i]->E*truthshowers[i]->z;
	}
	z_truth/=Etruth_tot;
	//Etot_truth->Fill(Etruth_tot);

	// Compare to thrown values
	double Etot_thrown=0.0;
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		Etot_thrown += mcthrowns[i]->energy();
		for(unsigned int j=0; j<locBCALShowers.size(); j++){
			double z = locBCALShowers[j]->z;
			//Erec_over_Ethrown_vs_z->Fill(z, locBCALShowers[j]->E/mcthrowns[i]->energy());

			double E = locBCALShowers[j]->E*(1.106+(z-208.4)*(z-208.4)*6.851E-6);
			//E_over_Erec_vs_z->Fill(z, mcthrowns[i]->energy()/E);
		}
	}
	
	//Ereconstructed_vs_Ethrown->Fill(Etot_thrown, Etot_reconstructed);
	//Etruth_over_Ethrown_vs_z->Fill(z_truth, Etruth_tot/Etot_thrown);

	// Single thrown particle
	if(mcthrowns.size()==1){
		const DMCThrown* mcthrown = mcthrowns[0];
		if(mcthrown->momentum().Theta()>0.0001){
			double z = mcthrown->position().Z() + 65.0/tan(mcthrown->momentum().Theta());
			double Ethrown = 1.0; // for some reason, mcthrown->E is zero right now.
			// I fudge this for now since I know all of the events thrw 1.0GeV
			//Edeposited_over_Ethrown_vs_z->Fill(z, Ehit_tot/Ethrown);
		}
	}

	UnlockState();	


	/*
	//Optional: Save event to output REST file. Use this to create skims.
	dEventWriterREST->Write_RESTEvent(locEventLoop, "BCAL_Shower"); //string is part of output file name
	*/

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_BCAL_Shower::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_BCAL_Shower::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}


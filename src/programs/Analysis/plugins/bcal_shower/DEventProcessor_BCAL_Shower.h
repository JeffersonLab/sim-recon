// $Id$
//
//    File: DEventProcessor_BCAL_Shower.h
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DEventProcessor_BCAL_Shower_
#define _DEventProcessor_BCAL_Shower_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DEventWriterROOT.h>
#include <HDDM/DEventWriterREST.h>
#include <ANALYSIS/DHistogramActions.h>
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"

#include "DLorentzVector.h"


using namespace jana;
using namespace std;

class DEventProcessor_BCAL_Shower : public jana::JEventProcessor
{
	public:
		DEventProcessor_BCAL_Shower(){};
		~DEventProcessor_BCAL_Shower(){};
		const char* className(void){return "DEventProcessor_BCAL_Shower";}
		DVector3 Calc_CrudeVertex(const deque<const DKinematicData*>& locParticles) const;
		TTree *BCALPoint_Charged_neg;
		TTree *BCALPoint_Charged_pos;
		TTree *BCALPoint_Neutral;
		TTree *BCAL_Neutrals;
		TTree *FCAL_Neutrals;
		TTree *FCALClusterNeutrals;
		TTree *Triple_FCAL_Neutrals;
		TTree *Split_Gamma_Neutrals;
		TTree *Split_Gamma_Neutrals_raw;
		//uint32_t shower_energy;
		Float_t energy_shower;
		Float_t energy_raw_shower;
		Float_t energy_point;
		Float_t track_momentum;
		Float_t layer1_energysum;
		Float_t layer2_energysum;
		Float_t layer3_energysum;
		Float_t layer4_energysum;
		//uint32_t hit_energy;
		uint32_t module;
		uint32_t layer;
		uint32_t sector;
		uint32_t eventnum;
		Float_t theta;
		uint32_t channel;
		Float_t phi;
		Float_t r;
		Float_t z;
		Float_t L1_pathlength;
		Float_t L2_pathlength;
		Float_t L3_pathlength;
		Float_t L4_pathlength;
		Float_t E1;
		Float_t E1_raw;
		Float_t E2;
		Float_t E2_raw;
		Float_t E3;
		Float_t t1;
		Float_t t2;
		Float_t t3;
		Float_t z1;
		Float_t z2;
		Float_t z3;
		Float_t x1;
		Float_t x2;
		Float_t x3;
		Float_t y1;
		Float_t y2;
		Float_t y3;
		Float_t phi1;
		Float_t phi2;
		Float_t vertexZ;
		Float_t vertexX;
		Float_t vertexY;
		uint32_t BCALShowers_per_event;
		uint32_t FCALShowers_per_event;
		uint32_t FCALClusters_per_event;
		Float_t p1_mag;
		Float_t p1_raw_mag;
		Float_t p2_mag;
		Float_t p2_raw_mag;
		Float_t p3_mag;
		Float_t inv_mass;
		Float_t inv_mass_raw;
		Float_t bcal_E;
		Float_t bcal_x;
		Float_t bcal_y;
		Float_t bcal_z;
		Float_t bcal_phi;
		Float_t bcal_t;
		Float_t bcal_p;
		Float_t fcal_E;
		Float_t fcal_x;
		Float_t fcal_y;
		Float_t fcal_z;
		Float_t fcal_t;
		Float_t fcal_p;
		
		//Float_t charge;
		//uint32_t end;


	private:
		const DAnalysisUtilities* dAnalysisUtilities;
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, int locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		TH1F* two_gamma_mass;
		TH1F* bcal_diphoton_mass;
		TH1F* bcal_diphoton_mass_highE;
		//TH1F* two_gamma_mass_corr,
		//TH1F* two_gamma_mass_cut;
		TH1F* bcal_fcal_two_gamma_mass;
		//TH1F* bcal_fcal_two_gamma_mass_cut;
	//	TH2F* xy_shower;
		TH1F* z_shower;
		TH1F* E_shower;
		TH1F* Neutral_E;
		TH1F* Assoc_E;
		TH1F* two_fcal_gamma_mass;
		TH1F* mass_summed;
		TH1F* Theta_Hist_Both_BCAL;
		TH1F* Phi_Hist_Both_BCAL;
		TH1F* Psi_Hist_Both_BCAL;
		TH1F* Theta_Hist_Split_Gammas;
		TH1F* Phi_Hist_Split_Gammas;
		TH1F* Psi_Hist_Split_Gammas;
		TH1F* Theta_Hist_Both_FCAL;
		TH1F* Phi_Hist_Both_FCAL;
		TH1F* Psi_Hist_Both_FCAL;
		TH1F* two_fcal_gamma_mass_test;
		TH1F* VertexZ;
		TH1F* goodVertexZ;
		TH1F* Thrown_Gamma_Theta;
		TH1F* Thrown_Inv_Mass;
		TH1F* Thrown_Vertex;
		TH1F* Thrown_Vertex_Frequency;
		TH2F* E1_Vs_E2CosTheta;
		TH1F* Pi0_velocity;
		TH1F* Cluster_Energy;
		TH1F* Point_Energy;
		TH1F* Point_Module;
		TH1F* Point_z;
		TH2F* Layer1_Energy_vs_Channel;
		TH2F* Layer2_Energy_vs_Channel;
		TH2F* Layer3_Energy_vs_Channel;
		TH2F* Layer4_Energy_vs_Channel;
		TH1F* Layer1_Energy;
		TH1F* Layer2_Energy;
		TH1F* Layer3_Energy;
		TH1F* Layer4_Energy;
		TH1F* Layer1_Energy_v2;
		TH1F* Layer2_Energy_v2;
		TH1F* Layer3_Energy_v2;
		TH1F* Layer4_Energy_v2;
		TH1F* Layer1_Energy_pos;
		TH1F* Layer2_Energy_pos;
		TH1F* Layer3_Energy_pos;
		TH1F* Layer4_Energy_pos;
		TH1F* BCALShowerTrack_Energy;
		TH1F* All_Layer_Energy;
		TH1F* Time_Diff;
		TH1F* Eoverp_plus_cuts;
		TH1F* Eoverp_minus_cuts;
		TH1F* Eoverp_plus_nocuts;
		TH1F* Eoverp_minus_nocuts;
		TH2F* Evsp_plus;
		TH2F* Evsp_minus;
		TH2F* Eoverpvsp_plus;
		TH2F* Eoverpvsp_minus;
		TH2F* FCAL_Evsp;
		TH2F* FCAL_Eoverpvsp;
		TH1F* FCAL_Eoverp_cuts;
		TH1F* FCAL_Eoverp_nocuts;
		TH1F* Point_E_M1S1L1;
		TH1F* Point_E_M12S2L2;
		TH1F* Point_E_M25S3L3;
		TH1F* Point_E_M37S4L4;
	//	TH2F* Erec_over_Ethrown_vs_z;
	//	TH2F* Ecorr_over_Erec_vs_z;
	//	TH2F* Ereconstructed_vs_Ethrown;
	//	TH1F* Etot_truth,
		TH1F* Etot_hits;
	//	TH2F* Etruth_over_Ethrown_vs_z;
	//	TH2F *Edeposited_over_Ethrown_vs_z;

		double dTargetZCenter;

		//DTrackFinder *finder;

		const DEventWriterROOT* dEventWriterROOT;
		const DEventWriterREST* dEventWriterREST;
};

#endif // _DEventProcessor_BCAL_Shower_


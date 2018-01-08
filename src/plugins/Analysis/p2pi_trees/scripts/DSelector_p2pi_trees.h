#ifndef DSelector_p2pi_trees_h
#define DSelector_p2pi_trees_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

class DSelector_p2pi_trees : public DSelector
{
	public:

		DSelector_p2pi_trees(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_p2pi_trees(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dProtonWrapper;
		DChargedTrackHypothesis* dPiPlusWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;
		TH1I* dHist_pMomentumMeasured;
		TH1I* dHist_piplusMomentumMeasured;
		TH1I* dHist_piminusMomentumMeasured;
		TH2I* dHist_Proton_dEdx_P;
		TH1I* dHist_KinFitChiSq;
		TH1I* dHist_KinFitCL;
		TH1I* dHist_M2pi;
		TH1I* dHist_t;
		TH2I* dHist_CosTheta_Psi;;
		TH1I* dHist_phi;
		TH1I* dHist_Phi;
		TH1I* dHist_psi;
		TH1I* dHist_pDeltap;
		TH1I* dHist_pipDeltap;
		TH1I* dHist_pimDeltap;
		TH1I* dHist_pDeltap_Measured;
		TH1I* dHist_pipDeltap_Measured;
		TH1I* dHist_pimDeltap_Measured;

		// Cut parameters
		TF1* fMinProton_dEdx;
		TF1* fMaxPion_dEdx;
		Double_t dMinKinFitCL;
		Double_t dMaxKinFitChiSq;
		Double_t dMinBeamEnergy;
		Double_t dMaxBeamEnergy;
		Double_t dMin2piMass;
		Double_t dMax2piMass;
		Double_t dMinMissingMassSquared;
		Double_t dMaxMissingMassSquared;

	ClassDef(DSelector_p2pi_trees, 0);
};

void DSelector_p2pi_trees::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
}

#endif // DSelector_p2pi_trees_h

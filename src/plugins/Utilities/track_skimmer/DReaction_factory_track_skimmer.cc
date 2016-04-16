// $Id$
//
//    File: DReaction_factory_track_skimmer.cc
// Created: Tue Jan 13 11:08:16 EST 2015
// Creator: Paul (on Darwin Pauls-MacBook-Pro.local 14.0.0 i386)
//

#include "DReaction_factory_track_skimmer.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_track_skimmer::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = NULL; //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** skim_pi0 ****************************************************/

	locReaction = new DReaction("skim_pi0");

	//pi0 -> g, g
	DReactionStep* locReactionStep_Pi01 = new DReactionStep();
	locReactionStep_Pi01->Set_InitialParticleID(Pi0);
	locReactionStep_Pi01->Add_FinalParticleID(Gamma);
	locReactionStep_Pi01->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep_Pi01);
	dReactionStepPool.push_back(locReactionStep_Pi01); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_2pi0 ****************************************************/

	locReaction = new DReaction("skim_2pi0");

	//X -> 2pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g

	//pi0 -> g, g
	DReactionStep* locReactionStep_Pi02 = new DReactionStep();
	locReactionStep_Pi02->Set_InitialParticleID(Pi0);
	locReactionStep_Pi02->Add_FinalParticleID(Gamma);
	locReactionStep_Pi02->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep_Pi02);
	dReactionStepPool.push_back(locReactionStep_Pi02); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_3pi0 ****************************************************/

	locReaction = new DReaction("skim_3pi0");

	//X -> 3pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi02); //pi0 -> g, g

	//pi0 -> g, g
	DReactionStep* locReactionStep_Pi03 = new DReactionStep();
	locReactionStep_Pi03->Set_InitialParticleID(Pi0);
	locReactionStep_Pi03->Add_FinalParticleID(Gamma);
	locReactionStep_Pi03->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep_Pi03);
	dReactionStepPool.push_back(locReactionStep_Pi03); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_ks_2piq ****************************************************/

	locReaction = new DReaction("skim_ks_2piq");

	//ks -> pi+, pi-
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(KShort);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(KShort, 0.3, 0.7);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, KShort, false, 800, 0.3, 0.7, "KShort"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_ks_2pi0 ****************************************************/

	locReaction = new DReaction("skim_ks_2pi0");

	//ks -> 2pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(KShort);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi02); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(KShort, 0.3, 0.7);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, KShort, false, 800, 0.3, 0.7, "KShort"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_eta_2g ****************************************************/

	locReaction = new DReaction("skim_eta_2g");

	//eta -> g, g
	DReactionStep* locReactionStep_Eta2g = new DReactionStep();
	locReactionStep_Eta2g->Set_InitialParticleID(Eta);
	locReactionStep_Eta2g->Add_FinalParticleID(Gamma);
	locReactionStep_Eta2g->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep_Eta2g);
	dReactionStepPool.push_back(locReactionStep_Eta2g); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_eta_3pi0 ****************************************************/

	locReaction = new DReaction("skim_eta_3pi0");

	//eta -> 3pi0
	DReactionStep* locReactionStep_Eta3Pi0 = new DReactionStep();
	locReactionStep_Eta3Pi0->Set_InitialParticleID(Eta);
	locReactionStep_Eta3Pi0->Add_FinalParticleID(Pi0);
	locReactionStep_Eta3Pi0->Add_FinalParticleID(Pi0);
	locReactionStep_Eta3Pi0->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep_Eta3Pi0);
	dReactionStepPool.push_back(locReactionStep_Eta3Pi0); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi02); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi03); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_eta_3piq ****************************************************/

	locReaction = new DReaction("skim_eta_3piq");

	//eta -> pi+, pi-, pi0
	DReactionStep* locReactionStep_Eta3Piq = new DReactionStep();
	locReactionStep_Eta3Piq->Set_InitialParticleID(Eta);
	locReactionStep_Eta3Piq->Add_FinalParticleID(PiPlus);
	locReactionStep_Eta3Piq->Add_FinalParticleID(PiMinus);
	locReactionStep_Eta3Piq->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep_Eta3Piq);
	dReactionStepPool.push_back(locReactionStep_Eta3Piq); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_omega_3piq ****************************************************/

	locReaction = new DReaction("skim_omega_3piq");

	//omega -> pi+, pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(omega);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(omega, 0.4, 1.2);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 1600, 0.4, 1.2, "Omega"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_omega_pi0g ****************************************************/

	locReaction = new DReaction("skim_omega_pi0g");

	//omega -> pi0, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(omega);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(omega, 0.4, 1.2);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 1600, 0.4, 1.2, "Omega"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_etaprm_2piqeta_2g ****************************************************/

	locReaction = new DReaction("skim_etaprm_2piqeta_2g");

	//eta' -> pi+, pi-, eta
	DReactionStep* locReactionStep_EtaPrime2PiqEta = new DReactionStep();
	locReactionStep_EtaPrime2PiqEta->Set_InitialParticleID(EtaPrime);
	locReactionStep_EtaPrime2PiqEta->Add_FinalParticleID(PiPlus);
	locReactionStep_EtaPrime2PiqEta->Add_FinalParticleID(PiMinus);
	locReactionStep_EtaPrime2PiqEta->Add_FinalParticleID(Eta);
	locReaction->Add_ReactionStep(locReactionStep_EtaPrime2PiqEta);
	dReactionStepPool.push_back(locReactionStep_EtaPrime2PiqEta); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Eta2g); //eta -> g, g

	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Set_InvariantMassCut(EtaPrime, 0.6, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, EtaPrime, false, 1400, 0.6, 1.3, "EtaPrime"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_etaprm_2piqeta_3pi0 ****************************************************/

	locReaction = new DReaction("skim_etaprm_2piqeta_3pi0");

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_EtaPrime2PiqEta); //eta' -> pi+, pi-, eta
	locReaction->Add_ReactionStep(locReactionStep_Eta3Pi0); //eta -> 3pi0
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi02); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi03); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Set_InvariantMassCut(EtaPrime, 0.6, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, EtaPrime, false, 1400, 0.6, 1.3, "EtaPrime"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_etaprm_2piqeta_3piq ****************************************************/

	locReaction = new DReaction("skim_etaprm_2piqeta_3piq");

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_EtaPrime2PiqEta); //eta' -> pi+, pi-, eta
	locReaction->Add_ReactionStep(locReactionStep_Eta3Piq); //eta -> pi+, pi-, pi0
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Set_InvariantMassCut(EtaPrime, 0.6, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, EtaPrime, false, 1400, 0.6, 1.3, "EtaPrime"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_etaprm_2piqg ****************************************************/

	locReaction = new DReaction("skim_etaprm_2piqg");

	//eta' -> pi+, pi-, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(EtaPrime);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(EtaPrime, 0.6, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, EtaPrime, false, 1400, 0.6, 1.3, "EtaPrime"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_etaprm_2pi0eta_2g ****************************************************/

	locReaction = new DReaction("skim_etaprm_2pi0eta_2g");

	//eta' -> 2pi0, eta
	DReactionStep* locReactionStep_EtaPrime2Pi0Eta = new DReactionStep();
	locReactionStep_EtaPrime2Pi0Eta->Set_InitialParticleID(EtaPrime);
	locReactionStep_EtaPrime2Pi0Eta->Add_FinalParticleID(Pi0);
	locReactionStep_EtaPrime2Pi0Eta->Add_FinalParticleID(Pi0);
	locReactionStep_EtaPrime2Pi0Eta->Add_FinalParticleID(Eta);
	locReaction->Add_ReactionStep(locReactionStep_EtaPrime2Pi0Eta);
	dReactionStepPool.push_back(locReactionStep_EtaPrime2Pi0Eta); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Eta2g); //eta -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi02); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Set_InvariantMassCut(EtaPrime, 0.6, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, EtaPrime, false, 1400, 0.6, 1.3, "EtaPrime"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_etaprm_2pi0eta_3pi0 ****************************************************/

	locReaction = new DReaction("skim_etaprm_2pi0eta_3pi0");

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_EtaPrime2Pi0Eta); //eta' -> 2pi0, eta
	locReaction->Add_ReactionStep(locReactionStep_Eta3Pi0); //eta -> 3pi0
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi02); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi03); //pi0 -> g, g

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Set_InvariantMassCut(EtaPrime, 0.6, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, EtaPrime, false, 1400, 0.6, 1.3, "EtaPrime"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_etaprm_2pi0eta_3piq ****************************************************/

	locReaction = new DReaction("skim_etaprm_2pi0eta_3piq");

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_EtaPrime2Pi0Eta); //eta' -> 2pi0, eta
	locReaction->Add_ReactionStep(locReactionStep_Eta3Piq); //eta -> pi+, pi-, pi0
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi02); //pi0 -> g, g
	locReaction->Add_ReactionStep(locReactionStep_Pi03); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(Eta, 0.3, 0.8);
	locReaction->Set_InvariantMassCut(EtaPrime, 0.6, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Eta, false, 1000, 0.3, 0.8, "Eta"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, EtaPrime, false, 1400, 0.6, 1.3, "EtaPrime"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_phi_2kq ****************************************************/

	locReaction = new DReaction("skim_phi_2kq");

	//phi -> K+, K-
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(phiMeson);
	locReactionStep->Add_FinalParticleID(KPlus);
	locReactionStep->Add_FinalParticleID(KMinus);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(phiMeson, 0.8, 1.2);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, phiMeson, false, 800, 0.8, 1.2, "Phi"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_phi_3piq ****************************************************/

	locReaction = new DReaction("skim_phi_3piq");

	//phi -> pi+, pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(phiMeson);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(phiMeson, 0.8, 1.2);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, phiMeson, false, 800, 0.8, 1.2, "Phi"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_lambda ****************************************************/

	locReaction = new DReaction("skim_lambda");

	//lambda -> p, pi-
	DReactionStep* locReactionStep_Lambda = new DReactionStep();
	locReactionStep_Lambda->Set_InitialParticleID(Lambda);
	locReactionStep_Lambda->Add_FinalParticleID(Proton);
	locReactionStep_Lambda->Add_FinalParticleID(PiMinus);
	locReaction->Add_ReactionStep(locReactionStep_Lambda);
	dReactionStepPool.push_back(locReactionStep_Lambda); //register so will be deleted later: prevent memory leak

	locReaction->Set_InvariantMassCut(Lambda, 1.0, 1.2);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Lambda, false, 800, 1.0, 1.2, "Lambda"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_sigma0 ****************************************************/

	locReaction = new DReaction("skim_sigma0");

	//sigma0 -> lambda, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Sigma0);
	locReactionStep->Add_FinalParticleID(Lambda);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Lambda); //lambda -> p, pi-

	locReaction->Set_InvariantMassCut(Lambda, 1.0, 1.2);
	locReaction->Set_InvariantMassCut(Sigma0, 1.1, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Lambda, false, 800, 1.0, 1.2, "Lambda"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Sigma0, false, 800, 1.1, 1.3, "Sigma0"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_sigma+ ****************************************************/

	locReaction = new DReaction("skim_sigma+");

	//sigma+ -> p, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(SigmaPlus);
	locReactionStep->Add_FinalParticleID(Proton);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(SigmaPlus, 1.1, 1.3);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, SigmaPlus, false, 800, 1.1, 1.3, "Sigma+"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_xi- ****************************************************/

	locReaction = new DReaction("skim_xi-");

	//xi- -> lambda, pi-
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(XiMinus);
	locReactionStep->Add_FinalParticleID(Lambda);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Lambda); //lambda -> p, pi-

	locReaction->Set_InvariantMassCut(Lambda, 1.0, 1.2);
	locReaction->Set_InvariantMassCut(XiMinus, 1.1, 1.5);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Lambda, false, 800, 1.0, 1.2, "Lambda"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, XiMinus, false, 800, 1.1, 1.5, "Xi-"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** skim_xi0 ****************************************************/

	locReaction = new DReaction("skim_xi0");

	//xi0 -> lambda, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Xi0);
	locReactionStep->Add_FinalParticleID(Lambda);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Reuse steps: Save time & memory //However, don't reuse within the SAME DReaction!
	locReaction->Add_ReactionStep(locReactionStep_Lambda); //lambda -> p, pi-
	locReaction->Add_ReactionStep(locReactionStep_Pi01); //pi0 -> g, g

	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(Lambda, 1.0, 1.2);
	locReaction->Set_InvariantMassCut(Xi0, 1.1, 1.5);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 510, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Lambda, false, 800, 1.0, 1.2, "Lambda"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Xi0, false, 800, 1.1, 1.5, "Xi0"));

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_track_skimmer::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}


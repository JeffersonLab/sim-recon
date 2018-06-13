// $Id$
//
//    File: DReaction_factory_B3pi_eff_missgamma.cc
// Created: Fri Jun 30 00:38:22 EDT 2017
// Creator: jmhardin (on Linux ifarm1402.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//


#include "DReaction_factory_B3pi_eff_missgamma.h"
//#include "DCustomAction_dEdxCut.h"

//------------------
// brun
//------------------
jerror_t DReaction_factory_B3pi_eff_missgamma::brun(JEventLoop* locEventLoop, int32_t locRunNumber)
{
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_B3pi_eff_missgamma::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = NULL; //create with a unique name for each DReaction object. CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/************************************************** B3pi_eff_missgamma Reaction Definition *************************************************/

	locReaction = new DReaction("B3pi_eff_missgamma");
        double pull_hist_confidence_level = 0.05;

        //DKinFitType locKinFitType = d_P4Fit;
        //DKinFitType locKinFitType = d_VertexFit;
        DKinFitType locKinFitType = d_P4AndVertexFit;
        //Required: DReactionSteps to specify the channel and decay chain you want to study
                //Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

        std::deque<Particle_t> off_proton;
        off_proton.push_back(Proton);

/*
        // g, p -> omega, p                                                                   
        locReactionStep = new DReactionStep();
        locReactionStep->Set_InitialParticleID(Gamma);
        locReactionStep->Set_TargetParticleID(Proton);
        locReactionStep->Add_FinalParticleID(omega);
        locReactionStep->Add_FinalParticleID(Proton);
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep);
        //register so will be deleted later: prevent memory leak                                                                             

        // omega -> pi+, pi-, pi0                                                             
        locReactionStep = new DReactionStep();
        locReactionStep->Set_InitialParticleID(omega);
        locReactionStep->Set_KinFitConstrainInitMassFlag(false);
        locReactionStep->Add_FinalParticleID(PiPlus);
        locReactionStep->Add_FinalParticleID(PiMinus);
        locReactionStep->Add_FinalParticleID(Pi0);
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep);
        //register so will be deleted later: prevent memory leak                       
        // pi0 -> g, g                                                                        
        locReactionStep = new DReactionStep();
        locReactionStep->Set_InitialParticleID(Pi0);
        locReactionStep->Set_KinFitConstrainInitMassFlag(true);
        locReactionStep->Add_FinalParticleID(Gamma,true);
        locReactionStep->Add_FinalParticleID(Gamma);
        //locReactionStep->Set_KinFitConstrainInitMassFlag(false);                            
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep);
*/
        vector<Particle_t> step0_final;
        step0_final.push_back(omega);
        step0_final.push_back(Proton);
        vector<Particle_t> step1_final;
        step1_final.push_back(PiPlus);
        step1_final.push_back(PiMinus);
        step1_final.push_back(Pi0);
        vector<Particle_t> step2_final;
        step2_final.push_back(Gamma);

        locReactionStep = new DReactionStep(Gamma,Proton,step0_final);
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep);

        locReactionStep = new DReactionStep(omega,step1_final);
        //locReactionStep = new DReactionStep(omega,step1_final);       
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep);

        locReactionStep = new DReactionStep(Pi0,step2_final,Gamma);
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep);



        /**************************************************** pippimpi0_withmiss Control Settings ****************************************************/

        // Highly Recommended: Set EventStore skim query (use with "eventstore" source)
                // This will skip creating particle combos for events that aren't in the skims you list
                // Query should be comma-separated list of skims to boolean-AND together

        // Recommended: Type of kinematic fit to perform (default is d_NoFit)
                //fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DReaction.h
                //Options: d_NoFit (default), d_P4Fit, d_VertexFit, d_P4AndVertexFit
                //P4 fits automatically constrain decaying particle masses, unless they are manually disabled
        //locReaction->Set_KinFitType(d_P4AndVertexFit);
        locReaction->Set_KinFitType(locKinFitType);
        locReaction->Set_EventStoreSkims("q+, q-"); //boolean-AND of skims

        // Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
       //Dep  locReaction->Set_MaxPhotonRFDeltaT(3.5*dBeamBunchPeriod); //should be minimum cut value
       locReaction->Set_NumPlusMinusRFBunches(2); //should be minimum cut value

        // Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" (or user custom) factory.
                // Important: Keep cut large: Can have many ghost and accidental tracks that look "good"
        locReaction->Set_MaxExtraGoodTracks(4);
        // Highly Recommended: Enable ROOT TTree output for this DReaction
        // string is file name (must end in ".root"!!): doen't need to be unique, feel free to change
        locReaction->Enable_TTreeOutput("tree_pippimpi0_eff_missgamma.root", true); //true/false: do/don't save unused hypotheses

        /************************************************** pippimpi0_withmiss Pre-Combo Custom Cuts *************************************************/

        // Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
        // Example: pi0 -> g, g cut


        /**************************************************** B3pi Analysis Actions ****************************************************/

//        locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.06, 0.06));

        //Dep locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.06, 0.06));
        locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.06, 0.06));

        //Cris's PID
        // HISTOGRAM PID
        locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

        //locReaction->Add_AnalysisAction(new DCutAction_BeamEnergy(locReaction,false,3.0,12.0));
        locReaction->Add_AnalysisAction(new DCutAction_BeamEnergy(locReaction,false,8.4,9.05));

        // PID                                                                                                                                                 
        locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
        //locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut(locReaction, false)); //false: focus on keeping signal                                     

        locReaction->Add_AnalysisAction(new DCutAction_dEdx(locReaction)); //false: measured data                            
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_TOF)); //false: measured data                            
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL)); //false: measured data                                   
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_FCAL)); //false: measured data                                   
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_TOF)); //false: measured data                                    
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, PiPlus, SYS_BCAL)); //false: measured data                                   
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, PiPlus, SYS_FCAL)); //false: measured data                                   
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Gamma, SYS_BCAL)); //false: measured data                                    
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Gamma, SYS_FCAL)); //false: measured data                                    
        locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

        //KINEMATIC FIT
        locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, pull_hist_confidence_level, true));
        locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, -1.0)); //require kinematic fit to converge

        locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

        // Pi0                                                                                                                                                         
        locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 600, 0.0, 0.3, "Pi0")); //false: measured data  

        // Missing Mass Squared (Hist and Cut)                                                                                                                    
        locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 600, -0.06, 0.06));
        locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -1, 1, "FullRange"));



        // Omega Mass (Hist and Cut)                                                                                                                                   
//        locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 1000, 0., 2.0, "Omega"));  //false: measured data                
//        locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 1000, 0., 2.0, "Omega_Kinfit")); //true: kinfit                   
//        locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, omega, off_proton, false, 1000, 0., 2.0, "OffProt")); //true: kinfit                   
//        locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, omega, off_proton, true, 1000, 0., 2.0, "OffProt_Kinfit")); //true: kinfit                   
        //locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, true, 0.69, 0.88)); //~ +/- 3sigma-ish //true: kinfit        
        // Kinematics of final selection                                                                                                                               
        locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final"));
        //false: fill histograms with measured particle data                                                                                                                                        

        _data.push_back(locReaction); //Register the DReaction with the factory                                                                                        


	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_B3pi_eff_missgamma::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}


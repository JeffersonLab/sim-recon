#include "DSelector_p2pi_trees.h"

void DSelector_p2pi_trees::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "DSelector_p2pi_trees.root"; //"" for none
	dOutputTreeFileName = "tree_DSelector_p2pi_trees.root"; //"" for none
	//dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	//dFlatTreeName = ""; //if blank, default name will be chosen

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));//false: use measured data
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "KinFit")); //true: use kinfit data;
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 2.0, Unknown, SYS_BCAL));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.4, 9.05));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.24, 0.24);
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
	dHist_pMomentumMeasured = new TH1I("pMomentumMeasured", ";p Momentum Measured (GeV)", 100, 0.0, 2);
	dHist_piplusMomentumMeasured = new TH1I("piplusMomentumMeasured", ";#pi^{+} Momentum Measured (GeV)", 600, 0.0, 12);
	dHist_piminusMomentumMeasured = new TH1I("piminusMomentumMeasured", ";#pi^{-} Momentum Measured (GeV)", 600, 0.0, 12);

	dHist_Proton_dEdx_P = new TH2I("Proton_dEdx_P", " ;p_{proton} GeV/c; dE/dx (keV/cm)", 250, 0.0, 5.0, 250, 0.0, 25.);
	dHist_KinFitChiSq = new TH1I("KinFitChiSq", ";Kinematic Fit #chi^{2}/NDF", 250, 0., 25.);
	dHist_KinFitCL = new TH1I("KinFitCL", ";Kinematic Fit Confidence Level", 100, 0., 1.);

	dHist_M2pi = new TH1I("M2pi", ";M_{#pi^{+}#pi^{-}} (GeV/c^{2})", 600, 0, 1.2);
	dHist_t = new TH1I("t", ";|t| (GeV/c)^{2}", 100, 0.0, 2.0);
	dHist_CosTheta_Psi = new TH2I("CosTheta_Psi", "; #psi; cos#theta", 360, -180., 180, 200, -1., 1.);
	dHist_phi = new TH1I("phi", ";phi (degrees)", 360,-180,180);
	dHist_Phi = new TH1I("Phi", ";Phi (degrees)", 360,-180,180);
	dHist_psi = new TH1I("psi", ";psi (degrees)", 360,-180,180);

	dHist_pDeltap = new TH1I("pDeltap","; proton: Thrown Fit p - KinFit p/Thrown p",100,-0.2,0.2);
	dHist_pipDeltap = new TH1I("pipDeltap","; #pi^{+}: Thrown Fit p - KinFit p/Thrown p",100,-0.2,0.2);
	dHist_pimDeltap = new TH1I("pimDeltap","; #pi^{-}: Thrown Fit p - KinFit p/ Thrown p",100,-0.2,0.2);

	dHist_pDeltap_Measured = new TH1I("pDeltap_Measured","; proton: Thrown Fit p - Measured p/Thrown p",100,-0.2,0.2);
	dHist_pipDeltap_Measured = new TH1I("pipDeltap_Measured","; #pi^{+}: Thrown Fit p - Measured p/Thrown p",100,-0.2,0.2);
	dHist_pimDeltap_Measured = new TH1I("pimDeltap_Measured","; #pi^{-}: Thrown Fit p - Measured p/ Thrown p",100,-0.2,0.2);

	// EXAMPLE CUT PARAMETERS:
	fMinProton_dEdx = new TF1("fMinProton_dEdx", "exp(-1.*[0]*x + [1]) + [2]", 0., 10.);
	fMinProton_dEdx->SetParameters(4.0, 2.5, 1.25);
	fMaxPion_dEdx = new TF1("fMaxPion_dEdx", "exp(-1.*[0]*x + [1]) + [2]", 0., 10.);
	fMaxPion_dEdx->SetParameters(4.0, 2.0, 2.5);
	dMinKinFitCL = 5.73303e-7; //5.73303e-7;
	dMaxKinFitChiSq = 5.0;
	dMinBeamEnergy = 8.4;
	dMaxBeamEnergy = 9.0;
	dMin2piMass = 0.5;
	dMax2piMass = 1.1;
	// dMinMissingMassSquared = -0.04;
	// dMaxMissingMassSquared = 0.04;
	dMinMissingMassSquared = -1;
	dMaxMissingMassSquared = 1;



	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_p2pi_trees::Process(Long64_t locEntry)
{
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search
	set<Int_t> locUsedSoFar_Proton, locUsedSoFar_PiPlus, locUsedSoFar_PiMinus;

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass, locUsedSoFar_2pi, locUsedSoFar_Angles;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/


	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEbeam_Thrown = dThrownBeam->Get_P4().E();

	//Loop over throwns
	TLorentzVector locProtonP4_Thrown;
	TLorentzVector locPiPlusP4_Thrown;
	TLorentzVector locPiMinusP4_Thrown;

	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);
		Particle_t thrown_pid = dThrownWrapper->Get_PID();
		// cout << " loc_i=" << loc_i << " thrown_pid=" << thrown_pid << endl;
                TLorentzVector locP4_Thrown = dThrownWrapper->Get_P4();
		if (loc_i == 0) locProtonP4_Thrown = locP4_Thrown;    // assume order of particles as PID is zero at the moment
		if (loc_i == 1) locPiPlusP4_Thrown = locP4_Thrown;
		if (loc_i == 2) locPiMinusP4_Thrown = locP4_Thrown;
		
	}
	/*locProtonP4_Thrown.Print();
	locPiPlusP4_Thrown.Print();
	locPiMinusP4_Thrown.Print();
	cout << endl;
	*/



	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();
		Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locProtonP4_Measured + locPiPlusP4_Measured + locPiMinusP4_Measured;
		TLorentzVector loc2piP4 = locPiPlusP4 + locPiMinusP4;


	cout << "MM2 =" << locMissingP4_Measured.M2() << endl;



		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/


		/**************************************** EXAMPLE: PID CUT ACTION ************************************************/


		// Proton CDC dE/dx histogram and cut 
		double locProton_dEdx_CDC = dProtonWrapper->Get_dEdx_CDC()*1e6;
		if(locProton_dEdx_CDC < fMinProton_dEdx->Eval(locProtonP4.P())) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		
		// Pi+/- CDC dE/dx histogram cut (histograms in HistComboPID action)
		double locPiPlus_dEdx_CDC = dPiPlusWrapper->Get_dEdx_CDC()*1e6;
		double locPiMinus_dEdx_CDC = dPiMinusWrapper->Get_dEdx_CDC()*1e6;
		if(locPiPlus_dEdx_CDC > fMaxPion_dEdx->Eval(locPiPlusP4.P()) || locPiMinus_dEdx_CDC > fMaxPion_dEdx->Eval(locPiMinusP4.P())) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}


		/************************************ EXAMPLE: SELECTION CUTS AND HISTOGRAMS  ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinusTrackID);

		// beam energy cut for SDME
		if(locBeamP4.E() < dMinBeamEnergy || locBeamP4.E() > dMaxBeamEnergy) {
                        dComboWrapper->Set_IsComboCut(true);
                        continue;
                }

		// Cut on Missing mass
		if((locMissingMassSquared < dMinMissingMassSquared ) || (locMissingMassSquared >  dMaxMissingMassSquared)){
			dComboWrapper->Set_IsComboCut(true);    
			continue; 
		}

		// kinematic fit CL cut
		dHist_KinFitChiSq->Fill(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit());
		dHist_KinFitCL->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit());
		if(dComboWrapper->Get_ConfidenceLevel_KinFit() <= dMinKinFitCL) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

		// rho mass histogram and cut
		map<Particle_t, set<Int_t> > locUsedThisCombo_2piMass;
		locUsedThisCombo_2piMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_2piMass[PiMinus].insert(locPiMinusTrackID);
		if(loc2piP4.M() < dMin2piMass || loc2piP4.M() > dMax2piMass) {
			dComboWrapper->Set_IsComboCut(true);
                        continue;
		} 
		if(locUsedSoFar_2pi.find(locUsedThisCombo_2piMass) == locUsedSoFar_2pi.end())
		{
			dHist_M2pi->Fill(loc2piP4.M());
			locUsedSoFar_2pi.insert(locUsedThisCombo_2piMass);
		}


		/**************************************** EXAMPLE: HISTOGRAM energies and momenta *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E());
			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared->Fill(locMissingMassSquared);
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}


		if(locUsedSoFar_Proton.find(locProtonTrackID) == locUsedSoFar_Proton.end())
		{
		        dHist_pMomentumMeasured->Fill(locProtonP4.Vect().Mag());
			dHist_Proton_dEdx_P->Fill(locProtonP4.P(), locProton_dEdx_CDC);
		        dHist_pDeltap->Fill(locProtonP4_Thrown.Vect().Mag() > 0? (locProtonP4_Thrown.Vect().Mag()-locProtonP4.Vect().Mag())/locProtonP4_Thrown.Vect().Mag() : 0);
		        dHist_pDeltap_Measured->Fill(locProtonP4_Thrown.Vect().Mag() > 0? (locProtonP4_Thrown.Vect().Mag()-locProtonP4_Measured.Vect().Mag())/locProtonP4_Thrown.Vect().Mag() : 0);
			locUsedSoFar_Proton.insert(locProtonTrackID);
		}

		if(locUsedSoFar_PiPlus.find(locPiPlusTrackID) == locUsedSoFar_PiPlus.end())
		{
		        dHist_piplusMomentumMeasured->Fill(locPiPlusP4.Vect().Mag());
		        dHist_pipDeltap->Fill(locPiPlusP4_Thrown.Vect().Mag() > 0? (locPiPlusP4_Thrown.Vect().Mag()-locPiPlusP4.Vect().Mag())/locPiPlusP4_Thrown.Vect().Mag() : 0);
		        dHist_pipDeltap_Measured->Fill(locPiPlusP4_Thrown.Vect().Mag() > 0? (locPiPlusP4_Thrown.Vect().Mag()-locPiPlusP4_Measured.Vect().Mag())/locPiPlusP4_Thrown.Vect().Mag() : 0);
			locUsedSoFar_PiPlus.insert(locPiPlusTrackID);
		}

		if(locUsedSoFar_PiMinus.find(locPiMinusTrackID) == locUsedSoFar_PiMinus.end())
		{
		        dHist_piminusMomentumMeasured->Fill(locPiMinusP4.Vect().Mag());
		        dHist_pimDeltap->Fill(locPiMinusP4_Thrown.Vect().Mag() > 0? (locPiMinusP4_Thrown.Vect().Mag()-locPiMinusP4.Vect().Mag())/locPiMinusP4_Thrown.Vect().Mag() : 0);
		        dHist_pimDeltap_Measured->Fill(locPiMinusP4_Thrown.Vect().Mag() > 0? (locPiMinusP4_Thrown.Vect().Mag()-locPiMinusP4_Measured.Vect().Mag())/locPiMinusP4_Thrown.Vect().Mag() : 0);
			locUsedSoFar_PiMinus.insert(locPiMinusTrackID);
		}

		cout << "Ebeam Thrown P4.E=" << dThrownBeam->Get_P4().E()  << " Kinfit P4.E=" << locBeamP4.E() << " P4 Measured =" << locBeamP4_Measured.E() << endl;
		cout << "Proton Thrown P4.E=" << locProtonP4_Thrown.E() << " Kinfit P4.E=" << locProtonP4.E() << " P4 Measured =" << locProtonP4_Measured.E() << " CL=" << dComboWrapper->Get_ConfidenceLevel_KinFit() << endl;
		cout << "PiPlus Thrown P4.E=" << locPiPlusP4_Thrown.E() << " Kinfit P4.E=" << locPiPlusP4.E() << " P4 Measured =" << locPiPlusP4_Measured.E() << endl;
		cout << "PiMinus Thrown P4.E=" << locPiMinusP4_Thrown.E() << " Kinfit P4.E=" << locPiMinusP4.E() << " P4 Measured =" << locPiMinusP4_Measured.E() << endl << endl;

		// calculate kinematic and angular variables
		double t = (locProtonP4 - dTargetP4).M2();
		TLorentzRotation resonanceBoost( -loc2piP4.BoostVector() );   // boost into 2pi frame
		TLorentzVector beam_res = resonanceBoost * locBeamP4;
		TLorentzVector recoil_res = resonanceBoost * locProtonP4;
		TLorentzVector p1_res = resonanceBoost * locPiPlusP4;
		TLorentzVector p2_res = resonanceBoost * locPiMinusP4;


		// normal to the production plane
		TVector3 y = (locBeamP4.Vect().Unit().Cross(-locProtonP4.Vect().Unit())).Unit();
  
		// choose helicity frame: z-axis opposite recoil proton in rho rest frame
		TVector3 z = -1. * recoil_res.Vect().Unit();
		TVector3 x = y.Cross(z).Unit();
		TVector3 angles(   (p1_res.Vect()).Dot(x),
                     (p1_res.Vect()).Dot(y),
                     (p1_res.Vect()).Dot(z) );

		double CosTheta = angles.CosTheta();
  
		double phi = angles.Phi();
  
		TVector3 eps(1.0, 0.0, 0.0); // beam polarization vector
		double Phi = atan2(y.Dot(eps), locBeamP4.Vect().Unit().Dot(eps.Cross(y)));

		double psi = phi - Phi;
		if(psi < -3.14159) psi += 2*3.14159;
		  if(psi > 3.14159) psi -= 2*3.14159;

		// cout << " phi=" << phi << " Phi=" << Phi << " psi=" << psi << endl;

		map<Particle_t, set<Int_t> > locUsedThisCombo_Angles;
		locUsedThisCombo_Angles[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_Angles[Proton].insert(locProtonTrackID);
		locUsedThisCombo_Angles[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_Angles[PiMinus].insert(locPiMinusTrackID);
		if(locUsedSoFar_Angles.find(locUsedThisCombo_Angles) == locUsedSoFar_Angles.end())
		{
			dHist_t->Fill(fabs(t));
			dHist_CosTheta_Psi->Fill(psi*180./3.14159, CosTheta);
			dHist_phi->Fill(phi*180./3.14159);
			dHist_Phi->Fill(Phi*180./3.14159);
			dHist_psi->Fill(psi*180./3.14159);
			locUsedSoFar_Angles.insert(locUsedThisCombo_Angles);
		}
	}

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

		/*
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		{
			dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		}
	       

		//FILL FLAT TREE
		//Fill_FlatTree(); //for the active combo
		} // end of combo loop
		*/

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/

	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();


	return kTRUE;
}

void DSelector_p2pi_trees::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}

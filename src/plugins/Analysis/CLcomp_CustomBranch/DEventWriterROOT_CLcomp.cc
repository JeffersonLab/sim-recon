// $Id$
//
//    File: DEventWriterROOT_CLcomp.cc
// Created: Thu Feb  1 13:48:42 EST 2018
// Creator: aebarnes (on Linux egbert 2.6.32-696.13.2.el6.x86_64 x86_64)
//

#include "DEventWriterROOT_CLcomp.h"

//GLUEX TTREE DOCUMENTATION: https://halldweb.jlab.org/wiki/index.php/Analysis_TTreeFormat

void DEventWriterROOT_CLcomp::Create_CustomBranches_DataTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop, const DReaction* locReaction, bool locIsMCDataFlag) const
{
/*
	//EXAMPLES: Create a branch for a single object (e.g. Int_t, Float_t, TVector3):
	//If filling for a specific particle, the branch name should match the particle branch name
	locBranchRegister.Register_Single<UInt_t>("DummyUInt");
	locBranchRegister.Register_Single<Float_t>("PiPlus__DummyFloat");
	locBranchRegister.Register_Single<TVector3>("Dummy3Vector");
	locBranchRegister.Register_Single<TLorentzVector>("PiPlus__Dummy4Vector");
*/

/*
	//EXAMPLES: Create a branch to hold an array of fundamental type:
	//If filling for a specific particle, the branch name should match the particle branch name
	//locArraySizeString is the name of the branch whose variable that contains the size of the array for that tree entry
		//To match the default TTree branches, use either: 'NumThrown', 'NumBeam', 'NumChargedHypos', 'NumNeutralHypos', or 'NumCombos', as appropriate
	unsigned int locInitArraySize = 10; //if too small, will auto-increase as needed, but requires new calls //if too large, uses more memory than needed
	locBranchRegister.Register_Single<UInt_t>( "DummyArraySize"); //you must store the size of the fundamental array for each entry!!
	locBranchRegister.Register_FundamentalArray<Int_t>("PiPlus__DummyIntArray", "DummyArraySize", locInitArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>("DummyFloatArray", "DummyArraySize", locInitArraySize);
*/
	unsigned int locInitArraySize = 10;

        locBranchRegister.Register_Single<UInt_t>( "ComboArraySize" );
        locBranchRegister.Register_FundamentalArray<double>("kpkm_CL__double", "ComboArraySize", locInitArraySize );
        locBranchRegister.Register_FundamentalArray<double>("p2pi_CL__double", "ComboArraySize", locInitArraySize );
        locBranchRegister.Register_FundamentalArray<double>("IM__double", "ComboArraySize", locInitArraySize );

/*
	//EXAMPLES: Create a branch to hold a TClonesArray of TObject type:
	//If filling for a specific particle, the branch name should match the particle branch name
	unsigned int locInitObjectArraySize = 10; //if too small, will auto-increase as needed, but requires new calls //if too large, uses more memory than needed
	locBranchRegister.Register_ClonesArray<TLorentzVector>("PiPlus__Dummy4VectorArray", locInitObjectArraySize);
	locBranchRegister.Register_ClonesArray<TVector3>("Dummy3VectorArray", locInitObjectArraySize);
*/
}

void DEventWriterROOT_CLcomp::Create_CustomBranches_ThrownTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop) const
{
	//EXAMPLES: See Create_CustomBranches_DataTree
}

void DEventWriterROOT_CLcomp::Fill_CustomBranches_DataTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DReaction* locReaction, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
	const DMCThrownMatching* locMCThrownMatching, const DDetectorMatches* locDetectorMatches,
	const vector<const DBeamPhoton*>& locBeamPhotons, const vector<const DChargedTrackHypothesis*>& locChargedHypos,
	const vector<const DNeutralParticleHypothesis*>& locNeutralHypos, const deque<const DParticleCombo*>& locParticleCombos) const
{
	//The array indices of the particles/combos in the main TTree branches match the vectors of objects passed into this function
		//So if you want to add custom data for each (e.g.) charged track, the correspondence to the main arrays is 1 <--> 1

/*
	//EXAMPLES: Fill a branch for a fundamental data type (e.g. Int_t, Float_t):
	//!!!!! YOU MUST BE SURE THAT template type matches the type you used to create the branch
	locTreeFillData->Fill_Single<UInt_t>("DummyUInt", 14); //14: dummy value
	locTreeFillData->Fill_Single<Float_t>("PiPlus__DummyFloat", TMath::Pi()); //pi: dummy value
*/

	DAnalysisUtilities* dAnalysisUtils = new DAnalysisUtilities(locEventLoop);
        DKinFitUtils_GlueX* dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
        DKinFitter* dKinFitter = new DKinFitter(dKinFitUtils);

        // Set up uniqueness tracking
		// map of a pair of candidateIDs with their corresponding CL or invariant mass
	map<pair<oid_t, oid_t>, double > locKK_map;
        map<pair<oid_t, oid_t>, double > locPiPi_map;
        map<pair<oid_t, oid_t>, double > locIM_map;

	double locConfidenceLevel = 0;
        double locConfidenceLevel2pi = 0;
        double locIM = 0;

	for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
        {
                const DParticleCombo* locParticleCombo = locParticleCombos[loc_i];
                const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
		if (locKinFitResults == NULL)
                        locConfidenceLevel = -1.0;
                else
                        locConfidenceLevel = locKinFitResults->Get_ConfidenceLevel();
		dKinFitter->Reset_NewEvent();

                //CREATE DKINFITPARTICLE OBJECTS FOR EACH PARTICLE
                const DBeamPhoton* locBeamPhoton = static_cast<const DBeamPhoton*>(locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured());
                shared_ptr<DKinFitParticle> locKinFitParticle_BeamPhoton = dKinFitUtils->Make_BeamParticle(locBeamPhoton);
                shared_ptr<DKinFitParticle> locKinFitParticle_Target = dKinFitUtils->Make_TargetParticle(Proton);

                // Get the Pi+, pi- chargedtrackhypotheses and make particles out of that
                shared_ptr<DKinFitParticle> locKinFitParticle_PiPlus;
                shared_ptr<DKinFitParticle> locKinFitParticle_PiMinus;
                shared_ptr<DKinFitParticle> locKinFitParticle_KPlus = dKinFitUtils->Make_DetectedParticle(locParticleCombo->Get_ParticleComboStep(0)->Get_FinalParticle_Measured(0));
                shared_ptr<DKinFitParticle> locKinFitParticle_KMinus = dKinFitUtils->Make_DetectedParticle(locParticleCombo->Get_ParticleComboStep(0)->Get_FinalParticle_Measured(1));
                shared_ptr<DKinFitParticle> locKinFitParticle_Proton = dKinFitUtils->Make_DetectedParticle(locParticleCombo->Get_ParticleComboStep(0)->Get_FinalParticle_Measured(2));

                // Get invariant mass of KK
                locIM = (locKinFitParticle_KPlus->Get_P4() + locKinFitParticle_KMinus->Get_P4()).M();

                vector<const JObject*> FinalObjs = locParticleCombo->Get_FinalParticle_SourceObjects();
                oid_t pos_trackID;
                oid_t neg_trackID;
                for(size_t loc_i = 0; loc_i < FinalObjs.size() - 1; ++loc_i) // size = 3, same order as reaction (K+, K-, p)
                {
                        auto particle = dynamic_cast<const DChargedTrack*>(FinalObjs[loc_i]);
                        if(particle == nullptr) continue;

                        if(particle->Get_Charge() > 0)
                        {
                                const DKinematicData *locKinData_PiPlus = dynamic_cast<const DKinematicData*>(particle->Get_Hypothesis(PiPlus));
				if(locKinData_PiPlus == nullptr) continue;

                                locKinFitParticle_PiPlus = dKinFitUtils->Make_DetectedParticle(locKinData_PiPlus);
                                pos_trackID = particle->candidateid;

                        }
                        else if(particle->Get_Charge() < 0)
                        {
                                const DKinematicData *locKinData_PiMinus = dynamic_cast<const DKinematicData*>(particle->Get_Hypothesis(PiMinus));
                                if(locKinData_PiMinus == nullptr) continue;
                                locKinFitParticle_PiMinus = dKinFitUtils->Make_DetectedParticle(locKinData_PiMinus);
                                neg_trackID = particle->candidateid;
                        }
                }

		pair<oid_t, oid_t> locUsedThisCombo = make_pair(pos_trackID, neg_trackID);

		if (locKinFitParticle_PiPlus == nullptr || locKinFitParticle_PiMinus == nullptr)
                {
                        locConfidenceLevel2pi = -1.0;
                }
		else
		{

                	// SETUP THE CONSTRAINTS
                	dKinFitter->Reset_NewFit(); //discards everything from the previous kinematic fit

                	// P4 Constraint
                	set<shared_ptr<DKinFitParticle>> locInitialParticles, locFinalParticles;
                	locInitialParticles.insert(locKinFitParticle_BeamPhoton);  locInitialParticles.insert(locKinFitParticle_Target);
                	locFinalParticles.insert(locKinFitParticle_PiPlus);  locFinalParticles.insert(locKinFitParticle_PiMinus); locFinalParticles.insert(locKinFitParticle_Proton);

                	shared_ptr<DKinFitConstraint_P4> locP4Constraint = dKinFitUtils->Make_P4Constraint(locInitialParticles, locFinalParticles);
                	TVector3 locMomentum = locKinFitParticle_BeamPhoton->Get_Momentum() + locKinFitParticle_Target->Get_Momentum();
                	locMomentum -= (locKinFitParticle_PiPlus->Get_Momentum() + locKinFitParticle_PiMinus->Get_Momentum() + locKinFitParticle_Proton->Get_Momentum());
                	locP4Constraint->Set_InitP3Guess(locMomentum);
                	dKinFitter->Add_Constraint(locP4Constraint);

                	// Vertex constraint
                	set<shared_ptr<DKinFitParticle>> locFullConstrainParticles, locNoConstrainParticles;
                	locFullConstrainParticles.insert(locKinFitParticle_PiPlus);
                	locFullConstrainParticles.insert(locKinFitParticle_PiMinus);
                	locFullConstrainParticles.insert(locKinFitParticle_Proton);
                	locNoConstrainParticles.insert(locKinFitParticle_BeamPhoton); locNoConstrainParticles.insert(locKinFitParticle_Target);

                	vector<shared_ptr<DKinFitParticle>> locVectFinalParticles;
                	locVectFinalParticles.push_back(locKinFitParticle_PiPlus);
                	locVectFinalParticles.push_back(locKinFitParticle_PiMinus);
                	locVectFinalParticles.push_back(locKinFitParticle_Proton);
                	TVector3 locVertexGuess = dAnalysisUtils->Calc_CrudeVertex(locVectFinalParticles);
                	shared_ptr<DKinFitConstraint_Vertex> locDecayVertexConstraint = dKinFitUtils->Make_VertexConstraint(locFullConstrainParticles, locNoConstrainParticles, locVertexGuess);
                	dKinFitter->Add_Constraint(locDecayVertexConstraint);

                	// PERFORM THE KINEMATIC FIT
        		dKinFitter->Set_DebugLevel(0); //e.g. 500 for everything
			dKinFitter->Fit_Reaction();

			locConfidenceLevel2pi = dKinFitter->Get_ConfidenceLevel();
		}

		// Check for uniqueness. Keep information based on best KK confidence level for a given pair of track IDs.
			// This might not be the best approach. Investigate this further.
                map<pair<oid_t, oid_t>, double >::iterator it = locKK_map.find(locUsedThisCombo);
                if (it != locKK_map.end() && locConfidenceLevel > locKK_map.at(locUsedThisCombo))
                {
                        locKK_map.at(locUsedThisCombo) = locConfidenceLevel;
                        locPiPi_map.at(locUsedThisCombo) = locConfidenceLevel2pi;
                }
                else
                {
                        locKK_map.insert( make_pair(locUsedThisCombo, locConfidenceLevel) );
                        locPiPi_map.insert( make_pair(locUsedThisCombo, locConfidenceLevel2pi) );
                        locIM_map.insert( make_pair(locUsedThisCombo, locIM) );
                }

	}

        locTreeFillData->Fill_Single<UInt_t>("ComboArraySize", locKK_map.size());
        size_t loc_index = 0;
        for (auto const& x : locKK_map)
        {
                locTreeFillData->Fill_Array<double>("kpkm_CL__double", x.second, loc_index);
                loc_index++;
        }
        loc_index = 0;
        for (auto const& x : locPiPi_map)
        {
                locTreeFillData->Fill_Array<double>("p2pi_CL__double", x.second, loc_index);
                loc_index++;
        }
        loc_index = 0;
        for (auto const& x : locIM_map)
        {
                locTreeFillData->Fill_Array<double>("IM__double", x.second, loc_index);
                loc_index++;
        }
/*
	//EXAMPLES: Fill a branch for a TObject data type (e.g. TVector3, TLorentzVector):
	//!!!!! YOU MUST BE SURE THAT template type matches the type you used to create the branch
	TVector3 locPosition(0.0, 0.0, 65.0);
	locTreeFillData->Fill_Single<TVector3>("Dummy3Vector", locPosition);
	TLorentzVector locP4(1.0, 2.0, 3.0, 4.0);
	locTreeFillData->Fill_Single<TLorentzVector>("PiPlus__Dummy4Vector", locP4);
*/

/*
	//EXAMPLES: Fill a branch with an array of fundamental type:
	//!!!!! YOU MUST BE SURE THAT template type matches the type you used to create the branch
	Int_t locOutputArraySize = 7;
	locTreeFillData->Fill_Single<UInt_t>("DummyArraySize", locOutputArraySize);
	for(Int_t loc_i = 0; loc_i < locOutputArraySize; ++loc_i)
	{
		Int_t locValue = loc_i - 14; //dummy number
		locTreeFillData->Fill_Array<Int_t>("PiPlus__DummyIntArray", locValue, loc_i);
		locTreeFillData->Fill_Array<Float_t>("DummyFloatArray", locValue + 9999, loc_i);
	}
*/

/*
	//EXAMPLES: Fill a branch with a TClonesArray of TObject type:
	//!!!!! YOU MUST BE SURE THAT template type matches the type you used to create the branch
	for(Int_t loc_i = 0; loc_i < 15; ++loc_i)
	{
		TLorentzVector locNewP4(99.0, 2.0, 3.0, 4.0);
		locTreeFillData->Fill_Array<TLorentzVector>("PiPlus__Dummy4VectorArray", locNewP4, loc_i);
		TVector3 locNewPosition(100.0, 0.0, 65.0);
		locTreeFillData->Fill_Array<TVector3>("Dummy3VectorArray", locNewPosition, loc_i);
	}
*/
}

void DEventWriterROOT_CLcomp::Fill_CustomBranches_ThrownTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns) const
{
	//EXAMPLES: See Fill_CustomBranches_DataTree
}


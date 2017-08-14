#include "ANALYSIS/DSourceComboP4Handler.h"
#include "ANALYSIS/DSourceComboer.h"

/*************************************************** PHOTON P4 RECONSTRUCTION **************************************************
*
* The exact photon momentum is unknown until its production vertex is known.
* However, that vertex is combo-dependent. We'd like to make cuts on the pi0 mass globally in advance, rather than on a combo-by-combo basis.
* This would be a huge savings in time and memory.
*
* The momentum of the hypotheses calculated in DNeutralParticle is based on the DVertex (class) position.
* This vertex is determined from all of the "good" charged tracks in the event, typically by doing a kinematic fit.
*
* However, there a several potential problems with using this vertex:
* 1) There may have been extra (junk, accidental) tracks reconstructed in the event. These will throw off the vertex position.
*    And, if there are only 2 tracks, it can throw it off considerably.
* 2) The track position & momentum are different for each hypothesis, and it's not clear in advance which hypothesis should be used.
* 3) The photons may have come from a detached vertex, which can be 10+ cm downstream of the main vertex.
*
* So while the DVertex is OK for someone looking for a quick estimate, it should not be used in an actual analysis.
*
* Now, as we said before, the true photon momentum is combo-dependent, and we want to do loose mass cuts in a combo-independent way.
* So, we can compute all of the p4's at a specific vertex position (e.g. center of the target), rather than separately for each combo.
* But, how much of an impact will a given error in the vertex position have on the calculated 2-photon invariant mass?
*
* The calculation below determines that the maximum 2-photon-mass error occurs when both photons are at 45 degrees near the eta peak (less near pi0).
* Specifically: z_err = 5cm yields a mass error of ~20 MeV, 10cm -> ~40 MeV, 15cm -> ~60 MeV, etc.
*
* So, what is the maximum delta_m we can tolerate with our loose cuts?
* The idea is that the loose cuts should be wide enough around the signal peak to provide enough statistics for sideband subtraction.
* So, no matter what we choose, it won't affect the signal peak.  But we also don't want to affect the sidebands too much.
* E.g. pi0 mass peak is from ~110 -> ~160 MeV, loose cut is 50 -> 220 MeV at the moment
* Therefore, you probably want to keep the maximum delta_m at around the 20 MeV level.
* This means that the max z_error should be about 5 cm
*
* Therefore, every 10 cm, from 5cm upstream of the target to 15 cm downstream of the target (detached vertex) (5 bins):
* Compute p4s at the centers of these vertex-z bins and do loose mass cuts
*
*******************************************************************************************************************************/

/************************************************ 2-PHOTON MASS ERROR DERIVATION ***********************************************
*
* For this error estimate, consider the decay pi0 -> 2g  (or eta -> 2g).
* The equation for the invariant mass squared of the 2 photons is:
* m^2 = (E1 + E2)^2 - (px1 + px2)^2 - (py1 + py2)^2 - (pz1 + pz2)^2
*
* The difference in mass squared due to the error is: (using "_g" prefix for guess and "_t" prefix for "true")
* delta(m^2) = m_g^2 - m^2_t
* delta(m^2) = (px1 + px2)^2_t + (py1 + py2)^2_t + (pz1 + pz2)^2_t - (px1 + px2)^2_g - (py1 + py2)^2_g - (pz1 + pz2)^2_g
*
* Simplifying, we will choose the 2 photons to have opposite phi's at 0 & 180 degrees
* Thus, py1 = py2 = 0 for guess & true momentum.
* delta(m^2) = (px1 + px2)^2_t + (pz1 + pz2)^2_t - (px1 + px2)^2_g - (pz1 + pz2)^2_g
* Also, px1_unit = -px2_unit, but we'll use this a little later.
*
* Now, assume that both photons have the same E & theta: Whatever the worst-case is, it will be the same for both photons.
* Since from before px1_unit = -px2_unit, Now px1_t = -px2_t and px1_g = -px2_g. Also, pz1_t = pz2_t = pz_t, & pz1_g = pz2_g = pz_g
* delta(m^2) = 4*pz_t^2 - 4*pz_g^2
*
* Now for photons, pz = E*dz/d, where d is the distance between the shower and the vertex, and dz is the z-component of d.
* Plugging this in gives:
* delta(m^2) = 4*E^2*[(dz_t/d_t)^2 - (dz_g/d_g)^2]
* Using dz_g/d_g = cos(theta_g) gives:
* delta(m^2) = 4*E^2*[(dz_t/d_t)^2 - cos^2(theta_g)]
*
* Now, m_g^2 = (E1 + E2)^2 - (px1 + px2)^2_g - (py1 + py2)^2_g - (pz1 + pz2)^2_g
* However, we've defined our photon guesses pretty narrowly: py = 0, px1 = -px2, pz1 = pz2, E1 = E2
* Thus, m_g^2 = (2E)^2 - (2pz_g)^2
* Plugging in pz_g = E*dz_g/d_g yields:
* Thus, m_g^2 = 4*E^2*[1 - (dz_g/d_g)^2]
* Using dz_g/d_g = cos(theta_g)
* Thus, m_g^2 = 4*E^2*sin^2(theta_g)
* Rearranging gives: 4*E^2 = m_g^2/sin^2(theta_g)
*
* Plugging the above into delta(m^2)
* delta(m^2) = m_g^2*[(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g)
*
* But, we want delta_m, not delta(m^2)
* delta_m = m_g - m_t, m_t = sqrt(m_g^2 - delta(m^2))
* delta_m = m_g - sqrt(m_g^2 - delta(m^2))
* delta_m = m_g - sqrt(m_g^2 - m_g^2*[(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g))
*
* Rearrange and cancel terms
* delta_m = m_g - m_g*sqrt(1 - [(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g))
* delta_m = m_g - m_g*sqrt([sin^2(theta_g) - (dz_t/d_t)^2 + cos^2(theta_g)]/sin^2(theta_g))
* delta_m = m_g - m_g*sqrt[1 - (dz_t/d_t)^2]/sin(theta_g)
* delta_m = m_g - m_g*sqrt[(d_t^2 - dz_t^2)/d_t^2]/sin(theta_g)
*
* Note that d_t^2 - dz_t^2 = dx^2 (since dy is zero)
* delta_m = m_g - m_g*sqrt[dx^2/d_t^2]/sin(theta_g)
* delta_m = m_g - m_g*dx/(d_t*sin(theta_g))
*
* Getting the true dz_t & d_t in terms of some z_error:
* d_t = sqrt(dx^2 + dz_t^2), dz_t = dz_g + z_error
* delta_m = m_g - m_g*dx/(sin(theta_g)*sqrt(dx^2 + (dz_g + z_error)^2))
* And dz_g = dx/tan(theta_g)
* delta_m = m_g - m_g*dx/(sin(theta_g)*sqrt(dx^2 + (dx/tan(theta_g) + z_error)^2))
* delta_m = m_g - m_g/(sin(theta_g)*sqrt(1 + (1/tan(theta_g) + z_error/dx)^2))
*
* For BCAL, dx = 65.
* For the FCAL, dx = dz*tan(theta), dz = 650 - z_error (approx 650):
* delta_m = m_g - m_g/(sin(theta_g)*sqrt(1 + (1 + z_error/(650 - z_error))^2/tan^2(theta_g)))
* delta_m = m_g - m_g/(cos(theta_g)*sqrt(tan^2(theta_g) + (1 + z_error/(650 - z_error))^2))
*
* delta_m Is larger at higher m_g, max at 45 degrees (and is thus small for FCAL)
* In fact, for the FCAL, delta_m is ~25 MeV for the eta mass when the z_error is 30cm (max: center of target + detached vertex)
* Therefore, if the center of the target is used, the error is negligible compared to the width of the mass cut.
*
* For the BCAL:
* With m_g near eta mass, z_error = 15: delta_m_max = ~60 MeV
* With m_g near eta mass, z_error = 10: delta_m_max = ~40 MeV
* With m_g near eta mass, z_error = 5: delta_m_max = ~20 MeV
* With m_g near eta mass, z_error = 3: delta_m_max = ~15 MeV
* With m_g near eta mass, z_error = 2: delta_m_max = ~9 MeV
* Instead of the above, you can of course plot the delta_m for real data, and get something similar
* So, choose a z_error of 5: compute at center of 10-cm-wide bins.
*
* OK, but what about eta -> 3pi0?  Or Lambda -> pi0, n?
* eta -> 3pi0: Max error for a given pi0 is ~small, not bad when combined with 3 others: it's fine as long as cut is wide.
* pi0, n: Neutron is likely slow, similar to charged tracks: Error is much larger, cannot combo massive neutrals without exact vertex position
* Well, OK, we can COMBO them, but we can't place mass cuts.
*
*******************************************************************************************************************************/


namespace DAnalysis
{

DSourceComboP4Handler::DSourceComboP4Handler(DSourceComboer* locSourceComboer, bool locCreateHistsFlag) : dSourceComboer(locSourceComboer)
{
	gPARMS->SetDefaultParameter("COMBO:DEBUG_LEVEL", dDebugLevel);

	//INVARIANT MASS CUTS: MESONS
	dInvariantMassCuts.emplace(Pi0, std::make_pair(0.08, 0.19));
	dInvariantMassCuts.emplace(KShort, std::make_pair(0.3, 0.7));
	dInvariantMassCuts.emplace(Eta, std::make_pair(0.35, 0.75));
	dInvariantMassCuts.emplace(omega, std::make_pair(0.4, 1.2));
	dInvariantMassCuts.emplace(EtaPrime, std::make_pair(0.6, 1.3));
	dInvariantMassCuts.emplace(phiMeson, std::make_pair(0.8, 1.2));
	dInvariantMassCuts.emplace(D0, std::make_pair(1.8, 1.92));
//	dInvariantMassCuts.emplace(Jpsi, std::make_pair(2.7, 3.5));

	//INVARIANT MASS CUTS: BARYONS
	dInvariantMassCuts.emplace(Lambda, std::make_pair(1.0, 1.2));
	dInvariantMassCuts.emplace(AntiLambda, dInvariantMassCuts[Lambda]);
	dInvariantMassCuts.emplace(Sigma0, std::make_pair(1.1, 1.3));
	dInvariantMassCuts.emplace(SigmaPlus, dInvariantMassCuts[Sigma0]);
	dInvariantMassCuts.emplace(SigmaMinus, dInvariantMassCuts[Sigma0]);
	dInvariantMassCuts.emplace(XiMinus, std::make_pair(1.1, 1.5));
	dInvariantMassCuts.emplace(Xi0, dInvariantMassCuts[XiMinus]);
	dInvariantMassCuts.emplace(OmegaMinus, std::make_pair(1.32, 2.22));
	dInvariantMassCuts.emplace(Lambda_c, std::make_pair(2.0, 2.6));

	//MISSING MASS CUTS & MASS HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!! //I have no idea why this is needed, but without it it crashes.  Sigh. 
	{
		//Cuts are vs. beam energy
		//None missing & photon
		dMissingMassSquaredCuts[Unknown].first = new TF1("df_MissingMassSquaredCut_NoneLow", "[0]", 0.0, 12.0);
		dMissingMassSquaredCuts[Unknown].first->SetParameter(0, -0.1);
		dMissingMassSquaredCuts[Unknown].second = new TF1("df_MissingMassSquaredCut_NoneHigh", "[0]", 0.0, 12.0);
		dMissingMassSquaredCuts[Unknown].second->SetParameter(0, 0.1);
		dMissingMassSquaredCuts.emplace(Gamma, dMissingMassSquaredCuts[Unknown]);

		//Pi+/-
		dMissingMassSquaredCuts[PiPlus].first = new TF1("df_MissingMassSquaredCut_PiPlusLow", "[0]", 0.0, 12.0);
		dMissingMassSquaredCuts[PiPlus].first->SetParameter(0, -1.0);
		dMissingMassSquaredCuts[PiPlus].second = new TF1("df_MissingMassSquaredCut_PiPlusHigh", "[0]", 0.0, 12.0);
		dMissingMassSquaredCuts[PiPlus].second->SetParameter(0, 1.0);
		dMissingMassSquaredCuts.emplace(PiMinus, dMissingMassSquaredCuts[PiPlus]);
		dMissingMassSquaredCuts.emplace(Pi0, dMissingMassSquaredCuts[PiPlus]);

		//Proton/Neutron
		dMissingMassSquaredCuts[Proton].first = new TF1("df_MissingMassSquaredCut_ProtonLow", "[0]", 0.0, 12.0);
		dMissingMassSquaredCuts[Proton].first->SetParameter(0, -0.5);
		dMissingMassSquaredCuts[Proton].second = new TF1("df_MissingMassSquaredCut_ProtonHigh", "[0]", 0.0, 12.0);
		dMissingMassSquaredCuts[Proton].second->SetParameter(0, 4.41);
		dMissingMassSquaredCuts.emplace(Neutron, dMissingMassSquaredCuts[Proton]);

/*
		//e+/-
		dMissingMassSquaredCuts[Electron].first = new TF1("df_MissingMassSquaredCut_PiPlusLow", "[0]", 0.0, 12.0);
		dMissingMassSquaredCuts[Electron].first->SetParameter(0, -1.0);
		dMissingMassSquaredCuts[Electron].second = new TF1("df_MissingMassSquaredCut_PiPlusHigh", "[0]", 0.0, 12.0);
		dMissingMassSquaredCuts[Electron].second->SetParameter(0, 1.0);
		dMissingMassSquaredCuts.emplace(Positron, dMissingMassSquaredCuts[Electron]);
*/
		//Missing E cuts
		dMissingECuts = std::pair<TF1*, TF1*>(nullptr, nullptr); //COMPARE:
/*
		dMissingECuts.first = new TF1("df_MissingECut_NoneLow", "[0]", 0.0, 12.0);
		dMissingECuts.first->SetParameter(0, -3.0);
		dMissingECuts.second = new TF1("df_MissingECut_NoneHigh", "[0]", 0.0, 12.0);
		dMissingECuts.second->SetParameter(0, 3.0);
*/
		if(!locCreateHistsFlag)
		{
			japp->RootUnLock(); //RELEASE ROOT LOCK!!
			return;
		}

		//HISTOGRAMS
		//get and change to the base (file/global) directory
		TDirectory* locCurrentDir = gDirectory;

		string locDirName = "Independent";
		TDirectoryFile* locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		locDirName = "Combo_Construction";
		locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		locDirName = "Invariant_Mass";
		locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		//INVARIANT MASS HISTOGRAMS: 2GAMMA
		{
			string locHistName = "InvariantMass_2Gamma_BCAL";
			auto locHist = gDirectory->Get(locHistName.c_str());
			dHistMap_2GammaMass[SYS_BCAL] = (locHist != nullptr) ? static_cast<TH1*>(locHist) : new TH1I(locHistName.c_str(), ";BCAL 2#gamma Invariant Mass (GeV/c^{2})", 2000, 0.0, 2.0);

			locHistName = "InvariantMass_2Gamma_FCAL";
			locHist = gDirectory->Get(locHistName.c_str());
			dHistMap_2GammaMass[SYS_FCAL] = (locHist != nullptr) ? static_cast<TH1*>(locHist) : new TH1I(locHistName.c_str(), ";FCAL 2#gamma Invariant Mass (GeV/c^{2})", 2000, 0.0, 2.0);

			locHistName = "InvariantMass_2Gamma_BCALFCAL";
			locHist = gDirectory->Get(locHistName.c_str());
			dHistMap_2GammaMass[SYS_NULL] = (locHist != nullptr) ? static_cast<TH1*>(locHist) : new TH1I(locHistName.c_str(), ";BCAL/FCAL 2#gamma Invariant Mass (GeV/c^{2})", 2000, 0.0, 2.0);
		}

		//INVARIANT MASS HISTOGRAMS: PIDS
		for(const auto& locPIDPair : dInvariantMassCuts)
		{
			auto locPID = locPIDPair.first;
			auto& locMassPair = locPIDPair.second;
			string locHistName = string("InvariantMass_") + ParticleType(locPID);
			auto locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = string("From All Decay Products;") + string(ParticleName_ROOT(locPID)) + string(" Invariant Mass (GeV/c^{2})");
				auto locMinMass = locMassPair.first - 0.2;
				if(locMinMass < 0.0)
					locMinMass = 0.0;
				auto locMaxMass = locMassPair.second + 0.2;
				auto locNumBins = 1000.0*(locMaxMass - locMinMass);
				dHistMap_InvariantMass[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), locNumBins, locMinMass, locMaxMass);
			}
			else
				dHistMap_InvariantMass[locPID] = static_cast<TH1*>(locHist);
		}
		gDirectory->cd("..");

		locDirName = "Missing_Mass";
		locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		//MISSING MASS HISTOGRAMS
		for(const auto& locPIDPair : dMissingMassSquaredCuts)
		{
			auto locPID = locPIDPair.first;
			auto& locMassPair = locPIDPair.second;
			string locHistName = string("MissingMassVsBeamEnergy_") + ((locPID != Unknown) ? ParticleType(locPID) : "None");
			auto locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = string("From All Production Mechanisms;Beam Energy (GeV);") + string((locPID != Unknown) ? ParticleName_ROOT(locPID) : "None") + string(" Missing Mass Squared (GeV/c^{2})^{2}");
				auto locMinMass = locMassPair.first->Eval(12.0) - 0.2; //assume widest at highest energy
				auto locMaxMass = locMassPair.second->Eval(12.0) + 0.2;
				auto locNumBins = 1000.0*(locMaxMass - locMinMass);
				dHistMap_MissingMassSquaredVsBeamEnergy[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 600, 0.0, 12.0, locNumBins, locMinMass, locMaxMass);
			}
			else
				dHistMap_MissingMassSquaredVsBeamEnergy[locPID] = static_cast<TH2*>(locHist);
		}

		//NONE MISSING, MISSING E:
		{
			string locHistName = "NoneMissing_MissingEVsBeamEnergy_PreMissMassSqCut";
			auto locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = string("None Missing: From All Production Mechanisms;Beam Energy (GeV); Missing E (GeV)");
				dHist_NoneMissing_MissingEVsBeamEnergy_PreMissMassSqCut = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 600, 0.0, 12.0, 1200, -6.0, 6.0);
			}
			else
				dHist_NoneMissing_MissingEVsBeamEnergy_PreMissMassSqCut = static_cast<TH2*>(locHist);

			locHistName = "NoneMissing_MissingEVsBeamEnergy_PostMissMassSqCut";
			locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = string("None Missing: From All Production Mechanisms;Beam Energy (GeV); Missing E (GeV)");
				dHist_NoneMissing_MissingEVsBeamEnergy_PostMissMassSqCut = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 600, 0.0, 12.0, 1200, -6.0, 6.0);
			}
			else
				dHist_NoneMissing_MissingEVsBeamEnergy_PostMissMassSqCut = static_cast<TH2*>(locHist);
		}

		//NONE MISSING, MISSING PT:
		{
			string locHistName = "NoneMissing_MissingPtVsMissingE_PreMissMassSqCut";
			auto locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = string("None Missing: From All Production Mechanisms; Missing E (GeV); Missing Transverse Momentum (GeV/c)");
				dHist_NoneMissing_MissingPtVsMissingE_PreMissMassSqCut = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 1200, -6.0, 6.0, 800, 0.0, 4.0);
			}
			else
				dHist_NoneMissing_MissingPtVsMissingE_PreMissMassSqCut = static_cast<TH2*>(locHist);

			locHistName = "NoneMissing_MissingPtVsMissingE_PostMissMassSqCut";
			locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = string("None Missing: From All Production Mechanisms; Missing E (GeV); Missing Transverse Momentum (GeV/c)");
				dHist_NoneMissing_MissingPtVsMissingE_PostMissMassSqCut = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 1200, -6.0, 6.0, 600, 0.0, 3.0);
			}
			else
				dHist_NoneMissing_MissingPtVsMissingE_PostMissMassSqCut = static_cast<TH2*>(locHist);
		}

		locCurrentDir->cd();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

DLorentzVector DSourceComboP4Handler::Get_P4_NotMassiveNeutral(Particle_t locPID, const JObject* locObject, const DVector3& locVertex, bool locAccuratePhotonsFlag) const
{
	if(ParticleCharge(locPID) != 0)
		return static_cast<const DChargedTrack*>(locObject)->Get_Hypothesis(locPID)->lorentzMomentum();

	//assume is photon!
	auto locNeutralShower = static_cast<const DNeutralShower*>(locObject);
	if(!locAccuratePhotonsFlag)
	{
		signed char locVertexZBin = DSourceComboInfo::Get_VertexZIndex_ZIndependent();
		if(locNeutralShower->dDetectorSystem == SYS_BCAL)
			locVertexZBin = dSourceComboTimeHandler->Get_PhotonVertexZBin(locVertex.Z());
		if(locVertexZBin == DSourceComboInfo::Get_VertexZIndex_OutOfRange())
			locVertexZBin = dSourceComboTimeHandler->Get_VertexZBin_TargetCenter(); //we need a zbin for BCAL showers, but it is unknown: we must pick something: center of target
		auto& locKinematicData = dPhotonKinematics.find(locVertexZBin)->second.find(locNeutralShower)->second;
		return locKinematicData->lorentzMomentum();
	}

	//vertex is known: calculate exactly
	auto locMomentum = locNeutralShower->dEnergy*((locNeutralShower->dSpacetimeVertex.Vect() - locVertex).Unit());
	return DLorentzVector(locMomentum, locNeutralShower->dEnergy);
}

DLorentzVector DSourceComboP4Handler::Calc_MassiveNeutralP4(const DNeutralShower* locNeutralShower, Particle_t locPID, const DVector3& locVertex, double locRFVertexTime) const
{
	//locRFVertexTime: the RF time propagated to the vertex, through any decaying particles if necessary
	auto locMass = ParticleMass(locPID);
	auto locPath = locNeutralShower->dSpacetimeVertex.Vect() - locVertex;
	double locDeltaT = locNeutralShower->dSpacetimeVertex.T() - locRFVertexTime;

	double locBeta = locPath.Mag()/(locDeltaT*29.9792458);
	if(locBeta >= 1.0)
		locBeta = dMaxMassiveNeutralBeta;
	if(locBeta < 0.0)
		locBeta = 0.0;

	auto locGamma = 1.0/sqrt(1.0 - locBeta*locBeta);
	auto locPMag = locGamma*locBeta*locMass;
	locPath.SetMag(locPMag); //is now the momentum!

	auto locEnergy = sqrt(locPMag*locPMag + locMass*locMass);
	return DLorentzVector(locPath, locEnergy);
}

DLorentzVector DSourceComboP4Handler::Calc_P4_NoMassiveNeutrals(const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DVector3& locVertex, signed char locVertexZBin, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag)
{
	//locVertexZBin MUST be a separate argument, to signal special bins!!!
	if(!locVertexCombo->Get_IsComboingZIndependent() && ((locVertexZBin == DSourceComboInfo::Get_VertexZIndex_Unknown()) || (locVertexZBin == DSourceComboInfo::Get_VertexZIndex_OutOfRange())))
		locVertexZBin = dSourceComboTimeHandler->Get_VertexZBin_TargetCenter(); //we need a zbin for BCAL showers, but it is unknown: we must pick something: center of target

	auto locHasPhotons = !DAnalysis::Get_SourceParticles(locVertexCombo->Get_SourceParticles(true, d_Neutral), Gamma).empty();
	auto locIterator = (locHasPhotons && locAccuratePhotonsFlag) ? dFinalStateP4ByCombo.end() : dFinalStateP4ByCombo.find(std::make_pair(locVertexCombo, locVertexZBin));
	if(locIterator != dFinalStateP4ByCombo.end())
	{
		if(dDebugLevel >= 10)
			cout << "Read Calc_P4_NoMassiveNeutrals: Combo " << locVertexCombo << " P4: " << locIterator->second.Px() << ", " << locIterator->second.Py() << ", " << locIterator->second.Pz() << ", " << locIterator->second.E() << endl;
		return locIterator->second;
	}

	//loop over particles
	//vertex-z bin may be different for decay products! (detached vertex)
	//save/retrieve masses by combo instead
	DLorentzVector locTotalP4 = Calc_P4_SourceParticles(locVertexCombo, locVertex, 0.0, locAccuratePhotonsFlag);

	//loop over decays
	auto locFurtherDecayCombos = locVertexCombo->Get_FurtherDecayCombos();
	for(const auto& locCombosByUsePair : locFurtherDecayCombos)
	{
		auto locDetachedVertexFlag = IsDetachedVertex(std::get<0>(locCombosByUsePair.first));
		auto locDecayVertexZBin = std::get<1>(locCombosByUsePair.first);
		for(const auto& locCombo : locCombosByUsePair.second)
		{
			auto locIsVertexKnown = dSourceComboVertexer->Get_IsVertexKnown(false, locReactionCombo, locCombo, locBeamParticle); //false: will be detached if needed
			auto locDecayVertex = (locDetachedVertexFlag && locIsVertexKnown) ? dSourceComboVertexer->Get_Vertex(false, locReactionCombo, locCombo, locBeamParticle) : locVertex;
			locTotalP4 += Calc_P4_NoMassiveNeutrals(locReactionCombo, locCombo, locDecayVertex, locDecayVertexZBin, locBeamParticle, locAccuratePhotonsFlag);
		}
	}

	//save the results and return
	if(!locAccuratePhotonsFlag || !locHasPhotons)
		dFinalStateP4ByCombo.emplace(std::make_pair(locVertexCombo, locVertexZBin), locTotalP4);

	if(!locAccuratePhotonsFlag)
	{
		auto locSourceParticles = locVertexCombo->Get_SourceParticles();
		if((locSourceParticles.size() == 2) && (locSourceParticles[0].first == Gamma) && (locSourceParticles[1].first == Gamma))
		{
			auto locSystem1 = static_cast<const DNeutralShower*>(locSourceParticles[0].second)->dDetectorSystem;
			auto locSystem2 = static_cast<const DNeutralShower*>(locSourceParticles[1].second)->dDetectorSystem;
			auto locSystem = (locSystem1 != locSystem2) ? SYS_NULL : locSystem1;
			d2GammaInvariantMasses[locSystem].push_back(locTotalP4.M());
		}
	}

	if(dDebugLevel >= 10)
		cout << "Save Calc_P4_NoMassiveNeutrals: Combo " << locVertexCombo << " P4: " << locTotalP4.Px() << ", " << locTotalP4.Py() << ", " << locTotalP4.Pz() << ", " << locTotalP4.E() << endl;
	return locTotalP4;
}

DLorentzVector DSourceComboP4Handler::Calc_P4_SourceParticles(const DSourceCombo* locVertexCombo, const DVector3& locVertex, double locRFVertexTime, bool locAccuratePhotonsFlag)
{
	//z-bin is kept separate from locVertex because it may indicate special values, or the vertex may not be known yet
	DLorentzVector locTotalP4(0.0, 0.0, 0.0, 0.0);
	auto locSourceParticles = locVertexCombo->Get_SourceParticles(false); //false: NOT the whole chain
	for(const auto& locParticlePair : locSourceParticles)
	{
		auto locPID = locParticlePair.first;
		if((ParticleCharge(locPID) == 0) && (ParticleMass(locPID) > 0.0))
		{
			auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);
			auto locParticleP4 = Calc_MassiveNeutralP4(locNeutralShower, locPID, locVertex, locRFVertexTime);
			if(dDebugLevel >= 20)
				cout << "pid, pointer, pxyzE = " << locPID << ", " << locParticlePair.second << ", " << locParticleP4.Px() << ", " << locParticleP4.Py() << ", " << locParticleP4.Pz() << ", " << locParticleP4.E() << endl;
			locTotalP4 += locParticleP4;
		}
		else
		{
			auto locParticleP4 = Get_P4_NotMassiveNeutral(locParticlePair.first, locParticlePair.second, locVertex, locAccuratePhotonsFlag);
			if(dDebugLevel >= 20)
				cout << "pid, pointer, pxyzE = " << locPID << ", " << locParticlePair.second << ", " << locParticleP4.Px() << ", " << locParticleP4.Py() << ", " << locParticleP4.Pz() << ", " << locParticleP4.E() << endl;
			locTotalP4 += locParticleP4;
		}
	}

	if(dDebugLevel >= 10)
		cout << "Calc_P4_SourceParticles: Combo " << locVertexCombo << " P4: " << locTotalP4.Px() << ", " << locTotalP4.Py() << ", " << locTotalP4.Pz() << ", " << locTotalP4.E() << endl;

	return locTotalP4;
}

bool DSourceComboP4Handler::Calc_P4_Decay(bool locIsProductionVertex, bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceComboUse& locDecayUse, const DSourceCombo* locDecayCombo, DVector3 locVertex, int locRFBunch, double locRFVertexTime, DLorentzVector& locDecayP4, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag)
{
	auto locDecayPID = std::get<0>(locDecayUse);
	auto locDecayVertexZBin = std::get<1>(locDecayUse);
	auto locHasMassiveNeutrals = dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locDecayUse));

	locDecayP4 = DLorentzVector(0.0, 0.0, 0.0, 0.0);
	if(!IsDetachedVertex(locDecayPID))
	{
		if(!locHasMassiveNeutrals)
			locDecayP4 = Calc_P4_NoMassiveNeutrals(locReactionFullCombo, locDecayCombo, locVertex, locDecayVertexZBin, locBeamParticle, locAccuratePhotonsFlag);
		else if(!Calc_P4_HasMassiveNeutrals(locIsProductionVertex, locIsPrimaryProductionVertex, locReactionFullCombo, locDecayCombo, locVertex, locRFBunch, locRFVertexTime, DSourceComboUse(Unknown, 0, nullptr, false, Unknown), locDecayP4, locBeamParticle, locAccuratePhotonsFlag))
			return false;
		if(dDebugLevel >= 10)
			cout << "Calc_P4_Decay: Combo " << locDecayCombo << " P4: " << locDecayP4.Px() << ", " << locDecayP4.Py() << ", " << locDecayP4.Pz() << ", " << locDecayP4.E() << endl;
		return true;
	}

	//detached vertex!
	if(!locHasMassiveNeutrals)
	{
		if(locAccuratePhotonsFlag)
			locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locReactionFullCombo, locDecayCombo, locBeamParticle);
		else
			locVertex.SetXYZ(0.0, 0.0, dSourceComboTimeHandler->Get_PhotonVertexZBinCenter(locDecayVertexZBin));
		locDecayP4 = Calc_P4_NoMassiveNeutrals(locReactionFullCombo, locDecayCombo, locVertex, locDecayVertexZBin, locBeamParticle, locAccuratePhotonsFlag);
	}
	else if(!locAccuratePhotonsFlag)
	{
		//p4 better have already been calculated: if not, it means it was uncalcable (e.g. vertex unknown)
		auto locDecayP4LookupTuple = std::make_tuple(locIsProductionVertex, locReactionFullCombo, locDecayCombo, locRFBunch, locBeamParticle);
		auto locDecayIterator = dFinalStateP4ByCombo_HasMassiveNeutrals.find(locDecayP4LookupTuple);
		if(locDecayIterator == dFinalStateP4ByCombo_HasMassiveNeutrals.end())
		{
			if(locBeamParticle == nullptr)
				return false; //failed!

			//try without beam particle
			locDecayP4LookupTuple = std::make_tuple(locIsProductionVertex, locReactionFullCombo, locDecayCombo, locRFBunch, (const DKinematicData*)nullptr);
			locDecayIterator = dFinalStateP4ByCombo_HasMassiveNeutrals.find(locDecayP4LookupTuple);
			if(locDecayIterator == dFinalStateP4ByCombo_HasMassiveNeutrals.end())
				return false; //failed!
		}
		locDecayP4 = locDecayIterator->second;
	}
	else //recalculate regardless
	{
		locVertex = dSourceComboVertexer->Get_Vertex(false, locReactionFullCombo, locDecayCombo, locBeamParticle);
		auto locPrimaryVertex = dSourceComboVertexer->Get_Vertex(true, locReactionFullCombo, locReactionFullCombo, locBeamParticle);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsPrimaryProductionVertex, locReactionFullCombo, locDecayCombo, locBeamParticle);
		locRFVertexTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertex.Z(), locRFBunch, locTimeOffset);
		if(!Calc_P4_HasMassiveNeutrals(false, locIsPrimaryProductionVertex, locReactionFullCombo, locDecayCombo, locVertex, locRFBunch, locRFVertexTime, DSourceComboUse(Unknown, 0, nullptr, false, Unknown), locDecayP4, locBeamParticle, locAccuratePhotonsFlag))
			return false;
	}

	if(dDebugLevel >= 10)
		cout << "Calc_P4_Decay: Combo " << locDecayCombo << " P4: " << locDecayP4.Px() << ", " << locDecayP4.Py() << ", " << locDecayP4.Pz() << ", " << locDecayP4.E() << endl;
	return true;
}

bool DSourceComboP4Handler::Calc_P4_HasMassiveNeutrals(bool locIsProductionVertex, bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, DVector3 locVertex, int locRFBunch, double locRFVertexTime, const DSourceComboUse& locToExcludeUse, DLorentzVector& locTotalP4, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag)
{
	auto locP4LookupTuple = std::make_tuple(locIsProductionVertex, locReactionFullCombo, locVertexCombo, locRFBunch, locBeamParticle);
	auto locHasPhotons = !DAnalysis::Get_SourceParticles(locVertexCombo->Get_SourceParticles(true, d_Neutral), Gamma).empty();

	if(!locHasPhotons || !locAccuratePhotonsFlag)
	{
		auto locIterator = dFinalStateP4ByCombo_HasMassiveNeutrals.find(locP4LookupTuple);
		if(locIterator != dFinalStateP4ByCombo_HasMassiveNeutrals.end())
		{
			locTotalP4 = locIterator->second;
			if(dDebugLevel >= 10)
				cout << "Read Calc_P4_HasMassiveNeutrals: Combo " << locVertexCombo << " P4: " << locTotalP4.Px() << ", " << locTotalP4.Py() << ", " << locTotalP4.Pz() << ", " << locTotalP4.E() << endl;
			return true;
		}
	}

	//final state particles
	locTotalP4 = DLorentzVector(0.0, 0.0, 0.0, 0.0);
	auto locParticlesP4 = Calc_P4_SourceParticles(locVertexCombo, locVertex, locRFVertexTime, locAccuratePhotonsFlag);
	if(dDebugLevel >= 20)
		cout << "combo, source total pxyzE = " << locVertexCombo << ", " << locParticlesP4.Px() << ", " << locParticlesP4.Py() << ", " << locParticlesP4.Pz() << ", " << locParticlesP4.E() << endl;
	locTotalP4 += locParticlesP4;

	//now loop over decays
	auto locFurtherDecayCombos = locVertexCombo->Get_FurtherDecayCombos();
	for(const auto& locCombosByUsePair : locFurtherDecayCombos)
	{
		if(locCombosByUsePair.first == locToExcludeUse)
			continue; //e.g. when calc'ing missing p4

		for(const auto& locDecayCombo : locCombosByUsePair.second)
		{
			DLorentzVector locDecayP4;
			if(!Calc_P4_Decay(locIsProductionVertex, locIsPrimaryProductionVertex, locReactionFullCombo, locCombosByUsePair.first, locDecayCombo, locVertex, locRFBunch, locRFVertexTime, locDecayP4, locBeamParticle, locAccuratePhotonsFlag))
				return false;
			if(dDebugLevel >= 20)
				cout << "combo, add decayp4 = " << locDecayCombo << ", " << locDecayP4.Px() << ", " << locDecayP4.Py() << ", " << locDecayP4.Pz() << ", " << locDecayP4.E() << endl;
			locTotalP4 += locDecayP4;
		}
	}

	if(!locAccuratePhotonsFlag || !locHasPhotons)
		dFinalStateP4ByCombo_HasMassiveNeutrals.emplace(locP4LookupTuple, locTotalP4);
	if(dDebugLevel >= 10)
		cout << "Save Calc_P4_HasMassiveNeutrals: Combo " << locVertexCombo << " P4: " << locTotalP4.Px() << ", " << locTotalP4.Py() << ", " << locTotalP4.Pz() << ", " << locTotalP4.E() << endl;
	return true;
}

bool DSourceComboP4Handler::Get_InvariantMassCut(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, bool locAccuratePhotonsFlag, pair<float, float>& locMinMaxMassCuts_GeV) const
{
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return false; //no cut to place!!
	locMinMaxMassCuts_GeV = locCutIterator->second;

	if(locAccuratePhotonsFlag)
		return true;

	auto locNumPhotons = DAnalysis::Get_SourceParticles(locSourceCombo->Get_SourceParticles(true, d_Neutral), Gamma).size();
	if(locNumPhotons == 0)
		return true;
	auto locMassError = (locNumPhotons > 2) ? 0.5*double(locNumPhotons)*d2PhotonInvariantMassCutError : d2PhotonInvariantMassCutError;
	locMinMaxMassCuts_GeV.first -= locMassError;
	locMinMaxMassCuts_GeV.second += locMassError;
	return true;
}

bool DSourceComboP4Handler::Cut_InvariantMass_NoMassiveNeutrals(const DSourceCombo* locVertexCombo, Particle_t locDecayPID, const DVector3& locVertex, signed char locVertexZBin, bool locAccuratePhotonsFlag)
{
	//Z-bin necessary to signal the special negative bins!!
	//Don't call if it contains massive neutrals! Call the other cut function instead!!
	pair<float, float> locMassCuts;
	if(!Get_InvariantMassCut(locVertexCombo, locDecayPID, locAccuratePhotonsFlag, locMassCuts))
		return true; //no cut to place!!

	if(!locVertexCombo->Get_IsComboingZIndependent() && (locVertexZBin == DSourceComboInfo::Get_VertexZIndex_Unknown()) && !locAccuratePhotonsFlag)
		return true; //don't cut yet: will cut later when vertex known or accurate photons
	auto locInvariantMass = Calc_P4_NoMassiveNeutrals(nullptr, locVertexCombo, locVertex, locVertexZBin, nullptr, locAccuratePhotonsFlag).M();

	//save and cut
	auto locCutResult = ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));

	if(!locAccuratePhotonsFlag) //else has already been saved during inaccurate stage
	{
		auto locSaveTuple = std::make_tuple(locDecayPID, locVertexCombo);
		if(dInvariantMassFilledSet.find(locSaveTuple) == dInvariantMassFilledSet.end())
		{
			if(dDebugLevel >= 20)
				cout << "fill no-massive hist" << endl;
			dInvariantMasses[locDecayPID].emplace_back(locInvariantMass);
			dInvariantMassFilledSet.insert(locSaveTuple);
		}
	}
	if(dDebugLevel >= 10)
		cout << "accurate flag, decay pid, z, zbin, mass, cut min/max, pass flag: " << locAccuratePhotonsFlag << ", " << locDecayPID << ", " << locVertex.Z() << ", " << int(locVertexZBin) << ", " << locInvariantMass << ", " << locMassCuts.first << ", " << locMassCuts.second << ", " << locCutResult << endl;
	return locCutResult;
}

bool DSourceComboP4Handler::Cut_InvariantMass_HasMassiveNeutral(bool locIsProductionVertex, bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, Particle_t locDecayPID, double locPrimaryVertexZ, const DVector3& locVertex, double locTimeOffset, vector<int>& locValidRFBunches, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag)
{
	if(locValidRFBunches.empty())
		return true; //massive neutral p4 not defined, can't cut

	pair<float, float> locMassCuts;
	if(!Get_InvariantMassCut(locVertexCombo, locDecayPID, locAccuratePhotonsFlag, locMassCuts))
		return true; //no cut to place!!

	auto locVertexZBin = dSourceComboTimeHandler->Get_PhotonVertexZBin(locVertex.Z());
	if(!locVertexCombo->Get_IsComboingZIndependent() && (locVertexZBin == DSourceComboInfo::Get_VertexZIndex_Unknown()) && !locAccuratePhotonsFlag)
		return true; //don't cut yet: will cut later when vertex known or accurate photons

	//function for calculating and cutting the invariant mass for each rf bunch
	auto CalcAndCut_InvariantMass = [&](int locRFBunch) -> bool
	{
		auto locRFVertexTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, locTimeOffset);

		DLorentzVector locTotalP4(0.0, 0.0, 0.0, 0.0);
		if(!Calc_P4_HasMassiveNeutrals(locIsProductionVertex, locIsPrimaryProductionVertex, locReactionFullCombo, locVertexCombo, locVertex, locRFBunch, locRFVertexTime, DSourceComboUse(Unknown, 0, nullptr, false, Unknown), locTotalP4, locBeamParticle, locAccuratePhotonsFlag))
			return true; //can't cut it yet!

		auto locInvariantMass = locTotalP4.M();
		auto locPassCutFlag = ((locInvariantMass < locMassCuts.first) || (locInvariantMass > locMassCuts.second));
		if(dDebugLevel >= 10)
			cout << "has-mass neutral: accurate flag, decay pid, z, mass, cut min/max, pass flag: " << locAccuratePhotonsFlag << ", " << locDecayPID << ", " << locVertex.Z() << ", " << locTotalP4.M() << ", " << locMassCuts.first << ", " << locMassCuts.second << ", " << locPassCutFlag << endl;

		//save and cut
		if(!locAccuratePhotonsFlag) //else has already been saved during inaccurate stage
		{
			auto locSaveTuple = std::make_tuple(locDecayPID, locIsProductionVertex, locReactionFullCombo, locVertexCombo, locRFBunch);
			if(dInvariantMassFilledSet_MassiveNeutral.find(locSaveTuple) == dInvariantMassFilledSet_MassiveNeutral.end())
			{
				if(dDebugLevel >= 20)
					cout << "fill hist: massive neutrals" << endl;
				dInvariantMasses[locDecayPID].emplace_back(locInvariantMass);
				dInvariantMassFilledSet_MassiveNeutral.insert(locSaveTuple);
			}
		}
		return locPassCutFlag;
	};

	//apply the function
	locValidRFBunches.erase(std::remove_if(locValidRFBunches.begin(), locValidRFBunches.end(), CalcAndCut_InvariantMass), locValidRFBunches.end());
	return !locValidRFBunches.empty();
}

bool DSourceComboP4Handler::Cut_InvariantMass_HasMassiveNeutral_OrPhotonVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, vector<int>& locValidRFBunches)
{
	//do 2 things at once (where vertex is known) (hence the really long function name):
		//calc & cut invariant mass: when massive neutral present
		//calc & cut invariant mass: when vertex-z was unknown with only charged tracks, but is known now, and contains BCAL photons (won't happen very often)
	if(dDebugLevel >= 10)
		cout << "Cut_InvariantMass_HasMassiveNeutral_OrPhotonVertex()" << endl;
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locReactionFullCombo, nullptr).Z();
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();

	//loop over vertices in reverse step order //dependency for calculating invariant mass
	auto locStepVertexInfos = DAnalysis::Get_StepVertexInfos_ReverseOrderByStep(locReactionVertexInfo);
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //unknown position!

		//see if vertex position has been found yet
		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		if(!dSourceComboVertexer->Get_IsVertexKnown_NoBeam(locIsProductionVertex, locVertexPrimaryFullCombo))
			continue; //vertex not found yet!

		auto locVertexComboUse = dSourceComboer->Get_SourceComboUse(locStepVertexInfo);
		auto locHasMassiveNeutrals = dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locVertexComboUse));

		auto locVertexConstrainedByChargedFlag = dSourceComboVertexer->Get_VertexDeterminableWithCharged(locStepVertexInfo);
		if(dDebugLevel >= 10)
			cout << "determinable, has-flag: " << locVertexConstrainedByChargedFlag << ", " << locHasMassiveNeutrals << endl;
		if(locVertexConstrainedByChargedFlag && !locHasMassiveNeutrals)
			continue; //no massive neutrals, and vertex found with charged: mass cuts already placed

		if(locVertexPrimaryFullCombo->Get_IsComboingZIndependent() && !locHasMassiveNeutrals)
			continue; //nothing new (e.g. new vertex found with neutrals, but all of the neutrals in the combo are FCAL photons anyway: nothing to do

		//get vertex position and time offset
		auto locVertex = dSourceComboVertexer->Get_Vertex_NoBeam(locIsProductionVertex, locVertexPrimaryFullCombo);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsPrimaryProductionVertex, locReactionFullCombo, locVertexPrimaryFullCombo, nullptr);

		//get all combos at this vertex, with their uses, in reverse dependency order
		vector<pair<DSourceComboUse, vector<const DSourceCombo*>>> locSourceCombosAndUses_ThisVertex = DAnalysis::Get_SourceCombosAndUses_ThisVertex(locVertexPrimaryFullCombo);
		locSourceCombosAndUses_ThisVertex.emplace(locSourceCombosAndUses_ThisVertex.begin(), locVertexComboUse, vector<const DSourceCombo*>{locVertexPrimaryFullCombo}); //not ideal ...

		//loop over combos & uses at this vertex in dependency order (in reverse)
		for(auto locIterator = locSourceCombosAndUses_ThisVertex.rbegin(); locIterator != locSourceCombosAndUses_ThisVertex.rend(); ++locIterator)
		{
			auto locDecayComboUse = locIterator->first;
			auto locDecayPID = std::get<0>(locDecayComboUse);
			if((locDecayPID == Unknown) || std::get<3>(locDecayComboUse))
				continue; //no mass cut to place!

			auto locDecayHasMassiveNeutrals = dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locDecayComboUse));

			//loop over combos
			for(const auto& locDecayCombo : locIterator->second)
			{
				if(locDecayCombo->Get_IsComboingZIndependent() && !locDecayHasMassiveNeutrals)
					continue; //no massive neutrals, no BCAL showers: already cut!
				if(!Cut_InvariantMass_HasMassiveNeutral(locIsProductionVertex, locIsPrimaryProductionVertex, locReactionFullCombo, locDecayCombo, locDecayPID, locPrimaryVertexZ, locVertex, locTimeOffset, locValidRFBunches, nullptr, false))
					return false; //failed mass cut for all possible rf bunches!
			}
		}
	}

	return true;
}

bool DSourceComboP4Handler::Cut_InvariantMass_MissingMassVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch)
{
	//calc & cut invariant mass: when vertex-z was unknown without the beam energy
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locReactionFullCombo, locBeamParticle).Z();
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();

	//loop over vertices in reverse step order //dependency for calculating invariant mass
	auto locStepVertexInfos = DAnalysis::Get_StepVertexInfos_ReverseOrderByStep(locReactionVertexInfo);
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //unknown position!

		//see if vertex position has been found yet
		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();

		auto locVertexComboUse = dSourceComboer->Get_SourceComboUse(locStepVertexInfo);
		auto locHasMassiveNeutrals = dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locVertexComboUse));

		if(dSourceComboVertexer->Get_VertexDeterminableWithCharged(locStepVertexInfo))
			continue; //vertex found with charged: mass cuts already placed

		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		if(dSourceComboVertexer->Get_VertexDeterminableWithPhotons(locStepVertexInfo))
			continue; //vertex found with photons: mass cuts already placed

		if(locVertexPrimaryFullCombo->Get_IsComboingZIndependent() && !locHasMassiveNeutrals)
			continue; //nothing new (e.g. new vertex found with beam, but all of the neutrals in the combo are FCAL photons anyway: nothing to do

		//get vertex position and time offset
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsPrimaryProductionVertex, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle);

		//get all combos at this vertex, with their uses, in reverse dependency order
		vector<pair<DSourceComboUse, vector<const DSourceCombo*>>> locSourceCombosAndUses_ThisVertex = DAnalysis::Get_SourceCombosAndUses_ThisVertex(locVertexPrimaryFullCombo);
		locSourceCombosAndUses_ThisVertex.emplace(locSourceCombosAndUses_ThisVertex.begin(), locVertexComboUse, vector<const DSourceCombo*>{locVertexPrimaryFullCombo}); //not ideal ...

		//loop over combos & uses at this vertex in dependency order (in reverse)
		for(auto locIterator = locSourceCombosAndUses_ThisVertex.rbegin(); locIterator != locSourceCombosAndUses_ThisVertex.rend(); ++locIterator)
		{
			auto locDecayComboUse = locIterator->first;
			auto locDecayPID = std::get<0>(locDecayComboUse);
			if((locDecayPID == Unknown) || std::get<3>(locDecayComboUse))
				continue; //no mass cut to place!
			auto locDecayHasMassiveNeutrals = dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locDecayComboUse));

			//loop over combos
			for(const auto& locDecayCombo : locIterator->second)
			{
				if(locDecayCombo->Get_IsComboingZIndependent() && !locDecayHasMassiveNeutrals)
					continue; //no massive neutrals, no BCAL showers: already cut!
				vector<int> locRFBunches{locRFBunch};
				if(!Cut_InvariantMass_HasMassiveNeutral(locIsProductionVertex, locIsPrimaryProductionVertex, locReactionFullCombo, locDecayCombo, locDecayPID, locPrimaryVertexZ, locVertex, locTimeOffset, locRFBunches, nullptr, false))
					return false; //failed mass cut for all possible rf bunches!
			}
		}
	}

	return true;
}

bool DSourceComboP4Handler::Cut_MissingMass(const DReaction* locReaction, const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch)
{
	if(dDebugLevel >= 5)
		cout << "DSourceComboP4Handler::Cut_MissingMass()" << endl;

	//see if can evaluate
	if(locReaction->Get_IsInclusiveFlag() || (locBeamParticle == nullptr))
		return true; //cannot cut
	auto locStepVertexInfo = locReactionVertexInfo->Get_StepVertexInfo(0);
	if(!locStepVertexInfo->Get_ProductionVertexFlag())
		return true; //initial decaying particle: can't compute missing mass
	if(locStepVertexInfo->Get_DanglingVertexFlag())
		return true; //unknown position!

	//get mass cuts
	auto locMissingPIDs = locReaction->Get_MissingPIDs();
	if(locMissingPIDs.size() > 1)
		return true; //cut has no meaning
	auto locMissingPID = locMissingPIDs.empty() ? Unknown : locMissingPIDs[0];
	auto locCutIterator = dMissingMassSquaredCuts.find(locMissingPID);
	if(locCutIterator == dMissingMassSquaredCuts.end())
		return true; //no cut defined
	auto locCutPair = locCutIterator->second;

	if(dDebugLevel >= 5)
		cout << "Cuts retrieved for missing particle: " << locMissingPID << endl;

	//get primary vertex info
	auto locPrimaryVertex = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locReactionFullCombo, locBeamParticle);
	auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(true, locReactionFullCombo, locReactionFullCombo, locBeamParticle);
	auto locRFVertexTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertex.Z(), locRFBunch, locTimeOffset);

	//compute final-state p4
	DLorentzVector locFinalStateP4(0.0, 0.0, 0.0, 0.0);
	//it may not have massive neutrals, but this function is the worst-case (hardest, most-complete) scenario
	if(!Calc_P4_HasMassiveNeutrals(true, true, locReactionFullCombo, locReactionFullCombo, locPrimaryVertex, locRFBunch, locRFVertexTime, DSourceComboUse(Unknown, 0, nullptr, false, Unknown), locFinalStateP4, locBeamParticle, true))
		return true; //can't cut it yet!

	//compute missing p4
	//ASSUMES FIXED TARGET EXPERIMENT!
	auto locTargetPID = locReaction->Get_ReactionStep(0)->Get_TargetPID();
	auto locMissingP4 = locBeamParticle->lorentzMomentum() + DLorentzVector(0.0, 0.0, 0.0, ParticleMass(locTargetPID)) - locFinalStateP4;

	auto locBeamEnergy = locBeamParticle->lorentzMomentum().E();
	auto locMissingMassSquared_Min = locCutPair.first->Eval(locBeamEnergy);
	auto locMissingMassSquared_Max = locCutPair.second->Eval(locBeamEnergy);
	if(dDebugLevel >= 5)
	{
		cout << "missing pxyzE = " << locMissingP4.Px() << ", " << locMissingP4.Py() << ", " << locMissingP4.Pz() << ", " << locMissingP4.E() << endl;
		cout << "beam E, missing mass^2 min/max/measured = " << locBeamEnergy << ", " << locMissingMassSquared_Min << ", " << locMissingMassSquared_Max << ", " << locMissingP4.M2() << endl;
	}

	//save and cut
	dMissingMassPairs[locMissingPID].emplace_back(locBeamEnergy, locMissingP4.M2());
	auto locPassCutFlag = ((locMissingP4.M2() >= locMissingMassSquared_Min) && (locMissingP4.M2() <= locMissingMassSquared_Max));
	if(locMissingPID == Unknown)
	{
		dMissingEVsBeamEnergy_PreMissMassSqCut.emplace_back(locBeamEnergy, locMissingP4.E());
		dMissingPtVsMissingE_PreMissMassSqCut.emplace_back(locMissingP4.E(), locMissingP4.Perp());
		if(locPassCutFlag)
		{
			if((dMissingECuts.first != nullptr) && (dMissingECuts.second != nullptr)) //cut missing E!!
				locPassCutFlag = ((locMissingP4.E() >= dMissingECuts.first->Eval(locBeamEnergy)) && (locMissingP4.E() <= dMissingECuts.second->Eval(locBeamEnergy)));
			dMissingEVsBeamEnergy_PostMissMassSqCut.emplace_back(locBeamEnergy, locMissingP4.E());
			dMissingPtVsMissingE_PostMissMassSqCut.emplace_back(locMissingP4.E(), locMissingP4.Perp());
		}
	}
	return locPassCutFlag;
}

bool DSourceComboP4Handler::Cut_InvariantMass_AccuratePhotonKinematics(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch)
{
	if(dDebugLevel >= 5)
		cout << "DSourceComboP4Handler::Cut_InvariantMass_CorrectPhotonKinematics()" << endl;

	if(DAnalysis::Get_SourceParticles(locReactionFullCombo->Get_SourceParticles(true, d_Neutral), Gamma).empty())
		return true; //nothing has changed!

	//if photons in combo, cut again using correct p4 (rather than estimated)
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locReactionFullCombo, nullptr).Z();
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();

	//loop over vertices in reverse step order //dependency for calculating invariant mass
	auto locStepVertexInfos = DAnalysis::Get_StepVertexInfos_ReverseOrderByStep(locReactionVertexInfo);
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //unknown position!

		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		if(DAnalysis::Get_SourceParticles(locVertexPrimaryFullCombo->Get_SourceParticles(true, d_Neutral), Gamma).empty())
			continue; //nothing has changed

		//get vertex position and time offset
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsPrimaryProductionVertex, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle);

		//get all combos at this vertex, with their uses, in reverse dependency order
		auto locVertexComboUse = dSourceComboer->Get_SourceComboUse(locStepVertexInfo);
		vector<pair<DSourceComboUse, vector<const DSourceCombo*>>> locSourceCombosAndUses_ThisVertex = DAnalysis::Get_SourceCombosAndUses_ThisVertex(locVertexPrimaryFullCombo);
		locSourceCombosAndUses_ThisVertex.emplace(locSourceCombosAndUses_ThisVertex.begin(), locVertexComboUse, vector<const DSourceCombo*>{locVertexPrimaryFullCombo}); //not ideal ...

		//loop over combos & uses at this vertex in dependency order (in reverse)
		for(auto locIterator = locSourceCombosAndUses_ThisVertex.rbegin(); locIterator != locSourceCombosAndUses_ThisVertex.rend(); ++locIterator)
		{
			auto locDecayComboUse = locIterator->first;
			auto locDecayPID = std::get<0>(locDecayComboUse);
			if((locDecayPID == Unknown) || std::get<3>(locDecayComboUse))
				continue; //no mass cut to place!
			auto locDecayHasMassiveNeutrals = dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locDecayComboUse));

			//loop over combos
			for(const auto& locDecayCombo : locIterator->second)
			{
				if(DAnalysis::Get_SourceParticles(locDecayCombo->Get_SourceParticles(true, d_Neutral), Gamma).empty())
					continue; //nothing has changed
				if(locDecayHasMassiveNeutrals)
				{
					vector<int> locBeamBunches{locRFBunch};
					if(!Cut_InvariantMass_HasMassiveNeutral(locIsProductionVertex, locIsPrimaryProductionVertex, locReactionFullCombo, locDecayCombo, locDecayPID, locPrimaryVertexZ, locVertex, locTimeOffset, locBeamBunches, locBeamParticle, true))
						return false;
				}
				else
				{
					if(!Cut_InvariantMass_NoMassiveNeutrals(locDecayCombo, locDecayPID, locVertex, std::get<1>(locDecayComboUse), true))
						return false;
				}
			}
		}
	}
	return true;
}

void DSourceComboP4Handler::Fill_Histograms(void)
{
	japp->WriteLock("DSourceComboP4Handler");
	{
		for(auto& locPIDPair : dInvariantMasses)
		{
			for(auto& locMass : locPIDPair.second)
				dHistMap_InvariantMass[locPIDPair.first]->Fill(locMass);
		}
		for(auto& locPIDPair : dMissingMassPairs)
		{
			for(auto& locMassPair : locPIDPair.second)
				dHistMap_MissingMassSquaredVsBeamEnergy[locPIDPair.first]->Fill(locMassPair.first, locMassPair.second);
		}
		for(auto& locSystemPair : d2GammaInvariantMasses)
		{
			for(auto& locMass : locSystemPair.second)
				dHistMap_2GammaMass[locSystemPair.first]->Fill(locMass);
		}

		for(auto& locMissingPair : dMissingEVsBeamEnergy_PreMissMassSqCut)
			dHist_NoneMissing_MissingEVsBeamEnergy_PreMissMassSqCut->Fill(locMissingPair.first, locMissingPair.second);
		for(auto& locMissingPair : dMissingEVsBeamEnergy_PostMissMassSqCut)
			dHist_NoneMissing_MissingEVsBeamEnergy_PostMissMassSqCut->Fill(locMissingPair.first, locMissingPair.second);
		for(auto& locMissingPair : dMissingPtVsMissingE_PreMissMassSqCut)
			dHist_NoneMissing_MissingPtVsMissingE_PreMissMassSqCut->Fill(locMissingPair.first, locMissingPair.second);
		for(auto& locMissingPair : dMissingPtVsMissingE_PostMissMassSqCut)
			dHist_NoneMissing_MissingPtVsMissingE_PostMissMassSqCut->Fill(locMissingPair.first, locMissingPair.second);
	}
	japp->Unlock("DSourceComboP4Handler");

	vector<pair<float, float>>().swap(dMissingEVsBeamEnergy_PreMissMassSqCut);
	vector<pair<float, float>>().swap(dMissingEVsBeamEnergy_PostMissMassSqCut);
	vector<pair<float, float>>().swap(dMissingPtVsMissingE_PreMissMassSqCut);
	vector<pair<float, float>>().swap(dMissingPtVsMissingE_PostMissMassSqCut);

	for(auto& locPIDPair : dInvariantMasses)
		decltype(locPIDPair.second)().swap(locPIDPair.second);
	for(auto& locPIDPair : dMissingMassPairs)
		decltype(locPIDPair.second)().swap(locPIDPair.second);
	for(auto& locSystemPair : d2GammaInvariantMasses)
		decltype(locSystemPair.second)().swap(locSystemPair.second);
}

} //end namespace


// $Id$
//
//    File: DMCThrownMatching_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DMCThrownMatching_factory.h"
#include "TMath.h"

DMCThrownMatching_factory::DMCThrownMatching_factory(void)
{
	dMCThrownComparisonPIDs.push_back(Gamma);
	dMCThrownComparisonPIDs.push_back(Positron);
	dMCThrownComparisonPIDs.push_back(Electron);
	dMCThrownComparisonPIDs.push_back(MuonPlus);
	dMCThrownComparisonPIDs.push_back(MuonMinus);
	dMCThrownComparisonPIDs.push_back(Neutrino);
	dMCThrownComparisonPIDs.push_back(PiPlus);
	dMCThrownComparisonPIDs.push_back(PiMinus);
	dMCThrownComparisonPIDs.push_back(KPlus);
	dMCThrownComparisonPIDs.push_back(KMinus);
	dMCThrownComparisonPIDs.push_back(Neutron);
	dMCThrownComparisonPIDs.push_back(Proton);
	dMCThrownComparisonPIDs.push_back(AntiProton);
}

//------------------
// init
//------------------
jerror_t DMCThrownMatching_factory::init(void)
{
	dMinimumMatchFOM = -1.0;
	dDebugLevel = 0;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DMCThrownMatching_factory::brun(jana::JEventLoop* locEventLoop, int runnumber)
{
	gPARMS->SetDefaultParameter("MCMATCH:MINMATCHFOM", dMinimumMatchFOM);
	gPARMS->SetDefaultParameter("MCMATCH:DEBUGLEVEL", dDebugLevel);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCThrownMatching_factory::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
 	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);
 	vector<const DMCThrown*> locOriginalMCThrowns_Charged;
 	vector<const DMCThrown*> locOriginalMCThrowns_Neutral;

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		if(!Check_IsValidMCComparisonPID(locMCThrowns, locMCThrowns[loc_i]))
			continue;
		if(ParticleCharge((Particle_t)locMCThrowns[loc_i]->type) == 0)
			locOriginalMCThrowns_Neutral.push_back(locMCThrowns[loc_i]);
		else
			locOriginalMCThrowns_Charged.push_back(locMCThrowns[loc_i]);
	}
	if(dDebugLevel > 0)
		cout << "input #thrown, ok charged # thrown, ok neutral # thrown = " << locMCThrowns.size() << ", " << locOriginalMCThrowns_Charged.size() << ", " << locOriginalMCThrowns_Neutral.size() << endl;

	if(locMCThrowns.empty())
		return NOERROR;

 	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

 	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);

 	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses);

 	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses_Reaction;
	locEventLoop->Get(locChargedTrackHypotheses_Reaction, "Reaction");
	locChargedTrackHypotheses.insert(locChargedTrackHypotheses.end(), locChargedTrackHypotheses_Reaction.begin(), locChargedTrackHypotheses_Reaction.end());

 	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses);

	DMCThrownMatching* locMCThrownMatching = new DMCThrownMatching();

	Find_GenReconMatches_ChargedHypo(locOriginalMCThrowns_Charged, locChargedTrackHypotheses, locMCThrownMatching);
	Find_GenReconMatches_ChargedTrack(locOriginalMCThrowns_Charged, locChargedTracks, locMCThrownMatching);

	Find_GenReconMatches_NeutralHypo(locOriginalMCThrowns_Neutral, locNeutralParticleHypotheses, locMCThrownMatching);
	Find_GenReconMatches_NeutralParticle(locOriginalMCThrowns_Neutral, locNeutralParticles, locMCThrownMatching);

	if(dDebugLevel > 0)
	{
		cout << "Charged Track Matching Summary:" << endl;
		for(size_t loc_i = 0; loc_i < locChargedTrackHypotheses.size(); ++loc_i)
		{
			double locP = locChargedTrackHypotheses[loc_i]->momentum().Mag();
			double locTheta = locChargedTrackHypotheses[loc_i]->momentum().Theta()*180.0/TMath::Pi();
			double locPhi = locChargedTrackHypotheses[loc_i]->momentum().Phi()*180.0/TMath::Pi();
			Particle_t locPID = locChargedTrackHypotheses[loc_i]->PID();
			cout << "charged info: " << locChargedTrackHypotheses[loc_i]->candidateid << ", " << ParticleType(locPID) << ", " << locP << ", " << locTheta << ", " << locPhi << endl;
			cout << "matched info: ";
			const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypotheses[loc_i]);
			if(locMCThrown == NULL)
			{
				cout << "NO MATCH." << endl;
				continue;
			}
			locP = locMCThrown->momentum().Mag();
			locTheta = locMCThrown->momentum().Theta()*180.0/TMath::Pi();
			locPhi = locMCThrown->momentum().Phi()*180.0/TMath::Pi();
			locPID = locMCThrown->PID();
			cout << ParticleType(locPID) << ", " << locP << ", " << locTheta << ", " << locPhi << endl;
		}
		cout << "Unmatched Charged DMCThrowns:" << endl;
		for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		{
			if(ParticleCharge(locMCThrowns[loc_i]->PID()) == 0)
				continue;
			if(locMCThrownMatching->Get_MatchingChargedTrack(locMCThrowns[loc_i]) != NULL)
				continue;
			double locP = locMCThrowns[loc_i]->momentum().Mag();
			double locTheta = locMCThrowns[loc_i]->momentum().Theta()*180.0/TMath::Pi();
			double locPhi = locMCThrowns[loc_i]->momentum().Phi()*180.0/TMath::Pi();
			Particle_t locPID = locMCThrowns[loc_i]->PID();
			cout << "thrown info: " << ParticleType(locPID) << ", " << locP << ", " << locTheta << ", " << locPhi << endl;
		}
	}

	_data.push_back(locMCThrownMatching);

	return NOERROR;
}

bool DMCThrownMatching_factory::Check_IsValidMCComparisonPID(const vector<const DMCThrown*>& locAllMCThrowns, const DMCThrown* locMCThrown) const
{
	Particle_t locPID = (Particle_t)locMCThrown->type;
	Particle_t locParentPID = Unknown;
	int locParentID = locMCThrown->parentid;
	if(locParentID != 0)
	{
		for(size_t loc_i = 0; loc_i < locAllMCThrowns.size(); ++loc_i)
		{
			if(locAllMCThrowns[loc_i]->myid != locParentID)
				continue;
			locParentPID = (Particle_t)locAllMCThrowns[loc_i]->type;
			break;
		}
		if(locParentPID == Unknown)
			return false; //parent not listed or totally unknown
	}

	bool locValidCompareTypeFlag = false;
	for(size_t loc_j = 0; loc_j < dMCThrownComparisonPIDs.size(); ++loc_j)
	{
		if(dMCThrownComparisonPIDs[loc_j] == locParentPID) //e.g. on mu+, and pi+ is parent
			return false;
		if(dMCThrownComparisonPIDs[loc_j] == locPID)
			locValidCompareTypeFlag = true;
	}
	return locValidCompareTypeFlag;
}

void DMCThrownMatching_factory::Find_GenReconMatches_ChargedTrack(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DChargedTrack*>& locChargedTracks, DMCThrownMatching* locMCThrownMatching) const
{
	//assumes Find_GenReconMatches_ChargedHypo has been called first!
	map<const DChargedTrackHypothesis*, const DMCThrown*> locChargedHypoToThrownMap;
	map<const DMCThrown*, deque<const DChargedTrackHypothesis*> > locThrownToChargedHypoMap;
	locMCThrownMatching->Get_ChargedHypoToThrownMap(locChargedHypoToThrownMap);
	locMCThrownMatching->Get_ThrownToChargedHypoMap(locThrownToChargedHypoMap);

	map<const DChargedTrack*, const DMCThrown*> locChargedToThrownMap;
	map<const DMCThrown*, const DChargedTrack*> locThrownToChargedMap;

	map<const DMCThrown*, deque<const DChargedTrackHypothesis*> >::iterator locIterator;
	for(locIterator = locThrownToChargedHypoMap.begin(); locIterator != locThrownToChargedHypoMap.end(); ++locIterator)
	{
		for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j)
		{
			if(locChargedTracks[loc_j]->Get_BestFOM()->candidateid == (locIterator->second)[0]->candidateid) //match
			{
				locChargedToThrownMap[locChargedTracks[loc_j]] = locIterator->first;
				locThrownToChargedMap[locIterator->first] = locChargedTracks[loc_j];
			}
		}
	}

	locMCThrownMatching->Set_ChargedToThrownMap(locChargedToThrownMap);
	locMCThrownMatching->Set_ThrownToChargedMap(locThrownToChargedMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_ChargedHypo(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DChargedTrackHypothesis*>& locInputChargedTrackHypothesisVector, DMCThrownMatching* locMCThrownMatching) const
{
	map<const DChargedTrackHypothesis*, const DMCThrown*> locChargedToThrownMap;
	map<const DMCThrown*, deque<const DChargedTrackHypothesis*> > locThrownToChargedMap;

	const DChargedTrackHypothesis* locChargedTrackHypothesis;
	const DMCThrown* locMCThrown;
	size_t locBestChargedTrackHypothesisIndex = 0, locBestMCThrownIndex = 0;
	double locMatchFOM, locBestMatchFOM;
	vector<const DMCThrown*> locMCThrownVector = locInputMCThrownVector;
	vector<const DChargedTrackHypothesis*> locChargedTrackHypothesisVector = locInputChargedTrackHypothesisVector;
	if(dDebugLevel > 0)
		cout << "START IT!" << endl;
	while((!locMCThrownVector.empty()) && (!locChargedTrackHypothesisVector.empty()))
	{
		if(dDebugLevel > 0)
			cout << "Begin loop!" << endl;
		bool locMatchFoundFlag = false;
		locBestMatchFOM = dMinimumMatchFOM;
		for(size_t loc_i = 0; loc_i < locMCThrownVector.size(); ++loc_i)
		{
			locMCThrown = locMCThrownVector[loc_i];
			for(size_t loc_j = 0; loc_j < locChargedTrackHypothesisVector.size(); ++loc_j)
			{
				locChargedTrackHypothesis = locChargedTrackHypothesisVector[loc_j];
				if(ParticleCharge(locChargedTrackHypothesis->PID()) != ParticleCharge((Particle_t)(locMCThrown->type)))
					continue; //wrong charge
				locMatchFOM = Calc_MatchFOM(locMCThrown->momentum(), locChargedTrackHypothesis->momentum());

				if(dDebugLevel > 0)
				{
					cout << "MATCHING: MCTHROWN: ";
					cout << ParticleType((Particle_t)(locMCThrown->type)) << ", " << locMCThrown->momentum().Mag() << ", " << locMCThrown->momentum().Theta()*180.0/TMath::Pi() << ", " << locMCThrown->momentum().Phi()*180.0/TMath::Pi() << endl;
					cout << "MATCHING: CHARGEDHYPO: ";
					cout << ParticleType(locChargedTrackHypothesis->PID()) << ", " << locChargedTrackHypothesis->momentum().Mag() << ", " << locChargedTrackHypothesis->momentum().Theta()*180.0/TMath::Pi() << ", " << locChargedTrackHypothesis->momentum().Phi()*180.0/TMath::Pi() << endl;
					cout << "MATCHING: FOM, candidate id: " << locMatchFOM << ", " << locChargedTrackHypothesis->candidateid << endl;
				}

				if(locMatchFOM >= locBestMatchFOM)
				{
					locMatchFoundFlag = true;
					locBestMatchFOM = locMatchFOM;
					locBestMCThrownIndex = loc_i;
					locBestChargedTrackHypothesisIndex = loc_j;
				}
			}
		}

		if(!locMatchFoundFlag) //no more good matches!1
			break;

		locMCThrown = locMCThrownVector[locBestMCThrownIndex];
		locChargedTrackHypothesis = locChargedTrackHypothesisVector[locBestChargedTrackHypothesisIndex];

		locChargedToThrownMap[locChargedTrackHypothesis] = locMCThrown;
		locChargedTrackHypothesisVector.erase(locChargedTrackHypothesisVector.begin() + locBestChargedTrackHypothesisIndex);
		locMCThrownVector.erase(locMCThrownVector.begin() + locBestMCThrownIndex);

		//automatically add all other DChargedTrackHypothesis objects from the same DChargedTrack to this match.
		deque<const DChargedTrackHypothesis*> locMatchedChargedHypos(1, locChargedTrackHypothesis);
		for(int loc_i = locChargedTrackHypothesisVector.size() - 1; loc_i >= 0; --loc_i)
		{
			if(dDebugLevel > 0)
			{
				cout << "CHECKING: CHARGEDHYPO: ";
				cout << ParticleType(locChargedTrackHypothesisVector[loc_i]->PID()) << ", " << locChargedTrackHypothesisVector[loc_i]->momentum().Mag() << ", " << locChargedTrackHypothesisVector[loc_i]->momentum().Theta()*180.0/TMath::Pi() << ", " << locChargedTrackHypothesisVector[loc_i]->momentum().Phi()*180.0/TMath::Pi() << endl;
				cout << "best id, test id = " << locChargedTrackHypothesis->candidateid << ", " << locChargedTrackHypothesisVector[loc_i]->candidateid << endl;
			}
			if(locChargedTrackHypothesisVector[loc_i]->candidateid == locChargedTrackHypothesis->candidateid)
			{
				if(dDebugLevel > 0)
					cout << "save!" << endl;
				locChargedToThrownMap[locChargedTrackHypothesisVector[loc_i]] = locMCThrown;
				locMatchedChargedHypos.push_back(locChargedTrackHypothesisVector[loc_i]);
				locChargedTrackHypothesisVector.erase(locChargedTrackHypothesisVector.begin() + loc_i);
			}
		}
		locThrownToChargedMap[locMCThrown] = locMatchedChargedHypos;
	}

	locMCThrownMatching->Set_ChargedHypoToThrownMap(locChargedToThrownMap);
	locMCThrownMatching->Set_ThrownToChargedHypoMap(locThrownToChargedMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_NeutralParticle(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DNeutralParticle*>& locNeutralParticles, DMCThrownMatching* locMCThrownMatching) const
{
	//assumes Find_GenReconMatches_NeutralHypo has been called first!
	map<const DNeutralParticleHypothesis*, const DMCThrown*> locNeutralHypoToThrownMap;
	map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> > locThrownToNeutralHypoMap;
	locMCThrownMatching->Get_NeutralHypoToThrownMap(locNeutralHypoToThrownMap);
	locMCThrownMatching->Get_ThrownToNeutralHypoMap(locThrownToNeutralHypoMap);

	vector<const DNeutralShower*> locNeutralShowerVector_Matched;
	vector<const DNeutralShower*> locNeutralShowerVector_Check;

	map<const DNeutralParticle*, const DMCThrown*> locNeutralToThrownMap;
	map<const DMCThrown*, const DNeutralParticle*> locThrownToNeutralMap;

	map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> >::iterator locIterator;
	for(locIterator = locThrownToNeutralHypoMap.begin(); locIterator != locThrownToNeutralHypoMap.end(); ++locIterator)
	{
		(locIterator->second)[0]->GetT(locNeutralShowerVector_Matched);
		for(size_t loc_j = 0; loc_j < locNeutralParticles.size(); ++loc_j)
		{
			locNeutralParticles[loc_j]->GetT(locNeutralShowerVector_Check);
			if(locNeutralShowerVector_Check[0] == locNeutralShowerVector_Matched[0])
			{
				locNeutralToThrownMap[locNeutralParticles[loc_j]] = locIterator->first;
				locThrownToNeutralMap[locIterator->first] = locNeutralParticles[loc_j];
			}
		}
	}

	locMCThrownMatching->Set_NeutralToThrownMap(locNeutralToThrownMap);
	locMCThrownMatching->Set_ThrownToNeutralMap(locThrownToNeutralMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_NeutralHypo(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DNeutralParticleHypothesis*>& locInputNeutralParticleHypothesisVector, DMCThrownMatching* locMCThrownMatching) const
{
	map<const DNeutralParticleHypothesis*, const DMCThrown*> locNeutralToThrownMap;
	map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> > locThrownToNeutralMap;

	const DNeutralParticleHypothesis* locNeutralParticleHypothesis;
	const DMCThrown* locMCThrown;
	size_t locBestNeutralParticleHypothesisIndex = 0, locBestMCThrownIndex = 0;
	double locMatchFOM, locBestMatchFOM;
	vector<const DMCThrown*> locMCThrownVector = locInputMCThrownVector;
	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypothesisVector = locInputNeutralParticleHypothesisVector;

	while((!locMCThrownVector.empty()) && (!locNeutralParticleHypothesisVector.empty()))
	{
		bool locMatchFoundFlag = false;
		locBestMatchFOM = dMinimumMatchFOM;
		for(size_t loc_i = 0; loc_i < locMCThrownVector.size(); ++loc_i)
		{
			locMCThrown = locMCThrownVector[loc_i];
			for(size_t loc_j = 0; loc_j < locNeutralParticleHypothesisVector.size(); ++loc_j)
			{
				locNeutralParticleHypothesis = locNeutralParticleHypothesisVector[loc_j];
				locMatchFOM = Calc_MatchFOM(locMCThrown->momentum(), locNeutralParticleHypothesis->momentum());
				if(locMatchFOM >= locBestMatchFOM)
				{
					locMatchFoundFlag = true;
					locBestMatchFOM = locMatchFOM;
					locBestMCThrownIndex = loc_i;
					locBestNeutralParticleHypothesisIndex = loc_j;
				}
			}
		}

		if(!locMatchFoundFlag) //no more good matches!1
			break;

		locMCThrown = locMCThrownVector[locBestMCThrownIndex];
		locNeutralParticleHypothesis = locNeutralParticleHypothesisVector[locBestNeutralParticleHypothesisIndex];

		locNeutralToThrownMap[locNeutralParticleHypothesis] = locMCThrown;
		locNeutralParticleHypothesisVector.erase(locNeutralParticleHypothesisVector.begin() + locBestNeutralParticleHypothesisIndex);
		locMCThrownVector.erase(locMCThrownVector.begin() + locBestMCThrownIndex);

		//automatically add all other DNeutralParticleHypothesis objects from the same DNeutralShower to this match.
		deque<const DNeutralParticleHypothesis*> locMatchedNeutralHypos(1, locNeutralParticleHypothesis);
		vector<const DNeutralShower*> locNeutralShowerVector_Matched;
		vector<const DNeutralShower*> locNeutralShowerVector_Check;
		locNeutralParticleHypothesis->GetT(locNeutralShowerVector_Matched);
		for(int loc_i = locNeutralParticleHypothesisVector.size() - 1; loc_i >= 0; --loc_i)
		{
			locNeutralParticleHypothesisVector[loc_i]->GetT(locNeutralShowerVector_Check);
			if(locNeutralShowerVector_Check[0] == locNeutralShowerVector_Matched[0])
			{
				locNeutralToThrownMap[locNeutralParticleHypothesisVector[loc_i]] = locMCThrown;
				locMatchedNeutralHypos.push_back(locNeutralParticleHypothesisVector[loc_i]);
				locNeutralParticleHypothesisVector.erase(locNeutralParticleHypothesisVector.begin() + loc_i);
			}
		}
		locThrownToNeutralMap[locMCThrown] = locMatchedNeutralHypos;
	}

	locMCThrownMatching->Set_NeutralHypoToThrownMap(locNeutralToThrownMap);
	locMCThrownMatching->Set_ThrownToNeutralHypoMap(locThrownToNeutralMap);
}

double DMCThrownMatching_factory::Calc_MatchFOM(const DVector3& locMomentum_Thrown, const DVector3& locMomentum_Detected) const
{
	double locDeltaPOverP = (locMomentum_Thrown.Mag() > 0.0) ? (locMomentum_Thrown.Mag() - locMomentum_Detected.Mag())/locMomentum_Thrown.Mag() : 9.9E9; //mean is zero
	double locDeltaPOverPSigma = 0.03;
	double locDeltaPOverPChiSq = locDeltaPOverP*locDeltaPOverP/(locDeltaPOverPSigma*locDeltaPOverPSigma);

	double locDeltaTheta = locMomentum_Thrown.Angle(locMomentum_Detected); //mean is zero
	double locDeltaThetaSigma = 0.5*TMath::Pi()/180.0; //0.5 degrees: phi dominates (theta is ~0.03 degrees)
	double locDeltaThetaChiSq = locDeltaTheta*locDeltaTheta/(locDeltaThetaSigma*locDeltaThetaSigma);

	double locFOM = TMath::Prob(locDeltaPOverPChiSq + locDeltaThetaChiSq, 2);
	if(dDebugLevel > 10)
		cout << "delta p over p, delta theta, total FOM = " << locDeltaPOverP << ", " << locDeltaTheta << ", " << locFOM << endl;
	return locFOM;
}

//------------------
// erun
//------------------
jerror_t DMCThrownMatching_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DMCThrownMatching_factory::fini(void)
{
	return NOERROR;
}



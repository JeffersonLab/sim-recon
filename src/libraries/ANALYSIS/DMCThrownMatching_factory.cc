// $Id$
//
//    File: DMCThrownMatching_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DMCThrownMatching_factory.h"
#include "TMath.h"

//------------------
// init
//------------------
jerror_t DMCThrownMatching_factory::init(void)
{
	dMinimumMatchFOM = -1.0;
	dDebugLevel = 0;
	dMaximumTOFMatchDistance = 10.0; //cm
	dMaximumFCALMatchDistance = 10.0; //cm
	dMaximumBCALMatchAngleDegrees = 5.0;
	dTargetCenter = 65.0; //set me from geometry!
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DMCThrownMatching_factory::brun(jana::JEventLoop* locEventLoop, int runnumber)
{
	gPARMS->SetDefaultParameter("MCMATCH:MINMATCHFOM", dMinimumMatchFOM);
	gPARMS->SetDefaultParameter("MCMATCH:DEBUGLEVEL", dDebugLevel);
	gPARMS->SetDefaultParameter("MCMATCH:MAXMATCHTOFDISTANCE", dMaximumTOFMatchDistance);
	gPARMS->SetDefaultParameter("MCMATCH:MAXMATCHFCALDISTANCE", dMaximumFCALMatchDistance);
	gPARMS->SetDefaultParameter("MCMATCH:MAXMATCHBCALANGLE", dMaximumBCALMatchAngleDegrees);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCThrownMatching_factory::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
 	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns, "FinalState");
 	vector<const DMCThrown*> locOriginalMCThrowns_Charged;
 	vector<const DMCThrown*> locOriginalMCThrowns_Neutral;

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
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

 	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses);

	DMCThrownMatching* locMCThrownMatching = new DMCThrownMatching();

	Find_GenReconMatches_ChargedHypo(locOriginalMCThrowns_Charged, locChargedTrackHypotheses, locMCThrownMatching);
	Find_GenReconMatches_ChargedTrack(locOriginalMCThrowns_Charged, locChargedTracks, locMCThrownMatching);

	Find_GenReconMatches_NeutralHypo(locOriginalMCThrowns_Neutral, locNeutralParticleHypotheses, locMCThrownMatching);
	Find_GenReconMatches_NeutralParticle(locOriginalMCThrowns_Neutral, locNeutralParticles, locMCThrownMatching);

	Find_GenReconMatches_TOFPoints(locEventLoop, locMCThrownMatching);
	Find_GenReconMatches_BCALShowers(locEventLoop, locMCThrownMatching);
	Find_GenReconMatches_FCALShowers(locEventLoop, locMCThrownMatching);

	if(dDebugLevel > 0)
	{
		cout << "Charged Track Matching Summary:" << endl;
		double locMatchFOM = 0.0;
		for(size_t loc_i = 0; loc_i < locChargedTrackHypotheses.size(); ++loc_i)
		{
			double locP = locChargedTrackHypotheses[loc_i]->momentum().Mag();
			double locTheta = locChargedTrackHypotheses[loc_i]->momentum().Theta()*180.0/TMath::Pi();
			double locPhi = locChargedTrackHypotheses[loc_i]->momentum().Phi()*180.0/TMath::Pi();
			Particle_t locPID = locChargedTrackHypotheses[loc_i]->PID();
			cout << "charged info: " << locChargedTrackHypotheses[loc_i]->candidateid << ", " << ParticleType(locPID) << ", " << locP << ", " << locTheta << ", " << locPhi << endl;
			cout << "matched info: ";
			const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypotheses[loc_i], locMatchFOM);
			if(locMCThrown == NULL)
			{
				cout << "NO MATCH." << endl;
				continue;
			}
			locP = locMCThrown->momentum().Mag();
			locTheta = locMCThrown->momentum().Theta()*180.0/TMath::Pi();
			locPhi = locMCThrown->momentum().Phi()*180.0/TMath::Pi();
			locPID = locMCThrown->PID();
			cout << ParticleType(locPID) << ", " << locP << ", " << locTheta << ", " << locPhi << ", " << locMatchFOM << endl;
		}
		cout << "Unmatched Charged DMCThrowns:" << endl;
		for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		{
			if(ParticleCharge(locMCThrowns[loc_i]->PID()) == 0)
				continue;
			if(locMCThrownMatching->Get_MatchingChargedTrack(locMCThrowns[loc_i], locMatchFOM) != NULL)
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

void DMCThrownMatching_factory::Find_GenReconMatches_FCALShowers(JEventLoop* locEventLoop, DMCThrownMatching* locMCThrownMatching) const
{
	vector<const DFCALTruthShower*> locFCALTruthShowers;
	const DFCALTruthShower* locFCALTruthShower;
	locEventLoop->Get(locFCALTruthShowers);

	vector<const DFCALShower*> locFCALShowers;
	const DFCALShower* locFCALShower;
	locEventLoop->Get(locFCALShowers);

	map<const DFCALShower*, pair<const DFCALTruthShower*, double> > locFCALShowerToTruthMap;
	map<const DFCALTruthShower*, pair<const DFCALShower*, double> > locFCALTruthToShowerMap;

	size_t locBestFCALShowerIndex = 0, locBestFCALTruthShowerIndex = 0;
	double locBestMatchDistance, locMatchDistance;
	while((!locFCALShowers.empty()) && (!locFCALTruthShowers.empty()))
	{
		bool locMatchFoundFlag = false;
		double locLargestEnergy = 0.0;
		for(size_t loc_i = 0; loc_i < locFCALTruthShowers.size(); ++loc_i)
		{
			locFCALTruthShower = locFCALTruthShowers[loc_i];
			for(size_t loc_j = 0; loc_j < locFCALShowers.size(); ++loc_j)
			{
				locFCALShower = locFCALShowers[loc_j];

				DVector3 locFCALShowerPosition = locFCALShower->getPosition();

				//propagate truth shower from the FCAL face to the depth of the measured shower
				double locDeltaZ = locFCALShowerPosition.Z() - locFCALTruthShower->z();
				DVector3 locTrueMomentum(locFCALTruthShower->px(), locFCALTruthShower->py(), locFCALTruthShower->pz());
				double locDeltaPathLength = locDeltaZ/cos(locTrueMomentum.Theta());
				double locTrueX = locFCALTruthShower->x() + locDeltaPathLength*sin(locTrueMomentum.Theta())*cos(locTrueMomentum.Phi());
				double locTrueY = locFCALTruthShower->y() + locDeltaPathLength*sin(locTrueMomentum.Theta())*sin(locTrueMomentum.Phi());
				DVector3 locFCALTruthShowerPosition(locTrueX, locTrueY, locFCALShowerPosition.Z());

				locMatchDistance = (locFCALShowerPosition - locFCALTruthShowerPosition).Mag();
				if(locMatchDistance > dMaximumFCALMatchDistance)
					continue;

				//keep the shower with the largest energy
				if(locFCALShower->getEnergy() > locLargestEnergy)
				{
					locBestMatchDistance = locMatchDistance;
					locMatchFoundFlag = true;
					locLargestEnergy = locFCALShower->getEnergy();
					locBestFCALTruthShowerIndex = loc_i;
					locBestFCALShowerIndex = loc_j;
				}
			}
		}

		if(!locMatchFoundFlag) //no more good matches!1
			break;

		locFCALTruthShower = locFCALTruthShowers[locBestFCALTruthShowerIndex];
		locFCALShower = locFCALShowers[locBestFCALShowerIndex];

		locFCALShowerToTruthMap[locFCALShower] = pair<const DFCALTruthShower*, double>(locFCALTruthShower, locBestMatchDistance);
		locFCALTruthToShowerMap[locFCALTruthShower] = pair<const DFCALShower*, double>(locFCALShower, locBestMatchDistance);

		locFCALShowers.erase(locFCALShowers.begin() + locBestFCALShowerIndex);
		locFCALTruthShowers.erase(locFCALTruthShowers.begin() + locBestFCALTruthShowerIndex);
	}

	locMCThrownMatching->Set_FCALShowerToTruthMap(locFCALShowerToTruthMap);
	locMCThrownMatching->Set_FCALTruthToShowerMap(locFCALTruthToShowerMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_BCALShowers(JEventLoop* locEventLoop, DMCThrownMatching* locMCThrownMatching) const
{
	vector<const DBCALTruthShower*> locBCALTruthShowers;
	const DBCALTruthShower* locBCALTruthShower;
	locEventLoop->Get(locBCALTruthShowers);

	vector<const DBCALShower*> locBCALShowers;
	const DBCALShower* locBCALShower;
	locEventLoop->Get(locBCALShowers);

	map<const DBCALShower*, pair<const DBCALTruthShower*, double> > locBCALShowerToTruthMap;
	map<const DBCALTruthShower*, pair<const DBCALShower*, double> > locBCALTruthToShowerMap;

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	//map of truth shower to dmcthrown
	map<const DBCALTruthShower*, const DMCThrown*> locBCALTruthToMCThrownMap;
	for(size_t loc_i = 0; loc_i < locBCALTruthShowers.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locMCThrowns.size(); ++loc_j)
		{
			if(locBCALTruthShowers[loc_i]->track != locMCThrowns[loc_j]->myid)
				continue;
			locBCALTruthToMCThrownMap[locBCALTruthShowers[loc_i]] = locMCThrowns[loc_j];
			break;
		}
	}

	size_t locBestBCALShowerIndex = 0, locBestBCALTruthShowerIndex = 0;
	while((!locBCALShowers.empty()) && (!locBCALTruthShowers.empty()))
	{
		bool locMatchFoundFlag = false;
		double locLargestEnergy = 0.0;
		double locBestUnitCircleArcLength = 0.0;
		for(size_t loc_i = 0; loc_i < locBCALTruthShowers.size(); ++loc_i)
		{
			locBCALTruthShower = locBCALTruthShowers[loc_i];
			for(size_t loc_j = 0; loc_j < locBCALShowers.size(); ++loc_j)
			{
				locBCALShower = locBCALShowers[loc_j];

				DVector3 locBCALShowerPosition(locBCALShower->x, locBCALShower->y, locBCALShower->z);
				DVector3 locBCALTruthShowerPosition(locBCALTruthShower->r*cos(locBCALTruthShower->phi), locBCALTruthShower->r*sin(locBCALTruthShower->phi), locBCALTruthShower->z);

				//compare theta & phi
				map<const DBCALTruthShower*, const DMCThrown*>::const_iterator locIterator = locBCALTruthToMCThrownMap.find(locBCALTruthShower);
				DVector3 locProductionVertex(0.0, 0.0, dTargetCenter); //target center
				if(locIterator != locBCALTruthToMCThrownMap.end())
					locProductionVertex = locIterator->second->position();
				double locTrueTheta = (locBCALTruthShowerPosition - locProductionVertex).Theta();
				double locReconstructedTheta = (locBCALShowerPosition - locProductionVertex).Theta();

				double locDeltaPhi = locBCALShowerPosition.Phi() - locBCALTruthShower->phi;
				while(locDeltaPhi > TMath::Pi())
					locDeltaPhi -= 2.0*TMath::Pi();
				while(locDeltaPhi < -1.0*TMath::Pi())
					locDeltaPhi += 2.0*TMath::Pi();

				double locUnitCircleArcLength = acos(sin(locTrueTheta)*sin(locReconstructedTheta)*cos(locDeltaPhi) + cos(locTrueTheta)*cos(locReconstructedTheta));
//x = r*cos(phi)*sin(theta)
//y = r*sin(phi)*sin(theta)
//z = r*cos(theta)

//unit arclength = acos(unit dot product)
//unit dot product
	//cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2) + sin(phi1)*sin(theta1)*sin(phi2)*sin(theta2) + cos(theta1)*cos(theta2)
	//sin(theta1)*sin(theta2)*(cos(phi1)*cos(phi2) + sin(phi1)*sin(phi2)) + cos(theta1)*cos(theta2)
//unit arclength = acos(sin(theta1)*sin(theta2)*cos(phi1 - phi2) + cos(theta1)*cos(theta2));

				if((locUnitCircleArcLength*180.0/TMath::Pi()) > dMaximumBCALMatchAngleDegrees)
					continue;

				//keep the shower with the largest energy
				if(locBCALShower->E > locLargestEnergy)
				{
					locBestUnitCircleArcLength = locUnitCircleArcLength;
					locMatchFoundFlag = true;
					locLargestEnergy = locBCALShower->E;
					locBestBCALTruthShowerIndex = loc_i;
					locBestBCALShowerIndex = loc_j;
				}
			}
		}

		if(!locMatchFoundFlag) //no more good matches!1
			break;

		locBCALTruthShower = locBCALTruthShowers[locBestBCALTruthShowerIndex];
		locBCALShower = locBCALShowers[locBestBCALShowerIndex];

		locBCALShowerToTruthMap[locBCALShower] = pair<const DBCALTruthShower*, double>(locBCALTruthShower, locBestUnitCircleArcLength);
		locBCALTruthToShowerMap[locBCALTruthShower] = pair<const DBCALShower*, double>(locBCALShower, locBestUnitCircleArcLength);

		locBCALShowers.erase(locBCALShowers.begin() + locBestBCALShowerIndex);
		locBCALTruthShowers.erase(locBCALTruthShowers.begin() + locBestBCALTruthShowerIndex);
	}

	locMCThrownMatching->Set_BCALShowerToTruthMap(locBCALShowerToTruthMap);
	locMCThrownMatching->Set_BCALTruthToShowerMap(locBCALTruthToShowerMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_TOFPoints(JEventLoop* locEventLoop, DMCThrownMatching* locMCThrownMatching) const
{
	vector<const DTOFTruth*> locTOFTruths;
	const DTOFTruth* locTOFTruth;
	locEventLoop->Get(locTOFTruths);

	vector<const DTOFPoint*> locTOFPoints;
	const DTOFPoint* locTOFPoint;
	locEventLoop->Get(locTOFPoints);

	map<const DTOFPoint*, pair<const DTOFTruth*, double> > locTOFPointToTruthMap;
	map<const DTOFTruth*, pair<const DTOFPoint*, double> > locTOFTruthToPointMap;

	size_t locBestTOFPointIndex = 0, locBestTOFTruthIndex = 0;
	double locBestMatchDistance, locMatchDistance;
	while((!locTOFPoints.empty()) && (!locTOFTruths.empty()))
	{
		bool locMatchFoundFlag = false;
		double locLargestEnergy = 0.0;
		for(size_t loc_i = 0; loc_i < locTOFTruths.size(); ++loc_i)
		{
			locTOFTruth = locTOFTruths[loc_i];
			for(size_t loc_j = 0; loc_j < locTOFPoints.size(); ++loc_j)
			{
				locTOFPoint = locTOFPoints[loc_j];
				DVector3 locTOFPointPosition = locTOFPoint->pos;

				//DTOFPoint and DTOFTruth reported at different z's (I think center vs. detector face): propagate truth information to the reconstructed z
				double locDeltaZ = locTOFPointPosition.Z() - locTOFTruth->z;
				DVector3 locTrueMomentum(locTOFTruth->px, locTOFTruth->py, locTOFTruth->pz);
				double locDeltaPathLength = locDeltaZ/cos(locTrueMomentum.Theta());
				double locTrueX = locTOFTruth->x + locDeltaPathLength*sin(locTrueMomentum.Theta())*cos(locTrueMomentum.Phi());
				double locTrueY = locTOFTruth->y + locDeltaPathLength*sin(locTrueMomentum.Theta())*sin(locTrueMomentum.Phi());
				DVector3 locTOFTruthPosition(locTrueX, locTrueY, locTOFPointPosition.Z());

				locMatchDistance = (locTOFTruthPosition - locTOFPointPosition).Mag();
				if(locMatchDistance > dMaximumTOFMatchDistance)
					continue;

				//keep the hit with the largest energy
				if(locTOFTruth->E > locLargestEnergy)
				{
					locBestMatchDistance = locMatchDistance;
					locMatchFoundFlag = true;
					locLargestEnergy = locTOFTruth->E;
					locBestTOFTruthIndex = loc_i;
					locBestTOFPointIndex = loc_j;
				}
			}
		}

		if(!locMatchFoundFlag) //no more good matches!
			break;

		locTOFTruth = locTOFTruths[locBestTOFTruthIndex];
		locTOFPoint = locTOFPoints[locBestTOFPointIndex];

		locTOFPointToTruthMap[locTOFPoint] = pair<const DTOFTruth*, double>(locTOFTruth, locBestMatchDistance);
		locTOFTruthToPointMap[locTOFTruth] = pair<const DTOFPoint*, double>(locTOFPoint, locBestMatchDistance);

		locTOFPoints.erase(locTOFPoints.begin() + locBestTOFPointIndex);
		locTOFTruths.erase(locTOFTruths.begin() + locBestTOFTruthIndex);
	}

	locMCThrownMatching->Set_TOFPointToTruthMap(locTOFPointToTruthMap);
	locMCThrownMatching->Set_TOFTruthToPointMap(locTOFTruthToPointMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_ChargedTrack(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DChargedTrack*>& locChargedTracks, DMCThrownMatching* locMCThrownMatching) const
{
	//assumes Find_GenReconMatches_ChargedHypo has been called first!
	map<const DChargedTrackHypothesis*, pair<const DMCThrown*, double> > locChargedHypoToThrownMap;
	map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> > locThrownToChargedHypoMap;
	locMCThrownMatching->Get_ChargedHypoToThrownMap(locChargedHypoToThrownMap);
	locMCThrownMatching->Get_ThrownToChargedHypoMap(locThrownToChargedHypoMap);

	map<const DChargedTrack*, pair<const DMCThrown*, double> > locChargedToThrownMap;
	map<const DMCThrown*, pair<const DChargedTrack*, double> > locThrownToChargedMap;

	map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> >::iterator locIterator;
	for(locIterator = locThrownToChargedHypoMap.begin(); locIterator != locThrownToChargedHypoMap.end(); ++locIterator)
	{
		for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j)
		{
			if(locChargedTracks[loc_j]->Get_BestFOM()->candidateid == (locIterator->second.first)[0]->candidateid) //match
			{
				locChargedToThrownMap[locChargedTracks[loc_j]] = pair<const DMCThrown*, double>(locIterator->first, locIterator->second.second);
				locThrownToChargedMap[locIterator->first] = pair<const DChargedTrack*, double>(locChargedTracks[loc_j], locIterator->second.second);
			}
		}
	}

	locMCThrownMatching->Set_ChargedToThrownMap(locChargedToThrownMap);
	locMCThrownMatching->Set_ThrownToChargedMap(locThrownToChargedMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_ChargedHypo(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DChargedTrackHypothesis*>& locInputChargedTrackHypothesisVector, DMCThrownMatching* locMCThrownMatching) const
{
	map<const DChargedTrackHypothesis*, pair<const DMCThrown*, double> > locChargedToThrownMap;
	map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> > locThrownToChargedMap;

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
				locMatchFOM = Calc_MatchFOM(locMCThrown->momentum(), locChargedTrackHypothesis->momentum(), locChargedTrackHypothesis->errorMatrix());

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

		locChargedToThrownMap[locChargedTrackHypothesis] = pair<const DMCThrown*, double>(locMCThrown, locBestMatchFOM);
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
				locChargedToThrownMap[locChargedTrackHypothesisVector[loc_i]] = pair<const DMCThrown*, double>(locMCThrown, locBestMatchFOM);
				locMatchedChargedHypos.push_back(locChargedTrackHypothesisVector[loc_i]);
				locChargedTrackHypothesisVector.erase(locChargedTrackHypothesisVector.begin() + loc_i);
			}
		}
		locThrownToChargedMap[locMCThrown] = pair<deque<const DChargedTrackHypothesis*>, double>(locMatchedChargedHypos, locBestMatchFOM);
	}

	locMCThrownMatching->Set_ChargedHypoToThrownMap(locChargedToThrownMap);
	locMCThrownMatching->Set_ThrownToChargedHypoMap(locThrownToChargedMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_NeutralParticle(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DNeutralParticle*>& locNeutralParticles, DMCThrownMatching* locMCThrownMatching) const
{
	//assumes Find_GenReconMatches_NeutralHypo has been called first!
	map<const DNeutralParticleHypothesis*, pair<const DMCThrown*, double> > locNeutralHypoToThrownMap;
	map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> > locThrownToNeutralHypoMap;
	locMCThrownMatching->Get_NeutralHypoToThrownMap(locNeutralHypoToThrownMap);
	locMCThrownMatching->Get_ThrownToNeutralHypoMap(locThrownToNeutralHypoMap);

	vector<const DNeutralShower*> locNeutralShowerVector_Matched;
	vector<const DNeutralShower*> locNeutralShowerVector_Check;

	map<const DNeutralParticle*, pair<const DMCThrown*, double> > locNeutralToThrownMap;
	map<const DMCThrown*, pair<const DNeutralParticle*, double> > locThrownToNeutralMap;

	map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> >::iterator locIterator;
	for(locIterator = locThrownToNeutralHypoMap.begin(); locIterator != locThrownToNeutralHypoMap.end(); ++locIterator)
	{
		(locIterator->second.first)[0]->GetT(locNeutralShowerVector_Matched);
		for(size_t loc_j = 0; loc_j < locNeutralParticles.size(); ++loc_j)
		{
			locNeutralParticles[loc_j]->GetT(locNeutralShowerVector_Check);
			if(locNeutralShowerVector_Check[0] == locNeutralShowerVector_Matched[0])
			{
				locNeutralToThrownMap[locNeutralParticles[loc_j]] = pair<const DMCThrown*, double>(locIterator->first, locIterator->second.second);
				locThrownToNeutralMap[locIterator->first] = pair<const DNeutralParticle*, double>(locNeutralParticles[loc_j], locIterator->second.second);
			}
		}
	}

	locMCThrownMatching->Set_NeutralToThrownMap(locNeutralToThrownMap);
	locMCThrownMatching->Set_ThrownToNeutralMap(locThrownToNeutralMap);
}

void DMCThrownMatching_factory::Find_GenReconMatches_NeutralHypo(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DNeutralParticleHypothesis*>& locInputNeutralParticleHypothesisVector, DMCThrownMatching* locMCThrownMatching) const
{
	map<const DNeutralParticleHypothesis*, pair<const DMCThrown*, double> > locNeutralToThrownMap;
	map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> > locThrownToNeutralMap;

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
				locMatchFOM = Calc_MatchFOM(locMCThrown->momentum(), locNeutralParticleHypothesis->momentum(), locNeutralParticleHypothesis->errorMatrix());
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

		locNeutralToThrownMap[locNeutralParticleHypothesis] = pair<const DMCThrown*, double>(locMCThrown, locBestMatchFOM);
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
				locNeutralToThrownMap[locNeutralParticleHypothesisVector[loc_i]] = pair<const DMCThrown*, double>(locMCThrown, locBestMatchFOM);
				locMatchedNeutralHypos.push_back(locNeutralParticleHypothesisVector[loc_i]);
				locNeutralParticleHypothesisVector.erase(locNeutralParticleHypothesisVector.begin() + loc_i);
			}
		}
		locThrownToNeutralMap[locMCThrown] = pair<deque<const DNeutralParticleHypothesis*>, double>(locMatchedNeutralHypos, locBestMatchFOM);
	}

	locMCThrownMatching->Set_NeutralHypoToThrownMap(locNeutralToThrownMap);
	locMCThrownMatching->Set_ThrownToNeutralHypoMap(locThrownToNeutralMap);
}

double DMCThrownMatching_factory::Calc_MatchFOM(const DVector3& locMomentum_Thrown, const DVector3& locMomentum_Detected, const DMatrixDSym& locInputCovarianceMatrix) const
{
	//Detailed Calculation:
	DVector3 locDeltaP3 = locMomentum_Detected - locMomentum_Thrown;
	double locTotalDeltaSq = locDeltaP3.Mag2();
	double locTotalDelta = sqrt(locTotalDeltaSq);

	//dx = detected_x - thrown_x
	//total_delta = sqrt(dx^2 + dy^2 + dz^2)
	//partial = 1/2*(1/total_delta)*2*dx*1 = dx/total_delta
	DMatrixDSym locCovarianceMatrix = locInputCovarianceMatrix;
	DMatrix locJacobian(1, 7);
	locJacobian.Zero();
	locJacobian(0, 0) = locDeltaP3.Px()/locTotalDelta;
	locJacobian(0, 1) = locDeltaP3.Py()/locTotalDelta;
	locJacobian(0, 2) = locDeltaP3.Pz()/locTotalDelta;
	double locVariance = (locCovarianceMatrix.Similarity(locJacobian))(0, 0);

	double locChiSq = locTotalDeltaSq/locVariance;
	double locFOM = TMath::Prob(locChiSq, 3);
	if(dDebugLevel > 10)
		cout << "delta pxyz, var pxyz, total delsq, total var, total FOM = " << locDeltaP3.Px() << ", " << locDeltaP3.Py() << ", " << locDeltaP3.Pz() << ", " << locInputCovarianceMatrix(0, 0) << ", " << locInputCovarianceMatrix(1, 1) << ", " << locInputCovarianceMatrix(2, 2) << ", " << locTotalDeltaSq << ", " << locVariance << ", " << locFOM << endl;

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



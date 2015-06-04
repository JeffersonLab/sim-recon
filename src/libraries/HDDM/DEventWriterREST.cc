#include "DEventWriterREST.h"

//file-scope so shared amongst threads; only accessed below via locks
size_t gRESTNumOutputThreads = 0;
map<string, pair<ofstream*, hddm_r::ostream*> >* gRESTOutputFilePointers = NULL;

DEventWriterREST::DEventWriterREST(JEventLoop* locEventLoop, string locOutputFileBaseName) : dOutputFileBaseName(locOutputFileBaseName)
{
	japp->WriteLock("RESTWriter");
	{
		++gRESTNumOutputThreads;
		if(gRESTOutputFilePointers == NULL)
			gRESTOutputFilePointers = new map<string, pair<ofstream*, hddm_r::ostream*> >();
	}
	japp->Unlock("RESTWriter");
	
	HDDM_USE_COMPRESSION = true;
	string locCompressionString = "Turn on/off compression of the output HDDM stream. Set to \"0\" to turn off (it's on by default)";
	gPARMS->SetDefaultParameter("HDDM:USE_COMPRESSION", HDDM_USE_COMPRESSION, locCompressionString);

	HDDM_USE_INTEGRITY_CHECKS = true;
	string locIntegrityString = "Turn on/off automatic integrity checking on the output HDDM stream. Set to \"0\" to turn off (it's on by default)";
	gPARMS->SetDefaultParameter("HDDM:USE_INTEGRITY_CHECKS", HDDM_USE_INTEGRITY_CHECKS, locIntegrityString);
}

bool DEventWriterREST::Write_RESTEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const
{
	std::vector<const DMCReaction*> reactions;
	locEventLoop->Get(reactions);

	std::vector<const DRFTime*> rftimes;
	locEventLoop->Get(rftimes);

	std::vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	std::vector<const DFCALShower*> fcalshowers;
	locEventLoop->Get(fcalshowers);

	std::vector<const DBCALShower*> bcalshowers;
	locEventLoop->Get(bcalshowers);

	std::vector<const DTOFPoint*> tofpoints;
	locEventLoop->Get(tofpoints);

	std::vector<const DSCHit*> starthits;
	locEventLoop->Get(starthits);

	std::vector<const DTrackTimeBased*> tracks;
	locEventLoop->Get(tracks);

	std::vector<const DDetectorMatches*> locDetectorMatches;
	locEventLoop->Get(locDetectorMatches);

	//Check to see if there are any objects to write out.  If so, don't write out an empty event
	bool locOutputDataPresentFlag = false;
	if((!reactions.empty()) || (!rftimes.empty()) || (!locBeamPhotons.empty()) || (!tracks.empty()))
		locOutputDataPresentFlag = true;
	else if((!fcalshowers.empty()) || (!bcalshowers.empty()) || (!tofpoints.empty()) || (!starthits.empty()))
		locOutputDataPresentFlag = true;
	//don't need to check detector matches: no matches if none of the above objects
	if(!locOutputDataPresentFlag)
		return true; //had correct response to data

	string locOutputFileName = Get_OutputFileName(locOutputFileNameSubString);

	hddm_r::HDDM locRecord;
	hddm_r::ReconstructedPhysicsEventList res = locRecord.addReconstructedPhysicsEvents(1);

	// load the run and event numbers
	JEvent& event = locEventLoop->GetJEvent();
	res().setRunNo(event.GetRunNumber());
	res().setEventNo(event.GetEventNumber());

	// push any DMCReaction objects to the output record
	for (size_t i=0; i < reactions.size(); i++)
	{
		hddm_r::ReactionList rea = res().addReactions(1);
		rea().setType(reactions[i]->type);
		rea().setWeight(reactions[i]->weight);
		rea().setEbeam(reactions[i]->beam.energy());
		rea().setTargetType(reactions[i]->target.PID());

		if(i != 0)
			break;

		std::vector<const DMCThrown*> throwns;
		locEventLoop->Get(throwns);
		hddm_r::VertexList ver = rea().getVertices();
		DLorentzVector locPreviousX4(-9.9E9, -9.9E9, -9.9E9, -9.9E9);
		for(size_t it=0; it < throwns.size(); ++it)
		{
			DLorentzVector locThrownX4(throwns[it]->position(), throwns[it]->time());
			if((locThrownX4.T() != locPreviousX4.T()) || (locThrownX4.Vect() != locPreviousX4.Vect()))
			{
				//new vertex
				ver = rea().addVertices(1);
				hddm_r::OriginList ori = ver().addOrigins(1);
				ori().setT(locThrownX4.T());
				ori().setVx(locThrownX4.X());
				ori().setVy(locThrownX4.Y());
				ori().setVz(locThrownX4.Z());
				locPreviousX4 = locThrownX4;
			}

			hddm_r::ProductList pro = ver().addProducts(1);
			pro().setId(throwns[it]->myid);
			pro().setParentId(throwns[it]->parentid);
			int pdgtype = throwns[it]->pdgtype;
			if (pdgtype == 0)
				pdgtype = PDGtype((Particle_t)throwns[it]->type);
			pro().setPdgtype(pdgtype);
			hddm_r::MomentumList mom = pro().addMomenta(1);
			mom().setE(throwns[it]->energy());
			mom().setPx(throwns[it]->px());
			mom().setPy(throwns[it]->py());
			mom().setPz(throwns[it]->pz());
		}
	}

	// push any DRFTime objects to the output record
	for (size_t i=0; i < rftimes.size(); i++)
	{
		hddm_r::RFtimeList rf = res().addRFtimes(1);
		rf().setTsync(rftimes[i]->dTime);
	}

	// push any DBeamPhoton objects to the output record
	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
	{
		if(locBeamPhotons[loc_i]->t0_detector() == SYS_TAGM)
		{
			hddm_r::TagmBeamPhotonList locTagmBeamPhotonList = res().addTagmBeamPhotons(1);
			locTagmBeamPhotonList().setT(locBeamPhotons[loc_i]->time());
			locTagmBeamPhotonList().setE(locBeamPhotons[loc_i]->energy());
		}
		else if(locBeamPhotons[loc_i]->t0_detector() == SYS_TAGH)
		{
			hddm_r::TaghBeamPhotonList locTaghBeamPhotonList = res().addTaghBeamPhotons(1);
			locTaghBeamPhotonList().setT(locBeamPhotons[loc_i]->time());
			locTaghBeamPhotonList().setE(locBeamPhotons[loc_i]->energy());
		}
	}

	// push any DFCALShower objects to the output record
	for (size_t i=0; i < fcalshowers.size(); i++)
	{
		hddm_r::FcalShowerList fcal = res().addFcalShowers(1);
		DVector3 pos = fcalshowers[i]->getPosition();
		DVector3 poserr = fcalshowers[i]->getPositionError();
		fcal().setX(pos(0));
		fcal().setY(pos(1));
		fcal().setZ(pos(2));
		fcal().setT(fcalshowers[i]->getTime());
		fcal().setE(fcalshowers[i]->getEnergy());
		fcal().setXerr(poserr(0));
		fcal().setYerr(poserr(1));
		fcal().setZerr(poserr(2));
		fcal().setTerr(0);
		fcal().setEerr(0);
		fcal().setXycorr(0);
		fcal().setXzcorr(0);
		fcal().setYzcorr(0);
		fcal().setEzcorr(0);
		fcal().setTzcorr(0);
	}

	// push any DBCALShower objects to the output record
	for (size_t i=0; i < bcalshowers.size(); i++)
	{
		hddm_r::BcalShowerList bcal = res().addBcalShowers(1);
		DVector3 pos(bcalshowers[i]->x,bcalshowers[i]->y,bcalshowers[i]->z);
		bcal().setX(bcalshowers[i]->x);
		bcal().setY(bcalshowers[i]->y);
		bcal().setZ(bcalshowers[i]->z);
		bcal().setT(bcalshowers[i]->t);
		bcal().setE(bcalshowers[i]->E);
		bcal().setXerr(bcalshowers[i]->xErr);
		bcal().setYerr(bcalshowers[i]->yErr);
		bcal().setZerr(bcalshowers[i]->zErr);
		bcal().setTerr(bcalshowers[i]->tErr);
		bcal().setEerr(0);
		bcal().setXycorr(0);
		bcal().setXzcorr(0);
		bcal().setYzcorr(0);
		bcal().setEzcorr(0);
		bcal().setTzcorr(0);

		//N_cell
                hddm_r::BcalClusterList bcalcluster = bcal().addBcalClusters(1);
                bcalcluster().setNcell(bcalshowers[i]->N_cell);
	}

	// push any DTOFPoint objects to the output record
	for (size_t i=0; i < tofpoints.size(); i++)
	{
		hddm_r::TofPointList tof = res().addTofPoints(1);
		tof().setX(tofpoints[i]->pos(0));
		tof().setY(tofpoints[i]->pos(1));
		tof().setZ(tofpoints[i]->pos(2));
		tof().setT(tofpoints[i]->t);
		tof().setDE(tofpoints[i]->dE);

		//Status //Assume compiler optimizes multiplication
		hddm_r::TofStatusList tofstatus = tof().addTofStatuses(1);
		int locStatus = tofpoints[i]->dHorizontalBar + 45*tofpoints[i]->dVerticalBar;
		locStatus += 45*4*tofpoints[i]->dHorizontalBarStatus + 45*4*4*tofpoints[i]->dVerticalBarStatus;
		tofstatus().setStatus(locStatus);
	}

	// push any DSCHit objects to the output record
	for (size_t i=0; i < starthits.size(); i++)
	{
		hddm_r::StartHitList hit = res().addStartHits(1);
		hit().setSector(starthits[i]->sector);
		hit().setT(starthits[i]->t);
		hit().setDE(starthits[i]->dE);
	}

	// push any DTrackTimeBased objects to the output record
	for (size_t i=0; i < tracks.size(); ++i)
	{
		hddm_r::ChargedTrackList tra = res().addChargedTracks(1);
		tra().setCandidateId(tracks[i]->candidateid);
		tra().setPtype(tracks[i]->PID());

		hddm_r::TrackFitList fit = tra().addTrackFits(1);
		fit().setNdof(tracks[i]->Ndof);
		fit().setChisq(tracks[i]->chisq);
		fit().setX0(tracks[i]->x());
		fit().setY0(tracks[i]->y());
		fit().setZ0(tracks[i]->z());
		fit().setPx(tracks[i]->px());
		fit().setPy(tracks[i]->py());
		fit().setPz(tracks[i]->pz());
		fit().setT0(tracks[i]->t0());
		fit().setT0err(tracks[i]->t0_err());
		fit().setT0det(tracks[i]->t0_detector());

		DMatrixDSym errors = tracks[i]->TrackingErrorMatrix();
		fit().setE11(errors(0,0));
		fit().setE12(errors(0,1));
		fit().setE13(errors(0,2));
		fit().setE14(errors(0,3));
		fit().setE15(errors(0,4));
		fit().setE22(errors(1,1));
		fit().setE23(errors(1,2));
		fit().setE24(errors(1,3));
		fit().setE25(errors(1,4));
		fit().setE33(errors(2,2));
		fit().setE34(errors(2,3));
		fit().setE35(errors(2,4));
		fit().setE44(errors(3,3));
		fit().setE45(errors(3,4));
		fit().setE55(errors(4,4));

		hddm_r::HitlayersList locHitLayers = tra().addHitlayerses(1);
		locHitLayers().setCDCrings(tracks[i]->dCDCRings);
		locHitLayers().setFDCplanes(tracks[i]->dFDCPlanes);

		hddm_r::McmatchList locMCMatches = tra().addMcmatchs(1);
		locMCMatches().setIthrown(tracks[i]->dMCThrownMatchMyID);
		locMCMatches().setNumhitsmatch(tracks[i]->dNumHitsMatchedToThrown);

		if (tracks[i]->dNumHitsUsedFordEdx_FDC + tracks[i]->dNumHitsUsedFordEdx_CDC > 0)
		{
			hddm_r::DEdxDCList elo = tra().addDEdxDCs(1);
			elo().setNsampleFDC(tracks[i]->dNumHitsUsedFordEdx_FDC);
			elo().setNsampleCDC(tracks[i]->dNumHitsUsedFordEdx_CDC);
			elo().setDxFDC(tracks[i]->ddx_FDC);
			elo().setDxCDC(tracks[i]->ddx_CDC);
			elo().setDEdxFDC(tracks[i]->ddEdx_FDC);
			elo().setDEdxCDC(tracks[i]->ddEdx_CDC);
		}
	}

	// push any DMCTrigger objects to the output record
	std::vector<const DMCTrigger*> triggers;
	locEventLoop->Get(triggers);
	for (size_t i=0; i < triggers.size(); i++)
	{
		hddm_r::TriggerList trigger = res().addTriggers(1);
		trigger().setL1a(triggers[i]->L1a_fired);
		trigger().setL1b(triggers[i]->L1b_fired);
		trigger().setL1c(triggers[i]->L1c_fired);
	}

	// push any DDetectorMatches objects to the output record
	for(size_t loc_i = 0; loc_i < locDetectorMatches.size(); ++loc_i)
	{
		hddm_r::DetectorMatchesList matches = res().addDetectorMatcheses(1);
		for(size_t loc_j = 0; loc_j < tracks.size(); ++loc_j)
		{
			vector<DBCALShowerMatchParams> locBCALShowerMatchParamsVector;
			locDetectorMatches[loc_i]->Get_BCALMatchParams(tracks[loc_j], locBCALShowerMatchParamsVector);
			for(size_t loc_k = 0; loc_k < locBCALShowerMatchParamsVector.size(); ++loc_k)
			{
				hddm_r::BcalMatchParams_v2List bcalList = matches().addBcalMatchParams_v2s(1);
				bcalList().setTrack(loc_j);

				const DBCALShower* locBCALShower = locBCALShowerMatchParamsVector[loc_k].dBCALShower;
				size_t locBCALindex = 0;
				for(; locBCALindex < bcalshowers.size(); ++locBCALindex)
				{
					if(bcalshowers[locBCALindex] == locBCALShower)
						break;
				}
				bcalList().setShower(locBCALindex);

				bcalList().setDeltaphi(locBCALShowerMatchParamsVector[loc_k].dDeltaPhiToShower);
				bcalList().setDeltaz(locBCALShowerMatchParamsVector[loc_k].dDeltaZToShower);
				bcalList().setDx(locBCALShowerMatchParamsVector[loc_k].dx);
				bcalList().setPathlength(locBCALShowerMatchParamsVector[loc_k].dPathLength);
				bcalList().setTflight(locBCALShowerMatchParamsVector[loc_k].dFlightTime);
				bcalList().setTflightvar(locBCALShowerMatchParamsVector[loc_k].dFlightTimeVariance);
			}

			vector<DFCALShowerMatchParams> locFCALShowerMatchParamsVector;
			locDetectorMatches[loc_i]->Get_FCALMatchParams(tracks[loc_j], locFCALShowerMatchParamsVector);
			for (size_t loc_k = 0; loc_k < locFCALShowerMatchParamsVector.size(); ++loc_k)
			{
				hddm_r::FcalMatchParamsList fcalList = matches().addFcalMatchParamses(1);
				fcalList().setTrack(loc_j);

				const DFCALShower* locFCALShower = locFCALShowerMatchParamsVector[loc_k].dFCALShower;
				size_t locFCALindex = 0;
				for(; locFCALindex < fcalshowers.size(); ++locFCALindex)
				{
					if(fcalshowers[locFCALindex] == locFCALShower)
						break;
				}
				fcalList().setShower(locFCALindex);

				fcalList().setDoca(locFCALShowerMatchParamsVector[loc_k].dDOCAToShower);
				fcalList().setDx(locFCALShowerMatchParamsVector[loc_k].dx);
				fcalList().setPathlength(locFCALShowerMatchParamsVector[loc_k].dPathLength);
				fcalList().setTflight(locFCALShowerMatchParamsVector[loc_k].dFlightTime);
				fcalList().setTflightvar(locFCALShowerMatchParamsVector[loc_k].dFlightTimeVariance);
			}

			vector<DTOFHitMatchParams> locTOFHitMatchParamsVector;
			locDetectorMatches[loc_i]->Get_TOFMatchParams(tracks[loc_j], locTOFHitMatchParamsVector);
			for(size_t loc_k = 0; loc_k < locTOFHitMatchParamsVector.size(); ++loc_k)
			{
				hddm_r::TofMatchParams_v2List tofList = matches().addTofMatchParams_v2s(1);
				tofList().setTrack(loc_j);

				size_t locTOFindex = 0;
				for(; locTOFindex < tofpoints.size(); ++locTOFindex)
				{
					if(tofpoints[locTOFindex] == locTOFHitMatchParamsVector[loc_k].dTOFPoint)
						break;
				}
				tofList().setHit(locTOFindex);

				tofList().setThit(locTOFHitMatchParamsVector[loc_k].dHitTime);
				tofList().setThitvar(locTOFHitMatchParamsVector[loc_k].dHitTimeVariance);
				tofList().setEhit(locTOFHitMatchParamsVector[loc_k].dHitEnergy);

				tofList().setDEdx(locTOFHitMatchParamsVector[loc_k].dEdx);
				tofList().setPathlength(locTOFHitMatchParamsVector[loc_k].dPathLength);
				tofList().setTflight(locTOFHitMatchParamsVector[loc_k].dFlightTime);
				tofList().setTflightvar(locTOFHitMatchParamsVector[loc_k].dFlightTimeVariance);

				tofList().setDeltax(locTOFHitMatchParamsVector[loc_k].dDeltaXToHit);
				tofList().setDeltay(locTOFHitMatchParamsVector[loc_k].dDeltaYToHit);
			}

			vector<DSCHitMatchParams> locSCHitMatchParamsVector;
			locDetectorMatches[loc_i]->Get_SCMatchParams(tracks[loc_j], locSCHitMatchParamsVector);
			for(size_t loc_k = 0; loc_k < locSCHitMatchParamsVector.size(); ++loc_k)
			{
				hddm_r::ScMatchParamsList scList = matches().addScMatchParamses(1);
				scList().setTrack(loc_j);

				size_t locSCindex = 0;
				for(; locSCindex < starthits.size(); ++locSCindex)
				{
					if(starthits[locSCindex] == locSCHitMatchParamsVector[loc_k].dSCHit)
						break;
				}
				scList().setHit(locSCindex);

				scList().setDEdx(locSCHitMatchParamsVector[loc_k].dEdx);
				scList().setDeltaphi(locSCHitMatchParamsVector[loc_k].dDeltaPhiToHit);
				scList().setEhit(locSCHitMatchParamsVector[loc_k].dHitEnergy);
				scList().setPathlength(locSCHitMatchParamsVector[loc_k].dPathLength);
				scList().setTflight(locSCHitMatchParamsVector[loc_k].dFlightTime);
				scList().setTflightvar(locSCHitMatchParamsVector[loc_k].dFlightTimeVariance);
				scList().setThit(locSCHitMatchParamsVector[loc_k].dHitTime);
				scList().setThitvar(locSCHitMatchParamsVector[loc_k].dHitTimeVariance);
			}

			double locFlightTimePCorrelation = 0.0;
			if(locDetectorMatches[loc_i]->Get_FlightTimePCorrelation(tracks[loc_j], SYS_BCAL, locFlightTimePCorrelation))
			{
				hddm_r::TflightPCorrelationList correlationList = matches().addTflightPCorrelations(1);
				correlationList().setTrack(loc_j);
				correlationList().setSystem(SYS_BCAL);
				correlationList().setCorrelation(locFlightTimePCorrelation);
			}
			if(locDetectorMatches[loc_i]->Get_FlightTimePCorrelation(tracks[loc_j], SYS_FCAL, locFlightTimePCorrelation))
			{
				hddm_r::TflightPCorrelationList correlationList = matches().addTflightPCorrelations(1);
				correlationList().setTrack(loc_j);
				correlationList().setSystem(SYS_FCAL);
				correlationList().setCorrelation(locFlightTimePCorrelation);
			}
			if(locDetectorMatches[loc_i]->Get_FlightTimePCorrelation(tracks[loc_j], SYS_TOF, locFlightTimePCorrelation))
			{
				hddm_r::TflightPCorrelationList correlationList = matches().addTflightPCorrelations(1);
				correlationList().setTrack(loc_j);
				correlationList().setSystem(SYS_TOF);
				correlationList().setCorrelation(locFlightTimePCorrelation);
			}
			if(locDetectorMatches[loc_i]->Get_FlightTimePCorrelation(tracks[loc_j], SYS_START, locFlightTimePCorrelation))
			{
				hddm_r::TflightPCorrelationList correlationList = matches().addTflightPCorrelations(1);
				correlationList().setTrack(loc_j);
				correlationList().setSystem(SYS_START);
				correlationList().setCorrelation(locFlightTimePCorrelation);
			}
		}

		for(size_t loc_j = 0; loc_j < bcalshowers.size(); ++loc_j)
		{
			double locDeltaPhi = 0.0, locDeltaZ = 0.0;
			if(!locDetectorMatches[loc_i]->Get_DistanceToNearestTrack(bcalshowers[loc_j], locDeltaPhi, locDeltaZ))
				continue;

			hddm_r::BcalDOCAtoTrack_v2List bcalDocaList = matches().addBcalDOCAtoTrack_v2s(1);
			bcalDocaList().setShower(loc_j);
			bcalDocaList().setDeltaphi(locDeltaPhi);
			bcalDocaList().setDeltaz(locDeltaZ);
		}

		for(size_t loc_j = 0; loc_j < fcalshowers.size(); ++loc_j)
		{
			double locDistance = 0.0;
			if(!locDetectorMatches[loc_i]->Get_DistanceToNearestTrack(fcalshowers[loc_j], locDistance))
				continue;

			hddm_r::FcalDOCAtoTrackList fcalDocaList = matches().addFcalDOCAtoTracks(1);
			fcalDocaList().setShower(loc_j);
			fcalDocaList().setDoca(locDistance);
		}
	}

	// write the resulting record to the output stream
	bool locWriteStatus = Write_RESTEvent(locOutputFileName, locRecord);
	locRecord.clear();
	return locWriteStatus;
}

string DEventWriterREST::Get_OutputFileName(string locOutputFileNameSubString) const
{
	string locOutputFileName = dOutputFileBaseName;
	if (locOutputFileNameSubString != "")
		locOutputFileName += string("_") + locOutputFileNameSubString;
	return (locOutputFileName + string(".hddm"));
}

bool DEventWriterREST::Write_RESTEvent(string locOutputFileName, hddm_r::HDDM& locRecord) const
{
	japp->WriteLock("RESTWriter");
	{
		//check to see if the REST file is open
		if(gRESTOutputFilePointers->find(locOutputFileName) != gRESTOutputFilePointers->end())
		{
			//open: get pointer, write event
			hddm_r::ostream* locOutputRESTFileStream = (*gRESTOutputFilePointers)[locOutputFileName].second;
			*(locOutputRESTFileStream) << locRecord;
			japp->Unlock("RESTWriter");
			return true;
		}

		//not open: open it
		pair<ofstream*, hddm_r::ostream*> locRESTFilePointers(NULL, NULL);
		locRESTFilePointers.first = new ofstream(locOutputFileName.c_str());
		if(!locRESTFilePointers.first->is_open())
		{
			//failed to open
			delete locRESTFilePointers.first;
			japp->Unlock("RESTWriter");
			return false;
		}
		locRESTFilePointers.second = new hddm_r::ostream(*locRESTFilePointers.first);

		// enable on-the-fly bzip2 compression on output stream
		if(HDDM_USE_COMPRESSION)
		{
			jout << " Enabling bz2 compression of output HDDM file stream" << std::endl;
			locRESTFilePointers.second->setCompression(hddm_r::k_bz2_compression);
		}
		else
			jout << " HDDM compression disabled" << std::endl;

		// enable a CRC data integrity check at the end of each event record
		if(HDDM_USE_INTEGRITY_CHECKS)
		{
			jout << " Enabling CRC data integrity check in output HDDM file stream" << std::endl;
			locRESTFilePointers.second->setIntegrityChecks(hddm_r::k_crc32_integrity);
		}
		else
			jout << " HDDM integrity checks disabled" << std::endl;

		// write a comment record at the head of the file
		hddm_r::HDDM locCommentRecord;
		hddm_r::ReconstructedPhysicsEventList res = locCommentRecord.addReconstructedPhysicsEvents(1);
		hddm_r::CommentList comment = res().addComments();
		comment().setText("this is a REST event stream, yadda yadda");
		*(locRESTFilePointers.second) << locCommentRecord;
		locCommentRecord.clear();

		//write the event
		*(locRESTFilePointers.second) << locRecord;

		//store the stream pointers
		(*gRESTOutputFilePointers)[locOutputFileName] = locRESTFilePointers;
	}
	japp->Unlock("RESTWriter");

	return true;
}

DEventWriterREST::~DEventWriterREST(void)
{
	japp->WriteLock("RESTWriter");
	{
		--gRESTNumOutputThreads;
		if(gRESTNumOutputThreads > 0)
		{
			japp->Unlock("RESTWriter");
			return; //not the last thread writing to REST files
		}

		if(gRESTOutputFilePointers == NULL)
		{
			japp->Unlock("RESTWriter");
			return; //not the last thread writing to REST files
		}
		
		//last thread writing to REST files: close all files and free all memory
		map<string, pair<ofstream*, hddm_r::ostream*> >::iterator locIterator;
		for(locIterator = gRESTOutputFilePointers->begin(); locIterator != gRESTOutputFilePointers->end(); ++locIterator)
		{
			string locOutputFileName = locIterator->first;
			if (locIterator->second.second != NULL)
				delete locIterator->second.second;
			if (locIterator->second.first != NULL)
				delete locIterator->second.first;
			std::cout << "Closed REST file " << locOutputFileName << std::endl;
		}
		delete gRESTOutputFilePointers;
		gRESTOutputFilePointers = NULL;
	}
	japp->Unlock("RESTWriter");
}


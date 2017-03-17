#include "DEventWriterHDDM.h"

#include <DANA/DApplication.h>
#include <JANA/JCalibration.h>

int& DEventWriterHDDM::Get_NumEventWriterThreads(void) const
{
	// must be read/used entirely in "HDDMWriter" lock
	static int locNumEventWriterThreads = 0;
	return locNumEventWriterThreads;
}

map<string, pair<ofstream*, hddm_s::ostream*> >& DEventWriterHDDM::Get_HDDMOutputFilePointers(void) const
{
	// must be read/used entirely in "HDDMWriter" lock
	// cannot do individual file locks, because the map itself can be modified
	static map<string, pair<ofstream*, hddm_s::ostream*> > locHDDMOutputFilePointers;
	return locHDDMOutputFilePointers;
}

DEventWriterHDDM::DEventWriterHDDM(JEventLoop* locEventLoop, string locOutputFileBaseName) : dOutputFileBaseName(locOutputFileBaseName)
{
	japp->WriteLock("HDDMWriter");
	{
		++Get_NumEventWriterThreads();
	}
	japp->Unlock("HDDMWriter");
	
	HDDM_USE_COMPRESSION = true;
	string locCompressionString = "Turn on/off compression of the output HDDM stream. Set to \"0\" to turn off (it's on by default)";
	gPARMS->SetDefaultParameter("HDDM:USE_COMPRESSION", HDDM_USE_COMPRESSION, locCompressionString);

	HDDM_USE_INTEGRITY_CHECKS = true;
	string locIntegrityString = "Turn on/off automatic integrity checking on the output HDDM stream. Set to \"0\" to turn off (it's on by default)";
	gPARMS->SetDefaultParameter("HDDM:USE_INTEGRITY_CHECKS", HDDM_USE_INTEGRITY_CHECKS, locIntegrityString);

    HDDM_DATA_VERSION_STRING = "";
    gPARMS->SetDefaultParameter("HDDM:DATAVERSIONSTRING", HDDM_DATA_VERSION_STRING, "");

    CCDB_CONTEXT_STRING = "";
    // if we can get the calibration context from the DANA interface, then save this as well
    DApplication *dapp = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
    if (dapp) {
        JEvent& event = locEventLoop->GetJEvent();
        JCalibration *jcalib = dapp->GetJCalibration(event.GetRunNumber());
        if (jcalib) {
            CCDB_CONTEXT_STRING = jcalib->GetContext();
        }
    }

}

bool DEventWriterHDDM::Write_HDDMEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const
{
  
	vector<const DCDCHit*> CDCHits;
	vector<const DTOFHit*> TOFHits;
	vector<const DFCALHit*> FCALHits;
	vector<const DSCHit*> SCHits;
	vector<const DBCALDigiHit*> BCALDigiHits;
	vector<const DBCALTDCDigiHit*> BCALTDCDigiHits;
	vector<const DPSHit*> PSHits;
	vector<const DPSCHit*> PSCHits;
	vector<const DFDCHit*> FDCHits;
	vector<const DTAGHHit*> TAGHHits;
	vector<const DTAGMHit*> TAGMHits;
	vector<const DTPOLHit*> TPOLHits;

	locEventLoop->Get(CDCHits);
	locEventLoop->Get(TOFHits);
	locEventLoop->Get(FCALHits);
	locEventLoop->Get(BCALDigiHits);
	locEventLoop->Get(BCALTDCDigiHits);
	locEventLoop->Get(SCHits);
	locEventLoop->Get(PSHits);
	locEventLoop->Get(PSCHits);
	locEventLoop->Get(FDCHits);
	locEventLoop->Get(TAGHHits);
	locEventLoop->Get(TAGMHits);
	locEventLoop->Get(TPOLHits);


	if(CDCHits.size()== uint(0) && TOFHits.size()==uint(0) && FCALHits.size()==uint(0) && BCALDigiHits.size()==uint(0) && BCALTDCDigiHits.size()==uint(0) && SCHits.size()==uint(0) && PSHits.size()==uint(0) && PSCHits.size()==uint(0) && FDCHits.size()==uint(0) && TAGHHits.size()==uint(0) && TAGMHits.size()==uint(0) && TPOLHits.size()==uint(0))
	{
		return false;
	}

	//create an HDDM record to store the Events' hits and set the Event/Run Number
	hddm_s::HDDM* record = new hddm_s::HDDM;
	record->addPhysicsEvents();
	hddm_s::PhysicsEvent* pe = &record->getPhysicsEvent();
	pe->setEventNo(locEventLoop->GetJEvent().GetEventNumber());
	pe->setRunNo(locEventLoop->GetJEvent().GetRunNumber());
	//add a HitView which is necessary to create and save all of the data
	pe->addHitViews();
	hddm_s::HitView* hitv = &pe->getHitView();

	
	//Because HDDM groups hits by sub-unit of each detector we loop through the hits in each sub-dector and need to see if a sub-unit exists
	//in hddm so as to avoid duplication of sub-units in HDDM
	//====================================TPOL================================================
	for(uint i=0;i<TPOLHits.size();++i)
	{
		if(i==0)//if the TPOL has hits then on the first hit we need to build the TPOL
		{
			hitv->addTripletPolarimeters();
		}

		bool found=false;

		//get sectors.  We get the iterator outside the for loop so when the loop is broken the iterator is at the sector that
		//has been hit, making adding hits to it trivial
		hddm_s::TpolSectorList* TPOL_SectorList = &hitv->getTripletPolarimeter().getTpolSectors();
		hddm_s::TpolSectorList::iterator sectorIterator = TPOL_SectorList->begin();

		//loop over the unique sectors already created in the HDDM framework
		for(sectorIterator = TPOL_SectorList->begin(); sectorIterator != TPOL_SectorList->end(); sectorIterator++)
		{
		  //look to see if the same sub-unit of the sub-detector is hit again.  At worst this TPOLHit is new and all sectors already hit must be looped over
			if(int(TPOLHits[i]->ring) == sectorIterator->getRing() && TPOLHits[i]->sector == sectorIterator->getSector() )
			{
			        found=true;//we found it!
				break;//you can stop looping now to shorten the time and leave the iterator at the right spot
			}
		}

		//this part hasn't been hit yet.  We need to create it to hold the hit(s)
		if(found==false)
		{
		        hitv->getTripletPolarimeter().addTpolSectors();//currently the iterator is at "end", which is exactly where we want to add it
			sectorIterator = TPOL_SectorList->end()-1;//this might be a little redundant but the newly made sector is at end-1 so we set the iterator here
			sectorIterator->setRing(TPOLHits[i]->ring);//set the unique identifiers for the new sector
			sectorIterator->setSector(TPOLHits[i]->sector);
		}
		//Now that we either created the sub-unit or found the sub-unit we add a hit to it
		sectorIterator->addTpolHits();//either the iterator was at the sector that this hit belongs to or it is at the newly created sector
		hddm_s::TpolHitList* TPOL_HitList = &sectorIterator->getTpolHits();
		hddm_s::TpolHitList::iterator TPOL_HitIterator = TPOL_HitList->end()-1;//as above we ensure we are at the last element
		TPOL_HitIterator->setT(TPOLHits[i]->t);//and set the the proper values for the newly formed hit
		TPOL_HitIterator->setDE(TPOLHits[i]->dE);

	}


  //========================================TAGGER===========================================================

	if(TAGHHits.size() != uint(0) || TAGMHits.size() != uint(0) )//If either the TAGH or TAGM have hits create the Tagger
	{
		hitv->addTaggers();
	}

	//hodoscope
	for(uint i=0; i<TAGHHits.size(); ++i)
	{
		if(i==0)//if the hodoscope has at least one hit then on the first one make the set of hodoscope channels
		{
			hitv->getTagger().addHodoChannels();
		}

		bool found=false;
		//look to see if the same sub-unit of the sub-detector is hit again
		hddm_s::HodoChannelList* TAGH_ChannelList = &hitv->getTagger().getHodoChannels();
		hddm_s::HodoChannelList::iterator TAGH_ChannelIterator = TAGH_ChannelList->begin();

		//look to see if the same sub-unit of the sub-detector is hit again
		for(TAGH_ChannelIterator=TAGH_ChannelList->begin(); TAGH_ChannelIterator != TAGH_ChannelList->end(); TAGH_ChannelIterator++)
		{

			if(int(TAGHHits[i]->counter_id) == TAGH_ChannelIterator->getCounterId() )
			{
				found=true;
				break;
			}
		}

		if(found==false)
		{

			hitv->getTagger().addHodoChannels();
			TAGH_ChannelIterator = TAGH_ChannelList->end()-1;
			TAGH_ChannelIterator->setCounterId(TAGHHits[i]->counter_id);
			TAGH_ChannelIterator->setE(TAGHHits[i]->E);
		}

		//Add hits
		TAGH_ChannelIterator->addTaggerHits();
		hddm_s::TaggerHitList* TAGGER_HitList = &TAGH_ChannelIterator->getTaggerHits();
		hddm_s::TaggerHitList::iterator TAGGER_HitIterator = TAGGER_HitList->end()-1;
		TAGGER_HitIterator->setNpe(TAGHHits[i]->npe_fadc);
		TAGGER_HitIterator->setT(TAGHHits[i]->t);
		TAGGER_HitIterator->setTADC(TAGHHits[i]->time_fadc);

	}

	//microscope
	for(uint i=0;i<TAGMHits.size();++i)
	{
		if(i==0)
		{
			hitv->getTagger().addMicroChannels();
		}
		bool found=false;

		//look to see if the same sub-unit of the sub-detector is hit again
		hddm_s::MicroChannelList* TAGM_ChannelList= &hitv->getTagger().getMicroChannels();
		hddm_s::MicroChannelList::iterator TAGM_ChannelIterator=TAGM_ChannelList->begin();

		for(TAGM_ChannelIterator=TAGM_ChannelList->begin(); TAGM_ChannelIterator != TAGM_ChannelList->end(); TAGM_ChannelIterator++)
		{
			if(int(TAGMHits[i]->column) == TAGM_ChannelIterator->getColumn() && int(TAGMHits[i]->row) == TAGM_ChannelIterator->getRow() )
			{
				found = true;
				break;
			}
		}

		if(found == false)
		{

			hitv->getTagger().addMicroChannels();
			TAGM_ChannelIterator = TAGM_ChannelList->end()-1;
			TAGM_ChannelIterator->setColumn(TAGMHits[i]->column);
			TAGM_ChannelIterator->setRow(TAGMHits[i]->row);
			TAGM_ChannelIterator->setE(TAGMHits[i]->E);
		}

		TAGM_ChannelIterator->addTaggerHits();
		hddm_s::TaggerHitList* TAGGER_HitList = &TAGM_ChannelIterator->getTaggerHits();
		hddm_s::TaggerHitList::iterator TAGGER_HitIterator = TAGGER_HitList->end()-1;
		TAGGER_HitIterator->setNpe(TAGMHits[i]->npix_fadc);
		TAGGER_HitIterator->setT(TAGMHits[i]->t);
		TAGGER_HitIterator->setTADC(TAGMHits[i]->time_fadc);

	}



	//====================================FDC==========================================

	for(uint i=0;i<FDCHits.size();++i)
	{
		//if there is at least 1 FDC hit we need an FDC
		if(i==0)
		{
			hitv->addForwardDCs();
		}
		//look for the chamber by layer/module  (before searching for same strip/wire we must see if the chamber exists already
		bool foundChamber=false;
		hddm_s::FdcChamberList* FDC_ChamberList = &hitv->getForwardDC().getFdcChambers();
		hddm_s::FdcChamberList::iterator FDC_ChamberIterator = FDC_ChamberList->begin();

		for(FDC_ChamberIterator = FDC_ChamberList->begin(); FDC_ChamberIterator != FDC_ChamberList->end(); FDC_ChamberIterator++)
		{
			//it is important to note that in HDDM the Layer that is saved represents the clocking of the chamber.  1=0 degrees 2=60 3=120
			//in evio Layer is 1-3 where 1/3 are the cathode strips and 2 is the anode wire plane
			//the clocking can be calculated by accessing the glayer and using a little modular arithmatic
			if(((FDCHits[i]->gLayer-1)%3)+1 == FDC_ChamberIterator->getLayer() && FDCHits[i]->module == FDC_ChamberIterator->getModule() )
			{
				foundChamber = true;
				break;
			}
		}

		if(foundChamber == false)
		{

			hitv->getForwardDC().addFdcChambers();
			FDC_ChamberIterator = FDC_ChamberList->end()-1;
			FDC_ChamberIterator->setLayer(((FDCHits[i]->gLayer-1)%3)+1); //This will be dumped as the clocking and is not expected to match the Layer dumped from evio
			FDC_ChamberIterator->setModule(FDCHits[i]->module);

		}
		//now that we found the chamber (or created a new one) we need to see if it was a cathode or anode hit and see if it too is new
		if(FDCHits[i]->type == DFDCHit::AnodeWire)
		{

			//search either the wires or strips depending on the above
			bool found=false;

			hddm_s::FdcAnodeWireList* FDC_AnodeWireList= &FDC_ChamberIterator->getFdcAnodeWires();
			hddm_s::FdcAnodeWireList::iterator FDC_AnodeWireIterator = FDC_AnodeWireList->begin();

			for(FDC_AnodeWireIterator = FDC_AnodeWireList->begin(); FDC_AnodeWireIterator != FDC_AnodeWireList->end(); FDC_AnodeWireIterator++)
			{
				if(int(FDCHits[i]->element) == FDC_AnodeWireIterator->getWire() )
				{
					found=true;
					break;
				}
			}

			if(found == false)
			{
				FDC_ChamberIterator->addFdcAnodeWires();
				FDC_AnodeWireIterator = FDC_AnodeWireList->end()-1;
				FDC_AnodeWireIterator->setWire(FDCHits[i]->element);
			}

			FDC_AnodeWireIterator->addFdcAnodeHits();
			hddm_s::FdcAnodeHitList* FDC_AnodeWireHitList = &FDC_AnodeWireIterator->getFdcAnodeHits();
			hddm_s::FdcAnodeHitList::iterator FDC_AnodeWireHitIterator = FDC_AnodeWireHitList->end()-1;
			FDC_AnodeWireHitIterator->setT(FDCHits[i]->t);

		}
		else//this is a cathode hit.  Agnostic about being half length or not
		{
			bool found=false;

			hddm_s::FdcCathodeStripList* FDC_CathodeStripList = &FDC_ChamberIterator->getFdcCathodeStrips();
			hddm_s::FdcCathodeStripList::iterator FDC_CathodeStripIterator = FDC_CathodeStripList->begin();
			//look in the chamber to see if the cathode has already been hit
			for(FDC_CathodeStripIterator = FDC_CathodeStripList->begin(); FDC_CathodeStripIterator != FDC_CathodeStripList->end(); FDC_CathodeStripIterator++)
			{
				if(int(FDCHits[i]->element) == FDC_CathodeStripIterator->getStrip() && FDCHits[i]->plane == FDC_CathodeStripIterator->getPlane() )
				{
					found = true;
					break;
				}
			}

			if(found == false)
			{

				FDC_ChamberIterator->addFdcCathodeStrips();
				FDC_CathodeStripIterator=FDC_CathodeStripList->end()-1;
				FDC_CathodeStripIterator->setStrip(FDCHits[i]->element);
				FDC_CathodeStripIterator->setPlane(FDCHits[i]->plane);

			}

			FDC_CathodeStripIterator->addFdcCathodeHits();
			hddm_s::FdcCathodeHitList* FDC_CathodeStripHitList = &FDC_CathodeStripIterator->getFdcCathodeHits();
			hddm_s::FdcCathodeHitList::iterator FDC_CathodeStripHitIterator = FDC_CathodeStripHitList->end()-1;
			FDC_CathodeStripHitIterator->setT(FDCHits[i]->t);
			FDC_CathodeStripHitIterator->setQ(FDCHits[i]->q);

		}

	}



	//===============================================PS=============================================
	//in hddm the PS is in two parts PSFine and PSCoarse
	for(uint i=0; i<PSHits.size(); ++i)//Do the fine hits first
	{
		if(i == 0)
		{
			hitv->addPairSpectrometerFines(); //only create PS Fine Tiles if there is a hit in the
		}
		bool found=false;

		hddm_s::PsTileList* PS_TileList = &hitv->getPairSpectrometerFine().getPsTiles();
		hddm_s::PsTileList::iterator PS_TileIterator = PS_TileList->begin();

		//check for the same PS tile is hit
		for(PS_TileIterator = PS_TileList->begin(); PS_TileIterator != PS_TileList->end(); PS_TileIterator++)
		{
			if(int(PSHits[i]->arm) == PS_TileIterator->getArm() && PSHits[i]->column == PS_TileIterator->getColumn() )
			{
				found=true;
				break;
			}
		}

		if(found == false)
		{

			hitv->getPairSpectrometerFine().addPsTiles();
			PS_TileIterator = PS_TileList->end()-1;
			PS_TileIterator->setArm(PSHits[i]->arm);
			PS_TileIterator->setColumn(PSHits[i]->column);
		}

		PS_TileIterator->addPsHits();
		hddm_s::PsHitList* PS_HitList = &PS_TileIterator->getPsHits();
		hddm_s::PsHitList::iterator PS_HitIterator = PS_HitList->end()-1;
		PS_HitIterator->setT(PSHits[i]->t);
		PS_HitIterator->setDE(PSHits[i]->E);

	}
	//--------------------COARSE------------------------------------------
	for(uint i=0;i<PSCHits.size();++i)//repeat for the coarse hits
	{
		if(i==0)
		{
			hitv->addPairSpectrometerCoarses();
		}

		bool found=false;

		hddm_s::PscPaddleList* PS_PaddleList = &hitv->getPairSpectrometerCoarse().getPscPaddles();
		hddm_s::PscPaddleList::iterator PS_PaddleIterator = PS_PaddleList->begin();

		for(PS_PaddleIterator = PS_PaddleList->begin(); PS_PaddleIterator != PS_PaddleList->end(); PS_PaddleIterator++)
		{
			if(int(PSCHits[i]->arm) == PS_PaddleIterator->getArm() && PSCHits[i]->module == PS_PaddleIterator->getModule() )
			{
				found = true;
				break;
			}
		}

		if(found == false)
		{
			hitv->getPairSpectrometerCoarse().addPscPaddles();
			PS_PaddleIterator = PS_PaddleList->end()-1;
			PS_PaddleIterator->setArm(PSCHits[i]->arm);
			PS_PaddleIterator->setModule(PSCHits[i]->module);
		}

		PS_PaddleIterator->addPscHits();
		hddm_s::PscHitList* PSC_HitList = &PS_PaddleIterator->getPscHits();
		hddm_s::PscHitList::iterator pschitit = PSC_HitList->end()-1;
		pschitit->setT(PSCHits[i]->t);

	}



	//================================================Start Counter======================================

	for(uint i=0; i<SCHits.size(); ++i)
	{
		if(i == 0)
		{
			hitv->addStartCntrs();//if we have a hit add a start counter
		}
		bool found=false;

		hddm_s::StcPaddleList* SC_CounterList = &hitv->getStartCntr().getStcPaddles();
		hddm_s::StcPaddleList::iterator SC_CounterIterator = SC_CounterList->begin();

		//see if the sector has already been hit
		for(SC_CounterIterator = SC_CounterList->begin(); SC_CounterIterator != SC_CounterList->end(); SC_CounterIterator++)
		{
			if(SCHits[i]->sector == SC_CounterIterator->getSector() )
			{
				found = true;
				break;
			}
		}

		if(found == false)
		{
			hitv->getStartCntr().addStcPaddles();
			SC_CounterIterator=SC_CounterList->end()-1;
			SC_CounterIterator->setSector(SCHits[i]->sector);
		}

		SC_CounterIterator->addStcHits();
		hddm_s::StcHitList* schitl = &SC_CounterIterator->getStcHits();
		hddm_s::StcHitList::iterator schitit = schitl->end()-1;
		schitit->setT(SCHits[i]->t);
		schitit->setDE(SCHits[i]->dE);
	}


	//============================================BCAL=========================================

	//The BCAL is unique in that it needs the DigiHits put in HDDM ADC/TDC done separately
	if(BCALDigiHits.size() != uint(0) || BCALTDCDigiHits.size() != uint(0) )
	{
		hitv->addBarrelEMcals(); //still only need one BCAL if we have a hit of some kind
	}
	  //-------------------------------ADC--------------------------
	for(uint i=0; i<BCALDigiHits.size(); ++i)
	{
		bool found = false;

		hddm_s::BcalCellList* BCAL_CellList= &hitv->getBarrelEMcal().getBcalCells();
		hddm_s::BcalCellList::iterator BCAL_CellIterator=BCAL_CellList->begin();

		//note: these cells being searched are the same ones searched for TDCs
		for(BCAL_CellIterator = BCAL_CellList->begin(); BCAL_CellIterator != BCAL_CellList->end(); BCAL_CellIterator++)
		{
			if(BCALDigiHits[i]->sector == BCAL_CellIterator->getSector() && BCALDigiHits[i]->layer == BCAL_CellIterator->getLayer() && BCALDigiHits[i]->module == BCAL_CellIterator->getModule() )
			{
				found = true;
				break;
			}
		}

		if(found == false)
		{
			hitv->getBarrelEMcal().addBcalCells();
			BCAL_CellIterator = BCAL_CellList->end()-1;
			BCAL_CellIterator->setLayer(BCALDigiHits[i]->layer);
			BCAL_CellIterator->setSector(BCALDigiHits[i]->sector);
			BCAL_CellIterator->setModule(BCALDigiHits[i]->module);
		}

		BCAL_CellIterator->addBcalfADCDigiHits();
		hddm_s::BcalfADCDigiHitList* BCAL_FADCDigiHitList = &BCAL_CellIterator->getBcalfADCDigiHits();
		hddm_s::BcalfADCDigiHitList::iterator BCAL_FADCDigiHitIterator = BCAL_FADCDigiHitList->end()-1;
		BCAL_FADCDigiHitIterator->setEnd(BCALDigiHits[i]->end);
		BCAL_FADCDigiHitIterator->setPulse_time(BCALDigiHits[i]->pulse_time);
		BCAL_FADCDigiHitIterator->setPulse_integral(BCALDigiHits[i]->pulse_integral);
	}

	//------------------------TDC-----------------------------
	for(uint i=0; i<BCALTDCDigiHits.size(); ++i)
	{
		bool found = false;

		hddm_s::BcalCellList* BCAL_CellList = &hitv->getBarrelEMcal().getBcalCells();
		hddm_s::BcalCellList::iterator BCAL_CellIterator = BCAL_CellList->begin();

		for(BCAL_CellIterator = BCAL_CellList->begin(); BCAL_CellIterator != BCAL_CellList->end(); BCAL_CellIterator++)
		{
			if(BCALTDCDigiHits[i]->sector == uint(BCAL_CellIterator->getSector()) && BCALTDCDigiHits[i]->layer == uint(BCAL_CellIterator->getLayer()) && BCALTDCDigiHits[i]->module == uint(BCAL_CellIterator->getModule()) )
			{
				found = true;
				break;
			}
		}

		if(found == false)
		{
			hitv->getBarrelEMcal().addBcalCells();
			BCAL_CellIterator = BCAL_CellList->end()-1;
			BCAL_CellIterator->setLayer(BCALTDCDigiHits[i]->layer);
			BCAL_CellIterator->setSector(BCALTDCDigiHits[i]->sector);
			BCAL_CellIterator->setModule(BCALTDCDigiHits[i]->module);

		}

		BCAL_CellIterator->addBcalTDCDigiHits();
		hddm_s::BcalTDCDigiHitList* BCAL_TDCDigiHitList = &BCAL_CellIterator->getBcalTDCDigiHits();
		hddm_s::BcalTDCDigiHitList::iterator BCAL_TDCDigiHitIterator = BCAL_TDCDigiHitList->end()-1;
		BCAL_TDCDigiHitIterator->setEnd(BCALTDCDigiHits[i]->end);
		BCAL_TDCDigiHitIterator->setTime(BCALTDCDigiHits[i]->time);
	}


	//========================================FCAL=========================================================

	for(uint i=0; i<FCALHits.size(); ++i)
	{
		if(i == 0)
		{
			hitv->addForwardEMcals();
		}
		bool found = false;
		//FCAL only has one hit per block per event so we need not search
		hddm_s::FcalBlockList* FCAL_BlockList = &hitv->getForwardEMcal().getFcalBlocks();
		hddm_s::FcalBlockList::iterator FCAL_BlockIterator = FCAL_BlockList->begin();

		for(FCAL_BlockIterator = FCAL_BlockList->begin(); FCAL_BlockIterator != FCAL_BlockList->end(); FCAL_BlockIterator++)
				{
					if(FCALHits[i]->row==FCAL_BlockIterator->getRow() && FCALHits[i]->column==FCAL_BlockIterator->getColumn())
					{
						found=true;
						break;
					}
				}

		if(found==false)
		{
			hitv->getForwardEMcal().addFcalBlocks();
			FCAL_BlockIterator=FCAL_BlockList->end()-1;
			FCAL_BlockIterator->setColumn(FCALHits[i]->column);
			FCAL_BlockIterator->setRow(FCALHits[i]->row);
		}


		FCAL_BlockIterator->addFcalHits();
		hddm_s::FcalHitList* FCAL_HitList = &FCAL_BlockIterator->getFcalHits();
		hddm_s::FcalHitList::iterator FCAL_HitIterator = FCAL_HitList->end()-1;
		FCAL_HitIterator->setT(FCALHits[i]->t);
		FCAL_HitIterator->setE(FCALHits[i]->E);
	}


  //=============================================TOF=====================================================

	for(uint i=0; i<TOFHits.size(); ++i)
	{

		if(i == 0)
		{
			hitv->addForwardTOFs();//if we have a hit add the TOF
		}

		bool found = false;
		//gotta look for the same plane/bar
		hddm_s::FtofCounterList* TOF_CounterList = &hitv->getForwardTOF().getFtofCounters();
		hddm_s::FtofCounterList::iterator TOF_CounterIterator = TOF_CounterList->begin();

		for(TOF_CounterIterator = TOF_CounterList->begin(); TOF_CounterIterator != TOF_CounterList->end(); TOF_CounterIterator++)
		{
			if(TOFHits[i]->bar==TOF_CounterIterator->getBar() && TOFHits[i]->plane==TOF_CounterIterator->getPlane())
			{
				found=true;
				break;
			}
		}

		if(found==false)
		{
			hitv->getForwardTOF().addFtofCounters();
			TOF_CounterIterator=TOF_CounterList->end()-1;
			TOF_CounterIterator->setPlane(TOFHits[i]->plane);
			TOF_CounterIterator->setBar(TOFHits[i]->bar);

		}

		TOF_CounterIterator->addFtofHits();
		hddm_s::FtofHitList* ftofhitl=&TOF_CounterIterator->getFtofHits();
		hddm_s::FtofHitList::iterator ftofhitit=ftofhitl->end()-1;
		ftofhitit->setEnd(TOFHits[i]->end);
		ftofhitit->setT(TOFHits[i]->t);//walk corrected time
		ftofhitit->setDE(TOFHits[i]->dE);
	}


	//============================CDC=============================================

	for(uint i=0; i<CDCHits.size(); ++i)
	{
		if(i == 0)
		{
			hitv->addCentralDCs(); //if we have a hit then add the CDC
		}
		//CDC only has one hit per block per event so we need not search
		hitv->getCentralDC().addCdcStraws();

		hddm_s::CdcStrawList* CDC_StrawList = &hitv->getCentralDC().getCdcStraws();
		hddm_s::CdcStrawList::iterator CDC_StrawIterator = CDC_StrawList->end()-1;
		CDC_StrawIterator->setRing(CDCHits[i]->ring);
		CDC_StrawIterator->setStraw(CDCHits[i]->straw);

		CDC_StrawIterator->addCdcStrawHits();
		hddm_s::CdcStrawHitList* strawhitl = &CDC_StrawIterator->getCdcStrawHits();
		hddm_s::CdcStrawHitList::iterator cdcstrawhitit = strawhitl->end()-1;
		cdcstrawhitit->setQ(CDCHits[i]->q);
		cdcstrawhitit->setT(CDCHits[i]->t);
	}


	//*fout << *record; //stream the new record into the file

	// write the resulting record to the output stream
	string locOutputFileName = Get_OutputFileName(locOutputFileNameSubString);
	bool locWriteStatus = Write_HDDMEvent(locOutputFileName, *record);
	delete record;
	return locWriteStatus;
}

string DEventWriterHDDM::Get_OutputFileName(string locOutputFileNameSubString) const
{
	string locOutputFileName = dOutputFileBaseName;
	if (locOutputFileNameSubString != "")
		locOutputFileName += string("_") + locOutputFileNameSubString;
	return (locOutputFileName + string(".hddm"));
}

bool DEventWriterHDDM::Write_HDDMEvent(string locOutputFileName, hddm_s::HDDM& locRecord) const
{
	japp->WriteLock("HDDMWriter");
	{
		//check to see if the HDDM file is open
		if(Get_HDDMOutputFilePointers().find(locOutputFileName) != Get_HDDMOutputFilePointers().end())
		{
			//open: get pointer, write event
			hddm_s::ostream* locOutputHDDMFileStream = Get_HDDMOutputFilePointers()[locOutputFileName].second;
			*(locOutputHDDMFileStream) << locRecord;
			japp->Unlock("HDDMWriter");
			return true;
		}

		//not open: open it
		pair<ofstream*, hddm_s::ostream*> locHDDMFilePointers(NULL, NULL);
		locHDDMFilePointers.first = new ofstream(locOutputFileName.c_str());
		if(!locHDDMFilePointers.first->is_open())
		{
			//failed to open
			delete locHDDMFilePointers.first;
			japp->Unlock("HDDMWriter");
			return false;
		}
		locHDDMFilePointers.second = new hddm_s::ostream(*locHDDMFilePointers.first);

		// enable on-the-fly bzip2 compression on output stream
		if(HDDM_USE_COMPRESSION)
		{
			jout << " Enabling bz2 compression of output HDDM file stream" << std::endl;
			locHDDMFilePointers.second->setCompression(hddm_s::k_bz2_compression);
		}
		else
			jout << " HDDM compression disabled" << std::endl;

		// enable a CRC data integrity check at the end of each event record
		if(HDDM_USE_INTEGRITY_CHECKS)
		{
			jout << " Enabling CRC data integrity check in output HDDM file stream" << std::endl;
			locHDDMFilePointers.second->setIntegrityChecks(hddm_s::k_crc32_integrity);
		}
		else
			jout << " HDDM integrity checks disabled" << std::endl;


		//write the event
		*(locHDDMFilePointers.second) << locRecord;

		//store the stream pointers
		Get_HDDMOutputFilePointers()[locOutputFileName] = locHDDMFilePointers;
	}
	japp->Unlock("HDDMWriter");

	return true;
}

DEventWriterHDDM::~DEventWriterHDDM(void)
{
	japp->WriteLock("HDDMWriter");
	{
		--Get_NumEventWriterThreads();
		if(Get_NumEventWriterThreads() > 0)
		{
			japp->Unlock("HDDMWriter");
			return; //not the last thread writing to HDDM files
		}

		//last thread writing to HDDM files: close all files and free all memory
		map<string, pair<ofstream*, hddm_s::ostream*> >::iterator locIterator;
		for(locIterator = Get_HDDMOutputFilePointers().begin(); locIterator != Get_HDDMOutputFilePointers().end(); ++locIterator)
		{
			string locOutputFileName = locIterator->first;
			if (locIterator->second.second != NULL)
				delete locIterator->second.second;
			if (locIterator->second.first != NULL)
				delete locIterator->second.first;
			std::cout << "Closed HDDM file " << locOutputFileName << std::endl;
		}
		Get_HDDMOutputFilePointers().clear();
	}
	japp->Unlock("HDDMWriter");
}

int32_t DEventWriterHDDM::Convert_UnsignedIntToSigned(uint32_t locUnsignedInt) const
{
	//Convert uint32_t to int32_t
	//Scheme:
		//If bit 32 is zero, then the int32_t is the same as the uint32_t: Positive or zero
		//If bit 32 is one, and at least one other bit is 1, then the int32_t is -1 * uint32_t (after stripping the top bit)
		//If bit 32 is one, and all other bits are zero, then the int32_t is the minimum int: -(2^31)
	if((locUnsignedInt & 0x80000000) == 0)
		return int32_t(locUnsignedInt); //bit 32 is zero: positive or zero

	//bit 32 is 1. see if there is another bit set
	int32_t locTopBitStripped = int32_t(locUnsignedInt & uint32_t(0x7FFFFFFF)); //strip the top bit
	if(locTopBitStripped == 0)
		return numeric_limits<int32_t>::min(); //no other bit is set: minimum int
	return -1*locTopBitStripped; //return the negative
}

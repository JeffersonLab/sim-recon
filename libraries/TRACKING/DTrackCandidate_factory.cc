// $Id$
//
//    File: DTrackCandidate_factory.cc
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//


#include "DTrackCandidate_factory.h"

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory::evnt(JEventLoop *loop, int eventnumber)
{

	vector<const DTrackCandidate*> cdctrackcandidates;
	vector<const DTrackCandidate*> fdctrackcandidates;

	loop->Get(cdctrackcandidates, "CDC");
	loop->Get(fdctrackcandidates, "FDCCathodes");
	
	for(unsigned int i=0; i<cdctrackcandidates.size(); i++){
		DTrackCandidate *can = new DTrackCandidate;
		const DTrackCandidate *srccan = cdctrackcandidates[i];
		
		can->setMass(srccan->mass());
		can->setMomentum(srccan->momentum());
		can->setPosition(srccan->position());
		can->setCharge(srccan->charge());
		
		_data.push_back(can);
	}

	for(unsigned int i=0; i<fdctrackcandidates.size(); i++){
		DTrackCandidate *can = new DTrackCandidate;
		const DTrackCandidate *srccan = fdctrackcandidates[i];
		
		can->setMass(srccan->mass());
		can->setMomentum(srccan->momentum());
		can->setPosition(srccan->position());
		can->setCharge(srccan->charge());
		
		_data.push_back(can);
	}

	return NOERROR;
}


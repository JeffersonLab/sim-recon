// $Id$
//
//		File: DTOFPoint_factory.cc
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//
// Modified: Wed Feb 12 13:23:42 EST 2014 B. Zihlamnn
//					 use new TOF geometry with narrow long paddles
//					 and short paddles #22 and #23 for both north and south
//					 

#include <cassert>
#include <cmath>
using namespace std;

#include "DTOFPoint_factory.h"

bool Compare_TOFSpacetimeHitMatches_Distance(const DTOFPoint_factory::tof_spacetimehitmatch_t& locTOFSpacetimeHitMatch1, const DTOFPoint_factory::tof_spacetimehitmatch_t& locTOFSpacetimeHitMatch2)
{
	if(locTOFSpacetimeHitMatch2.dPositionWellDefinedFlag != locTOFSpacetimeHitMatch1.dPositionWellDefinedFlag)
		return locTOFSpacetimeHitMatch1.dPositionWellDefinedFlag; //one hit position is well defined and the other is not
	return (locTOFSpacetimeHitMatch1.delta_r < locTOFSpacetimeHitMatch2.delta_r);
};

//------------------
// brun
//------------------
jerror_t DTOFPoint_factory::brun(JEventLoop *loop, int runnumber)
{

	map<string, double> tofparms;
 	if( !loop->GetCalib("TOF/tof_parms", tofparms))
	{
		//cout<<"DTOFPoint_factory: loading values from TOF data base"<<endl;
		VELOCITY = tofparms["TOF_C_EFFECTIVE"];
		HALFPADDLE = tofparms["TOF_HALFPADDLE"];
		BARWIDTH = tofparms["TOF_PADDLEWIDTH"];
		E_THRESHOLD = tofparms["TOF_E_THRESHOLD"];
		ATTEN_LENGTH = tofparms["TOF_ATTEN_LENGTH"];

		E_THRESHOLD = 0.0;
	}
	else
	{
		cout << "DTOFPoint_factory: Error loading values from TOF data base" <<endl;
		VELOCITY = 15.; // set to some reasonable value
		HALFPADDLE = 126; // set to some reasonable value
		BARWIDTH = 6.;
		E_THRESHOLD = 0.0005;
		ATTEN_LENGTH = 400.;
	}

	if(eventLoop->GetCalib("TOF/propagation_speed", propagation_speed))
		jout << "Error loading /TOF/propagation_speed !" << endl;

	loop->Get(TOFGeom);

	//	dPositionMatchCut_DoubleEnded = 0.5*BARWIDTH + 1.5*3.0; //max if perfect precision (1/2 bar width) + ~3 sigma
	dPositionMatchCut_DoubleEnded = 1.5*BARWIDTH;

	return NOERROR;
}

DTOFPoint_factory::tof_spacetimehit_t* DTOFPoint_factory::Get_TOFSpacetimeHitResource(void)
{
	tof_spacetimehit_t* locTOFSpacetimeHit;
	if(dTOFSpacetimeHitPool_Available.empty())
	{
		locTOFSpacetimeHit = new tof_spacetimehit_t;
		dTOFSpacetimeHitPool_All.push_back(locTOFSpacetimeHit);
	}
	else
	{
		locTOFSpacetimeHit = dTOFSpacetimeHitPool_Available.back();
		dTOFSpacetimeHitPool_Available.pop_back();
	}
	return locTOFSpacetimeHit;
}

//------------------
// evnt
//------------------
jerror_t DTOFPoint_factory::evnt(JEventLoop *loop, int eventnumber)
{	
	// delete pool size if too large, preventing memory-leakage-like behavor.
	if(dTOFSpacetimeHitPool_All.size() > MAX_TOFSpacetimeHitPoolSize)
	{
		for(size_t loc_i = MAX_TOFSpacetimeHitPoolSize; loc_i < dTOFSpacetimeHitPool_All.size(); ++loc_i)
			delete dTOFSpacetimeHitPool_All[loc_i];
		dTOFSpacetimeHitPool_All.resize(MAX_TOFSpacetimeHitPoolSize);
	}
	dTOFSpacetimeHitPool_Available = dTOFSpacetimeHitPool_All;

	vector<const DTOFPaddleHit*> locTOFHitVector;
	loop->Get(locTOFHitVector);

	// create the hit spacetime information
	deque<tof_spacetimehit_t*> locTOFSpacetimeHits_Horizontal, locTOFSpacetimeHits_Vertical;
	set<tof_spacetimehit_t*> locUnusedTOFSpacetimeHits;
	for(size_t loc_i = 0; loc_i < locTOFHitVector.size(); ++loc_i)
	{
		const DTOFPaddleHit* locTOFHit = locTOFHitVector[loc_i];
		if(!((locTOFHit->E_north > E_THRESHOLD) || (locTOFHit->E_south > E_THRESHOLD)))
			continue;

		if(locTOFHit->orientation) //horizontal
		{
			tof_spacetimehit_t* locSpacetimeHit = Build_TOFSpacetimeHit_Horizontal(locTOFHit);
			locTOFSpacetimeHits_Horizontal.push_back(locSpacetimeHit);
			locUnusedTOFSpacetimeHits.insert(locSpacetimeHit);
		}
		else //vertical
		{
			tof_spacetimehit_t* locSpacetimeHit = Build_TOFSpacetimeHit_Vertical(locTOFHit);
			locTOFSpacetimeHits_Vertical.push_back(locSpacetimeHit);
			locUnusedTOFSpacetimeHits.insert(locSpacetimeHit);
		}
	}


	//find matches between planes and sort them by delta-r
	deque<tof_spacetimehitmatch_t> locTOFSpacetimeHitMatches; //use list for sorting, vector for resource pool
	for(size_t loc_i = 0; loc_i < locTOFSpacetimeHits_Horizontal.size(); ++loc_i)
	{
		tof_spacetimehit_t* locTOFSpacetimeHit_Horizontal = locTOFSpacetimeHits_Horizontal[loc_i];
		for(size_t loc_j = 0; loc_j < locTOFSpacetimeHits_Vertical.size(); ++loc_j)
		{
			tof_spacetimehit_t* locTOFSpacetimeHit_Vertical = locTOFSpacetimeHits_Vertical[loc_j];

			//If the position along BOTH paddles is not well defined, it's not clear that these hits belong together
				//If the hit multiplicity was low, could just generate all matches
				//However, the hit multiplicity is generally very high due to hits near the beamline: would have TOF hits everywhere
				//Therefore, don't keep as match: register separately, let track / tof matching salvage the situation
			if((!locTOFSpacetimeHit_Horizontal->dPositionWellDefinedFlag) && (!locTOFSpacetimeHit_Vertical->dPositionWellDefinedFlag))
				continue; //have no idea whether these hits go together: assume they don't

			float locDeltaX = locTOFSpacetimeHit_Horizontal->x - locTOFSpacetimeHit_Vertical->x;
			if(fabs(locDeltaX) > locTOFSpacetimeHit_Horizontal->pos_cut)
				continue;
			float locDeltaY = locTOFSpacetimeHit_Horizontal->y - locTOFSpacetimeHit_Vertical->y;
			if(fabs(locDeltaY) > locTOFSpacetimeHit_Vertical->pos_cut)
				continue;

			float locDeltaT = locTOFSpacetimeHit_Horizontal->t - locTOFSpacetimeHit_Vertical->t;
			float locTimeCut = (locTOFSpacetimeHit_Horizontal->t_cut > locTOFSpacetimeHit_Vertical->t_cut) ? locTOFSpacetimeHit_Horizontal->t_cut : locTOFSpacetimeHit_Vertical->t_cut;
			if(fabs(locDeltaT) > locTimeCut)
				continue;

			tof_spacetimehitmatch_t locTOFSpacetimeHitMatch;
			locTOFSpacetimeHitMatch.delta_t = locDeltaT;
			locTOFSpacetimeHitMatch.delta_r = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
			locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Horizontal = locTOFSpacetimeHit_Horizontal;
			locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Vertical = locTOFSpacetimeHit_Vertical;
			locTOFSpacetimeHitMatch.dPositionWellDefinedFlag = (locTOFSpacetimeHit_Horizontal->dPositionWellDefinedFlag == locTOFSpacetimeHit_Vertical->dPositionWellDefinedFlag);
			locTOFSpacetimeHitMatches.push_back(locTOFSpacetimeHitMatch);
		}
	}
	std::sort(locTOFSpacetimeHitMatches.begin(), locTOFSpacetimeHitMatches.end(), Compare_TOFSpacetimeHitMatches_Distance); //sort matches by delta_r


	// create DTOFPoints, in order of best matches (by delta_r)
	for(size_t loc_i = 0; loc_i < locTOFSpacetimeHitMatches.size(); ++loc_i)
	{
		tof_spacetimehit_t* locTOFSpacetimeHit_Horizontal = locTOFSpacetimeHitMatches[loc_i].dTOFSpacetimeHit_Horizontal;
		if(locUnusedTOFSpacetimeHits.find(locTOFSpacetimeHit_Horizontal) == locUnusedTOFSpacetimeHits.end())
			continue; //hit used in a previous successful match

		tof_spacetimehit_t* locTOFSpacetimeHit_Vertical = locTOFSpacetimeHitMatches[loc_i].dTOFSpacetimeHit_Vertical;
		if(locUnusedTOFSpacetimeHits.find(locTOFSpacetimeHit_Vertical) == locUnusedTOFSpacetimeHits.end())
			continue; //hit used in a previous successful match

		Create_TOFPoint(locTOFSpacetimeHit_Horizontal, locTOFSpacetimeHit_Vertical);

		//remove used hits from the unused list
		locUnusedTOFSpacetimeHits.erase(locTOFSpacetimeHit_Horizontal);
		locUnusedTOFSpacetimeHits.erase(locTOFSpacetimeHit_Vertical);
	}

	// Loop over unused/unmatched TOF Spacetime hits, and create separate DTOFPoint's for them: 
		// If the position is NOT well defined, set as NaN
	set<tof_spacetimehit_t*>::iterator locSetIterator = locUnusedTOFSpacetimeHits.begin();
	for(; locSetIterator != locUnusedTOFSpacetimeHits.end(); ++locSetIterator)
	{
		const tof_spacetimehit_t* locTOFSpacetimeHit = *locSetIterator;
		float locPointZ = locTOFSpacetimeHit->TOFHit->orientation ? TOFGeom[0]->CenterHPlane : TOFGeom[0]->CenterVPlane; //true/false = horizontal/vertical
		if(locTOFSpacetimeHit->dPositionWellDefinedFlag)
		{
			DTOFPoint* locTOFPoint = new DTOFPoint;
			locTOFPoint->AddAssociatedObject(locTOFSpacetimeHit->TOFHit);

			locTOFPoint->pos.SetXYZ(locTOFSpacetimeHit->x, locTOFSpacetimeHit->y, locPointZ);
			locTOFPoint->t = locTOFSpacetimeHit->t;
			locTOFPoint->dE = locTOFSpacetimeHit->TOFHit->dE;

			_data.push_back(locTOFPoint);
		}
		else //position not well defined: save, but report position as NaN
		{
continue;
			DTOFPoint* locTOFPoint = new DTOFPoint;
			locTOFPoint->AddAssociatedObject(locTOFSpacetimeHit->TOFHit);

			float locPointX = locTOFSpacetimeHit->TOFHit->orientation ? numeric_limits<float>::quiet_NaN() : locTOFSpacetimeHit->x; //true/false = horizontal/vertical
			float locPointY = locTOFSpacetimeHit->TOFHit->orientation ? locTOFSpacetimeHit->y : numeric_limits<float>::quiet_NaN(); //true/false = horizontal/vertical
			locTOFPoint->pos.SetXYZ(locPointX, locPointY, locPointZ);
			locTOFPoint->t = locTOFSpacetimeHit->t;

			float locEnergy = numeric_limits<float>::quiet_NaN();
			if(locTOFSpacetimeHit->TOFHit->E_north > E_THRESHOLD)
				locEnergy = locTOFSpacetimeHit->TOFHit->E_north;
			else if(locTOFSpacetimeHit->TOFHit->E_south > E_THRESHOLD)
				locEnergy = -1.0*locTOFSpacetimeHit->TOFHit->E_south;
			locTOFPoint->dE = locEnergy;

			_data.push_back(locTOFPoint);
		}
	}

	return NOERROR;
}

DTOFPoint_factory::tof_spacetimehit_t* DTOFPoint_factory::Build_TOFSpacetimeHit_Horizontal(const DTOFPaddleHit* locTOFHit)
{
	tof_spacetimehit_t* locTOFSpacetimeHit = Get_TOFSpacetimeHitResource();
	locTOFSpacetimeHit->TOFHit = locTOFHit;

	int bar = locTOFHit->bar;
	int id = 44*locTOFHit->orientation + locTOFHit->bar - 1;
	double v = propagation_speed[id];

	if((locTOFHit->bar < TOFGeom[0]->FirstShortBar) || (locTOFHit->bar > TOFGeom[0]->LastShortBar)) //double-ended bars
	{
		locTOFSpacetimeHit->dIsDoubleEndedBar = true;
		locTOFSpacetimeHit->y = TOFGeom[0]->bar2y(bar);
		if(locTOFHit->meantime != locTOFHit->meantime)
		{
			//NaN: only one energy hit above threshold on the double-ended bar
			locTOFSpacetimeHit->dPositionWellDefinedFlag = false;
			locTOFSpacetimeHit->x = 0.0;
			if(locTOFHit->E_north > E_THRESHOLD)
				locTOFSpacetimeHit->t = locTOFHit->t_north - HALFPADDLE/v;
			else
				locTOFSpacetimeHit->t = locTOFHit->t_south - HALFPADDLE/v;

			locTOFSpacetimeHit->pos_cut = 252.0;
			locTOFSpacetimeHit->t_cut = 10.0;
		}
		else
		{
			locTOFSpacetimeHit->dPositionWellDefinedFlag = true;
			locTOFSpacetimeHit->x = locTOFHit->pos;
			locTOFSpacetimeHit->t = locTOFHit->meantime;
			locTOFSpacetimeHit->pos_cut = dPositionMatchCut_DoubleEnded;
			// locTOFSpacetimeHit->t_cut = 1.0;
			locTOFSpacetimeHit->t_cut=10.0;
		}

		//printf("h: x %f y %f\n",locTOFSpacetimeHit->x,locTOFSpacetimeHit->y);
		return locTOFSpacetimeHit;
	}

	//single-ended bars
	locTOFSpacetimeHit->dIsDoubleEndedBar = false;
	locTOFSpacetimeHit->dPositionWellDefinedFlag = false;
	locTOFSpacetimeHit->pos_cut = 120.;
	locTOFSpacetimeHit->t_cut = 10.;

	if(locTOFHit->t_south != 0.)
	{
		locTOFSpacetimeHit->y = TOFGeom[0]->bar2y(bar,1);
		locTOFSpacetimeHit->x = -66.; //paddle extends from -6 -> -126
		locTOFSpacetimeHit->t = locTOFHit->t_south - 60.0/v; //60 is HALFPADDLE for single-ended bars
	}
	else
	{
		locTOFSpacetimeHit->y = TOFGeom[0]->bar2y(bar,0);
		locTOFSpacetimeHit->x = 66.; //paddle extends from 6 -> 126
		locTOFSpacetimeHit->t = locTOFHit->t_north - 60.0/v;
	}

	//printf("h: x %f y %f\n",locTOFSpacetimeHit->x,locTOFSpacetimeHit->y);
	return locTOFSpacetimeHit;
}

DTOFPoint_factory::tof_spacetimehit_t* DTOFPoint_factory::Build_TOFSpacetimeHit_Vertical(const DTOFPaddleHit* locTOFHit)
{
	tof_spacetimehit_t* locTOFSpacetimeHit = Get_TOFSpacetimeHitResource();
	locTOFSpacetimeHit->TOFHit = locTOFHit;

	int bar = locTOFHit->bar;
	int id = 44*locTOFHit->orientation + locTOFHit->bar - 1;
	double v = propagation_speed[id];

	if((locTOFHit->bar < TOFGeom[0]->FirstShortBar) || (locTOFHit->bar > TOFGeom[0]->LastShortBar))
	{
		//double-ended bars
		locTOFSpacetimeHit->dIsDoubleEndedBar = true;
		locTOFSpacetimeHit->x = TOFGeom[0]->bar2y(bar);
		if(locTOFHit->meantime != locTOFHit->meantime)
		{
			//NaN: only one energy hit above threshold on the double-ended bar
			locTOFSpacetimeHit->dPositionWellDefinedFlag = false;
			locTOFSpacetimeHit->y = 0.0;
			if(locTOFHit->E_north > E_THRESHOLD)
				locTOFSpacetimeHit->t = locTOFHit->t_north - HALFPADDLE/v;
			else
				locTOFSpacetimeHit->t = locTOFHit->t_south - HALFPADDLE/v;
			locTOFSpacetimeHit->pos_cut = 252.0;
			locTOFSpacetimeHit->t_cut = 10.0;
		}
		else
		{
			locTOFSpacetimeHit->dPositionWellDefinedFlag = true;
			locTOFSpacetimeHit->y = locTOFHit->pos;
			locTOFSpacetimeHit->t = locTOFHit->meantime;
			locTOFSpacetimeHit->pos_cut = dPositionMatchCut_DoubleEnded;
			//	locTOFSpacetimeHit->t_cut = 1.0;
			locTOFSpacetimeHit->t_cut=10.0;
		}

		//printf("h: x %f y %f\n",locTOFSpacetimeHit->x,locTOFSpacetimeHit->y);
		return locTOFSpacetimeHit;
	}

	//single-ended bars
	locTOFSpacetimeHit->dIsDoubleEndedBar = false;
	locTOFSpacetimeHit->dPositionWellDefinedFlag = false;
	locTOFSpacetimeHit->pos_cut = 120.;
	locTOFSpacetimeHit->t_cut = 10.;
	if(locTOFHit->t_south != 0.)
	{
		locTOFSpacetimeHit->x = TOFGeom[0]->bar2y(bar,0);
		locTOFSpacetimeHit->y = 66.; //paddle extends from 6 -> 126
		locTOFSpacetimeHit->t = locTOFHit->t_south - 60.0/v; //60 is HALFPADDLE for single-ended bars
	}
	else
	{
		locTOFSpacetimeHit->x = TOFGeom[0]->bar2y(bar,1);
		locTOFSpacetimeHit->y = -66.; //paddle extends from -6 -> -126
		locTOFSpacetimeHit->t = locTOFHit->t_north - 60.0/v;
	}

	//printf("h: x %f y %f\n",locTOFSpacetimeHit->x,locTOFSpacetimeHit->y);
	return locTOFSpacetimeHit;
}

void DTOFPoint_factory::Create_TOFPoint(const tof_spacetimehit_t* locTOFSpacetimeHit_Horizontal, const tof_spacetimehit_t* locTOFSpacetimeHit_Vertical)
{
	const DTOFPaddleHit* locTOFHit_Horizontal = locTOFSpacetimeHit_Horizontal->TOFHit;
	const DTOFPaddleHit* locTOFHit_Vertical = locTOFSpacetimeHit_Vertical->TOFHit;

	//reconstruct TOF hit information, using information from the best bar: one or both of the bars may have a PMT signal below threshold
	float locMatchX, locMatchY, locMatchZ, locMatchdE, locMatchT;
	if(locTOFSpacetimeHit_Horizontal->dPositionWellDefinedFlag && locTOFSpacetimeHit_Vertical->dPositionWellDefinedFlag)
	{
		//both bars have both PMT signals above threshold
		locMatchX = locTOFSpacetimeHit_Horizontal->x;
		locMatchY = locTOFSpacetimeHit_Vertical->y;
		locMatchZ = TOFGeom[0]->CenterMPlane; //z: midpoint between tof planes
		locMatchT = 0.5*(locTOFSpacetimeHit_Horizontal->t + locTOFSpacetimeHit_Vertical->t);
		locMatchdE = 0.5*(locTOFHit_Horizontal->dE + locTOFHit_Vertical->dE);
	}
	else if(locTOFSpacetimeHit_Horizontal->dPositionWellDefinedFlag)
	{
		//the vertical position from the vertical paddle is not well defined
		locMatchX = locTOFSpacetimeHit_Horizontal->x;
		locMatchY = locTOFSpacetimeHit_Horizontal->y;
		locMatchT = locTOFSpacetimeHit_Horizontal->t;
		locMatchZ = TOFGeom[0]->CenterHPlane; //z: center of horizontal plane
		locMatchdE = locTOFHit_Horizontal->dE;
	}
	else
	{
		//the horizontal position from the horizontal paddle is not well defined
		locMatchX = locTOFSpacetimeHit_Vertical->x;
		locMatchY = locTOFSpacetimeHit_Vertical->y;
		locMatchT = locTOFSpacetimeHit_Vertical->t;
		locMatchZ = TOFGeom[0]->CenterVPlane; //z: center of vertical plane
		locMatchdE = locTOFHit_Vertical->dE;
	}

	DTOFPoint* locTOFPoint = new DTOFPoint;
	locTOFPoint->AddAssociatedObject(locTOFHit_Horizontal);
	locTOFPoint->AddAssociatedObject(locTOFHit_Vertical);
	locTOFPoint->pos.SetXYZ(locMatchX, locMatchY, locMatchZ);
	locTOFPoint->t = locMatchT;
	locTOFPoint->dE = locMatchdE;

	_data.push_back(locTOFPoint);
}


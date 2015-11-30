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
	if(locTOFSpacetimeHitMatch2.dBothPositionsWellDefinedFlag != locTOFSpacetimeHitMatch1.dBothPositionsWellDefinedFlag)
		return locTOFSpacetimeHitMatch1.dBothPositionsWellDefinedFlag; //one hit position is well defined and the other is not
	return (locTOFSpacetimeHitMatch1.delta_r < locTOFSpacetimeHitMatch2.delta_r);
};

//------------------
// brun
//------------------
jerror_t DTOFPoint_factory::brun(JEventLoop *loop, int32_t runnumber)
{

	map<string, double> tofparms;
 	if( !loop->GetCalib("TOF/tof_parms", tofparms))
	{
		//cout<<"DTOFPoint_factory: loading values from TOF data base"<<endl;
		HALFPADDLE = tofparms["TOF_HALFPADDLE"];
		E_THRESHOLD = tofparms["TOF_E_THRESHOLD"];
		ATTEN_LENGTH = tofparms["TOF_ATTEN_LENGTH"];
	}
	else
	{
		cout << "DTOFPoint_factory: Error loading values from TOF data base" <<endl;
		HALFPADDLE = 126; // set to some reasonable value
		E_THRESHOLD = 0.0005;
		ATTEN_LENGTH = 400.;
	}

	if(eventLoop->GetCalib("TOF/propagation_speed", propagation_speed))
		jout << "Error loading /TOF/propagation_speed !" << endl;

	loop->GetSingle(dTOFGeometry);

	HALFPADDLE_ONESIDED = dTOFGeometry->SHORTBARLENGTH/2.0; //GET FROM GEOMETRY??
	double locBeamHoleWidth = dTOFGeometry->LONGBARLENGTH - 2.0*dTOFGeometry->SHORTBARLENGTH;
	ONESIDED_PADDLE_MIDPOINT_MAG = HALFPADDLE_ONESIDED + locBeamHoleWidth/2.0;

	dPositionMatchCut_DoubleEnded = 9.0; //1.5*BARWIDTH
//	dTimeMatchCut_PositionWellDefined = 1.0;
	dTimeMatchCut_PositionWellDefined = 10.0;
	dTimeMatchCut_PositionNotWellDefined = 10.0;

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
jerror_t DTOFPoint_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
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

			tof_spacetimehitmatch_t locTOFSpacetimeHitMatch;
			if(!Match_Hits(locTOFSpacetimeHit_Horizontal, locTOFSpacetimeHit_Vertical, locTOFSpacetimeHitMatch))
				continue; //not a match

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

		Create_MatchedTOFPoint(locTOFSpacetimeHit_Horizontal, locTOFSpacetimeHit_Vertical);

		//remove used hits from the unused list
		locUnusedTOFSpacetimeHits.erase(locTOFSpacetimeHit_Horizontal);
		locUnusedTOFSpacetimeHits.erase(locTOFSpacetimeHit_Vertical);
	}

	// Loop over unused/unmatched TOF Spacetime hits, and create separate DTOFPoint's for them
	set<tof_spacetimehit_t*>::iterator locSetIterator = locUnusedTOFSpacetimeHits.begin();
	for(; locSetIterator != locUnusedTOFSpacetimeHits.end(); ++locSetIterator)
		Create_UnMatchedTOFPoint(*locSetIterator);

	return NOERROR;
}

DTOFPoint_factory::tof_spacetimehit_t* DTOFPoint_factory::Build_TOFSpacetimeHit_Horizontal(const DTOFPaddleHit* locTOFHit)
{
	tof_spacetimehit_t* locTOFSpacetimeHit = Get_TOFSpacetimeHitResource();
	locTOFSpacetimeHit->TOFHit = locTOFHit;

	int bar = locTOFHit->bar;
	int id = 44 + locTOFHit->bar - 1;
	double v = propagation_speed[id];

	if((locTOFHit->bar < dTOFGeometry->FirstShortBar) || (locTOFHit->bar > dTOFGeometry->LastShortBar)) //double-ended bars
	{
		locTOFSpacetimeHit->dIsDoubleEndedBar = true;
		locTOFSpacetimeHit->dIsSingleEndedNorthPaddle = false;
		locTOFSpacetimeHit->y = dTOFGeometry->bar2y(bar);
		if(locTOFHit->meantime != locTOFHit->meantime)
		{
			//NaN: only one energy hit above threshold on the double-ended bar
			locTOFSpacetimeHit->dPositionWellDefinedFlag = false;
			locTOFSpacetimeHit->x = 0.0;
			if(locTOFHit->E_north > E_THRESHOLD)
				locTOFSpacetimeHit->t = locTOFHit->t_north - HALFPADDLE/v;
			else
				locTOFSpacetimeHit->t = locTOFHit->t_south - HALFPADDLE/v;

			locTOFSpacetimeHit->pos_cut = 1000.0;
			locTOFSpacetimeHit->t_cut = dTimeMatchCut_PositionNotWellDefined;
		}
		else
		{
			locTOFSpacetimeHit->dPositionWellDefinedFlag = true;
			locTOFSpacetimeHit->x = locTOFHit->pos;
			locTOFSpacetimeHit->t = locTOFHit->meantime;
			locTOFSpacetimeHit->pos_cut = dPositionMatchCut_DoubleEnded;
			locTOFSpacetimeHit->t_cut = dTimeMatchCut_PositionWellDefined;
		}

		//printf("h: x %f y %f\n",locTOFSpacetimeHit->x,locTOFSpacetimeHit->y);
		return locTOFSpacetimeHit;
	}

	//single-ended bars
	locTOFSpacetimeHit->dIsDoubleEndedBar = false;
	locTOFSpacetimeHit->dPositionWellDefinedFlag = false;
	locTOFSpacetimeHit->pos_cut = 1000.0;
	locTOFSpacetimeHit->t_cut = dTimeMatchCut_PositionNotWellDefined;

	if(locTOFHit->t_south != 0.)
	{
		locTOFSpacetimeHit->dIsSingleEndedNorthPaddle = false;
		locTOFSpacetimeHit->y = dTOFGeometry->bar2y(bar,1);
		locTOFSpacetimeHit->x = -1.0*ONESIDED_PADDLE_MIDPOINT_MAG;
		locTOFSpacetimeHit->t = locTOFHit->t_south - HALFPADDLE_ONESIDED/v;
	}
	else
	{
		locTOFSpacetimeHit->dIsSingleEndedNorthPaddle = true;
		locTOFSpacetimeHit->y = dTOFGeometry->bar2y(bar,0);
		locTOFSpacetimeHit->x = ONESIDED_PADDLE_MIDPOINT_MAG;
		locTOFSpacetimeHit->t = locTOFHit->t_north - HALFPADDLE_ONESIDED/v;
	}

	//printf("h: x %f y %f\n",locTOFSpacetimeHit->x,locTOFSpacetimeHit->y);
	return locTOFSpacetimeHit;
}

DTOFPoint_factory::tof_spacetimehit_t* DTOFPoint_factory::Build_TOFSpacetimeHit_Vertical(const DTOFPaddleHit* locTOFHit)
{
	tof_spacetimehit_t* locTOFSpacetimeHit = Get_TOFSpacetimeHitResource();
	locTOFSpacetimeHit->TOFHit = locTOFHit;

	int bar = locTOFHit->bar;
	int id = locTOFHit->bar - 1;
	double v = propagation_speed[id];

	if((locTOFHit->bar < dTOFGeometry->FirstShortBar) || (locTOFHit->bar > dTOFGeometry->LastShortBar))
	{
		//double-ended bars
		locTOFSpacetimeHit->dIsDoubleEndedBar = true;
		locTOFSpacetimeHit->dIsSingleEndedNorthPaddle = false;
		locTOFSpacetimeHit->x = dTOFGeometry->bar2y(bar);
		if(locTOFHit->meantime != locTOFHit->meantime)
		{
			//NaN: only one energy hit above threshold on the double-ended bar
			locTOFSpacetimeHit->dPositionWellDefinedFlag = false;
			locTOFSpacetimeHit->y = 0.0;
			if(locTOFHit->E_north > E_THRESHOLD)
				locTOFSpacetimeHit->t = locTOFHit->t_north - HALFPADDLE/v;
			else
				locTOFSpacetimeHit->t = locTOFHit->t_south - HALFPADDLE/v;
			locTOFSpacetimeHit->pos_cut = 1000.0;
			locTOFSpacetimeHit->t_cut = dTimeMatchCut_PositionNotWellDefined;
		}
		else
		{
			locTOFSpacetimeHit->dPositionWellDefinedFlag = true;
			locTOFSpacetimeHit->y = locTOFHit->pos;
			locTOFSpacetimeHit->t = locTOFHit->meantime;
			locTOFSpacetimeHit->pos_cut = dPositionMatchCut_DoubleEnded;
			locTOFSpacetimeHit->t_cut = dTimeMatchCut_PositionWellDefined;
		}

		//printf("h: x %f y %f\n",locTOFSpacetimeHit->x,locTOFSpacetimeHit->y);
		return locTOFSpacetimeHit;
	}

	//single-ended bars
	locTOFSpacetimeHit->dIsDoubleEndedBar = false;
	locTOFSpacetimeHit->dPositionWellDefinedFlag = false;
	locTOFSpacetimeHit->pos_cut = 1000.0;
	locTOFSpacetimeHit->t_cut = dTimeMatchCut_PositionNotWellDefined;
	if(locTOFHit->t_south != 0.)
	{
		locTOFSpacetimeHit->dIsSingleEndedNorthPaddle = false;
		locTOFSpacetimeHit->x = dTOFGeometry->bar2y(bar,0);
		locTOFSpacetimeHit->y = -1.0*ONESIDED_PADDLE_MIDPOINT_MAG;
		locTOFSpacetimeHit->t = locTOFHit->t_south - HALFPADDLE_ONESIDED/v;
	}
	else
	{
		locTOFSpacetimeHit->dIsSingleEndedNorthPaddle = true;
		locTOFSpacetimeHit->x = dTOFGeometry->bar2y(bar,1);
		locTOFSpacetimeHit->y = ONESIDED_PADDLE_MIDPOINT_MAG;
		locTOFSpacetimeHit->t = locTOFHit->t_north - HALFPADDLE_ONESIDED/v;
	}

	//printf("h: x %f y %f\n",locTOFSpacetimeHit->x,locTOFSpacetimeHit->y);
	return locTOFSpacetimeHit;
}

bool DTOFPoint_factory::Match_Hits(tof_spacetimehit_t* locTOFSpacetimeHit_Horizontal, tof_spacetimehit_t* locTOFSpacetimeHit_Vertical, tof_spacetimehitmatch_t& locTOFSpacetimeHitMatch)
{
	//make sure that single-ended paddles don't match each other
	if((!locTOFSpacetimeHit_Vertical->dIsDoubleEndedBar) && (!locTOFSpacetimeHit_Horizontal->dIsDoubleEndedBar))
		return false; //unphysical

	//make sure that (e.g. horizontal) single-ended paddles don't match (e.g. vertical) paddles on the opposite side
	if(!locTOFSpacetimeHit_Horizontal->dIsDoubleEndedBar)
	{
		//horizontal is single-ended
		if(locTOFSpacetimeHit_Horizontal->dIsSingleEndedNorthPaddle)
		{
			//horizontal is on north (+x) side
			if(locTOFSpacetimeHit_Vertical->TOFHit->bar < dTOFGeometry->FirstShortBar)
				return false; //vertical is on south (-x) side: CANNOT MATCH
		}
		else
		{
			//horizontal is on south (-x) side
			if(locTOFSpacetimeHit_Vertical->TOFHit->bar > dTOFGeometry->LastShortBar)
				return false; //vertical is on north (+x) side: CANNOT MATCH
		}
	}
	else if(!locTOFSpacetimeHit_Vertical->dIsDoubleEndedBar)
	{
		//vertical is single-ended
		if(locTOFSpacetimeHit_Vertical->dIsSingleEndedNorthPaddle)
		{
			//vertical is on north (+y) side
			if(locTOFSpacetimeHit_Horizontal->TOFHit->bar < dTOFGeometry->FirstShortBar)
				return false; //horizontal is on south (-y) side: CANNOT MATCH
		}
		else
		{
			//vertical is on south (-y) side
			if(locTOFSpacetimeHit_Horizontal->TOFHit->bar > dTOFGeometry->LastShortBar)
				return false; //horizontal is on north (+y) side: CANNOT MATCH
		}
	}

	//If the position along BOTH paddles is not well defined, cannot tell whether these hits match or not
		//If the hit multiplicity was low, could just generate all matches
		//However, the hit multiplicity is generally very high due to hits near the beamline: would have TOF hits everywhere
		//Therefore, don't keep as match: register separately, let track / tof matching salvage the situation
	if((!locTOFSpacetimeHit_Horizontal->dPositionWellDefinedFlag) && (!locTOFSpacetimeHit_Vertical->dPositionWellDefinedFlag))
		return false; //have no idea whether these hits go together: assume they don't

	float locDeltaX = locTOFSpacetimeHit_Horizontal->x - locTOFSpacetimeHit_Vertical->x;
	if(fabs(locDeltaX) > locTOFSpacetimeHit_Horizontal->pos_cut)
		return false;

	float locDeltaY = locTOFSpacetimeHit_Horizontal->y - locTOFSpacetimeHit_Vertical->y;
	if(fabs(locDeltaY) > locTOFSpacetimeHit_Vertical->pos_cut)
		return false;

	float locDeltaT = locTOFSpacetimeHit_Horizontal->t - locTOFSpacetimeHit_Vertical->t;
	float locTimeCut = (locTOFSpacetimeHit_Horizontal->t_cut > locTOFSpacetimeHit_Vertical->t_cut) ? locTOFSpacetimeHit_Horizontal->t_cut : locTOFSpacetimeHit_Vertical->t_cut;
	if(fabs(locDeltaT) > locTimeCut)
		return false;

	locTOFSpacetimeHitMatch.delta_t = locDeltaT;
	locTOFSpacetimeHitMatch.delta_r = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
	locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Horizontal = locTOFSpacetimeHit_Horizontal;
	locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Vertical = locTOFSpacetimeHit_Vertical;
	locTOFSpacetimeHitMatch.dBothPositionsWellDefinedFlag = (locTOFSpacetimeHit_Horizontal->dPositionWellDefinedFlag == locTOFSpacetimeHit_Vertical->dPositionWellDefinedFlag);

	return true;
}

void DTOFPoint_factory::Create_MatchedTOFPoint(const tof_spacetimehit_t* locTOFSpacetimeHit_Horizontal, const tof_spacetimehit_t* locTOFSpacetimeHit_Vertical)
{
	const DTOFPaddleHit* locTOFHit_Horizontal = locTOFSpacetimeHit_Horizontal->TOFHit;
	const DTOFPaddleHit* locTOFHit_Vertical = locTOFSpacetimeHit_Vertical->TOFHit;

	//reconstruct TOF hit information, using information from the best bar: one or both of the bars may have a PMT signal below threshold
	float locMatchX, locMatchY, locMatchZ, locMatchdE, locMatchT;
	if(locTOFSpacetimeHit_Horizontal->dPositionWellDefinedFlag && locTOFSpacetimeHit_Vertical->dPositionWellDefinedFlag)
	{
		//both bars have both PMT signals above threshold
			//is x/y resolution from energy calibration better than x/y resolution from paddle edges?
		locMatchX = locTOFSpacetimeHit_Horizontal->x;
		locMatchY = locTOFSpacetimeHit_Vertical->y;
		locMatchZ = dTOFGeometry->CenterMPlane; //z: midpoint between tof planes
		locMatchT = 0.5*(locTOFSpacetimeHit_Horizontal->t + locTOFSpacetimeHit_Vertical->t);
		locMatchdE = 0.5*(locTOFHit_Horizontal->dE + locTOFHit_Vertical->dE);
	}
	else if(locTOFSpacetimeHit_Horizontal->dPositionWellDefinedFlag)
	{
		//the vertical position from the vertical paddle is not well defined
		locMatchX = locTOFSpacetimeHit_Horizontal->x;
		locMatchY = locTOFSpacetimeHit_Horizontal->y;
		locMatchT = locTOFSpacetimeHit_Horizontal->t;
		locMatchZ = dTOFGeometry->CenterHPlane; //z: center of horizontal plane
		locMatchdE = locTOFHit_Horizontal->dE;
	}
	else
	{
		//the horizontal position from the horizontal paddle is not well defined
		locMatchX = locTOFSpacetimeHit_Vertical->x;
		locMatchY = locTOFSpacetimeHit_Vertical->y;
		locMatchT = locTOFSpacetimeHit_Vertical->t;
		locMatchZ = dTOFGeometry->CenterVPlane; //z: center of vertical plane
		locMatchdE = locTOFHit_Vertical->dE;
	}

	DTOFPoint* locTOFPoint = new DTOFPoint;
	locTOFPoint->AddAssociatedObject(locTOFHit_Horizontal);
	locTOFPoint->AddAssociatedObject(locTOFHit_Vertical);
	locTOFPoint->pos.SetXYZ(locMatchX, locMatchY, locMatchZ);
	locTOFPoint->t = locMatchT;
	locTOFPoint->tErr = 0.0; //SET ME
	locTOFPoint->dE = locMatchdE;

	locTOFPoint->dHorizontalBar = locTOFHit_Horizontal->bar;
	locTOFPoint->dVerticalBar = locTOFHit_Vertical->bar;

	//Status: 0 if no hit (or none above threshold), 1 if only North hit above threshold, 2 if only South hit above threshold, 3 if both hits above threshold
	locTOFPoint->dHorizontalBarStatus = int(locTOFHit_Horizontal->E_north > E_THRESHOLD) + 2*int(locTOFHit_Horizontal->E_south > E_THRESHOLD);
	locTOFPoint->dVerticalBarStatus = int(locTOFHit_Vertical->E_north > E_THRESHOLD) + 2*int(locTOFHit_Vertical->E_south > E_THRESHOLD);

	_data.push_back(locTOFPoint);
}

void DTOFPoint_factory::Create_UnMatchedTOFPoint(const tof_spacetimehit_t* locTOFSpacetimeHit)
{
	const DTOFPaddleHit* locPaddleHit = locTOFSpacetimeHit->TOFHit;
	bool locIsHorizontalBarFlag = (locPaddleHit->orientation == 1);
	float locPointZ = locIsHorizontalBarFlag ? dTOFGeometry->CenterHPlane : dTOFGeometry->CenterVPlane;
	if(locTOFSpacetimeHit->dPositionWellDefinedFlag)
	{
		//Position is well defined
		DTOFPoint* locTOFPoint = new DTOFPoint;
		locTOFPoint->AddAssociatedObject(locPaddleHit);

		locTOFPoint->pos.SetXYZ(locTOFSpacetimeHit->x, locTOFSpacetimeHit->y, locPointZ);
		locTOFPoint->t = locTOFSpacetimeHit->t;
		locTOFPoint->tErr = 0.0; //SET ME
		locTOFPoint->dE = locPaddleHit->dE;

		locTOFPoint->dHorizontalBar = locIsHorizontalBarFlag ? locPaddleHit->bar : 0;
		locTOFPoint->dVerticalBar = locIsHorizontalBarFlag ? 0 : locPaddleHit->bar;

		//Status: 0 if no hit (or none above threshold), 1 if only North hit above threshold, 2 if only South hit above threshold, 3 if both hits above threshold
		locTOFPoint->dHorizontalBarStatus = locIsHorizontalBarFlag ? 3 : 0;
		locTOFPoint->dVerticalBarStatus = locIsHorizontalBarFlag ? 0 : 3;

		_data.push_back(locTOFPoint);
	}
	else
	{
		//position not well defined: save anyway:
			//Will use track matching to define position in the other direction
			//Then, will update the hit energy and time based on that position

		DTOFPoint* locTOFPoint = new DTOFPoint;
		locTOFPoint->AddAssociatedObject(locPaddleHit);

		float locPointX = locTOFSpacetimeHit->x;
		float locPointY = locTOFSpacetimeHit->y;
		locTOFPoint->pos.SetXYZ(locPointX, locPointY, locPointZ);
		locTOFPoint->t = locTOFSpacetimeHit->t;
		locTOFPoint->tErr = 0.0; //SET ME

		bool locNorthAboveThresholdFlag = (locPaddleHit->E_north > E_THRESHOLD);

		locTOFPoint->dHorizontalBar = locIsHorizontalBarFlag ? locPaddleHit->bar : 0;
		locTOFPoint->dVerticalBar = locIsHorizontalBarFlag ? 0 : locPaddleHit->bar;

		int locBarStatus = locNorthAboveThresholdFlag ? 1 : 2;
		locTOFPoint->dHorizontalBarStatus = locIsHorizontalBarFlag ? locBarStatus : 0;
		locTOFPoint->dVerticalBarStatus = locIsHorizontalBarFlag ? 0 : locBarStatus;

		//Energy: Propagate to paddle mid-point
		double locDeltaXToMidPoint = locTOFSpacetimeHit->dIsDoubleEndedBar ? HALFPADDLE : HALFPADDLE_ONESIDED;
		float locEnergy = locNorthAboveThresholdFlag ? locPaddleHit->E_north : locPaddleHit->E_south;
		locEnergy *= exp(locDeltaXToMidPoint/ATTEN_LENGTH);
		locTOFPoint->dE = locEnergy;

		_data.push_back(locTOFPoint);
	}
}

jerror_t DTOFPoint_factory::fini(void)
{
	for(size_t loc_i = 0; loc_i < dTOFSpacetimeHitPool_All.size(); ++loc_i)
		delete dTOFSpacetimeHitPool_All[loc_i];
	return NOERROR;
}

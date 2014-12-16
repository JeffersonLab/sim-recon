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

#define MAXTOFHITS 50

bool Compare_TOFSpacetimeHitMatches_Distance(const DTOFPoint_factory::tof_spacetimehitmatch_t& locTOFSpacetimeHitMatch1, const DTOFPoint_factory::tof_spacetimehitmatch_t& locTOFSpacetimeHitMatch2)
{
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

//------------------
// evnt
//------------------
jerror_t DTOFPoint_factory::evnt(JEventLoop *loop, int eventnumber)
{	
	vector<const DTOFPaddleHit*> locTOFHitVector;
	loop->Get(locTOFHitVector);

	// create the hit spacetime information
	deque<tof_spacetimehit_t> locTOFSpacetimeHits_Horizontal, locTOFSpacetimeHits_Vertical;
	for(size_t loc_i = 0; loc_i < locTOFHitVector.size(); ++loc_i)
	{
		const DTOFPaddleHit* locTOFHit = locTOFHitVector[loc_i];
		if(!((locTOFHit->E_north > E_THRESHOLD) || (locTOFHit->E_south > E_THRESHOLD)))
			continue;

		if(locTOFHit->orientation) //horizontal
			locTOFSpacetimeHits_Horizontal.push_back(Build_TOFSpacetimeHit_Horizontal(locTOFHit));
		else //vertical
			locTOFSpacetimeHits_Vertical.push_back(Build_TOFSpacetimeHit_Vertical(locTOFHit));
	}
	

	//find matches between planes and sort them by delta-r
	deque<tof_spacetimehitmatch_t> locTOFSpacetimeHitMatches; //use list for sorting, vector for resource pool
	for(size_t loc_i = 0; loc_i < locTOFSpacetimeHits_Horizontal.size(); ++loc_i)
	{
		const tof_spacetimehit_t* locTOFSpacetimeHit_Horizontal = &(locTOFSpacetimeHits_Horizontal[loc_i]);
		for(size_t loc_j = 0; loc_j < locTOFSpacetimeHits_Vertical.size(); ++loc_j)
		{
			const tof_spacetimehit_t* locTOFSpacetimeHit_Vertical = &(locTOFSpacetimeHits_Vertical[loc_j]);
			
			float locDeltaX = locTOFSpacetimeHit_Horizontal->x - locTOFSpacetimeHit_Vertical->x;
			if(fabs(locDeltaX) >= locTOFSpacetimeHit_Horizontal->pos_cut)
				continue;
			float locDeltaY = locTOFSpacetimeHit_Horizontal->y - locTOFSpacetimeHit_Vertical->y;
			if(fabs(locDeltaY) >= locTOFSpacetimeHit_Vertical->pos_cut)
				continue;

			float locDeltaT = locTOFSpacetimeHit_Horizontal->t - locTOFSpacetimeHit_Vertical->t;
			float locTimeCut = (locTOFSpacetimeHit_Horizontal->t_cut > locTOFSpacetimeHit_Vertical->t_cut) ? locTOFSpacetimeHit_Horizontal->t_cut : locTOFSpacetimeHit_Vertical->t_cut;
			if(fabs(locDeltaT) >= locTimeCut)
				continue;
			
			tof_spacetimehitmatch_t locTOFSpacetimeHitMatch;
			locTOFSpacetimeHitMatch.delta_t = locDeltaT;
			locTOFSpacetimeHitMatch.delta_r = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
			locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Horizontal = locTOFSpacetimeHit_Horizontal;
			locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Vertical = locTOFSpacetimeHit_Vertical;
			locTOFSpacetimeHitMatches.push_back(locTOFSpacetimeHitMatch);
		}
	}
	std::sort(locTOFSpacetimeHitMatches.begin(), locTOFSpacetimeHitMatches.end(), Compare_TOFSpacetimeHitMatches_Distance); //sort matches by delta_r


	// create DTOFPoints, in order of best matches (by delta_r)
	while(!locTOFSpacetimeHitMatches.empty())
	{
		tof_spacetimehitmatch_t& locTOFSpacetimeHitMatch = *(locTOFSpacetimeHitMatches.begin());
		const tof_spacetimehit_t* locTOFSpacetimeHit_Horizontal = locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Horizontal;
		const tof_spacetimehit_t* locTOFSpacetimeHit_Vertical = locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Vertical;
		
		if(!Create_TOFPoint(locTOFSpacetimeHit_Horizontal, locTOFSpacetimeHit_Vertical))
		{
			//match rejected
			locTOFSpacetimeHitMatches.erase(locTOFSpacetimeHitMatches.begin());
			continue;
		}
		
		//remove matches that contain the hits that were just used
		deque<tof_spacetimehitmatch_t>::iterator locIterator = locTOFSpacetimeHitMatches.begin();
		while(locIterator != locTOFSpacetimeHitMatches.end())
		{
			tof_spacetimehitmatch_t& locTOFSpacetimeHitMatch = *locIterator;

			if(locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Horizontal == locTOFSpacetimeHit_Horizontal)
				locIterator = locTOFSpacetimeHitMatches.erase(locIterator); //erase advances the iterator
			else if(locTOFSpacetimeHitMatch.dTOFSpacetimeHit_Vertical == locTOFSpacetimeHit_Vertical)
				locIterator = locTOFSpacetimeHitMatches.erase(locIterator); //erase advances the iterator
			else
				++locIterator;
		}
	}
	
	return NOERROR;
}

DTOFPoint_factory::tof_spacetimehit_t DTOFPoint_factory::Build_TOFSpacetimeHit_Horizontal(const DTOFPaddleHit* locTOFHit)
{
	tof_spacetimehit_t locTOFSpacetimeHit;
	locTOFSpacetimeHit.TOFHit = locTOFHit;

	int bar = locTOFHit->bar;
	int id = 44*locTOFHit->orientation + locTOFHit->bar - 1;
	double v = propagation_speed[id];

	if((locTOFHit->bar < TOFGeom[0]->FirstShortBar) || (locTOFHit->bar > TOFGeom[0]->LastShortBar)) //double-ended bars
	{
		locTOFSpacetimeHit.y = TOFGeom[0]->bar2y(bar);
		if(locTOFHit->meantime != locTOFHit->meantime)
		{
			//NaN: only one energy hit above threshold on the double-ended bar
			locTOFSpacetimeHit.x = 0.0;
			if(locTOFHit->E_north > E_THRESHOLD)
				locTOFSpacetimeHit.t = locTOFHit->t_north - HALFPADDLE/v;
			else
				locTOFSpacetimeHit.t = locTOFHit->t_south - HALFPADDLE/v;

			locTOFSpacetimeHit.pos_cut = 252.0;
			locTOFSpacetimeHit.t_cut = 10.0;
		}
		else
		{
			locTOFSpacetimeHit.x = locTOFHit->pos;
			locTOFSpacetimeHit.t = locTOFHit->meantime;
			locTOFSpacetimeHit.pos_cut = dPositionMatchCut_DoubleEnded;
			// locTOFSpacetimeHit.t_cut = 1.0;
			locTOFSpacetimeHit.t_cut=10.0;
		}

		//printf("h: x %f y %f\n",locTOFSpacetimeHit.x,locTOFSpacetimeHit.y);
		return locTOFSpacetimeHit;
	}

	//single-ended bars
	locTOFSpacetimeHit.pos_cut = 120.;
	locTOFSpacetimeHit.t_cut = 10.;

	if(locTOFHit->t_south != 0.)
	{
		locTOFSpacetimeHit.y = TOFGeom[0]->bar2y(bar,1);
		locTOFSpacetimeHit.x = -66.; //paddle extends from -6 -> -126
		locTOFSpacetimeHit.t = locTOFHit->t_south - 60.0/v; //60 is HALFPADDLE for single-ended bars
	}
	else
	{
		locTOFSpacetimeHit.y = TOFGeom[0]->bar2y(bar,0);
		locTOFSpacetimeHit.x = 66.; //paddle extends from 6 -> 126
		locTOFSpacetimeHit.t = locTOFHit->t_north - 60.0/v;
	}

	//printf("h: x %f y %f\n",locTOFSpacetimeHit.x,locTOFSpacetimeHit.y);
	return locTOFSpacetimeHit;
}

DTOFPoint_factory::tof_spacetimehit_t DTOFPoint_factory::Build_TOFSpacetimeHit_Vertical(const DTOFPaddleHit* locTOFHit)
{
	tof_spacetimehit_t locTOFSpacetimeHit;
	locTOFSpacetimeHit.TOFHit = locTOFHit;

	int bar = locTOFHit->bar;
	int id = 44*locTOFHit->orientation + locTOFHit->bar - 1;
	double v = propagation_speed[id];

	if((locTOFHit->bar < TOFGeom[0]->FirstShortBar) || (locTOFHit->bar > TOFGeom[0]->LastShortBar))
	{
		//double-ended bars
		locTOFSpacetimeHit.x = TOFGeom[0]->bar2y(bar);
		if(locTOFHit->meantime != locTOFHit->meantime)
		{
			//NaN: only one energy hit above threshold on the double-ended bar
			locTOFSpacetimeHit.y = 0.0;
			if(locTOFHit->E_north > E_THRESHOLD)
				locTOFSpacetimeHit.t = locTOFHit->t_north - HALFPADDLE/v;
			else
				locTOFSpacetimeHit.t = locTOFHit->t_south - HALFPADDLE/v;
			locTOFSpacetimeHit.pos_cut = 252.0;
			locTOFSpacetimeHit.t_cut = 10.0;
		}
		else
		{
			locTOFSpacetimeHit.y = locTOFHit->pos;
			locTOFSpacetimeHit.t = locTOFHit->meantime;
			locTOFSpacetimeHit.pos_cut = dPositionMatchCut_DoubleEnded;
			//	locTOFSpacetimeHit.t_cut = 1.0;
			locTOFSpacetimeHit.t_cut=10.0;
		}

		//printf("h: x %f y %f\n",locTOFSpacetimeHit.x,locTOFSpacetimeHit.y);
		return locTOFSpacetimeHit;
	}

	//single-ended bars
	locTOFSpacetimeHit.pos_cut = 120.;
	locTOFSpacetimeHit.t_cut = 10.;
	if(locTOFHit->t_south != 0.)
	{
		locTOFSpacetimeHit.x = TOFGeom[0]->bar2y(bar,0);
		locTOFSpacetimeHit.y = 66.; //paddle extends from 6 -> 126
		locTOFSpacetimeHit.t = locTOFHit->t_south - 60.0/v; //60 is HALFPADDLE for single-ended bars
	}
	else
	{
		locTOFSpacetimeHit.x = TOFGeom[0]->bar2y(bar,1);
		locTOFSpacetimeHit.y = -66.; //paddle extends from -6 -> -126
		locTOFSpacetimeHit.t = locTOFHit->t_north - 60.0/v;
	}

	//printf("h: x %f y %f\n",locTOFSpacetimeHit.x,locTOFSpacetimeHit.y);
	return locTOFSpacetimeHit;
}

bool DTOFPoint_factory::Create_TOFPoint(const tof_spacetimehit_t* locTOFSpacetimeHit_Horizontal, const tof_spacetimehit_t* locTOFSpacetimeHit_Vertical)
{
	const DTOFPaddleHit* locTOFHit_Horizontal = locTOFSpacetimeHit_Horizontal->TOFHit;
	const DTOFPaddleHit* locTOFHit_Vertical = locTOFSpacetimeHit_Vertical->TOFHit;

	//see which bars had both PMT signals above threshold
	bool locOneSideBelowEThresholdFlag_Horizontal = false;
	if(locTOFHit_Horizontal->meantime != locTOFHit_Horizontal->meantime)
		locOneSideBelowEThresholdFlag_Horizontal = true; //NaN: only one PMT signal above threshold

	bool locOneSideBelowEThresholdFlag_Vertical = false;
	if(locTOFHit_Vertical->meantime != locTOFHit_Vertical->meantime)
		locOneSideBelowEThresholdFlag_Vertical = true; //NaN: only one PMT signal above threshold

	//Disable matching for certain scenarios until hit differentiation is improved
		//e.g. separate DTOFPoints due to energy deposited in adjacent paddles, PID routine for picking the correct DTOFPoint if adjacent, etc.
	if(locOneSideBelowEThresholdFlag_Horizontal && locOneSideBelowEThresholdFlag_Vertical)
		return false; //both bars have only one PMT signal above threshold
	else if(locOneSideBelowEThresholdFlag_Horizontal && (locTOFHit_Horizontal->bar < 22 || locTOFHit_Horizontal->bar>23))
		return false; //a horizontal double-sided bar has only one PMT signal above threshold
	else if(locOneSideBelowEThresholdFlag_Vertical && (locTOFHit_Vertical->bar < 22 || locTOFHit_Vertical->bar>23))
		return false; //a vertical double-sided bar has only one PMT signal above threshold

	//reconstruct TOF hit information, using information from the best bar: one or both of the bars may have a PMT signal below threshold
	float locMatchX, locMatchY, locMatchZ, locMatchdE, locMatchT;
	if(locOneSideBelowEThresholdFlag_Horizontal && locOneSideBelowEThresholdFlag_Vertical) //both bars have only one PMT signal above threshold
	{
		locMatchX = locTOFSpacetimeHit_Vertical->x;
		locMatchY = locTOFSpacetimeHit_Horizontal->y;
		locMatchZ = TOFGeom[0]->CenterMPlane; //z: midpoint between tof planes

		float locDepositedEnergy_Horizontal = 0.0;
		if(locTOFHit_Horizontal->E_north > E_THRESHOLD)
			locDepositedEnergy_Horizontal = locTOFHit_Horizontal->E_north * exp((HALFPADDLE - locMatchX)/ATTEN_LENGTH);
		else
			locDepositedEnergy_Horizontal = locTOFHit_Horizontal->E_south * exp((HALFPADDLE + locMatchX)/ATTEN_LENGTH);

		float locDepositedEnergy_Vertical = 0.0;
		if(locTOFHit_Vertical->E_north > E_THRESHOLD)
			locDepositedEnergy_Vertical = locTOFHit_Vertical->E_north * exp((HALFPADDLE - locMatchY)/ATTEN_LENGTH);
		else
			locDepositedEnergy_Vertical = locTOFHit_Vertical->E_south * exp((HALFPADDLE + locMatchY)/ATTEN_LENGTH);

		locMatchdE = 0.5*(locDepositedEnergy_Horizontal + locDepositedEnergy_Vertical);
		locMatchT = 0.5*(locTOFSpacetimeHit_Horizontal->t + locTOFSpacetimeHit_Vertical->t);
	}
	else if(locOneSideBelowEThresholdFlag_Horizontal) //the horizontal bar has only one PMT signal above threshold (and the vertical bar has both)
	{
		locMatchX = locTOFSpacetimeHit_Vertical->x;
		locMatchY = locTOFSpacetimeHit_Vertical->y;
		locMatchT = locTOFSpacetimeHit_Vertical->t;
		locMatchZ = TOFGeom[0]->CenterVPlane; //z: center of vertical plane
		locMatchdE = locTOFHit_Vertical->dE;
	}
	else if(locOneSideBelowEThresholdFlag_Vertical) //the vertical bar has only one PMT signal above threshold (and the horizontal bar has both)
	{
		locMatchX = locTOFSpacetimeHit_Horizontal->x;
		locMatchY = locTOFSpacetimeHit_Horizontal->y;
		locMatchT = locTOFSpacetimeHit_Horizontal->t;
		locMatchZ = TOFGeom[0]->CenterHPlane; //z: center of horizontal plane
		locMatchdE = locTOFHit_Horizontal->dE;
	}
	else //both bars have both PMT signals above threshold
	{
		locMatchX = locTOFSpacetimeHit_Horizontal->x;
		locMatchY = locTOFSpacetimeHit_Vertical->y;
		locMatchZ = TOFGeom[0]->CenterMPlane; //z: midpoint between tof planes
		locMatchT = 0.5*(locTOFSpacetimeHit_Horizontal->t + locTOFSpacetimeHit_Vertical->t);
		locMatchdE = 0.5*(locTOFHit_Horizontal->dE + locTOFHit_Vertical->dE);
	}

	DTOFPoint* locTOFPoint = new DTOFPoint;
	locTOFPoint->AddAssociatedObject(locTOFHit_Horizontal);
	locTOFPoint->AddAssociatedObject(locTOFHit_Vertical);
	locTOFPoint->pos.SetXYZ(locMatchX, locMatchY, locMatchZ);
	locTOFPoint->t = locMatchT;
	locTOFPoint->dE = locMatchdE;

	_data.push_back(locTOFPoint);

	return true;
}


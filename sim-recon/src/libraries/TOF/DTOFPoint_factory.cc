// $Id$
//
//    File: DTOFPoint_factory.cc
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>
#include <math.h>
using namespace std;

#include "DTOFPoint_factory.h"

#define MAXTOFHITS 50

bool Compare_TOFSpacetimeHitMatches_Distance(DTOFPoint_factory::tof_spacetimehitmatch_t *locTOFSpacetimeHitMatch1, DTOFPoint_factory::tof_spacetimehitmatch_t *locTOFSpacetimeHitMatch2){
	return (locTOFSpacetimeHitMatch1->delta_r < locTOFSpacetimeHitMatch2->delta_r);
};

//------------------
// brun
//------------------
jerror_t DTOFPoint_factory::brun(JEventLoop *loop, int runnumber)
{

  map<string, double> tofparms;
 
  if ( !loop->GetCalib("TOF/tof_parms", tofparms)){
    cout<<"DTOFPoint_factory: loading values from TOF data base"<<endl;
  } else {
    cout << "DTOFPoint_factory: Error loading values from TOF data base" <<endl;

    VELOCITY = 15.;  // set to some reasonable value
    HALFPADDLE = 126;   // set to some reasonable value
    BARWIDTH = 6.;
    return NOERROR;
  }

  VELOCITY     =    tofparms["TOF_C_EFFECTIVE"];
  HALFPADDLE   =    tofparms["TOF_HALFPADDLE"];
  BARWIDTH     =    tofparms["TOF_PADDLEWIDTH"];
  E_THRESHOLD  =    tofparms["TOF_E_THRESHOLD"];
  ATTEN_LENGTH =    tofparms["TOF_ATTEN_LENGTH"];

	MAX_TOFSpacetimeHits = 20;
	MAX_TOFSpacetimeHitMatches = 10;
	dPositionMatchCut_DoubleEnded = 0.5*BARWIDTH + 1.05*3.0; //max if perfect precision (1/2 bar width) + ~3 sigma
  return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t DTOFPoint_factory::evnt(JEventLoop *loop, int eventnumber)
{
	unsigned int loc_i, loc_j;
	tof_spacetimehit_t *locTOFSpacetimeHit;
	tof_spacetimehit_t *locTOFSpacetimeHit_Horizontal;
	tof_spacetimehit_t *locTOFSpacetimeHit_Vertical;
	const DTOFHit *locTOFHit, *locTOFHit_Horizontal, *locTOFHit_Vertical;
	float locTimeCut, locDeltaX, locDeltaY, locDeltaT;
	DTOFPoint *locTOFPoint;
	float locMatchX, locMatchY, locMatchZ, locMatchdE, locMatchT;

	vector<const DTOFHit*> locTOFHitVector;
	loop->Get(locTOFHitVector);
	deque<tof_spacetimehit_t*> locTOFSpacetimeHits_Horizontal;
	deque<tof_spacetimehit_t*> locTOFSpacetimeHits_Vertical;

	//increase pool size if necessary (used to minimize memory allocation time)
	for(loc_i = dTOFSpacetimeHitPool.size(); loc_i < locTOFHitVector.size(); loc_i++){
		locTOFSpacetimeHit = new tof_spacetimehit_t();
		dTOFSpacetimeHitPool.push_back(locTOFSpacetimeHit);
	}
	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if((locTOFHitVector.size() < MAX_TOFSpacetimeHits) && (dTOFSpacetimeHitPool.size() > MAX_TOFSpacetimeHits)){
		for(loc_i = MAX_TOFSpacetimeHits; loc_i < dTOFSpacetimeHitPool.size(); loc_i++) delete dTOFSpacetimeHitPool[loc_i];
		dTOFSpacetimeHitPool.resize(MAX_TOFSpacetimeHits);
	}

	// create the hit spacetime information, fill the vectors
	for (loc_i = 0; loc_i < locTOFHitVector.size(); loc_i++){
		locTOFHit = locTOFHitVector[loc_i];
		if (!((locTOFHit->E_north > E_THRESHOLD) || (locTOFHit->E_south > E_THRESHOLD)))
			continue; //

		locTOFSpacetimeHit = dTOFSpacetimeHitPool[loc_i];
		if (locTOFHit->orientation) { //horizontal
			if (locTOFHit->bar <= 40){ //double-ended bars
				locTOFSpacetimeHit->y = BARWIDTH*(locTOFHit->bar - 20.5) + ((locTOFHit->bar > 20) ? BARWIDTH : -1.0*BARWIDTH);
				if(!((locTOFHit->meantime >= 0.0) || (locTOFHit->meantime <= 0.0))){ //NaN: only one energy hit above threshold on the double-ended bar
					locTOFSpacetimeHit->x = 0.0;
					if (locTOFHit->E_north > E_THRESHOLD)
						locTOFSpacetimeHit->t = locTOFHit->t_north - HALFPADDLE/VELOCITY;
					else
						locTOFSpacetimeHit->t = locTOFHit->t_south - HALFPADDLE/VELOCITY;
					locTOFSpacetimeHit->pos_cut = 252.0;
					locTOFSpacetimeHit->t_cut = 10.0;
				}else{
					locTOFSpacetimeHit->x = locTOFHit->pos;
					locTOFSpacetimeHit->t = locTOFHit->meantime;
					locTOFSpacetimeHit->pos_cut = dPositionMatchCut_DoubleEnded;
					locTOFSpacetimeHit->t_cut = 1.0;
				}
			}else{ //single-ended bars
				locTOFSpacetimeHit->pos_cut = 120.;
				locTOFSpacetimeHit->t_cut = 10.;
				switch(locTOFHit->bar){
					case 41:
						locTOFSpacetimeHit->x = -66.; //paddle extends from -6 -> -126
						locTOFSpacetimeHit->y = -3.;
						locTOFSpacetimeHit->t = locTOFHit->t_south - 60.0/VELOCITY; //60 is HALFPADDLE for signle-ended bars
						break;
					case 43:
						locTOFSpacetimeHit->x = -66.; //paddle extends from -6 -> -126
						locTOFSpacetimeHit->y = 3.;
						locTOFSpacetimeHit->t = locTOFHit->t_south - 60.0/VELOCITY;
						break;
					case 42:
						locTOFSpacetimeHit->x = 66.; //paddle extends from 6 -> 126
						locTOFSpacetimeHit->y = -3.;
						locTOFSpacetimeHit->t = locTOFHit->t_north - 60.0/VELOCITY;
						break;
					case 44:
						locTOFSpacetimeHit->x = 66.; //paddle extends from 6 -> 126
						locTOFSpacetimeHit->y = 3.;
						locTOFSpacetimeHit->t = locTOFHit->t_north - 60.0/VELOCITY;
						break;
					default:
						break;
				}
			}
			locTOFSpacetimeHit->TOFHit = locTOFHit;
			locTOFSpacetimeHits_Horizontal.push_back(locTOFSpacetimeHit);
		} else {  //vertical

			if (locTOFHit->bar <= 40){ //double-ended bars
				locTOFSpacetimeHit->x = BARWIDTH*(locTOFHit->bar - 20.5) + ((locTOFHit->bar > 20) ? BARWIDTH : -1.0*BARWIDTH);
				if(!((locTOFHit->meantime >= 0.0) || (locTOFHit->meantime <= 0.0))){ //NaN: only one energy hit above threshold on the double-ended bar
					locTOFSpacetimeHit->y = 0.0;
					if (locTOFHit->E_north > E_THRESHOLD)
						locTOFSpacetimeHit->t = locTOFHit->t_north - HALFPADDLE/VELOCITY;
					else
						locTOFSpacetimeHit->t = locTOFHit->t_south - HALFPADDLE/VELOCITY;
					locTOFSpacetimeHit->pos_cut = 252.0;
					locTOFSpacetimeHit->t_cut = 10.0;
				}else{
					locTOFSpacetimeHit->y = locTOFHit->pos;
					locTOFSpacetimeHit->t = locTOFHit->meantime;
					locTOFSpacetimeHit->pos_cut = dPositionMatchCut_DoubleEnded;
					locTOFSpacetimeHit->t_cut = 1.0;
				}
			} else { //single-ended bars
				locTOFSpacetimeHit->pos_cut = 120.;
				locTOFSpacetimeHit->t_cut = 10.;
				switch(locTOFHit->bar){
					case 41:
						locTOFSpacetimeHit->x = -3.;
						locTOFSpacetimeHit->y = 66.; //paddle extends from 6 -> 126
						locTOFSpacetimeHit->t = locTOFHit->t_south - 60.0/VELOCITY;
						break;
					case 43:
						locTOFSpacetimeHit->x = 3.;
						locTOFSpacetimeHit->y = 66.; //paddle extends from 6 -> 126
						locTOFSpacetimeHit->t = locTOFHit->t_south - 60.0/VELOCITY;
						break;
					case 42:
						locTOFSpacetimeHit->x = -3.;
						locTOFSpacetimeHit->y = -66.; //paddle extends from -6 -> -126
						locTOFSpacetimeHit->t = locTOFHit->t_north - 60.0/VELOCITY;
						break;
					case 44:
						locTOFSpacetimeHit->x = 3.;
						locTOFSpacetimeHit->y = -66.; //paddle extends from -6 -> -126
						locTOFSpacetimeHit->t = locTOFHit->t_north - 60.0/VELOCITY;
						break;
					default:
						break;
				}
			}
			locTOFSpacetimeHit->TOFHit = locTOFHit;
			locTOFSpacetimeHits_Vertical.push_back(locTOFSpacetimeHit);
		} //end vertical hit
	} //end DTOFHit loop

	//find matches and sort them
	list<tof_spacetimehitmatch_t*> locTOFSpacetimeHitMatchList; //use list for sorting, vector for resource pool
	list<tof_spacetimehitmatch_t*>::iterator locListIterator;
	tof_spacetimehitmatch_t *locTOFSpacetimeHitMatch;
	for(loc_i = 0; loc_i < locTOFSpacetimeHits_Horizontal.size(); loc_i++){
		locTOFSpacetimeHit_Horizontal = locTOFSpacetimeHits_Horizontal[loc_i];
		for(loc_j = 0; loc_j < locTOFSpacetimeHits_Vertical.size(); loc_j++){
			locTOFSpacetimeHit_Vertical = locTOFSpacetimeHits_Vertical[loc_j];

			locDeltaX = locTOFSpacetimeHit_Horizontal->x - locTOFSpacetimeHit_Vertical->x;
			if (fabs(locDeltaX) >= locTOFSpacetimeHit_Horizontal->pos_cut)
				continue;
			locDeltaY = locTOFSpacetimeHit_Horizontal->y - locTOFSpacetimeHit_Vertical->y;
			if (fabs(locDeltaY) >= locTOFSpacetimeHit_Vertical->pos_cut)
				continue;
			locDeltaT = locTOFSpacetimeHit_Horizontal->t - locTOFSpacetimeHit_Vertical->t;
			locTimeCut = (locTOFSpacetimeHit_Horizontal->t_cut > locTOFSpacetimeHit_Vertical->t_cut) ? locTOFSpacetimeHit_Horizontal->t_cut : locTOFSpacetimeHit_Vertical->t_cut;
			if (fabs(locDeltaT) >= locTimeCut)
				continue;

			if(locTOFSpacetimeHitMatchList.size() == dTOFSpacetimeHitMatchPool.size()){ //increase pool size if necessary (used to minimize memory allocation time)
				locTOFSpacetimeHitMatch = new tof_spacetimehitmatch_t();
				dTOFSpacetimeHitMatchPool.push_back(locTOFSpacetimeHitMatch);
			}
			locTOFSpacetimeHitMatch = dTOFSpacetimeHitMatchPool[locTOFSpacetimeHitMatchList.size()];
			locTOFSpacetimeHitMatch->delta_t = locDeltaT;
			locTOFSpacetimeHitMatch->delta_r = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
			locTOFSpacetimeHitMatch->dTOFSpacetimeHit_Horizontal = locTOFSpacetimeHit_Horizontal;
			locTOFSpacetimeHitMatch->dTOFSpacetimeHit_Vertical = locTOFSpacetimeHit_Vertical;
			locTOFSpacetimeHitMatchList.push_back(locTOFSpacetimeHitMatch);
		}
	}
	locTOFSpacetimeHitMatchList.sort(Compare_TOFSpacetimeHitMatches_Distance); //sort matches by delta_r

	// delete pool size if too large, preventing memory-leakage-like behavor.
	if((locTOFSpacetimeHitMatchList.size() < MAX_TOFSpacetimeHitMatches) && (dTOFSpacetimeHitMatchPool.size() > MAX_TOFSpacetimeHitMatches)){
		for(loc_i = MAX_TOFSpacetimeHitMatches; loc_i < dTOFSpacetimeHitMatchPool.size(); loc_i++) delete dTOFSpacetimeHitMatchPool[loc_i];
		dTOFSpacetimeHitMatchPool.resize(MAX_TOFSpacetimeHitMatches);
	}

	// select best matches by delta_r
	while(locTOFSpacetimeHitMatchList.size() > 0){
		locTOFSpacetimeHitMatch = *(locTOFSpacetimeHitMatchList.begin());
		locTOFSpacetimeHit_Horizontal = locTOFSpacetimeHitMatch->dTOFSpacetimeHit_Horizontal;
		locTOFSpacetimeHit_Vertical = locTOFSpacetimeHitMatch->dTOFSpacetimeHit_Vertical;
		locTOFHit_Horizontal = locTOFSpacetimeHit_Horizontal->TOFHit;
		locTOFHit_Vertical = locTOFSpacetimeHit_Vertical->TOFHit;

		//see which bars had both PMT signals above threshold
		bool locOneSideBelowEThresholdFlag_Horizontal = false;
		bool locOneSideBelowEThresholdFlag_Vertical = false;
		if(!((locTOFHit_Horizontal->meantime >= 0.0) || (locTOFHit_Horizontal->meantime <= 0.0)))
			locOneSideBelowEThresholdFlag_Horizontal = true; //NaN: only one PMT signal above threshold
		if(!((locTOFHit_Vertical->meantime >= 0.0) || (locTOFHit_Vertical->meantime <= 0.0)))
			locOneSideBelowEThresholdFlag_Vertical = true; //NaN: only one PMT signal above threshold

		//Disable matching for certain scenarios until hit differentiation is improved
		if((locOneSideBelowEThresholdFlag_Horizontal == true) && (locOneSideBelowEThresholdFlag_Vertical == true)){ //both bars have only one PMT signal above threshold
			//currently, the software isn't optimized for this scenario, so disabled for now
				//e.g. separate DTOFPoints due to energy deposited in adjacent paddles, PID routine for picking the correct DTOFPoint if adjacent, etc.
			locTOFSpacetimeHitMatchList.erase(locTOFSpacetimeHitMatchList.begin());
			continue;
		}else if((locOneSideBelowEThresholdFlag_Horizontal == true) && (locTOFHit_Horizontal->bar < 41)){ //a horizontal double-sided bar has only one PMT signal above threshold
			//currently, the software isn't optimized for this scenario, so disabled for now
				//e.g. separate DTOFPoints due to energy deposited in adjacent paddles, PID routine for picking the correct DTOFPoint if adjacent, etc.
			locTOFSpacetimeHitMatchList.erase(locTOFSpacetimeHitMatchList.begin());
			continue;
		}else if((locOneSideBelowEThresholdFlag_Vertical == true) && (locTOFHit_Vertical->bar < 41)){ //a vertical double-sided bar has only one PMT signal above threshold
			//currently, the software isn't optimized for this scenario, so disabled for now
				//e.g. separate DTOFPoints due to energy deposited in adjacent paddles, PID routine for picking the correct DTOFPoint if adjacent, etc.
			locTOFSpacetimeHitMatchList.erase(locTOFSpacetimeHitMatchList.begin());
			continue;
		}

		locTOFPoint = new DTOFPoint;
		//reconstruct TOF hit information, using information from the best bar: one or both of the bars may have a PMT signal below threshold
		if((locOneSideBelowEThresholdFlag_Horizontal == true) && (locOneSideBelowEThresholdFlag_Vertical == true)){ //both bars have only one PMT signal above threshold
			locMatchX = locTOFSpacetimeHit_Vertical->x;
			locMatchY = locTOFSpacetimeHit_Horizontal->y;
			locMatchZ = 618.81; //z: midpoint between tof planes
			float locDepositedEnergy_Horizontal, locDepositedEnergy_Vertical;
			if(locTOFHit_Horizontal->E_north > E_THRESHOLD)
				locDepositedEnergy_Horizontal = locTOFHit_Horizontal->E_north * exp((HALFPADDLE - locMatchX)/ATTEN_LENGTH);
			else
				locDepositedEnergy_Horizontal = locTOFHit_Horizontal->E_south * exp((HALFPADDLE + locMatchX)/ATTEN_LENGTH);
			if(locTOFHit_Vertical->E_north > E_THRESHOLD)
				locDepositedEnergy_Vertical = locTOFHit_Vertical->E_north * exp((HALFPADDLE - locMatchY)/ATTEN_LENGTH);
			else
				locDepositedEnergy_Vertical = locTOFHit_Vertical->E_south * exp((HALFPADDLE + locMatchY)/ATTEN_LENGTH);
			locMatchdE = 0.5*(locDepositedEnergy_Horizontal + locDepositedEnergy_Vertical);
			locMatchT = 0.5*(locTOFSpacetimeHit_Horizontal->t + locTOFSpacetimeHit_Vertical->t);
		}else if(locOneSideBelowEThresholdFlag_Horizontal == true){ //the horizontal bar has only one PMT signal above threshold (and the vertical bar has both)
			locMatchX = locTOFSpacetimeHit_Vertical->x;
			locMatchY = locTOFSpacetimeHit_Vertical->y;
			locMatchT = locTOFSpacetimeHit_Vertical->t;
			locMatchZ = 617.52; //z: center of vertical plane
			locMatchdE = locTOFHit_Vertical->dE;
		}else if(locOneSideBelowEThresholdFlag_Vertical == true){ //the vertical bar has only one PMT signal above threshold (and the horizontal bar has both)
			locMatchX = locTOFSpacetimeHit_Horizontal->x;
			locMatchY = locTOFSpacetimeHit_Horizontal->y;
			locMatchT = locTOFSpacetimeHit_Horizontal->t;
			locMatchZ = 620.10; //z: center of horizontal plane
			locMatchdE = locTOFHit_Horizontal->dE;
		}else{ //both bars have both PMT signals above threshold
			locMatchX = locTOFSpacetimeHit_Horizontal->x;
			locMatchY = locTOFSpacetimeHit_Vertical->y;
			locMatchZ = 618.81; //z: midpoint between tof planes
			locMatchT = 0.5*(locTOFSpacetimeHit_Horizontal->t + locTOFSpacetimeHit_Vertical->t);
			locMatchdE = 0.5*(locTOFHit_Horizontal->dE + locTOFHit_Vertical->dE);
		}
		locTOFPoint->AddAssociatedObject(locTOFHit_Horizontal);
		locTOFPoint->AddAssociatedObject(locTOFHit_Vertical);
		locTOFPoint->pos.SetXYZ(locMatchX, locMatchY, locMatchZ);
		locTOFPoint->t = locMatchT;
		locTOFPoint->dE = locMatchdE;

		_data.push_back(locTOFPoint);
		//remove matches that contain the hits that were just used
		for(locListIterator = locTOFSpacetimeHitMatchList.begin(); locListIterator != locTOFSpacetimeHitMatchList.end(); locListIterator++){
			locTOFSpacetimeHitMatch = *locListIterator;
			if(locTOFSpacetimeHitMatch->dTOFSpacetimeHit_Horizontal == locTOFSpacetimeHit_Horizontal){
				locListIterator = locTOFSpacetimeHitMatchList.erase(locListIterator); //erase advances the iterator
				if(locTOFSpacetimeHitMatchList.size() == 0)
					break;
				locListIterator--; //for loop will advance iterator
			}else if(locTOFSpacetimeHitMatch->dTOFSpacetimeHit_Vertical == locTOFSpacetimeHit_Vertical){
				locListIterator = locTOFSpacetimeHitMatchList.erase(locListIterator); //erase advances the iterator
				if(locTOFSpacetimeHitMatchList.size() == 0)
					break; //is also the end, don't --, just break
				locListIterator--; //for loop will advance iterator
			}
		}
	}

	return NOERROR;
}



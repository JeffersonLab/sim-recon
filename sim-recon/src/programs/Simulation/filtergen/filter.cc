// $Id: smear.cc 2432 2007-02-06 04:19:48Z davidl $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>
#include "HDDM/hddm_s.h"


#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "

//-----------
// Filter
//-----------
bool Filter(s_HDDM_t *hddm_s)
{
	// Return "true" to keep event, "false" to throw it away
	
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return false;
	
	for(unsigned int i=0; i<PE->mult; i++){

		//------------- FCAL -------------
		double Efcal = 0.0;
      s_HitView_t *hits = PE->in[i].hitView;
      if (hits != HDDM_NULL &&
          hits->forwardEMcal != HDDM_NULL &&
          hits->forwardEMcal->fcalBlocks != HDDM_NULL){
      
			s_FcalBlocks_t *blocks = hits->forwardEMcal->fcalBlocks;
			for(unsigned int j=0; j<blocks->mult; j++){
				 s_FcalBlock_t *block = &blocks->in[j];
				 for(unsigned int k=0; k<block->fcalHits->mult; k++){
					  s_FcalHit_t *fcalhit = &block->fcalHits->in[k];
					  
					  Efcal += fcalhit->E;
				 } // k  (fcalhits)
			} // j  (blocks)
		} // if * != HDDM_NULL

//_DBG_<<"Efcal="<<Efcal<<endl;
		// There must be at least 0.5 GeV in the FCAL to pass the level1 trigger
		if(Efcal<0.5)return false;
		
		//------------- BCAL -------------
		double Ebcal = 0.0;
		if (hits != HDDM_NULL &&
			hits->barrelEMcal != HDDM_NULL &&
			hits->barrelEMcal->bcalCells != HDDM_NULL){
		
			// Loop over BCAL cells
			s_BcalCells_t *cells = hits->barrelEMcal->bcalCells;
			for(unsigned int j=0;j<cells->mult;j++){
				s_BcalCell_t *cell = &cells->in[j];
				if(cell->bcalHits != HDDM_NULL){
					for(unsigned int k=0; k<cell->bcalHits->mult; k++){
						s_BcalHit_t *hit = &cell->bcalHits->in[k];
						
						Ebcal += hit->E;
					}
				}
			} // j   (cells)
		} // if * != HDDM_NULL

//_DBG_<<"Ebcal="<<Ebcal<<endl;
		
		//------------- TOF -------------
		int Ntof_north = 0;
		int Ntof_south = 0;
		if (hits != HDDM_NULL &&
			hits->forwardTOF != HDDM_NULL &&
			hits->forwardTOF->ftofCounters != HDDM_NULL){
		
			s_FtofCounters_t* ftofCounters = hits->forwardTOF->ftofCounters;
		
			// Loop over counters
			s_FtofCounter_t *ftofCounter = ftofCounters->in;
			for(unsigned int j=0;j<ftofCounters->mult; j++, ftofCounter++){
				s_FtofNorthHit_t *ftofNorthHit = ftofCounter->ftofNorthHits->in;
				for(unsigned int k=0;k<ftofCounter->ftofNorthHits->mult; k++, ftofNorthHit++){
					// only count TOF hits in a reasonable range
					if(ftofNorthHit->t < 50.0 && ftofNorthHit->t > 0.0)Ntof_north++;
				}

				s_FtofSouthHit_t *ftofSouthHit = ftofCounter->ftofSouthHits->in;
				for(unsigned int k=0;k<ftofCounter->ftofSouthHits->mult; k++, ftofSouthHit++){
					// only count TOF hits in a reasonable range
					if(ftofSouthHit->t < 50.0 && ftofSouthHit->t > 0.0)Ntof_south++;
				}
			}
		}
		
		// We want the number of TOF coincidences which we'll estimate as the
		// lesser of the north and south hits
		int Ntof =  Ntof_north < Ntof_south ? Ntof_north:Ntof_south;
		
		// If there are no hits in the TOF or the BCAL has more energy, then
		// cut the event
//_DBG_<<"Ntof="<<Ntof<<endl;
		if(Ntof==0 || Ebcal>Efcal)return false;
		
		// Reject events with too many TOF hits
		if(Ntof>6)return false;
		
		return true;
	}
	
	return true;
}

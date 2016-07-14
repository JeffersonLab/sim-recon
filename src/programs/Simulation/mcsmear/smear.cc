// $Id$
//
// Created June 22, 2005  David Lawrence

#include "smear.h"

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <math.h>
#include "units.h"
#include <TF1.h>
#include <TH2.h>

#include "DRandom2.h"

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif


DCCALGeometry *ccalGeom = NULL;

//-----------
// Smear (constructor)
//-----------
Smear::Smear(mcsmear_config_t *in_config, JEventLoop *loop, string detectors_to_load) 
{   
	// create configuration classes
	config = in_config;
			
	// Create smearing classes 
	if(detectors_to_load == "all") {
		// default to loading all subdetectors
		smearers[SYS_BCAL]  = static_cast<Smearer*>(new BCALSmearer(loop,config));
		smearers[SYS_FCAL]  = static_cast<Smearer*>(new FCALSmearer(loop,config));
		smearers[SYS_CDC]   = static_cast<Smearer*>(new CDCSmearer(loop,config));
		smearers[SYS_FDC]   = static_cast<Smearer*>(new FDCSmearer(loop,config));
		smearers[SYS_TOF]   = static_cast<Smearer*>(new TOFSmearer(loop,config));
		smearers[SYS_START] = static_cast<Smearer*>(new SCSmearer(loop,config));
		smearers[SYS_TAGH]  = static_cast<Smearer*>(new TAGHSmearer(loop,config));
		smearers[SYS_TAGM]  = static_cast<Smearer*>(new TAGMSmearer(loop,config));
		smearers[SYS_PS]    = static_cast<Smearer*>(new PSSmearer(loop,config));
		smearers[SYS_PSC]   = static_cast<Smearer*>(new PSCSmearer(loop,config));
		smearers[SYS_DIRC]  = static_cast<Smearer*>(new FDIRCSmearer(loop,config));
		smearers[SYS_CCAL]  = static_cast<Smearer*>(new CCALSmearer(loop,config));
		smearers[SYS_FMWPC] = static_cast<Smearer*>(new FMWPCSmearer(loop,config));
	} else {
		// Parse string of system names
    	std::istringstream ss(detectors_to_load);
    	std::string token;
    	while(std::getline(ss, token, ',')) {
			DetectorSystem_t the_detector = static_cast<DetectorSystem_t>(atoi(token.c_str()));
			switch(the_detector) {
				case SYS_BCAL:   smearers[the_detector] = static_cast<Smearer*>(new BCALSmearer(loop,config));  break;
				case SYS_FCAL:   smearers[the_detector] = static_cast<Smearer*>(new FCALSmearer(loop,config));  break;
				case SYS_CDC:    smearers[the_detector] = static_cast<Smearer*>(new CDCSmearer(loop,config));  break;
				case SYS_FDC:    smearers[the_detector] = static_cast<Smearer*>(new FDCSmearer(loop,config));  break;
				case SYS_TOF:    smearers[the_detector] = static_cast<Smearer*>(new TOFSmearer(loop,config));  break;
				case SYS_START:  smearers[the_detector] = static_cast<Smearer*>(new SCSmearer(loop,config));  break;
				case SYS_TAGH:   smearers[the_detector] = static_cast<Smearer*>(new TAGHSmearer(loop,config));  break;
				case SYS_TAGM:   smearers[the_detector] = static_cast<Smearer*>(new TAGMSmearer(loop,config));  break;
				case SYS_PS:     smearers[the_detector] = static_cast<Smearer*>(new PSSmearer(loop,config));  break;
				case SYS_PSC:    smearers[the_detector] = static_cast<Smearer*>(new PSCSmearer(loop,config));  break;
				case SYS_DIRC:   smearers[the_detector] = static_cast<Smearer*>(new FDIRCSmearer(loop,config));  break;
				case SYS_CCAL:   smearers[the_detector] = static_cast<Smearer*>(new CCALSmearer(loop,config));  break;
				case SYS_FMWPC:  smearers[the_detector] = static_cast<Smearer*>(new FMWPCSmearer(loop,config));  break;
                default:  break;   // don't smear any other detectors
			}
		}
	}
}
		
//-----------
// Smear (destructor)
//-----------
Smear::~Smear() 
{   
	for(map<DetectorSystem_t, Smearer *>::iterator smearer_it = smearers.begin();
		smearer_it != smearers.end(); smearer_it++)
		delete smearer_it->second;
}
		
//-----------
// SmearEvent
//-----------
void Smear::SmearEvent(hddm_s::HDDM *record)
{
    GetAndSetSeeds(record);

	// Smear each detector system
	for(map<DetectorSystem_t, Smearer *>::iterator smearer_it = smearers.begin();
		smearer_it != smearers.end(); smearer_it++) {
        //cerr << "smearing " << SystemName(smearer_it->first) << endl;
		smearer_it->second->SmearEvent(record);
    }

}

//-----------
// SetSeeds
//-----------
void Smear::SetSeeds(const char *vals)
{
   /// This is called from the command line parser to
   /// set the initial seeds based on user input from
   /// the command line.
   //
   //
   stringstream ss(vals);
   Int_t seed1, seed2, seed3;
   ss >> seed1 >> seed2 >> seed3;
   UInt_t *useed1 = reinterpret_cast<UInt_t*>(&seed1);
   UInt_t *useed2 = reinterpret_cast<UInt_t*>(&seed2);
   UInt_t *useed3 = reinterpret_cast<UInt_t*>(&seed3);
   gDRandom.SetSeeds(*useed1, *useed2, *useed3);

   cout << "Seeds set from command line. Any random number" << endl;
   cout << "seeds found in the input file will be ignored!" << endl;
   config->IGNORE_SEEDS = true;
}

//-----------
// GetAndSetSeeds
//-----------
void Smear::GetAndSetSeeds(hddm_s::HDDM *record)
{
   // Check if non-zero seed values exist in the input HDDM file.
   // If so, use them to set the seeds for the random number
   // generator. Otherwise, make sure the seeds that are used
   // are stored in the output event.
   
   if (record == 0)
      return;
   else if (record->getReactions().size() == 0)
      return;

   hddm_s::ReactionList::iterator reiter = record->getReactions().begin();
   if (reiter->getRandoms().size() == 0) {
      // No seeds stored in event. Add them
      hddm_s::RandomList blank_rand = reiter->addRandoms();
      blank_rand().setSeed1(0);
      blank_rand().setSeed2(0);
      blank_rand().setSeed3(0);
      blank_rand().setSeed4(0);
   }

   UInt_t seed1, seed2, seed3;
   hddm_s::Random my_rand = reiter->getRandom();

   if (!config->IGNORE_SEEDS) {
      // Copy seeds from event record to local variables
      seed1 = my_rand.getSeed1();
      seed2 = my_rand.getSeed2();
      seed3 = my_rand.getSeed3();
      
      // If the seeds in the event are all zeros it means they
      // were not set. In this case, initialize seeds to constants
      // to guarantee the seeds are used if this input file were
      // smeared again with the same command a second time. These
      // are set here to the fractional part of the cube roots of
      // the first three primes, truncated to 9 digits.
      if ((seed1 == 0) || (seed2 == 0) || (seed3 == 0)){
         uint64_t eventNo = record->getPhysicsEvent().getEventNo();
         seed1 = 259921049 + eventNo;
         seed2 = 442249570 + eventNo;
         seed3 = 709975946 + eventNo;
      }
      
      // Set the seeds in the random generator.
      gDRandom.SetSeeds(seed1, seed2, seed3);
   }

   // Copy seeds from generator to local variables
   gDRandom.GetSeeds(seed1, seed2, seed3);

   // Copy seeds from local variables to event record
   my_rand.setSeed1(seed1);
   my_rand.setSeed2(seed2);
   my_rand.setSeed3(seed3);
}



#include "SCSmearer.h"

//-----------
// sc_config_t  (constructor)
//-----------
sc_config_t::sc_config_t(JEventLoop *loop) 
{
	// default values
	START_SIGMA           = 0.0; // 300ps
	START_PHOTONS_PERMEV  = 0.0; // used to be 8000 should be more like 200
	START_PADDLE_THRESHOLD  = 0.0;

	// Load data from CCDB
    cout << "Get START_COUNTER/start_parms parameters from CCDB..." << endl;
    map<string, double> startparms;
    if(loop->GetCalib("START_COUNTER/start_parms", startparms)) {
		jerr << "Problem loading START_COUNTER/start_parms from CCDB!" << endl;
	} else {
     	START_SIGMA = startparms["START_SIGMA"] ;
     	START_PHOTONS_PERMEV = startparms["START_PHOTONS_PERMEV"];
	}
	
	cout<<"get START_COUNTER/time_resol_paddle from calibDB"<<endl;
    vector <double> START_TIME_RESOLUTIONS_TEMP;
    if(loop->GetCalib("START_COUNTER/time_resol_paddle", START_TIME_RESOLUTIONS_TEMP)) {
    	jerr << "Problem loading START_COUNTER/time_resol_paddle from CCDB!" << endl;
    } else {
    	for (unsigned int i = 0; i < START_TIME_RESOLUTIONS_TEMP.size(); i++) {
       		START_TIME_RESOLUTIONS.push_back(START_TIME_RESOLUTIONS_TEMP.at(i));
    	}
    }
    
	cout<<"get START_COUNTER/paddle_mc_efficiency from calibDB"<<endl;
    if(loop->GetCalib("START_COUNTER/paddle_mc_efficiency", paddle_efficiencies)) {
    	jerr << "Problem loading START_COUNTER/paddle_mc_efficiency from CCDB!" << endl;
    }
}


//-----------
// SmearEvent
//-----------
void SCSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::StcPaddleList pads = record->getStcPaddles();
   hddm_s::StcPaddleList::iterator iter;
   for (iter = pads.begin(); iter != pads.end(); ++iter) {
      iter->deleteStcHits();
      hddm_s::StcTruthHitList thits = iter->getStcTruthHits();
      hddm_s::StcTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
		 if(!gDRandom.DecideToAcceptHit(sc_config->GetMCEfficiency(iter->getSector())))
		 	continue;

         // smear the time
         //double t = titer->getT() + gDRandom.SampleGaussian(sc_config->START_SIGMA);   // constant smearing
         double t = titer->getT() + gDRandom.SampleGaussian(sc_config->GetPaddleTimeResolution(iter->getSector()));
         // smear the energy
         double npe = titer->getDE() * 1000. *  sc_config->START_PHOTONS_PERMEV;
         npe = npe +  gDRandom.SampleGaussian(sqrt(npe));
         double NewE = npe/sc_config->START_PHOTONS_PERMEV/1000.;
         if (NewE > sc_config->START_PADDLE_THRESHOLD) {
            hddm_s::StcHitList hits = iter->addStcHits();
            hits().setT(t);
            hits().setDE(NewE);
         }
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteStcTruthHits();
   }
}

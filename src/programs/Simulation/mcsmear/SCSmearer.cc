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
         // smear the time
         double t = titer->getT() + SampleGaussian(sc_config->START_SIGMA);
         // smear the energy
         double npe = titer->getDE() * 1000. *  sc_config->START_PHOTONS_PERMEV;
         npe = npe +  SampleGaussian(sqrt(npe));
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

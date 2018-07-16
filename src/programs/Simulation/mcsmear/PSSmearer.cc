#include "PSSmearer.h"

//-----------
// ps_config_t  (constructor)
//-----------
ps_config_t::ps_config_t(JEventLoop *loop) 
{
	// default values
	PS_SIGMA = 0.200; // ns
	PS_NPIX_PER_GEV = 1.e5;
	PS_THRESHOLD          = 0.0;
}

	
//-----------
// SmearEvent
//-----------
void PSSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::PsTileList tiles = record->getPsTiles();
   hddm_s::PsTileList::iterator iter;
   for (iter = tiles.begin(); iter != tiles.end(); ++iter) {
      iter->deletePsHits();
      hddm_s::PsTruthHitList thits = iter->getPsTruthHits();
      hddm_s::PsTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time
         double t = titer->getT();
         double dE = titer->getDE();
         if(config->SMEAR_HITS) {
			 t += gDRandom.SampleGaussian(ps_config->PS_SIGMA);
         	 // convert energy deposition in number of fired pixels
         	 double npe = gDRandom.SamplePoisson( titer->getDE() * ps_config->PS_NPIX_PER_GEV);
			 dE = npe/ps_config->PS_NPIX_PER_GEV;
		 }
	 	 hddm_s::PsHitList hits = iter->addPsHits();
	 	 hits().setT(t);
	 	 hits().setDE(dE);
      }

      if (config->DROP_TRUTH_HITS)
         iter->deletePsTruthHits();
   }
}

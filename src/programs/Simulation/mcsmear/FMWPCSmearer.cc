#include "FMWPCSmearer.h"

//-----------
// fmwpc_config_t  (constructor)
//-----------
fmwpc_config_t::fmwpc_config_t(JEventLoop *loop) 
{
	// default values
	FMWPC_TSIGMA = 10.0;  // ns
 	FMWPC_ASIGMA = 0.5E-6;
 	FMWPC_THRESHOLD = 0.0;
}



//-----------
// SmearEvent
//-----------
void FMWPCSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::FmwpcChamberList chambers = record->getFmwpcChambers();
   hddm_s::FmwpcChamberList::iterator iter;
   for (iter = chambers.begin(); iter != chambers.end(); ++iter) {
      iter->deleteFmwpcHits();
      hddm_s::FmwpcTruthHitList thits = iter->getFmwpcTruthHits();
      hddm_s::FmwpcTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time and energy
         double t = titer->getT();
         double dE = titer->getDE();
         if(config->SMEAR_HITS) {
         	t += gDRandom.SampleGaussian(fmwpc_config->FMWPC_TSIGMA);
         	dE += gDRandom.SampleGaussian(fmwpc_config->FMWPC_ASIGMA);
		 }
         if (dE > fmwpc_config->FMWPC_THRESHOLD) {
            hddm_s::FmwpcHitList hits = iter->addFmwpcHits();
            hits().setT(t);
            hits().setDE(dE);
         }
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteFmwpcTruthHits();
   }
}

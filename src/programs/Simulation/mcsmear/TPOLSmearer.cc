#include "TPOLSmearer.h"

//-----------
// tpol_config_t  (constructor)
//-----------
tpol_config_t::tpol_config_t(JEventLoop *loop) 
{
	// default values
	TPOL_SIGMA_NS = 4.4;  // ns
	TPOL_SIGMA1_MEV = 0.03; // MeV
        TPOL_SIGMA2_MEV = 0.03; // MeV
	TPOL_THRESHOLD_MEV = 0.05; // MeV
}

	
//-----------
// SmearEvent
//-----------
void TPOLSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::TpolSectorList sectors = record->getTpolSectors();
   hddm_s::TpolSectorList::iterator iter;
   for (iter = sectors.begin(); iter != sectors.end(); ++iter) {
      iter->deleteTpolHits();
      hddm_s::TpolTruthHitList thits = iter->getTpolTruthHits();
      hddm_s::TpolTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time
         double t_ns = titer->getT();
         // smear the energy, convert to MeV
         double dE_MeV = titer->getDE() * 1e3;
         if(config->SMEAR_HITS) {
			 t_ns += gDRandom.SampleGaussian(tpol_config->TPOL_SIGMA_NS);
			 dE_MeV += (gDRandom.SampleGaussian(tpol_config->TPOL_SIGMA1_MEV) + gDRandom.SampleGaussian(tpol_config->TPOL_SIGMA2_MEV));
		 }
         // apply the threshold
         if (dE_MeV > tpol_config->TPOL_THRESHOLD_MEV) {
	        hddm_s::TpolHitList hits = iter->addTpolHits();
	        hits().setT(t_ns);
	        hits().setDE(dE_MeV);
         }
      }
      if (config->DROP_TRUTH_HITS)
         iter->deleteTpolTruthHits();
   }
}

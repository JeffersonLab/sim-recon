#include "PSCSmearer.h"

//-----------
// psc_config_t  (constructor)
//-----------
psc_config_t::psc_config_t(JEventLoop *loop)
{
	// default values
	PSC_SIGMA = 0.200; //ns
	PSC_PHOTONS_PERMEV = 5.e5;
	PSC_THRESHOLD         = 0.0;
}
	

//-----------
// SmearEvent
//-----------
void PSCSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::PscPaddleList paddles = record->getPscPaddles();
   hddm_s::PscPaddleList::iterator iter;
   for (iter = paddles.begin(); iter != paddles.end(); ++iter) {
      iter->deletePscHits();
      hddm_s::PscTruthHitList thits = iter->getPscTruthHits();
      hddm_s::PscTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time
         double t = titer->getT();
         double NewE = titer->getDE();
         if(config->SMEAR_HITS) {
         	t += gDRandom.SampleGaussian(psc_config->PSC_SIGMA);
         	double npe = titer->getDE() * 1000. *  psc_config->PSC_PHOTONS_PERMEV;
         	npe = npe +  gDRandom.SampleGaussian(sqrt(npe));
        	NewE = npe/psc_config->PSC_PHOTONS_PERMEV/1000.;
		 }
         if (NewE > psc_config->PSC_THRESHOLD) {
            hddm_s::PscHitList hits = iter->addPscHits();
            hits().setT(t);
            hits().setDE(NewE);
         }
      }

      if (config->DROP_TRUTH_HITS)
         iter->deletePscTruthHits();
   }
}
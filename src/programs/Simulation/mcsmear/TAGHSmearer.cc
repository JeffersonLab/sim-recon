#include "TAGHSmearer.h"


//-----------
// tagh_config_t  (constructor)
//-----------
tagh_config_t::tagh_config_t(JEventLoop *loop) 
{
	// default values
	TAGH_TSIGMA = 0.350;        // ns
	TAGH_FADC_TSIGMA = 0.450;   // ns
	TAGH_NPE_PER_GEV = 5.e5;

}


//-----------
// SmearEvent
//-----------
void TAGHSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::HodoChannelList taghs = record->getHodoChannels();
   hddm_s::HodoChannelList::iterator iter;
   for (iter = taghs.begin(); iter != taghs.end(); ++iter) {
      iter->deleteTaggerHits();
      hddm_s::TaggerTruthHitList thits = iter->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time
         double t = titer->getT();
         double tADC = titer->getT();
         double npe = titer->getDE() * tagh_config->TAGH_NPE_PER_GEV;

         if(config->SMEAR_HITS) {
        	t += gDRandom.SampleGaussian(tagh_config->TAGH_TSIGMA);
         	tADC += gDRandom.SampleGaussian(tagh_config->TAGH_FADC_TSIGMA);
         	npe = gDRandom.SamplePoisson(titer->getDE() * tagh_config->TAGH_NPE_PER_GEV);
		 }
         hddm_s::TaggerHitList hits = iter->addTaggerHits();
         hits().setT(t);
         hits().setTADC(tADC);
         hits().setNpe(npe);
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteTaggerTruthHits();
   }
}
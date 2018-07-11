#include "TAGMSmearer.h"

//-----------
// tagm_config_t  (constructor)
//-----------
tagm_config_t::tagm_config_t(JEventLoop *loop) {
	// default values
	TAGM_TSIGMA = 0.200;        // ns
	TAGM_FADC_TSIGMA = 0.350;   // ns
	TAGM_NPIX_PER_GEV = 1.e5;
}

//-----------
// SmearEvent
//-----------
void TAGMSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::MicroChannelList tagms = record->getMicroChannels();
   hddm_s::MicroChannelList::iterator iter;
   for (iter = tagms.begin(); iter != tagms.end(); ++iter) {
      iter->deleteTaggerHits();
      hddm_s::TaggerTruthHitList thits = iter->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time
         double t = titer->getT();
         double tADC = titer->getT();
         double npe = titer->getDE() * tagm_config->TAGM_NPIX_PER_GEV;
         if(config->SMEAR_HITS) {
         	t += gDRandom.SampleGaussian(tagm_config->TAGM_TSIGMA);
          	tADC += gDRandom.SampleGaussian(tagm_config->TAGM_FADC_TSIGMA);
          	npe = gDRandom.SamplePoisson(titer->getDE() * tagm_config->TAGM_NPIX_PER_GEV);
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

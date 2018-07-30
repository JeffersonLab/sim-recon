#include "CCALSmearer.h"

//-----------
// ccal_config_t  (constructor)
//-----------
ccal_config_t::ccal_config_t(JEventLoop *loop) {
  // default values
  // (This is just a rough estimate 11/30/2010 DL)
  CCAL_PHOT_STAT_COEF = 0.035/2.0;
  CCAL_BLOCK_THRESHOLD = 20.0*k_MeV;
  CCAL_SIGMA = 200.0e-3;
}



//-----------
// SmearEvemt
//-----------
void CCALSmearer::SmearEvent(hddm_s::HDDM *record){
  /// Smear the CCAL hits using the same procedure as the FCAL above.
  /// See those comments for details.
  
  //   if (!ccalGeom)
  //   ccalGeom = new DCCALGeometry();
  
  hddm_s::CcalBlockList blocks = record->getCcalBlocks();   
  hddm_s::CcalBlockList::iterator iter;
  for (iter = blocks.begin(); iter != blocks.end(); ++iter) {
    iter->deleteCcalHits();
    hddm_s::CcalTruthHitList thits = iter->getCcalTruthHits();   
    hddm_s::CcalTruthHitList::iterator titer;
    for (titer = thits.begin(); titer != thits.end(); ++titer) {
      // Simulation simulates a grid of blocks for simplicity. 
      // Do not bother smearing inactive blocks. They will be
      // discarded in DEventSourceHDDM.cc while being read in
      // anyway.
      
      if (!ccalGeom->isBlockActive(iter->getRow(), iter->getColumn()))
		continue;
      // Smear the energy and timing of the hit
      //      double sigma = ccal_config->CCAL_PHOT_STAT_COEF/sqrt(titer->getE()) ;
      
      // A.S.  new calibration of the CCAL
      double E = titer->getE();
      double t = titer->getT();

	  if(config->SMEAR_HITS) {
      	double nphav = E * 2.3e3; // per GeV  Corrections
      
      	if(nphav < 30)
			E *= gDRandom.SamplePoisson(nphav)/nphav;           //photostatistics
      	else 
			E *= 1.0 + gDRandom.SampleGaussian(1./sqrt(nphav)); //photostatistics
      
      	E *= 1.167 + gDRandom.SampleGaussian(0.006);          // calibration
      	t += gDRandom.SampleGaussian(ccal_config->CCAL_SIGMA);
      }
        
      // Apply a single block threshold. If the (smeared) energy is below this,
      // then set the energy and time to zero. 	 
      // A.S. 
      //         if (E > ccal_config->CCAL_BLOCK_THRESHOLD) {
      hddm_s::CcalHitList hits = iter->addCcalHits();
      hits().setE(E);
      hits().setT(t);
      //         }
    }
    
    if (config->DROP_TRUTH_HITS)
      iter->deleteCcalTruthHits();
  }
}

// $Id: smear.cc 7650 2011-03-29 22:52:30Z shepherd $
//
// Created June 22, 2005  David Lawrence
//
// Major revision March 6, 2012 David Lawrence

#include "BCALSmearer.h"

#include "DRandom2.h"

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

//TH1D *hNincident_particles = NULL;


// Defined in this file
//int32_t GetRunNumber(hddm_s::HDDM *record);

//-----------
// Smear
//-----------
void BCALSmearer::SmearEvent(hddm_s::HDDM *record)
{

   /// May 27, 2015: HDGeant now outputs BCAL hit data in terms of just
   /// an energy and a time, rather than using time histograms which track
   /// energy deposition for each cell over 400 ns or so.  This greatly
   /// decreases the processing time required for BCAL smearing.
   ///
   /// The data from HDGEANT contains attenuated energies and times.
   /// The job of mcsmear is to use that data to create the bcalfADCUpHit and
   /// bcalfADCDownHit structures. These are made by summing signals from multiple
   /// SiPMs and are what are used as the entry points for the reconstruction.
   ///
   /// The energies must be smeared due to sampling fluctuations, which are
   /// parameterized based on the total energy of the shower. The total
   /// energy is kept in the bcalTruthIncidentParticle structures in HDDM.
   ///
   /// In addition to the sampling fluctuations, Poisson statistics and
   /// dark pulses are applied.
   
   // n.b. This code is slightly more complex than it might otherwise be because
   // it uses sparsified lists as opposed to full data structures for every SiPM
   // and every summed cell. It is done this way for two reasons:
   //
   // 1.) Sparsified lists are quicker to loop through and make the code faster.
   //
   // 2.) Sparsified lists are easier to keep on the stack rather than the heap
   //     avoiding expensive, large memory allocations every event. Alternatively
   //     one could keep the large structures in global variables but that would 
   //     not allow for efficient multi-threading.
   
    // First, we extract the energies and times for hit cells
    map<bcal_index, CellHits> SiPMHits;
    vector<IncidentParticle_t> incident_particles;
    GetSiPMHits(record, SiPMHits, incident_particles);

    // Sampling fluctuations
	if(config->SMEAR_HITS)
    	ApplySamplingFluctuations(SiPMHits, incident_particles);
	
    // Merge hits associated with different incident particles
    MergeHits(SiPMHits, bcal_config->BCAL_TWO_HIT_RESO);

    // Poisson Statistics
	if(config->SMEAR_HITS) 
    	ApplyPoissonStatistics(SiPMHits);
   
    // Place all hit cells into list indexed by fADC ID
    map<int, SumHits> bcalfADC;
    SortSiPMHits(SiPMHits, bcalfADC, bcal_config->BCAL_TWO_HIT_RESO);

    // Electronic noise/Dark hits Smearing
	if(config->SMEAR_HITS) 
    	SimpleDarkHitsSmear(bcalfADC);
    
    // Apply energy threshold to dismiss low-energy hits
    map<int, fADCHitList> fADCHits;
    map<int, TDCHitList> TDCHits;
    FindHits(bcal_config->BCAL_ADC_THRESHOLD_MEV, bcalfADC, fADCHits, TDCHits);

    // Apply time smearing to emulate the fADC resolution
	if(config->SMEAR_HITS) 
    	ApplyTimeSmearing(bcal_config->BCAL_FADC_TIME_RESOLUTION, bcal_config->BCAL_TDC_TIME_RESOLUTION, fADCHits, TDCHits);
   
    // Copy hits into HDDM tree
    CopyBCALHitsToHDDM(fADCHits, TDCHits, record);
   
    bcalfADC.clear();
}

int inline BCALSmearer::GetCalibIndex(int module, int layer, int sector) {
   return bcal_config->BCAL_NUM_LAYERS*bcal_config->BCAL_NUM_SECTORS*(module-1) 
   		+ bcal_config->BCAL_NUM_SECTORS*(layer-1) + (sector-1);
}


//-----------
// GetSiPMHits
//-----------
void BCALSmearer::GetSiPMHits(hddm_s::HDDM *record,
                   			  map<bcal_index, CellHits> &SiPMHits,
                   			  vector<IncidentParticle_t> &incident_particles)
{
   /// Loop through input HDDM data and extract the energy and time info into
   /// CellHits objects.
   
   // Make sure HDDM stuctures exist.
   // In the case of no real BCAL hits, we may still want to emit
   // dark hit only events. In this case, we must create the BCAL 
   // tree here.
   hddm_s::BarrelEMcalList bcals = record->getBarrelEMcals();
   if (bcals.size() == 0){
      if(record->getHitViews().empty()){
		record->getPhysicsEvent().addHitViews();
	  }
     bcals = record->getHitViews().begin()->addBarrelEMcals();
   }

   // Loop over GEANT hits in BCAL
   hddm_s::BcalTruthHitList hits = record->getBcalTruthHits();
   hddm_s::BcalTruthHitList::iterator iter;
   for (iter = hits.begin(); iter != hits.end(); ++iter) {
      bcal_index idxup(iter->getModule(), iter->getLayer(),
                     iter->getSector(), 
                     iter->getIncident_id(),
                     bcal_index::kUp);
      bcal_index idxdn(iter->getModule(), iter->getLayer(),
                     iter->getSector(), 
                     iter->getIncident_id(),
                     bcal_index::kDown);

     double Z = iter->getZLocal();
     double dist_up = 390.0/2.0 + Z;
     double dist_dn = 390.0/2.0 - Z;
     
     int layer = 0;
     if (iter->getLayer() == 1){
     	layer = 1;
     } else if (iter->getLayer() == 2 || iter->getLayer() == 3){
        layer = 2;
     } else if (iter->getLayer() == 4 || iter->getLayer() == 5 || iter->getLayer() == 6) {
     	layer = 3;
     } else {
     	layer = 4;
     }
     int table_id = GetCalibIndex( iter->getModule(), layer, iter->getSector() );  // key the cell identification off of the upstream cell
     double cEff = bcal_config->GetEffectiveVelocity(table_id);
     //double attenuation_length = 0; // initialize variable
     //double attenuation_L1=-1., attenuation_L2=-1.;  // these parameters are ignored for now
     //bcal_config->GetAttenuationParameters(table_id, attenuation_length, attenuation_L1, attenuation_L2);
    
     // Get reference to existing CellHits, or create one if it doesn't exist
     CellHits &cellhitsup = SiPMHits[idxup];
     cellhitsup.Etruth = iter->getE(); // Energy deposited in the cell in GeV
     //cellhitsup.E = iter->getE()*exp(-dist_up/attenuation_length)*1000.; // in attenuated MeV
     cellhitsup.E = iter->getE()*exp(-dist_up/bcal_config->BCAL_ATTENUATION_LENGTH)*1000.; // in attenuated MeV
     cellhitsup.t = iter->getT() + dist_up/cEff; // in ns
     cellhitsup.end = CellHits::kUp; // Keep track of BCal end

     // Get reference to existing CellHits, or create one if it doesn't exist
     CellHits &cellhitsdn = SiPMHits[idxdn];
     cellhitsdn.Etruth = iter->getE(); // Energy deposited in the cell in GeV
     cellhitsdn.E = iter->getE()*exp(-dist_dn/bcal_config->BCAL_ATTENUATION_LENGTH)*1000.; // in attenuated MeV
     cellhitsdn.t = iter->getT() + dist_dn/cEff; // in ns
     cellhitsdn.end = CellHits::kDown; // Keep track of BCal end
   }

   // Loop over incident particle list
   hddm_s::BcalTruthIncidentParticleList iparts = 
                                    bcals().getBcalTruthIncidentParticles();
   hddm_s::BcalTruthIncidentParticleList::iterator piter;
   int pcount = 0;
   for (piter = iparts.begin(); piter != iparts.end(); ++piter) {
      incident_particles.push_back(IncidentParticle_t(*piter));
      if (piter->getId() != ++pcount) {
         // If this ever gets called, we'll need to implement a sort routine
         _DBG_ << "Incident particle order not preserved!" << endl;
         exit(-1);
      }
   }
   
   //if (hNincident_particles)
   //   hNincident_particles->Fill(incident_particles.size());
}

//-----------
// ApplySamplingFluctuations
//-----------
void BCALSmearer::ApplySamplingFluctuations(map<bcal_index, CellHits> &SiPMHits, vector<IncidentParticle_t> &incident_particles)
{
   /// Loop over the CellHits objects and apply sampling fluctuations.
   ///
   /// Here we apply a statistical error due to the sampling
   /// fluctuations. The total energy (Etruth) is integrated by hdgeant.
   /// We calculate a sigma based on the deposited energy only. In
   /// reality, the sampling fluctuations are also a function of the
   /// angle of the shower particles w.r.t. the fibers. We do not include
   /// any angular dependence at this time. To do so will require more
   /// information be passed from hdgeant.
   ///
   /// The error is applied by finding the ratio of the smeared
   /// cell energy to unsmeared cell energy and scaling the energy
   /// by it.
   
   if(bcal_config->NO_SAMPLING_FLUCTUATIONS)return;
   if(bcal_config->NO_SAMPLING_FLOOR_TERM)
   		bcal_config->BCAL_SAMPLINGCOEFB=0.0; // (redundant, yes, but located in more obvious place here)

   map<bcal_index, CellHits>::iterator iter=SiPMHits.begin();
   for(; iter!=SiPMHits.end(); iter++){
      CellHits &cellhits = iter->second;
      
      // Find fractional sampling sigma based on deposited energy (whole colorimeter, not just fibers)
      double Etruth = cellhits.Etruth;
      double sqrtterm = bcal_config->BCAL_SAMPLINGCOEFA / sqrt( Etruth );
      double linterm = bcal_config->BCAL_SAMPLINGCOEFB;
      double sigmaSamp = sqrt(sqrtterm*sqrtterm + linterm*linterm);   

      // Convert sigma into GeV
      sigmaSamp *= Etruth;

      // Randomly sample the fluctuation
      double Esmeared = gDRandom.Gaus(Etruth,sigmaSamp);

      // Calculate ratio of smeared to unsmeared
      double ratio = Esmeared/Etruth;

      // Scale attenuated energy
      cellhits.E *= ratio;
   }
}

//-----------
// MergeHits
//-----------
void BCALSmearer::MergeHits(map<bcal_index, CellHits> &SiPMHits, double Resolution)
{
   /// Combine all SiPM CellHits corresponding to the same
   /// cell but different incident particles into a single
   /// hit. This is done after the sampling fluctuations
   /// have been applied so there is no more dependence on
   /// the incident particle parameters.
   
   // Loop until no merges are made
   while(true){
      bool merge=false;
      bool merged=false;
      map<bcal_index, CellHits>::iterator iter1=SiPMHits.begin();
      for(;iter1!=SiPMHits.end(); iter1++){
         map<bcal_index, CellHits>::iterator iter2 = iter1;
         for(++iter2; iter2!=SiPMHits.end(); iter2++){
         
            // If hits are not from same module,layer,sector,end
            // then just continue the loop
            if(iter1->first.module != iter2->first.module)continue;
            if(iter1->first.layer  != iter2->first.layer )continue;
            if(iter1->first.sector != iter2->first.sector)continue;
            if(iter1->first.end    != iter2->first.end   )continue;

            // If hits are far enough apart in time, don't merge them
            if(fabs(iter1->second.t - iter2->second.t) < Resolution) merge = true;
            if(!merge)continue;

            // ----- Merge hits -----
            if(merge){
              // Get values
              double E1 = iter1->second.E;
              double t1 = iter1->second.t;
              double E2 = iter2->second.E;
              double t2 = iter2->second.t;
              // It may be possible that one or both of the hits we wish to merge
              // don't exist. Check for this and handle accordingly.
              if(E1!=0.0 && E2!=0.0){
                 iter1->second.E += E2;
                 if(t1 > t2) iter1->second.t = t2; // Keep the earlier of the two times
              }
              if(E1==0.0 && E2!=0.0){
                 iter1->second.E = E2;
                 iter1->second.t = t2;
                 iter2->second.E = 0.0;
                 iter2->second.t = 0.0;
              }
            }

            // Erase second one
            SiPMHits.erase(iter2);

            // Set flag that we did merge hits and break
            // the loops so we can try again.
            merged = true;
            break;
         }
         if(merged)break;
      }

      // When we make it through without merging any hits,
      // we're done so break out of the infinite while loop.
      if(!merged)break;
   }
}

//-----------
// ApplyPoissonStatistics
//-----------
void BCALSmearer::ApplyPoissonStatistics(map<bcal_index, CellHits> &SiPMHits)
{
   /// Loop over the CellHits objects and apply Poisson Statistics.
   ///
   /// Because the response of the SiPM is quantized in units of photo-electrons
   /// Poisson counting statistics should be applied. This will affect the
   /// smaller energy depositions more than the larger ones.
   ///
   /// We do this by converting the cell's attenuated energy into
   /// photo-electrons and then sampling from a Poisson distribution with that
   /// mean. The ratio of the quantized, sampled value to the unquantized
   /// integral (in PE) is used to scale the energy.

   if(bcal_config->NO_POISSON_STATISTICS) return;

   map<bcal_index, CellHits>::iterator iter=SiPMHits.begin();
   for(; iter!=SiPMHits.end(); iter++){
      CellHits &cellhits = iter->second;

      if(cellhits.E>0.0){
         // Convert to number of PE
         double mean_pe = cellhits.E/bcal_config->BCAL_mevPerPE;
         
         int Npe = gDRandom.Poisson(mean_pe);
         double ratio = (double)Npe/mean_pe;

         cellhits.E *= ratio;
      }
   }
}

//-----------
// SortSiPMHits
//-----------
void BCALSmearer::SortSiPMHits(map<bcal_index, CellHits> &SiPMHits, map<int, SumHits> &bcalfADC, double Resolution)
{
   /// Loop over the CellHits objects and copy pointers to them into SumHits objects.
   ///
   /// For the BCAL, multiple SiPMs are summed together. This routine gathers individual
   /// SiPM hits into single SumHits objects. Each SumHits represents a summed
   /// cell that is readout by an fADC channel. Since not every cell has signal in it, each
   /// SumHits object may not have as many input cells as SiPMs that will actually be
   /// contributing.
   
   // Loop over SiPMHits and copy a pointer to it to the correct SumHits
   // element in the bcalfADC container. The bcalfADC container is an STL map which is
   // convenient since it creates a new SumHits object for us if it doesn't exist,
   // but otherwise, returns a reference to the existing object.
   
   map<bcal_index, CellHits>::iterator iter = SiPMHits.begin();
   for(; iter!=SiPMHits.end(); iter++){
      
      // Get reference to SumHits object
      const bcal_index &idx = iter->first;
      int fADCId = dBCALGeom->fADCId( idx.module, idx.layer, idx.sector);
      SumHits &sumhits = bcalfADC[fADCId];
      
      // Add CellHits object to list in SumHits
      CellHits &cellhits = iter->second;
      sumhits.cellhits.push_back(&cellhits);
      
      // If this is the first cell added to the SumHits, assign its
      // values to the first elements of the data arrays. Otherwise,
      // test if the hit overlaps in time with an existing element.
      bool mergedUP = false;
      bool mergedDN = false;

      // Upstream
      if(cellhits.end == CellHits::kUp && cellhits.E != 0.0){
        if(sumhits.EUP.empty()){
          sumhits.EUP.push_back(cellhits.E);
          sumhits.tUP.push_back(cellhits.t);
        }else{
          for(int ii = 0; ii < (int)sumhits.EUP.size(); ii++){
            if(fabs(cellhits.t - sumhits.tUP[ii]) < Resolution){
              sumhits.EUP[ii] += cellhits.E;
              if(sumhits.tUP[ii] > cellhits.t) sumhits.tUP[ii] = cellhits.t; // Again, keep the earlier of the two times
              mergedUP = true;
              break;
            }
          }
          if (!mergedUP){
            sumhits.EUP.push_back(cellhits.E);
            sumhits.tUP.push_back(cellhits.t);
          }
        }
      }
      
      // Downstream
      if(cellhits.end == CellHits::kDown && cellhits.E != 0.0){
        if(sumhits.EDN.empty()){
          sumhits.EDN.push_back(cellhits.E);
          sumhits.tDN.push_back(cellhits.t);
        }else{
          for(int ii = 0; ii < (int)sumhits.EDN.size(); ii++){
            if(fabs(cellhits.t - sumhits.tDN[ii]) < Resolution){
              sumhits.EDN[ii] += cellhits.E;
              if(sumhits.tDN[ii] > cellhits.t) sumhits.tDN[ii] = cellhits.t; // Again, keep the earlier of the two times
              mergedDN = true;
              break;
            }
          }
          if (!mergedDN){
            sumhits.EDN.push_back(cellhits.E);
            sumhits.tDN.push_back(cellhits.t);
          }
        }
      }
   }
}

//-----------
// SimpleDarkHitsSmear
//-----------
void BCALSmearer::SimpleDarkHitsSmear(map<int, SumHits> &bcalfADC)
{
   /// Loop over the SumHits objects and add Electronic noise and
   /// Dark hits smearing.
   ///
   /// Take SumHits objects and add to their energy values a random
   /// energy as sampled from a Gaussian.  The Gaussian for each 
   /// BCAL layer is based on data taken in May of 2015.
   /// In future, data on a channel-by-channel basis will be implemented.

   if(bcal_config->NO_DARK_PULSES) return;
   
   double Esmeared = 0;
   double sigma = 0;
   
   double sigma1 = bcal_config->BCAL_LAYER1_SIGMA_SCALE*bcal_config->BCAL_MEV_PER_ADC_COUNT; 
   double sigma2 = bcal_config->BCAL_LAYER2_SIGMA_SCALE*bcal_config->BCAL_MEV_PER_ADC_COUNT; 
   double sigma3 = bcal_config->BCAL_LAYER3_SIGMA_SCALE*bcal_config->BCAL_MEV_PER_ADC_COUNT; 
   double sigma4 = bcal_config->BCAL_LAYER4_SIGMA_SCALE*bcal_config->BCAL_MEV_PER_ADC_COUNT; 

   // Loop over all fADC readout cells
   for(int imodule=1; imodule<=dBCALGeom->GetBCAL_Nmodules(); imodule++){

      int n_layers = dBCALGeom->GetBCAL_Nlayers();
      for(int fADC_lay=1; fADC_lay<=n_layers; fADC_lay++){
         if(fADC_lay == 1) 
         	sigma = sigma1;
         else if(fADC_lay == 2) 
         	sigma = sigma2;
         else if(fADC_lay == 3) 
         	sigma = sigma3;
         else if(fADC_lay == 4) 
         	sigma = sigma4;

		 // there's a constant number of sectors now...
         //int n_sectors = (fADC_lay <= dBCALGeom->GetBCAL_NInnerLayers())? dBCALGeom->GetBCAL_NInnerSectors() : dBCALGeom->GetBCAL_NOuterSectors();
		 int n_sectors = dBCALGeom->GetBCAL_Nsectors();
         for(int fADC_sec=1; fADC_sec<=n_sectors; fADC_sec++){

            // Use cellId(...) to convert fADC layer and sector into fADCId
            // (see dBCALGeom->fADCId)
            int fADCId = dBCALGeom->cellId(imodule, fADC_lay, fADC_sec);
            
            // Get SumHits object if it already exists or create new one 
            // if it doesn't.
            SumHits &sumhits = bcalfADC[fADCId];

            for(int ii = 0; ii < (int)sumhits.EUP.size(); ii++){
				Esmeared = gDRandom.Gaus(sumhits.EUP[ii],sigma);
				sumhits.EUP[ii] = Esmeared;
            }
            for(int ii = 0; ii < (int)sumhits.EDN.size(); ii++){
            	Esmeared = gDRandom.Gaus(sumhits.EDN[ii],sigma);
            	sumhits.EDN[ii] = Esmeared;
            }
         }
      }
   }
}

//-----------
// ApplyTimeSmearing
//-----------
void BCALSmearer::ApplyTimeSmearing(double sigma_ns, double sigma_ns_TDC, map<int, fADCHitList> &fADCHits, map<int, TDCHitList> &TDCHits)
{
   /// The fADC250 will extract a time from the samples by applying an algorithm
   /// to a few of the samples taken every 4ns. The perfect times from HDGeant
   /// must be smeared to reflect the timing resolution of the fADC250.
   /// The F1TDC250 does something similar, but with ???ns samples.

   if(bcal_config->NO_T_SMEAR) return;

   // This is hardwired. Perhaps a future warrior would like to make it somehow work with CCDB?
   double BCAL_TIMINGADCCOEFA = 0.055;
   double BCAL_TIMINGADCCOEFB = 0.000;

   map<int, fADCHitList>::iterator it = fADCHits.begin();
   for(; it!=fADCHits.end(); it++){
      fADCHitList &hitlist = it->second;
      
      // upstream
      for(unsigned int i=0; i<hitlist.uphits.size(); i++){
         double EGeV = hitlist.uphits[i].E/1000;
         double sqrtterm = BCAL_TIMINGADCCOEFA / sqrt(EGeV);
         double linterm = BCAL_TIMINGADCCOEFB;
         double sigma_ns_ADC = sqrt(sqrtterm*sqrtterm + linterm*linterm);
         hitlist.uphits[i].t += gDRandom.SampleGaussian(sigma_ns_ADC);
      }

      // downstream
      for(unsigned int i=0; i<hitlist.dnhits.size(); i++){
         double EGeV = hitlist.dnhits[i].E/1000;
         double sqrtterm = BCAL_TIMINGADCCOEFA / sqrt(EGeV);
         double linterm = BCAL_TIMINGADCCOEFB;
         double sigma_ns_ADC = sqrt(sqrtterm*sqrtterm + linterm*linterm);
         hitlist.dnhits[i].t += gDRandom.SampleGaussian(sigma_ns_ADC);
      }
   }

   map<int, TDCHitList>::iterator itTDC = TDCHits.begin();
   for(; itTDC!=TDCHits.end(); itTDC++){
      TDCHitList &TDChitlist = itTDC->second;
      
      // upstream
      for(unsigned int i=0; i<TDChitlist.uphits.size(); i++){
         TDChitlist.uphits[i] += gDRandom.SampleGaussian(sigma_ns_TDC);
      }

      // downstream
      for(unsigned int i=0; i<TDChitlist.dnhits.size(); i++){
         TDChitlist.dnhits[i] += gDRandom.SampleGaussian(sigma_ns_TDC);
      }
   }
}

//-----------
// FindHits
//-----------
void BCALSmearer::FindHits(double thresh_MeV, map<int, SumHits> &bcalfADC, map<int, fADCHitList> &fADCHits, map<int,TDCHitList> &TDCHits)
{
   /// Loop over Sumhits objects and find hits that cross the energy threshold (ADC)
   map<int, SumHits>::iterator iter = bcalfADC.begin();
   for(; iter!=bcalfADC.end(); iter++){
      
      int fADCId = iter->first;
      SumHits &sumhits = iter->second;

      vector<fADCHit> uphits;
      vector<fADCHit> dnhits;

      vector<double> uphitsTDC;
      vector<double> dnhitsTDC;

      // The histogram should have the signal size for the ADC, but the TDC
      // leg will actually have a larger size since the pre-amp gain will be
      // set differently. Scale the threshold down here to accomodate this.
      double preamp_gain_tdc = 5.0;
      double thresh_MeV_TDC = thresh_MeV/preamp_gain_tdc;
	  //the outermost layer of the detector is not equipped with TDCs, so don't generate any TDC hits
	  int layer = dBCALGeom->layer(fADCId);

      for(int ii = 0; ii < (int)sumhits.EUP.size(); ii++){
        // correct simulation efficiencies 
		if (config->APPLY_EFFICIENCY_CORRECTIONS
		 		&& !gDRandom.DecideToAcceptHit(bcal_config->GetEfficiencyCorrectionFactor(GetCalibIndex(dBCALGeom->module(fADCId),
		 																				  dBCALGeom->layer(fADCId),
		 																				  dBCALGeom->sector(fADCId)),
		 																				  DBCALGeometry::End::kUpstream)))
		 			continue;
		 			
        if(sumhits.EUP[ii] > thresh_MeV && sumhits.tUP[ii] < 2000) uphits.push_back(fADCHit(sumhits.EUP[ii],sumhits.tUP[ii])); // Fill uphits and dnhits with energies (in MeV)
        if(layer != 4 && sumhits.EUP[ii] > thresh_MeV_TDC && sumhits.tUP[ii] < 2000) uphitsTDC.push_back(sumhits.tUP[ii]);     // and times when they cross an energy threshold.
      }                                                                                                                        // Also fill TDC uphits and dnhits with times if
      for(int ii = 0; ii < (int)sumhits.EDN.size(); ii++){                                                                     // they are not layer 4 hits and cross threshold.
        // correct simulation efficiencies 
		if (config->APPLY_EFFICIENCY_CORRECTIONS
		 		&& !gDRandom.DecideToAcceptHit(bcal_config->GetEfficiencyCorrectionFactor(GetCalibIndex(dBCALGeom->module(fADCId),
		 																				  dBCALGeom->layer(fADCId),
		 																				  dBCALGeom->sector(fADCId)),
		 																				  DBCALGeometry::End::kDownstream)))
		 			continue;

        if(sumhits.EDN[ii] > thresh_MeV && sumhits.tDN[ii] < 2000) dnhits.push_back(fADCHit(sumhits.EDN[ii],sumhits.tDN[ii]));
        if(layer != 4 && sumhits.EDN[ii] > thresh_MeV_TDC && sumhits.tDN[ii] < 2000) dnhitsTDC.push_back(sumhits.tDN[ii]);
      }
      
      // If at least one ADC readout channel has a hit, add the readout cell to fADCHits
      if(uphits.size()>0 || dnhits.size()>0){
         fADCHitList &hitlist = fADCHits[fADCId];

         // The module, fADC layer, and fADC sector are encoded in fADCId
         // (n.b. yes, these are the same methods used for extracting 
         // similar quantities from the cellId.)
         hitlist.module = dBCALGeom->module(fADCId);
         hitlist.sumlayer = dBCALGeom->layer(fADCId);
         hitlist.sumsector = dBCALGeom->sector(fADCId);
         
         hitlist.uphits = uphits;
         hitlist.dnhits = dnhits;
      }
      
      // If at least one TDC readout channel has a hit, add the readout cell to TDCHits
      if(uphitsTDC.size()>0 || dnhitsTDC.size()>0){
         TDCHitList &hitlistTDC = TDCHits[fADCId];

         // The module, fADC layer, and fADC sector are encoded in fADCId
         // (n.b. yes, these are the same methods used for extracting 
         // similar quantities from the cellId.)
         hitlistTDC.module = dBCALGeom->module(fADCId);
         hitlistTDC.sumlayer = dBCALGeom->layer(fADCId);
         hitlistTDC.sumsector = dBCALGeom->sector(fADCId);
         
         hitlistTDC.uphits = uphitsTDC;
         hitlistTDC.dnhits = dnhitsTDC;
      }
   }
}

//-----------
// CopyBCALHitsToHDDM
//-----------
void BCALSmearer::CopyBCALHitsToHDDM(map<int, fADCHitList> &fADCHits,
						map<int, TDCHitList> &TDCHits,
                        hddm_s::HDDM *record)
{
   /// Loop over fADCHitList objects and copy the fADC hits into the HDDM tree.
   ///
   /// This will copy all of the hits found into the first physicsEvent found
   /// in the HDDM file. Note that the hits were formed from data that may
   /// have been combined from several physicsEvent structures in the HDDM
   /// event. No attempt is made to keep track of this so all hits are thrown
   /// into only a single physicsEvent.

   hddm_s::BarrelEMcalList bcals = record->getBarrelEMcals();
   if (bcals.size() == 0){
      if(record->getHitViews().empty()){
		record->getPhysicsEvent().addHitViews();
	  }
      bcals = record->getHitViews().begin()->addBarrelEMcals();
   }
   hddm_s::BcalCellList cells = bcals().getBcalCells();
   hddm_s::BcalCellList::iterator iter;
   for (iter = cells.begin(); iter != cells.end(); ++iter) {

      // Delete any existing bcalfADCDigiHit and bcalTDCDigiHit structures
       iter->deleteBcalfADCDigiHits();
       iter->deleteBcalTDCDigiHits();
   }

   // If we have no cells over threshold, then bail now.
   if (fADCHits.size() == 0 && TDCHits.size() == 0)
      return;
   
   // Create bcalfADCHit structures to hold our fADC hits
   map<int, fADCHitList>::iterator it;
   for (it = fADCHits.begin(); it != fADCHits.end(); ++it) {
      // Get pointer to our fADC cell information that needs to be copied to HDDM
      fADCHitList &hitlist = it->second;
      // Check if this cell is already present in the cells list
      cells = bcals().getBcalCells();
      for (iter = cells.begin(); iter != cells.end(); ++iter) {
         if (iter->getModule() == it->second.module &&
             iter->getSector() == it->second.sumsector &&
             iter->getLayer() == it->second.sumlayer)
         {
            break;
         }
      }
      if (iter == cells.end()) {
         iter = bcals().addBcalCells().begin();
         iter->setModule(hitlist.module);
         iter->setLayer(hitlist.sumlayer);
         iter->setSector(hitlist.sumsector);
      }
      
      // Copy hits into BcalfADCDigiHit HDDM structure.
      // Energies and times must be converted to units of ADC counts.
      // Because we use unsigned integers, times must be positive.  HDGEANT can output negative times,
      // so we can offset the times now to ensure they are positive before the conversion, then
      // fix the offset layer in the hit factories.  Also, any hit that still has a negative time
      // will be ignored.
      for (unsigned int i = 0; i < hitlist.uphits.size(); i++) {
      	int integer_time = round((hitlist.uphits[i].t-bcal_config->BCAL_BASE_TIME_OFFSET)/bcal_config->BCAL_NS_PER_ADC_COUNT);
      	if (integer_time >= 0){
            hddm_s::BcalfADCDigiHitList fadcs = iter->addBcalfADCDigiHits();
            fadcs().setEnd(bcal_index::kUp);
	    double integral = round(hitlist.uphits[i].E/bcal_config->BCAL_MEV_PER_ADC_COUNT);
	    
	    // fADC saturation based on waveforms from data
	    if(!bcal_config->NO_FADC_SATURATION) { 
		    if(integral > bcal_config->fADC_MinIntegral_Saturation[0][hitlist.sumlayer-1]) {
			    double y = integral; 
			    double a = bcal_config->fADC_Saturation_Linear[0][hitlist.sumlayer-1];
			    double b = bcal_config->fADC_Saturation_Quadratic[0][hitlist.sumlayer-1];
			    double c = bcal_config->fADC_MinIntegral_Saturation[0][hitlist.sumlayer-1];
			    // "invert" saturation correction for MC
			    integral = (1 - a*y + 2.*b*c*y - sqrt(1. - 2.*a*y + 4.*b*c*y + (a*a - 4.*b)*y*y))/(2.*b*y);
		    }
	    }
            fadcs().setPulse_integral(integral);
            fadcs().setPulse_time(integer_time);
        }
      }
      for (unsigned int i = 0; i < hitlist.dnhits.size(); i++) {
      	int integer_time = round((hitlist.dnhits[i].t-bcal_config->BCAL_BASE_TIME_OFFSET)/bcal_config->BCAL_NS_PER_ADC_COUNT);
      	if (integer_time >= 0){
            hddm_s::BcalfADCDigiHitList fadcs = iter->addBcalfADCDigiHits();
            fadcs().setEnd(bcal_index::kDown);
	    double integral = round(hitlist.dnhits[i].E/bcal_config->BCAL_MEV_PER_ADC_COUNT);
	    
	    // fADC saturation based on waveforms from data
	    if(!bcal_config->NO_FADC_SATURATION) { 
		    if(integral > bcal_config->fADC_MinIntegral_Saturation[1][hitlist.sumlayer-1]) {
			    double y = integral; 
			    double a = bcal_config->fADC_Saturation_Linear[1][hitlist.sumlayer-1];
			    double b = bcal_config->fADC_Saturation_Quadratic[1][hitlist.sumlayer-1];
			    double c = bcal_config->fADC_MinIntegral_Saturation[1][hitlist.sumlayer-1];
			    // "invert" saturation correction for MC
			    integral = (1 - a*y + 2.*b*c*y - sqrt(1. - 2.*a*y + 4.*b*c*y + (a*a - 4.*b)*y*y))/(2.*b*y);
                    }

	    }
            fadcs().setPulse_integral(integral);
            fadcs().setPulse_time(integer_time);
        } 
      }
   }
   
   // Create bcalTDCDigiHit structures to hold our F1TDC hits
   map<int, TDCHitList>::iterator ittdc;
   for (ittdc = TDCHits.begin(); ittdc != TDCHits.end(); ittdc++) {
      // Get pointer to our TDC hit information that needs to be copied to HDDM
      TDCHitList &hitlist = ittdc->second;
      // Check if this cell is already present in the cells list
      cells = bcals().getBcalCells();
      for (iter = cells.begin(); iter != cells.end(); ++iter) {
         if (iter->getModule() == ittdc->second.module &&
             iter->getSector() == ittdc->second.sumsector &&
             iter->getLayer() == ittdc->second.sumlayer)
         {
            break;
         }
      }
      if (iter == cells.end()) {
         iter = bcals().addBcalCells().begin();
         iter->setModule(hitlist.module);
         iter->setLayer(hitlist.sumlayer);
         iter->setSector(hitlist.sumsector);
      }
      
      // Copy hits into BcalTDCDigiHit HDDM structure.
      // Times must be converted to units of TDC counts.
      for (unsigned int i = 0; i < hitlist.uphits.size(); ++i) {
         int integer_time = round((hitlist.uphits[i]-bcal_config->BCAL_TDC_BASE_TIME_OFFSET)/bcal_config->BCAL_NS_PER_TDC_COUNT);
         if (integer_time >= 0){
            hddm_s::BcalTDCDigiHitList tdcs = iter->addBcalTDCDigiHits();
            tdcs().setEnd(bcal_index::kUp);
            tdcs().setTime(integer_time);
         }
      }
      for (unsigned int i = 0; i < hitlist.dnhits.size(); i++) {
         int integer_time = round((hitlist.dnhits[i]-bcal_config->BCAL_TDC_BASE_TIME_OFFSET)/bcal_config->BCAL_NS_PER_TDC_COUNT);
         if (integer_time >= 0){
            hddm_s::BcalTDCDigiHitList tdcs = iter->addBcalTDCDigiHits();
            tdcs().setEnd(bcal_index::kDown);
            tdcs().setTime(round((hitlist.dnhits[i]-bcal_config->BCAL_TDC_BASE_TIME_OFFSET)/bcal_config->BCAL_NS_PER_TDC_COUNT));
         }
      }
   }
}


//-----------
// bcal_config_t  (constructor)
//-----------
bcal_config_t::bcal_config_t(JEventLoop *loop) 
{
 	BCAL_SAMPLINGCOEFA        = 0.0; // 0.042 (from calibDB BCAL/bcal_parms)
 	BCAL_SAMPLINGCOEFB        = 0.0; // 0.013 (from calibDB BCAL/bcal_parms)
 	BCAL_TWO_HIT_RESO        = 0.0; // 50. (from calibDB BCAL/bcal_parms)
 	BCAL_mevPerPE             = 0.0; // Energy corresponding to one pixel firing in MeV - FIX 
	BCAL_C_EFFECTIVE          = 0.0; // constant effective velocity, assumed to be property of fibers

	BCAL_LAYER1_SIGMA_SCALE	  = 0.0; // Approximated from https://logbooks.jlab.org/entry/3339692 (10 degree, 1.4 V OB pedestal data)
	BCAL_LAYER2_SIGMA_SCALE	  = 0.0; // Approximated from https://logbooks.jlab.org/entry/3339692 (10 degree, 1.4 V OB pedestal data)
	BCAL_LAYER3_SIGMA_SCALE	  = 0.0; // Approximated from https://logbooks.jlab.org/entry/3339692 (10 degree, 1.4 V OB pedestal data)
	BCAL_LAYER4_SIGMA_SCALE	  = 0.0; // Approximated from https://logbooks.jlab.org/entry/3339692 (10 degree, 1.4 V OB pedestal data)
   // Values from logbook entry are in units of integrated ADC counts.  Average SiPM gain ~ 0.029 MeV per integrated ADC count. 

	// FIX - Pull from geometry?
 	BCAL_NUM_MODULES = 48;
 	BCAL_NUM_LAYERS = 4;
 	BCAL_NUM_SECTORS = 4;

 	BCAL_BASE_TIME_OFFSET     = 0; // -100.0 (from calibDB BCAL/base_time_offset)
 	BCAL_TDC_BASE_TIME_OFFSET = 0; // -100.0 (from calibDB BCAL/base_time_offset)

 	BCAL_ADC_THRESHOLD_MEV    = 0.0;  // MeV (To be updated/improved)
 	BCAL_FADC_TIME_RESOLUTION = 0.0;  // ns (To be updated/improved)
 	BCAL_TDC_TIME_RESOLUTION  = 0.0;  // ns (To be updated/improved)
 	BCAL_MEV_PER_ADC_COUNT    = 0.0;  // MeV per integrated ADC count (based on Spring 2015 calibrations)
 	BCAL_NS_PER_ADC_COUNT     = 0.0;  // 0.0625 ns per ADC count (from calibDB BCAL/digi_scales)
 	BCAL_NS_PER_TDC_COUNT     = 0.0;  // 0.0559 ns per TDC count (from calibDB BCAL/digi_scales)

	// BCAL flags
 	NO_T_SMEAR = false;
 	NO_DARK_PULSES = false;
 	NO_SAMPLING_FLUCTUATIONS = false;
 	NO_SAMPLING_FLOOR_TERM = false;
 	NO_POISSON_STATISTICS = false;
	NO_FADC_SATURATION = false;
		
	// Load parameters from CCDB
    cout << "get BCAL/bcal_smear_parms_v2 parameters from CCDB..." << endl;
    map<string, double> bcalparms;
    if(loop->GetCalib("BCAL/bcal_smear_parms_v2", bcalparms)) {
     	jerr << "Problem loading BCAL/bcal_smear_parms_v2 from CCDB!" << endl;
     } else {
     	BCAL_SAMPLINGCOEFA 		  = bcalparms["BCAL_SAMPLINGCOEFA"];
     	BCAL_SAMPLINGCOEFB 		  = bcalparms["BCAL_SAMPLINGCOEFB"];
     	BCAL_TWO_HIT_RESO 		  = bcalparms["BCAL_TWO_HIT_RESO"];
	BCAL_mevPerPE			  = bcalparms["BCAL_mevPerPE"];
	BCAL_C_EFFECTIVE		  = bcalparms["BCAL_C_EFFECTIVE"];
	BCAL_ATTENUATION_LENGTH		  = bcalparms["BCAL_ATTENUATION_LENGTH"];
	BCAL_ADC_THRESHOLD_MEV		  = bcalparms["BCAL_ADC_THRESHOLD_MEV"];
	BCAL_FADC_TIME_RESOLUTION	  = bcalparms["BCAL_FADC_TIME_RESOLUTION"]; 
	BCAL_TDC_TIME_RESOLUTION	  = bcalparms["BCAL_TDC_TIME_RESOLUTION"];
	BCAL_MEV_PER_ADC_COUNT 		  = bcalparms["BCAL_MEV_PER_ADC_COUNT"];
	BCAL_LAYER1_SIGMA_SCALE		  = bcalparms["BCAL_LAYER1_SIGMA_SCALE"];
	BCAL_LAYER2_SIGMA_SCALE		  = bcalparms["BCAL_LAYER2_SIGMA_SCALE"];
	BCAL_LAYER3_SIGMA_SCALE		  = bcalparms["BCAL_LAYER3_SIGMA_SCALE"];
	BCAL_LAYER4_SIGMA_SCALE		  = bcalparms["BCAL_LAYER4_SIGMA_SCALE"];

	}
	
    //cout << "Get BCAL/attenuation_parameters from CCDB..." <<endl;
    //vector< vector<double> > in_atten_parameters;
    //if(loop->GetCalib("BCAL/attenuation_parameters", in_atten_parameters)) {
    // 	jerr << "Problem loading BCAL/bcal_parms from CCDB!" << endl;
	//} else {
    // 	attenuation_parameters.clear();
    //
    // 	for (unsigned int i = 0; i < in_atten_parameters.size(); i++) {
    // 		attenuation_parameters.push_back( in_atten_parameters.at(i) );
    // 	}
	//}
		
     cout << "Get BCAL/digi_scales parameters from CCDB..." << endl;
     map<string, double> bcaldigiscales;
     if(loop->GetCalib("BCAL/digi_scales", bcaldigiscales)) {
     	jerr << "Problem loading BCAL/digi_scales from CCDB!" << endl;
     } else {
     	BCAL_NS_PER_ADC_COUNT = bcaldigiscales["BCAL_ADC_TSCALE"];
     	BCAL_NS_PER_TDC_COUNT = bcaldigiscales["BCAL_TDC_SCALE"];
   	}

    cout << "Get BCAL/base_time_offset parameters from CCDB..." << endl;
    map<string, double> bcaltimeoffsets;
    if(loop->GetCalib("BCAL/base_time_offset", bcaltimeoffsets)) {
     	jerr << "Problem loading BCAL/base_time_offset from CCDB!" << endl;
 	} else {
     	BCAL_BASE_TIME_OFFSET = bcaltimeoffsets["BCAL_BASE_TIME_OFFSET"];
    	BCAL_TDC_BASE_TIME_OFFSET = bcaltimeoffsets["BCAL_TDC_BASE_TIME_OFFSET"];
   	}
   	
   	// load per-channel efficiencies
	vector<double> raw_table;
	if(loop->GetCalib("BCAL/channel_mc_efficiency", raw_table)) {
    	jerr << "Problem loading BCAL/channel_mc_efficiency from CCDB!" << endl;
    } else {
   	    int channel = 0;

    	for (int module=1; module<=BCAL_NUM_MODULES; module++) {
        	for (int layer=1; layer<=BCAL_NUM_LAYERS; layer++) {
            	for (int sector=1; sector<=BCAL_NUM_SECTORS; sector++) {
	                channel_efficiencies.push_back( pair<double,double>(raw_table[channel],raw_table[channel+1]) );

                	channel += 2;
                }
            }
        }
        
    }

    std::vector<std::map<string,double> > saturation_ADC_pars;
    if(loop->GetCalib("/BCAL/ADC_saturation", saturation_ADC_pars))
	    jout << "Error loading /BCAL/ADC_saturation !" << endl;
    for (unsigned int i=0; i < saturation_ADC_pars.size(); i++) {
	    int end = (saturation_ADC_pars[i])["end"];
	    int layer = (saturation_ADC_pars[i])["layer"] - 1;
	    fADC_MinIntegral_Saturation[end][layer] = (saturation_ADC_pars[i])["par0"];
	    fADC_Saturation_Linear[end][layer] = (saturation_ADC_pars[i])["par1"];
	    fADC_Saturation_Quadratic[end][layer] = (saturation_ADC_pars[i])["par2"];
    } 
}


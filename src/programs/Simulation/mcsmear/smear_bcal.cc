// $Id: smear.cc 7650 2011-03-29 22:52:30Z shepherd $
//
// Created June 22, 2005  David Lawrence
//
// Major revision March 6, 2012 David Lawrence


#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <sstream>
#include <queue>
#include <cmath>
using namespace std;

#include <DHistogram.h>
#include <BCAL/DBCALGeometry.h>

#include "units.h"
#include <HDDM/hddm_s.hpp>
#include <TMath.h>
#include <TH2D.h>
#include <TSpline.h>
#include <TDirectory.h>

#include "DRandom2.h"

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

TH1D *hNincident_particles = NULL;


//..........................
// bcal_index is a utility class that encapsulates the
// module, layer, sector, and end in a single object that
// can be used as a key to index an STL map. 
//..........................
class bcal_index{
   public:
      enum EndType{
         kUp,
         kDown
      };
      
      bcal_index(unsigned int module, unsigned int layer,
                 unsigned int sector, unsigned int incident_id,
                 EndType end)
       : module(module),
         layer(layer),
         sector(sector),
         incident_id(incident_id),
         end(end)
      {}
   
      unsigned int module;
      unsigned int layer;
      unsigned int sector;
      unsigned int incident_id;
      EndType end;
      
      bool operator<(const bcal_index &idx) const {
         if (module < idx.module)
            return true;
         if (module > idx.module)
            return false;
         if (layer < idx.layer)
            return true;
         if (layer > idx.layer)
            return false;
         if (sector < idx.sector)
            return true;
         if (sector > idx.sector)
            return false;
         if (incident_id < idx.incident_id)
            return true;
         if (incident_id > idx.incident_id)
            return false;
         if ((end==kUp) && (idx.end==kDown))
            return true;
         return false;
      }
};

//..........................
// CellHits is a utility class that holds information
// regarding the energy and time of depostions in a cell
//..........................
class CellHits{
   public:
      enum EndType{
         kUp,
         kDown
      };

      CellHits() : E(0.0), t(0.0)
      {}
      
      double E;
      double t;
      double Etruth;
      EndType end;
};

//..........................
// SumHits is a utility class that is used to hold info
// from the SiPMs contributing to that readout channel.
// This includes a list of CellHits objects, but also
// the total number of SiPMs that should be in the sum
// and the total up/downstream energies and times.
//..........................
class SumHits{
   public:
      SumHits()
      {}
      
      vector<CellHits *> cellhits;
      vector<double> EUP;
      vector<double> tUP;
      vector<double> EDN;
      vector<double> tDN;
};

//..........................
// fADCHit is a utility class that is used to hold info
// for a single fADC hit. 
//..........................
class fADCHit{
   public:
      fADCHit(double E, double t) : E(E), t(t)
      {}
      
      double E;
      double t;
};

//..........................
// fADCHitList is a utility class that is used to hold info
// for a set of fADCHit objects. 
//..........................
class fADCHitList{
   public:
      fADCHitList()
      {}
      
      int module;
      int sumlayer;
      int sumsector;
      
      vector<fADCHit> uphits;
      vector<fADCHit> dnhits;
};

//..........................
// TDCHitList is a utility class that is used to hold info
// for a single F1TDC hit
//..........................
class TDCHitList{
   public:
      TDCHitList()
      {}
      
      int module;
      int sumlayer;
      int sumsector;
      
      vector<double> uphits;
      vector<double> dnhits;
};

//..........................
// IncidentParticle_t is a utility class for holding the
// parameters of particles recorded as incident on the 
// BCAL (shower causing)
//..........................
class IncidentParticle_t{
   public:
      IncidentParticle_t(hddm_s::BcalTruthIncidentParticle &ipart) {
         x = ipart.getX();
         y = ipart.getY();
         z = ipart.getZ();
         px = ipart.getPx();
         py = ipart.getPy();
         pz = ipart.getPz();
         ptype = ipart.getPtype();
      }

      float x,y,z;
      float px, py, pz;
      int ptype, track;
};


// Defined in this file
int32_t GetRunNumber(hddm_s::HDDM *record);
int GetCalibIndex(int module, int layer, int sector);
void GetAttenuationParameters(int id, double &attenuation_length, double &attenuation_L1, double &attenuation_L2);
double GetEffectiveVelocity(int id);
void GetSiPMHits(hddm_s::HDDM *record,
                    map<bcal_index, CellHits> &SiPMHits,
                    vector<IncidentParticle_t> &incident_particles);
void ApplySamplingFluctuations(map<bcal_index,
                               CellHits> &SiPMHits,
                               vector<IncidentParticle_t> &incident_particles);
void MergeHits(map<bcal_index, CellHits> &SiPMHits, double Resolution);
void ApplyPoissonStatistics(map<bcal_index, CellHits> &SiPMHits);
void SortSiPMHits(map<bcal_index,
                     CellHits> &SiPMHits,
                     map<int, SumHits> &bcalfADC, double Resolution);
void SimpleDarkHitsSmear(map<int, SumHits> &bcalfADC);
void ApplyTimeSmearing(double sigma_ns, double sigma_ns_TDC, map<int, fADCHitList> &fADCHits, map<int, TDCHitList> &TDCHits);
void FindHits(double thresh_MeV,
              map<int, SumHits> &bcalfADC,
              map<int, fADCHitList> &fADCHits,
              map<int, TDCHitList> &TDCHits);
void CopyBCALHitsToHDDM(map<int, fADCHitList> &fADCHits,
                        map<int, TDCHitList> &TDCHits,
                        hddm_s::HDDM *record);

// Mutex used to control accessing the ROOT global memory
extern pthread_mutex_t root_mutex;

// Flag used specifically for BCAL
extern bool SMEAR_BCAL;

// The following are all false by default, but can be
// set to true via command line parameters. Setting
// one of these to true will turn OFF the feature.
extern bool NO_E_SMEAR;
extern bool NO_T_SMEAR;
extern bool NO_DARK_PULSES;
extern bool NO_SAMPLING_FLUCTUATIONS;
extern bool NO_SAMPLING_FLOOR_TERM;
extern bool NO_POISSON_STATISTICS;
extern bool NO_THRESHOLD_CUT;

extern double BCAL_FADC_TIME_RESOLUTION; // ns
extern double BCAL_TDC_TIME_RESOLUTION; // ns
extern double BCAL_ADC_THRESHOLD_MEV; // MeV

// setup response parameters
extern double BCAL_SAMPLINGCOEFA;               // 0.042 (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLINGCOEFB;               // 0.013 (from calibDB BCAL/bcal_parms)
extern double BCAL_TWO_HIT_RESOL;               // 50. (from calibDB BCAL/bcal_parms)
extern double BCAL_mevPerPE;                    // (defined below)
extern double BCAL_NS_PER_ADC_COUNT;            // 0.0625 (defined in mcsmear.cc)
extern double BCAL_NS_PER_TDC_COUNT;            // 0.0559 (defined in mcsmear.cc)
extern double BCAL_MEV_PER_ADC_COUNT;           // 0.029 (defined in mcsmear.cc)

extern vector<vector<double> > attenuation_parameters; // Avg. of 525 (from calibDB BCAL/attenuation_parameters)
extern vector<double> effective_velocities;     // 16.75 (from calibDB BCAL/effective_velocities)

extern int BCAL_NUM_MODULES;
extern int BCAL_NUM_LAYERS;
extern int BCAL_NUM_SECTORS;

extern double BCAL_BASE_TIME_OFFSET;            // -100.0 (from calibDB BCAL/base_time_offset)
extern double BCAL_TDC_BASE_TIME_OFFSET;        // -100.0 (from calibDB BCAL/base_time_offset)

// The following are not currently in use
//extern double BCAL_DARKRATE_GHZ;                // 0.0176 (from calibDB BCAL/bcal_parms) for 4x4 array
//extern double BCAL_DEVICEPDE;                   // 0.21   (from calibDB BCAL/bcal_parms)
//extern double BCAL_SAMPLING_FRACT;              // 0.095  (from calibDB BCAL/bcal_parms)
//extern double BCAL_PHOTONSPERSIDEPERMEV_INFIBER;// 75 (from calibDB BCAL/bcal_parms) 
//extern double BCAL_SIGMA_SIG_RELATIVE;          // 0.105  (from calibDB BCAL/bcal_parms)
//extern double BCAL_SIGMA_PED_RELATIVE;          // 0.139  (from calibDB BCAL/bcal_parms)
//extern double BCAL_SIPM_GAIN_VARIATION;         // 0.04   (from calibDB BCAL/bcal_parms)
//extern double BCAL_INTWINDOW_NS;                // 100    (from calibDB BCAL/bcal_parms)
//extern double BCAL_AVG_DARK_DIGI_VALS_PER_EVENT;// 240 used to set thresholds
//extern double BCAL_TIMEDIFFCOEFA;               // 0.07 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
//extern double BCAL_TIMEDIFFCOEFB;               // 0.00 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)

//-----------
// SmearBCAL
//-----------
void SmearBCAL(hddm_s::HDDM *record)
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
    ApplySamplingFluctuations(SiPMHits, incident_particles);

    // Merge hits associated with different incident particles
    MergeHits(SiPMHits, BCAL_TWO_HIT_RESOL);

    // Poisson Statistics
    ApplyPoissonStatistics(SiPMHits);
   
    // Place all hit cells into list indexed by fADC ID
    map<int, SumHits> bcalfADC;
    SortSiPMHits(SiPMHits, bcalfADC, BCAL_TWO_HIT_RESOL);

    // Electronic noise/Dark hits Smearing
    SimpleDarkHitsSmear(bcalfADC);
   
    // Apply energy threshold to dismiss low-energy hits
    map<int, fADCHitList> fADCHits;
    map<int, TDCHitList> TDCHits;
    FindHits(BCAL_ADC_THRESHOLD_MEV, bcalfADC, fADCHits, TDCHits);

    // Apply time smearing to emulate the fADC resolution
    ApplyTimeSmearing(BCAL_FADC_TIME_RESOLUTION, BCAL_TDC_TIME_RESOLUTION, fADCHits, TDCHits);
   
    // Copy hits into HDDM tree
    CopyBCALHitsToHDDM(fADCHits, TDCHits, record);
   
    bcalfADC.clear();
}

int GetCalibIndex(int module, int layer, int sector)
{
   return BCAL_NUM_LAYERS*BCAL_NUM_SECTORS*(module-1) + BCAL_NUM_SECTORS*(layer-1) + (sector-1);
}

void GetAttenuationParameters(int id, double &attenuation_length, double &attenuation_L1, double &attenuation_L2)
{
   vector<double> &parms = attenuation_parameters.at(id);

   attenuation_length = parms[0];
   attenuation_L1 = parms[1];
   attenuation_L2 = parms[2];
}

double GetEffectiveVelocity(int id)
{
   return effective_velocities.at(id);
}

//-----------
// GetSiPMHits
//-----------
void GetSiPMHits(hddm_s::HDDM *record,
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
     }
     else if (iter->getLayer() == 2 || iter->getLayer() == 3){
       layer = 2;
     }
     else if (iter->getLayer() == 4 || iter->getLayer() == 5 || iter->getLayer() == 6){
       layer = 3;
     }
     else layer = 4;

     int table_id = GetCalibIndex( iter->getModule(), layer, iter->getSector() );  // key the cell identification off of the upstream cell

     double cEff = GetEffectiveVelocity(table_id);
     double attenuation_length = 0; // initialize variable
     double attenuation_L1=-1., attenuation_L2=-1.;  // these parameters are ignored for now
     GetAttenuationParameters(table_id, attenuation_length, attenuation_L1, attenuation_L2);

     // Get reference to existing CellHits, or create one if it doesn't exist
     CellHits &cellhitsup = SiPMHits[idxup];
     cellhitsup.Etruth = iter->getE(); // Energy deposited in the cell in GeV
     cellhitsup.E = iter->getE()*exp(-dist_up/attenuation_length)*1000.; // in attenuated MeV
     cellhitsup.t = iter->getT() + dist_up/cEff; // in ns
     cellhitsup.end = CellHits::kUp; // Keep track of BCal end

     // Get reference to existing CellHits, or create one if it doesn't exist
     CellHits &cellhitsdn = SiPMHits[idxdn];
     cellhitsdn.Etruth = iter->getE(); // Energy deposited in the cell in GeV
     cellhitsdn.E = iter->getE()*exp(-dist_dn/attenuation_length)*1000.; // in attenuated MeV
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
   
   if (hNincident_particles)
      hNincident_particles->Fill(incident_particles.size());
}

//-----------
// ApplySamplingFluctuations
//-----------
void ApplySamplingFluctuations(map<bcal_index, CellHits> &SiPMHits, vector<IncidentParticle_t> &incident_particles)
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
   
   if(NO_SAMPLING_FLUCTUATIONS)return;
   if(NO_SAMPLING_FLOOR_TERM)BCAL_SAMPLINGCOEFB=0.0; // (redundant, yes, but located in more obvious place here)

   map<bcal_index, CellHits>::iterator iter=SiPMHits.begin();
   for(; iter!=SiPMHits.end(); iter++){
      CellHits &cellhits = iter->second;
      
      // Find fractional sampling sigma based on deposited energy (whole colorimeter, not just fibers)
      double Etruth = cellhits.Etruth;
      double sigmaSamp = BCAL_SAMPLINGCOEFA / sqrt( Etruth ) + BCAL_SAMPLINGCOEFB;

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
void MergeHits(map<bcal_index, CellHits> &SiPMHits, double Resolution)
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
void ApplyPoissonStatistics(map<bcal_index, CellHits> &SiPMHits)
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
   // integral (in PE) is used to scale the energy.

   if(NO_POISSON_STATISTICS)return;

   map<bcal_index, CellHits>::iterator iter=SiPMHits.begin();
   for(; iter!=SiPMHits.end(); iter++){
      CellHits &cellhits = iter->second;

      if(cellhits.E>0.0){
         // Convert to number of PE
         double mean_pe = cellhits.E/BCAL_mevPerPE;
         
         int Npe = gDRandom.Poisson(mean_pe);
         double ratio = (double)Npe/mean_pe;

         cellhits.E *= ratio;
      }
   }
}

//-----------
// SortSiPMHits
//-----------
void SortSiPMHits(map<bcal_index, CellHits> &SiPMHits, map<int, SumHits> &bcalfADC, double Resolution)
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
      int fADCId = DBCALGeometry::fADCId( idx.module, idx.layer, idx.sector);
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
void SimpleDarkHitsSmear(map<int, SumHits> &bcalfADC)
{
   /// Loop over the SumHits objects and add Electronic noise and
   /// Dark hits smearing.
   ///
   /// Take SumHits objects and add to their energy values a random
   /// energy as sampled from a Gaussian.  The Gaussian for each 
   /// BCAL layer is based on data taken in May of 2015.
   /// In future, data on a channel-by-channel basis will be implemented.

   if(NO_DARK_PULSES)return;
   
   double Esmeared = 0;
   double sigma = 0;
   double sigma1 = 43.*BCAL_MEV_PER_ADC_COUNT;  // Approximated from https://logbooks.jlab.org/entry/3339692 (10 degree, 1.4 V OB pedestal data)
   double sigma2 = 46.*BCAL_MEV_PER_ADC_COUNT;  // Approximated from https://logbooks.jlab.org/entry/3339692 (10 degree, 1.4 V OB pedestal data)
   double sigma3 = 49.*BCAL_MEV_PER_ADC_COUNT;  // Approximated from https://logbooks.jlab.org/entry/3339692 (10 degree, 1.4 V OB pedestal data)
   double sigma4 = 52.*BCAL_MEV_PER_ADC_COUNT;  // Approximated from https://logbooks.jlab.org/entry/3339692 (10 degree, 1.4 V OB pedestal data)
   // Values from logbook entry are in units of integrated ADC counts.  Average SiPM gain ~ 0.029 MeV per integrated ADC count.

   // Loop over all fADC readout cells
   for(int imodule=1; imodule<=DBCALGeometry::NBCALMODS; imodule++){

      int n_layers = DBCALGeometry::NBCALLAYSIN + DBCALGeometry::NBCALLAYSOUT;
      for(int fADC_lay=1; fADC_lay<=n_layers; fADC_lay++){

         if(fADC_lay == 1) sigma = sigma1;
         else if(fADC_lay == 2) sigma = sigma2;
         else if(fADC_lay == 3) sigma = sigma3;
         else if(fADC_lay == 4) sigma = sigma4;

         int n_sectors = (fADC_lay <= DBCALGeometry::NBCALLAYSIN)? DBCALGeometry::NBCALSECSIN : DBCALGeometry::NBCALSECSOUT;
         for(int fADC_sec=1; fADC_sec<=n_sectors; fADC_sec++){

            // Use cellId(...) to convert fADC layer and sector into fADCId
            // (see DBCALGeometry::fADCId)
            int fADCId = DBCALGeometry::cellId(imodule, fADC_lay, fADC_sec);
            
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
void ApplyTimeSmearing(double sigma_ns, double sigma_ns_TDC, map<int, fADCHitList> &fADCHits, map<int, TDCHitList> &TDCHits)
{
   /// The fADC250 will extract a time from the samples by applying an algorithm
   /// to a few of the samples taken every 4ns. The perfect times from HDGeant
   /// must be smeared to reflect the timing resolution of the fADC250.
   /// The F1TDC250 does something similar, but with ???ns samples.

   if(NO_T_SMEAR) return;

   map<int, fADCHitList>::iterator it = fADCHits.begin();
   for(; it!=fADCHits.end(); it++){
      fADCHitList &hitlist = it->second;
      
      // upstream
      for(unsigned int i=0; i<hitlist.uphits.size(); i++){
         hitlist.uphits[i].t += gDRandom.Gaus(sigma_ns);
      }

      // downstream
      for(unsigned int i=0; i<hitlist.dnhits.size(); i++){
         hitlist.dnhits[i].t += gDRandom.Gaus(sigma_ns);
      }
   }

   map<int, TDCHitList>::iterator itTDC = TDCHits.begin();
   for(; itTDC!=TDCHits.end(); itTDC++){
      TDCHitList &TDChitlist = itTDC->second;
      
      // upstream
      for(unsigned int i=0; i<TDChitlist.uphits.size(); i++){
         TDChitlist.uphits[i] += gDRandom.Gaus(sigma_ns_TDC);
      }

      // downstream
      for(unsigned int i=0; i<TDChitlist.dnhits.size(); i++){
         TDChitlist.dnhits[i] += gDRandom.Gaus(sigma_ns_TDC);
      }
   }
}

//-----------
// FindHits
//-----------
void FindHits(double thresh_MeV, map<int, SumHits> &bcalfADC, map<int, fADCHitList> &fADCHits, map<int,TDCHitList> &TDCHits)
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

      // The histogram should have the signal size for the ADC, but the TADC
      // leg will actually have a larger size since the pre-amp gain will be
      // set differently. Scale the threshold down here to accomodate this.
      double preamp_gain_tdc = 5.0;
      double thresh_MeV_TDC = thresh_MeV/preamp_gain_tdc;
      
      //the outermost layer of the detector is not equipped with TDCs, so don't generate any TDC hits
      int layer = DBCALGeometry::layer(fADCId);

      for(int ii = 0; ii < (int)sumhits.EUP.size(); ii++){
        if(sumhits.EUP[ii] > thresh_MeV && sumhits.tUP[ii] < 2000) uphits.push_back(fADCHit(sumhits.EUP[ii],sumhits.tUP[ii])); // Fill uphits and dnhits with energies (in MeV)
        if(layer != 4 && sumhits.EUP[ii] > thresh_MeV_TDC && sumhits.tUP[ii] < 2000) uphitsTDC.push_back(sumhits.tUP[ii]);     // and times when they cross an energy threshold.
      }                                                                                                                        // Also fill TDC uphits and dnhits with times if
      for(int ii = 0; ii < (int)sumhits.EDN.size(); ii++){                                                                     // they are not layer 4 hits and cross threshold.
        if(sumhits.EDN[ii] > thresh_MeV && sumhits.tDN[ii] < 2000) dnhits.push_back(fADCHit(sumhits.EDN[ii],sumhits.tDN[ii]));
        if(layer != 4 && sumhits.EDN[ii] > thresh_MeV_TDC && sumhits.tDN[ii] < 2000) dnhitsTDC.push_back(sumhits.tDN[ii]);
      }
      
      // If at least one ADC readout channel has a hit, add the readout cell to fADCHits
      if(uphits.size()>0 || dnhits.size()>0){
         fADCHitList &hitlist = fADCHits[fADCId];

         // The module, fADC layer, and fADC sector are encoded in fADCId
         // (n.b. yes, these are the same methods used for extracting 
         // similar quantities from the cellId.)
         hitlist.module = DBCALGeometry::module(fADCId);
         hitlist.sumlayer = DBCALGeometry::layer(fADCId);
         hitlist.sumsector = DBCALGeometry::sector(fADCId);
         
         hitlist.uphits = uphits;
         hitlist.dnhits = dnhits;
      }
      
      // If at least one TDC readout channel has a hit, add the readout cell to TDCHits
      if(uphitsTDC.size()>0 || dnhitsTDC.size()>0){
         TDCHitList &hitlistTDC = TDCHits[fADCId];

         // The module, fADC layer, and fADC sector are encoded in fADCId
         // (n.b. yes, these are the same methods used for extracting 
         // similar quantities from the cellId.)
         hitlistTDC.module = DBCALGeometry::module(fADCId);
         hitlistTDC.sumlayer = DBCALGeometry::layer(fADCId);
         hitlistTDC.sumsector = DBCALGeometry::sector(fADCId);
         
         hitlistTDC.uphits = uphitsTDC;
         hitlistTDC.dnhits = dnhitsTDC;
      }
   }
}

//-----------
// CopyBCALHitsToHDDM
//-----------
void CopyBCALHitsToHDDM(map<int, fADCHitList> &fADCHits,
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
         int integer_time = round((hitlist.uphits[i].t-BCAL_BASE_TIME_OFFSET)/BCAL_NS_PER_ADC_COUNT);
         if (integer_time >= 0){
            hddm_s::BcalfADCDigiHitList fadcs = iter->addBcalfADCDigiHits();
            fadcs().setEnd(bcal_index::kUp);
            fadcs().setPulse_integral(round(hitlist.uphits[i].E/BCAL_MEV_PER_ADC_COUNT));
            fadcs().setPulse_time(integer_time);
         }
      }
      for (unsigned int i = 0; i < hitlist.dnhits.size(); i++) {
         int integer_time = round((hitlist.dnhits[i].t-BCAL_BASE_TIME_OFFSET)/BCAL_NS_PER_ADC_COUNT);
         if (integer_time >= 0){
            hddm_s::BcalfADCDigiHitList fadcs = iter->addBcalfADCDigiHits();
            fadcs().setEnd(bcal_index::kDown);
            fadcs().setPulse_integral(round(hitlist.dnhits[i].E/BCAL_MEV_PER_ADC_COUNT));
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
         int integer_time = round((hitlist.uphits[i]-BCAL_TDC_BASE_TIME_OFFSET)/BCAL_NS_PER_TDC_COUNT);
         if (integer_time >= 0){
            hddm_s::BcalTDCDigiHitList tdcs = iter->addBcalTDCDigiHits();
            tdcs().setEnd(bcal_index::kUp);
            tdcs().setTime(integer_time);
         }
      }
      for (unsigned int i = 0; i < hitlist.dnhits.size(); i++) {
         int integer_time = round((hitlist.dnhits[i]-BCAL_TDC_BASE_TIME_OFFSET)/BCAL_NS_PER_TDC_COUNT);
         if (integer_time >= 0){
            hddm_s::BcalTDCDigiHitList tdcs = iter->addBcalTDCDigiHits();
            tdcs().setEnd(bcal_index::kDown);
            tdcs().setTime(round((hitlist.dnhits[i]-BCAL_TDC_BASE_TIME_OFFSET)/BCAL_NS_PER_TDC_COUNT));
         }
      }
   }
}

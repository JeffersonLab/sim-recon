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
// CellSpectra is a utility class that holds pointers to
// histograms for each end of the cell and the total
// unsmeared, energy deposited in the cell
//..........................
class CellSpectra{
   public:
      CellSpectra() : hup(NULL), hdn(NULL), Etruth(0.0)
      {}
   
      DHistogram *hup;
      DHistogram *hdn;
      
      double Etruth;
};

//..........................
// SumSpectra is a utility class that is used to hold info
// from the SiPMs contributing to that readout channel.
// This includes a list of CellSpectra objects, but also
// the total number of SiPMs that should be in the sum
// so that dark hits for non-hit SiPMs can be added.
//..........................
class SumSpectra{
   public:
      SumSpectra() : NSiPMs(0), hup(NULL), hdn(NULL)
      {}
      
      vector<CellSpectra *> cellspectra;
      int NSiPMs;
   
      DHistogram *hup;
      DHistogram *hdn;
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
void bcalInit(void);
int32_t GetEventNumber(hddm_s::HDDM *record);
int32_t GetRunNumber(hddm_s::HDDM *record);
void GetSiPMSpectra(hddm_s::HDDM *record,
                    map<bcal_index, CellSpectra> &SiPMspectra,
                    vector<IncidentParticle_t> &incident_particles);
void ApplySamplingFluctuations(map<bcal_index,
                               CellSpectra> &SiPMspectra,
                               vector<IncidentParticle_t> &incident_particles);
void MergeSpectra(map<bcal_index, CellSpectra> &SiPMspectra);
void ApplyPoissonStatistics(map<bcal_index, CellSpectra> &SiPMspectra);
void ApplySiPMTimeJitter(map<bcal_index, CellSpectra> &SiPMspectra);
void AddDarkHits(map<bcal_index, CellSpectra> &SiPMspectra);
double AddDarkHitsToOne(DHistogram *h, int NSiPMs=1);
void SortSiPMSpectra(map<bcal_index,
                     CellSpectra> &SiPMspectra,
                     map<int, SumSpectra> &bcalfADC);
void AddDarkHitsForNonHitSiPMs(map<int, SumSpectra> &bcalfADC);
void ApplyElectronicPulseShape(map<int, SumSpectra> &bcalfADC);
void ApplyElectronicPulseShapeOneHisto(DHistogram *h);
void ApplyTimeSmearing(double sigma_ns, map<int, fADCHitList> &fADCHits);
void FindHits(double thresh_mV,
              map<int, SumSpectra> &bcalfADC,
              map<int, fADCHitList> &fADCHits);
void FindHitsOneHisto(double thresh_mV, DHistogram *h, vector<fADCHit> &hits);
void FindTDCHits(double thresh_mV,
                 map<int, SumSpectra> &bcalfADC,
                 map<int, TDCHitList> &F1TDCHits);
void FindTDCHitsOneHisto(double thresh_mV, DHistogram *h, vector<double> &hits);
void CopyBCALHitsToHDDM(map<int, fADCHitList> &fADCHits,
                        map<int, TDCHitList> &F1TDCHits,
                        hddm_s::HDDM *record);
DHistogram* GetHistoFromPool(void);
void ReturnHistoToPool(DHistogram*);
TSpline* MakeTSpline(void);
void SaveDebugHistos(int32_t eventNumber,
                     map<bcal_index, CellSpectra> &SiPMspectra,
                     const char *suffix="", const char *yunits="Amplitude");
void SaveDebugHistos(int32_t eventNumber,
                     map<int, SumSpectra> &bcalfADC, 
                     const char *suffix="");

// The following are set in bcalInit below
bool BCAL_INITIALIZED = false;
double BCAL_mevPerPE = 0.0;
DHistogram *BCAL_PULSE_SHAPE_HISTO = NULL;

// Mutexes and histogram pool
pthread_mutex_t bcal_init_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t histo_pool_mutex = PTHREAD_MUTEX_INITIALIZER;
queue<DHistogram*> histo_pool;

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
extern bool NO_TIME_JITTER;
extern bool NO_THRESHOLD_CUT;
extern bool BCAL_DEBUG_HISTS;

extern double BCAL_TDC_THRESHOLD; // mV
extern double BCAL_ADC_THRESHOLD; // mV
extern double BCAL_FADC_TIME_RESOLUTION; // ns

// setup response parameters
extern double BCAL_DARKRATE_GHZ;                // 0.0176 (from calibDB BCAL/bcal_parms) for 4x4 array
extern double BCAL_XTALK_FRACT;                 // 0.157  (from calibDB BCAL/bcal_parms)
extern double BCAL_DEVICEPDE;                   // 0.21   (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLING_FRACT;              // 0.095  (from calibDB BCAL/bcal_parms)
extern double BCAL_PHOTONSPERSIDEPERMEV_INFIBER;// 75 (from calibDB BCAL/bcal_parms) 
extern double BCAL_SAMPLINGCOEFA;               // 0.042 (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLINGCOEFB;               // 0.013 (from calibDB BCAL/bcal_parms)

// The following are not currently in use
//extern double BCAL_SIGMA_SIG_RELATIVE;          // 0.105  (from calibDB BCAL/bcal_parms)
//extern double BCAL_SIGMA_PED_RELATIVE;          // 0.139  (from calibDB BCAL/bcal_parms)
//extern double BCAL_SIPM_GAIN_VARIATION;         // 0.04   (from calibDB BCAL/bcal_parms)
//extern double BCAL_INTWINDOW_NS;                // 100    (from calibDB BCAL/bcal_parms)
//extern double BCAL_AVG_DARK_DIGI_VALS_PER_EVENT;// 240 used to set thresholds
//extern double BCAL_TIMEDIFFCOEFA;               // 0.07 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
//extern double BCAL_TIMEDIFFCOEFB;               // 0.00 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)


double BCAL_MIN_ENERGY_MEV = 10.0;
double BCAL_FADC_INTEGRATION_WINDOW_PRE  =  20.0; // time to integrate signal before threshold crossing in ns
double BCAL_FADC_INTEGRATION_WINDOW_POST = 180.0; // time to integrate signal after threshold crossing in ns

//-----------
// SmearBCAL
//-----------
void SmearBCAL(hddm_s::HDDM *record)
{

   /// The data from HDGEANT contains attenuated, energy-weighted, timing spectra.
   /// The job of mcsmear is to use that data to create the bcalfADCUpHit and
   /// bcalfADCDownHit structures. These are made by summing signals from multiple
   /// SiPMs and are what are used as the entry points for the reconstruction.
   ///
   /// The time spectra must be smeared due to sampling fluctuations, which are
   /// parameterized based on the total energy of the shower. The total
   /// energy is kept in the bcalTruthIncidentParticle structures in HDDM.
   ///
   /// In addition to the sampling fluctuations, Poisson statistics, timing jitter,
   /// dark pulses, and timewalk effects are all applied. The timewalk effect is obtained
   /// by convoluting the energy spectra with an electronic pulse shape and identifying
   /// the threshold crossing points.
   ///
   
   
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


   // Initialize BCAL globals on first call
   if (!BCAL_INITIALIZED)
      bcalInit();
   
   int32_t eventNumber = GetEventNumber(record);
   
   // First, we extract the time spectra for hit cells
   map<bcal_index, CellSpectra> SiPMspectra;
   vector<IncidentParticle_t> incident_particles;
   GetSiPMSpectra(record, SiPMspectra, incident_particles);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber, SiPMspectra, 
                      "_Raw", "Amplitude (attenuated-MeV)");

   // Sampling fluctuations
   ApplySamplingFluctuations(SiPMspectra, incident_particles);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber, SiPMspectra,
                      "_Sampling", "Amplitude (attenuated-MeV)");

   // Merge spectra associated with different incident particles
   MergeSpectra(SiPMspectra);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber, SiPMspectra,
                      "_Merged", "Amplitude (attenuated-MeV)");

   // Poisson Statistics
   ApplyPoissonStatistics(SiPMspectra);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber, SiPMspectra,
                      "_Poisson", "Amplitude (attenuated-MeV)");
   
   // Time Jitter
   ApplySiPMTimeJitter(SiPMspectra);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber, SiPMspectra,
                      "_TimeJitter", "Amplitude (attenuated-MeV)");
   
   // Add Dark Hits (for hit cells only at this point)
   AddDarkHits(SiPMspectra);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber, SiPMspectra,
                      "_DarkHits", "Amplitude (attenuated-MeV)");
   
   // Place all hit cell spectra into list indexed by fADC ID
   map<int, SumSpectra> bcalfADC;
   SortSiPMSpectra(SiPMspectra, bcalfADC);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber,bcalfADC, "_SummedCell");

   // Add Dark Hits (for cells without energy deposited)
   // (n.b. this may add elements to bcalfADC for summed cells with no signal)
   AddDarkHitsForNonHitSiPMs(bcalfADC);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber,bcalfADC, "_DarkHits_SummedCell");
   
   // Electronic Pulse shape
   ApplyElectronicPulseShape(bcalfADC);
   if (BCAL_DEBUG_HISTS)
      SaveDebugHistos(eventNumber,bcalfADC, "_Electronic");
   
   // Apply threshold to find ADC hits in summed spectra
   map<int, fADCHitList> fADCHits;
   FindHits(BCAL_ADC_THRESHOLD, bcalfADC, fADCHits);

   // Apply additional smearing for the time extracted from 4ns samples
   ApplyTimeSmearing(BCAL_FADC_TIME_RESOLUTION, fADCHits);
   
   // Apply threshold to find TDC hits in summed spectra
   map<int, TDCHitList> F1TDCHits;
   FindTDCHits(BCAL_TDC_THRESHOLD, bcalfADC, F1TDCHits);
   
   // Copy hits into HDDM tree
   CopyBCALHitsToHDDM(fADCHits, F1TDCHits, record);
   
   // Return histogram objects to the pool for recycling
   map<bcal_index, CellSpectra>::iterator iter = SiPMspectra.begin();
   for(; iter != SiPMspectra.end(); iter++){
      CellSpectra &cellspectra = iter->second;
      if (cellspectra.hup)
         ReturnHistoToPool(cellspectra.hup);
      if (cellspectra.hdn)
         ReturnHistoToPool(cellspectra.hdn);
   }
   SiPMspectra.clear();

   map<int, SumSpectra>::iterator biter = bcalfADC.begin();
   for(; biter != bcalfADC.end(); biter++){
      SumSpectra &sumspectra = biter->second;
      if (sumspectra.hup)
         ReturnHistoToPool(sumspectra.hup);
      if (sumspectra.hdn)
         ReturnHistoToPool(sumspectra.hdn);
   }
   bcalfADC.clear();
}

//-----------
// bcalInit
//-----------
void bcalInit(void)
{
   /// Pre-calculate some values that will be used every event

   // Enclose this entire routine in a mutex lock and double check 
   // the intialization flag to guarantee only one thread actually
   // runs this and that all other threads wait for it to finish
   // before continuing.
   pthread_mutex_lock(&bcal_init_mutex);

   // If we did not get here first, then initialization has already
   // been done. Return immediately.
   if (BCAL_INITIALIZED) {
      pthread_mutex_unlock(&bcal_init_mutex);
      return;
   }

   //================ Conversion factor: MeV per PE =======================
   BCAL_mevPerPE = 1.0/( BCAL_PHOTONSPERSIDEPERMEV_INFIBER *
                         BCAL_DEVICEPDE * BCAL_SAMPLING_FRACT );
   //======================================================================


   //================== Electronic Pulse Shape ============================
   // The following will initialize a table that will be used
   // for convoluting the electronic pulse_shape function. The
   // conversion is of the form of a matrix multiplication with
   // the input and output histograms being vectors and this table
   // being the transformation matrix. It is therefore a NBINSxNBINS
   // matrix with NBINS being the number of bins in the histogram
   // being convoluted.
   //
   // For simplicity, we use the same histogram definition for the 
   // input energy-weighted time histogram and the output electronic
   // pulse shape histogram. It's possible some speed optimization
   // could be done by making these different.
   
   cout << "Pre-evaluating BCAL pulse_shape function ..." << endl;

   // Get a histogram from the pool so it has the right dimensions.
   BCAL_PULSE_SHAPE_HISTO = GetHistoFromPool();

   int NBINS = BCAL_PULSE_SHAPE_HISTO->GetNbins();
   //BCAL_PULSE_SHAPE_MATRIX = new double[NBINS*NBINS];
   //BCAL_PULSE_SHAPE_MATRIX_FIRST_NONZERO_BIN = new int[NBINS];
   //BCAL_PULSE_SHAPE_MATRIX_LAST_NONZERO_BIN = new int[NBINS];
   
   // Get electronic pulse shape
   TSpline *pulse_shape = MakeTSpline();

   // Conversion factor between attenuated MeV and mV
   // (see slides from BCAL segmentation meeting 7/22/2011
   // https://halldweb1.jlab.org/wiki/images/3/3f/20110722_bcal.pdf )
   double QCD_counts_per_PE = 61.67;
   double LSB_pC_CAEN_V792 = 0.100;
   double gain_factor_for_test = 20.0;
   double SiPM_pulse_integral_pC = 1200.0;
   double SiPM_pulse_peak_mV = 2293.0;
   double ADC_gain_factor = 1.0; // (see note below)
   
   double mV_per_MeV = 1.0; // just a factor to start with so we can multiply/divide easily below
   mV_per_MeV *= QCD_counts_per_PE;      // QCD counts/PE
   mV_per_MeV *= LSB_pC_CAEN_V792;         // pC/PE
   mV_per_MeV /= BCAL_mevPerPE;         // pC/MeV
   mV_per_MeV /= gain_factor_for_test;      // convert to actual gain of pre-amp
   mV_per_MeV *= SiPM_pulse_peak_mV/SiPM_pulse_integral_pC; // mV/MeV
   mV_per_MeV *= ADC_gain_factor;          // account for gain in ADC leg of pre-amp

   // Pre-Amp Gains
   //-----------------
   // The ADC and TDC legs of the pre-amps will have different gains with
   // the TDC leg having a factor of 5 more gain than the ADC leg. These
   // numbers have floated around a bit over time, but this should be the
   // final design. Some of the above values were taken from test data where
   // a pre-amp gain of 20 was used. This is why it is divided out in the
   // calculation of mv_per_MeV above.
   //
   // The spectra are calculated for the ADC leg so if debugging histograms
   // are recorded, that is what you'll be looking at. For the TDC, the applied
   // threshold is scaled down to accomodate the higher gain.

   // Normalize pulse shape between t=15ns and t=100ns, excluding the
   // pre and post bumps. We scale by the amplitude to give the pulse
   // shape a height of 1 so that when multiplied by mV_per_MeV, the
   // amplitude will again be in mV.
   //double norm = pulse_shape->GetMinimum(15.0, 100.0);
   double norm = -1020.80; // mV not trivial to find minimum with TSpline3

   // time shift pulse shape just to bring it in frame better for zoomed in picture
   double toffset = 6.0;
   
   // Fill the BCAL_PULSE_SHAPE_HISTO histogram
   // n.b. DHistogram content arrays allocate an extra
   // element at the beginning so it can use bin=1 as the
   // first bin. Thus, our loops run from 1 to Nbins rather than from
   // 0 to Nbins-1 (that first bin may eventually
   // be used for underflow).
   float *psh = BCAL_PULSE_SHAPE_HISTO->GetContentPointer();
   int first_nonzero_bin = NBINS;
   int last_nonzero_bin = 1;
   for (int jbin=1; jbin<=NBINS; jbin++) {
      double t1 = BCAL_PULSE_SHAPE_HISTO->GetBinCenter(jbin);
      double t = t1+toffset;

      // The spline data starts at -7ns and ends at 180ns. Beyond those 
      // limits it tends to diverge. Impose a cut here to keep the signal
      // at zero when evaluating the spline beyond those limits.
      if (t<-7.0 || t>180.0) {
         psh[jbin] = 0.0;
      }
      else {
         psh[jbin] = (mV_per_MeV/norm)*pulse_shape->Eval(t);
      }

      // Keep track of first and last non-zero bins in the BCAL_PULSE_SHAPE_MATRIX
      // table. This is used to speed things up when it is finally applied later.
      if (psh[jbin] > 0.0) {
         last_nonzero_bin = jbin;
         if (jbin<first_nonzero_bin) {
            first_nonzero_bin = jbin;
         }
      }
   }
   
   // Trim empty bins off outsides of histogram to speed things
   // up at event time.
   int Nbins_trimmed = last_nonzero_bin - first_nonzero_bin + 1;
   float lo = BCAL_PULSE_SHAPE_HISTO->GetBinLowEdge(first_nonzero_bin);
   float hi = BCAL_PULSE_SHAPE_HISTO->GetBinLowEdge(last_nonzero_bin+1);
   DHistogram *h = new DHistogram(Nbins_trimmed, lo, hi);
   float *hptr = h->GetContentPointer();
   for(int ibin=1; ibin<=Nbins_trimmed; ibin++){
      hptr[ibin] = psh[first_nonzero_bin+ibin-1];
   }
   
   // Return original BCAL_PULSE_SHAPE_HISTO to pool and keep
   // trimmed histo in its place
   ReturnHistoToPool(BCAL_PULSE_SHAPE_HISTO);
   BCAL_PULSE_SHAPE_HISTO = h;
   
   // Optionally create ROOT 1-D histo of pulse shape
   if(BCAL_DEBUG_HISTS){
      TH1D *h_pulse_shape = BCAL_PULSE_SHAPE_HISTO->MakeTH1D("pulse_shape","SiPM Electronic Pulse Shape");
      h_pulse_shape->SetXTitle("time (ns)");
      h_pulse_shape->SetYTitle("Amplitude (mV/MeV(attenuated))");
   }   

   // Clean up
   delete pulse_shape; // Get rid of TSpline since we don't need it anymore
   
   // Create histos for debugging
   hNincident_particles = new TH1D("Nincident_particles","Number of particles marked as 'incident'", 201,-0.5, 200.5);
   
   cout << "Done." << endl;
   //======================================================================
   
   BCAL_INITIALIZED = true;
   pthread_mutex_unlock(&bcal_init_mutex);

}

//-----------
// GetEventNumber
//-----------
int32_t GetEventNumber(hddm_s::HDDM *record)
{
   return record->getPhysicsEvent().getEventNo();
}

//-----------
// GetRunNumber
//-----------
int32_t GetRunNumber(hddm_s::HDDM *record)
{
   return record->getPhysicsEvent().getRunNo();
}

//-----------
// GetSiPMSpectra
//-----------
void GetSiPMSpectra(hddm_s::HDDM *record, 
                   map<bcal_index, CellSpectra> &SiPMspectra,
                   vector<IncidentParticle_t> &incident_particles)
{
   /// Loop through input HDDM data and extract the timing spectrum info into
   /// CellSpectra objects.
   
   // Make sure HDDM stuctures exist.
   // In the case of no real BCAL hits, we may still want to emit
   // dark hit only events. In this case, we must create the BCAL 
   // tree here.
   hddm_s::BarrelEMcalList bcals = record->getBarrelEMcals();
   if (bcals.size() == 0)
      bcals = record->getHitViews().begin()->addBarrelEMcals();

   // Loop over GEANT hits in BCAL
   hddm_s::BcalSiPMSpectrumList specs = record->getBcalSiPMSpectrums();
   hddm_s::BcalSiPMSpectrumList::iterator iter;
   for (iter = specs.begin(); iter != specs.end(); ++iter) {
      bcal_index idx(iter->getModule(), iter->getLayer(),
                     iter->getSector(), 
                     iter->getBcalSiPMTruth().getIncident_id(),
                     (iter->getEnd() == 0)? bcal_index::kUp :
                                            bcal_index::kDown);

      // Get reference to existing CellSpectra, or create one if it doesn't exist
      CellSpectra &cellspectra = SiPMspectra[idx];

      DHistogram *h=NULL;
      if (idx.end == bcal_index::kUp) {
         if (!cellspectra.hup)
            cellspectra.hup = GetHistoFromPool();
         h = cellspectra.hup;
      }
      else {
         if (!cellspectra.hdn)
            cellspectra.hdn = GetHistoFromPool();
         h = cellspectra.hdn;
      }
         
      // The energy truth info for the cell is stored with the spectra for both
      // ends. (This is for when hdgeant cuts the spectrum for one end away
      // due to threshold). They should be the same if both spectra are here so
      // we overwrite the one data member.
      cellspectra.Etruth = iter->getBcalSiPMTruth().getE();

      double t = iter->getTstart();
      double bin_width = iter->getBin_width();
      stringstream ss(iter->getVals());

      // Extract values and use them to fill histo
      int i = 0;
      double dE;
      while (ss >> dE) {
         h->Fill(t, dE);
         t += bin_width;
         if (++i > 100000) {
            _DBG_ << "More than 100k iterations parsing BCAL time"
                  << " spectrum string! Something is wrong!" << endl;
            exit(-1);
         }
      }
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
void ApplySamplingFluctuations(map<bcal_index, CellSpectra> &SiPMspectra, vector<IncidentParticle_t> &incident_particles)
{
   /// Loop over the CellSpectra objects and apply sampling fluctuations.
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
   /// cell energy to unsmeared cell energy and scaling both Eup and Edn
   /// by it.
   
   if(NO_SAMPLING_FLUCTUATIONS)return;
   if(NO_SAMPLING_FLOOR_TERM)BCAL_SAMPLINGCOEFB=0.0; // (redundant, yes, but located in more obvious place here)

   map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
   for(; iter!=SiPMspectra.end(); iter++){
      CellSpectra &cellspectra = iter->second;
      
      // Find fractional sampling sigma based on deposited energy (whole colorimeter, not just fibers)
      double Etruth = cellspectra.Etruth;
      double sigmaSamp = BCAL_SAMPLINGCOEFA / sqrt( Etruth ) + BCAL_SAMPLINGCOEFB;

      // Convert sigma into GeV
      sigmaSamp *= Etruth;

      // Randomly sample the fluctuation
      double Esmeared = gDRandom.Gaus(Etruth,sigmaSamp);

      // Calculate ratio of smeared to unsmeared
      double ratio = Esmeared/Etruth;

      // Scale attenuated energy histos for each end
      if(cellspectra.hup)cellspectra.hup->Scale(ratio);
      if(cellspectra.hdn)cellspectra.hdn->Scale(ratio);
   }
}

//-----------
// MergeSpectra
//-----------
void MergeSpectra(map<bcal_index, CellSpectra> &SiPMspectra)
{
   /// Combine all SiPM CellSpectra corresponding to the same
   /// cell but different incident particles into a single
   /// spectrum. This is done after the sampling fluctuations
   /// have been applied so there is no more dependence on
   /// the incident particle parameters.
   
   // Loop until no merges are made
   while(true){
      bool merged=false;
      map<bcal_index, CellSpectra>::iterator iter1=SiPMspectra.begin();
      for(;iter1!=SiPMspectra.end(); iter1++){
         map<bcal_index, CellSpectra>::iterator iter2 = iter1;
         for(++iter2; iter2!=SiPMspectra.end(); iter2++){
         
            // If spectra are not from same module,layer,sector,end
            // then just continue the loop
            if(iter1->first.module != iter2->first.module)continue;
            if(iter1->first.layer  != iter2->first.layer )continue;
            if(iter1->first.sector != iter2->first.sector)continue;
            if(iter1->first.end    != iter2->first.end   )continue;

            // ----- Merge spectra -----
            // Get pointers to histograms
            DHistogram *hup1 = iter1->second.hup;
            DHistogram *hup2 = iter2->second.hup;
            DHistogram *hdn1 = iter1->second.hdn;
            DHistogram *hdn2 = iter2->second.hdn;
            
            // It's possible one or both of the histograms we wish to merge
            // don't exist. Check for this and handle accordingly.
            
            // Upstream
            if(hup1!=NULL && hup2!=NULL)hup1->Add(hup2);
            if(hup1==NULL && hup2!=NULL){
               iter1->second.hup=hup2;
               hup2=NULL;
            }

            // Downstream
            if(hdn1!=NULL && hdn2!=NULL)hdn1->Add(hdn2);
            if(hdn1==NULL && hdn2!=NULL){
               iter1->second.hdn=hdn2;
               hdn2=NULL;
            }

            // Erase second one
            if(hup2)ReturnHistoToPool(hup2);
            if(hdn2)ReturnHistoToPool(hdn2);
            SiPMspectra.erase(iter2);
            
            // Set flag that we did merge spectra and break
            // the loops so we can try again.
            merged = true;
            break;
         }
         if(merged)break;
      }
      
      // When we make it through without merging any spectra,
      // we're done so break out of the infinite while loop.
      if(!merged)break;
   }
}

//-----------
// ApplyPoissonStatistics
//-----------
void ApplyPoissonStatistics(map<bcal_index, CellSpectra> &SiPMspectra)
{
   /// Loop over the CellSpectra objects and apply Poisson Statistics.
   ///
   /// Because the response of the SiPM is quantized in units of photo-electrons
   /// Poisson counting statistics should be applied. This will affect the
   /// smaller energy depositions more than the larger ones.
   ///
   /// We do this by converting the integral of the attenuated energy into
   /// photo-electrons and then sampling from a Poisson distribution with that
   /// mean. The ratio of the quantized, sampled value to the unquantized
   // integral (in PE) is used to scale the spectrum.

   if(NO_POISSON_STATISTICS)return;

   map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
   for(; iter!=SiPMspectra.end(); iter++){
      CellSpectra &cellspectra = iter->second;
      
      // Integrate hup and hdn
      double Iup = (cellspectra.hup ? cellspectra.hup->Integral():0.0); // in attenuated MeV
      double Idn = (cellspectra.hdn ? cellspectra.hdn->Integral():0.0); // in attenuated MeV

      // Upstream end
      if(Iup>0.0){
         // Convert to number of PE
         double mean_pe_up = Iup/BCAL_mevPerPE;
         
         int Npe_up = gDRandom.Poisson(mean_pe_up);
         double ratio_up = (double)Npe_up/mean_pe_up;

         cellspectra.hup->Scale(ratio_up); // only get here if cellspectra.hup is non-NULL
      }

      // Downstream end
      if(Idn>0.0){
         // Convert to number of PE
         double mean_pe_dn = Idn/BCAL_mevPerPE;
         
         int Npe_dn = gDRandom.Poisson(mean_pe_dn);
         double ratio_dn = (double)Npe_dn/mean_pe_dn;
         
         cellspectra.hdn->Scale(ratio_dn); // only get here if cellspectra.hdn is non-NULL
      }
   }
}

//-----------
// ApplySiPMTimeJitter
//-----------
void ApplySiPMTimeJitter(map<bcal_index, CellSpectra> &SiPMspectra)
{
   /// Loop over the CellSpectra objects and apply SiPM time jitter.
   ///
   /// The photo detection device has an intrinsic timing resolution
   /// caused by jitter in the time between the input signal (light)
   /// and the generated electronic pulse. This is quoted in the 
   /// spec sheet as "Time Resolution (FWHM)" for a single photon.
   /// (see page 2 of the following link:
   /// http://sales.hamamatsu.com/assets/pdf/parts_S/s10362-33series_kapd1023e05.pdf)
   ///
   /// n.b. XP2020 PMT's have a very similar time jitter as the SiPMs
   ///
   /// To include this effect, each bin is converted into a number
   /// of photo-electrons and each of them is randomly displaced in time
   /// with the specified width. Because the photo-statistics is
   /// applied on a global scale, we do not re-quantize the bin contents
   /// into an integer number of photo-electrons. Rather, we split it
   /// into a number of pieces that is equal to the approximate number
   /// of photo-electrons coming from this bin and shift those pieces.
   /// Bins with trace amounts of energy are still kept and shifted
   /// as well to keep the energy integral constant.
   
   if(NO_TIME_JITTER)return;
   
   double fwhm_jitter = 0.600; // ns
   double sigma_jitter = fwhm_jitter/2.35; // ns
   
   // Borrow a histogram from the pool
   DHistogram *tmp = GetHistoFromPool();
   int Nbins = tmp->GetNbins();

   map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
   for(; iter!=SiPMspectra.end(); iter++){
      CellSpectra &cellspectra = iter->second;
      
      // Upstream
      if(cellspectra.hup){
         tmp->Reset();
         for(int ibin=1; ibin<=Nbins; ibin++){
            double E = cellspectra.hup->GetBinContent(ibin);
            if(E==0.0)continue;
            double tbin = cellspectra.hup->GetBinCenter(ibin);

            int Npe = 1 + (int)floor(E/BCAL_mevPerPE);
            double dE = E/(double)Npe;
            
            for(int i=0; i<Npe; i++){
               double t = gDRandom.Gaus(tbin, sigma_jitter);
               tmp->Fill(t, dE);
            }
         }
         *cellspectra.hup = *tmp; // overwrite histo with tmp
      }

      // Downstream
      if(cellspectra.hdn){
         tmp->Reset();
         for(int ibin=1; ibin<=Nbins; ibin++){
            double E = cellspectra.hdn->GetBinContent(ibin);
            if(E==0.0)continue;
            double tbin = cellspectra.hdn->GetBinCenter(ibin);

            int Npe = 1 + (int)floor(E/BCAL_mevPerPE);
            double dE = E/(double)Npe;
            
            for(int i=0; i<Npe; i++){
               double t = gDRandom.Gaus(tbin, sigma_jitter);
               tmp->Fill(t, dE);
            }
         }
         *cellspectra.hdn = *tmp; // overwrite histo with tmp
      }
   }
   
   // Return histogram to pool
   ReturnHistoToPool(tmp);
}

//-----------
// AddDarkHits
//-----------
void AddDarkHits(map<bcal_index, CellSpectra> &SiPMspectra)
{
   /// Loop over the CellSpectra objects and add Dark Hits.
   ///
   /// The SiPMs will randomly fire pixels at some rate, even when no
   /// light is present. These "dark hits" will add to the signal
   /// spectrum. This can affect the timing if the hits show up close
   /// to the leading edge of the real signal. They can also add up
   /// enough to occasionally cross threshold, even if no actual signal is 
   /// present. In this routine, dark hits are added to the SiPMs
   /// that already have signal. SiPM spectra with no signal will have
   /// dark hits added later in AddDarkHitsForNonHitSiPMs().
   ///
   
   if(NO_DARK_PULSES)return;

   map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
   for(; iter!=SiPMspectra.end(); iter++){
      CellSpectra &cellspectra = iter->second;
      
      if(cellspectra.hup)AddDarkHitsToOne(cellspectra.hup);
      if(cellspectra.hdn)AddDarkHitsToOne(cellspectra.hdn);
   }
}

//-----------
// AddDarkHitsToOne
//-----------
double AddDarkHitsToOne(DHistogram *h, int NSiPMs)
{
   /// Dark hits are added by using 1/rate as the mean time between
   /// hits. It is assumed that the timing distribution bewteen hits
   /// is exponential. The time differences are randomly sampled until
   /// a hit is found past the end of the histogram limit.
   ///
   /// The second parameter is used to scale the rate up for the case 
   /// when multiple SiPMs are summed.
   ///
   /// The return value is the total energy equivalent in MeV added
   /// due to dark hits. This includes primaries plus cross-talk.
   ///
   /// The parameters for the dark noise should come from the CCDB. For 
   /// a reference, one can look at GlueX-doc-1754. Slide 16 has the total
   /// dark rate (17.6MHz), and slide 17 has the cross-talk (15.7%)

   int Nbins = h->GetNbins();
   double low_edge = h->GetBinLowEdge(1);
   double high_edge = h->GetBinLowEdge(Nbins+1);
   double t = low_edge;

   double tau_ns = 1.0/((double)NSiPMs*BCAL_DARKRATE_GHZ);
   
   double dE_tot = 0;
   while(t < high_edge){ // actually break inside while loop so this is superfluous

      // Loop until we get a finite time difference
      // (both s=0 and s=1 are poorly behaved).
      double delta_t = 0.0;
      do{
         double s = gDRandom.Rndm();
         delta_t = tau_ns*log(1.0/(1.0-s)); // derived by setting s equal to integral fraction of exp(t/tau)
      }while(!isfinite(delta_t));
      
      t += delta_t;
      if(t > high_edge) break;
      
      double dE = BCAL_mevPerPE;  // single photo electron
      if(gDRandom.Rndm() <= BCAL_XTALK_FRACT){
         dE += BCAL_mevPerPE; // single cross-talk pixel
         if(gDRandom.Rndm() <= BCAL_XTALK_FRACT){
            dE += BCAL_mevPerPE; // second cross-talk pixel
         }
      }
      h->Fill(t, dE);
      dE_tot += dE;
   }

   return dE_tot;
}

//-----------
// SortSiPMSpectra
//-----------
void SortSiPMSpectra(map<bcal_index, CellSpectra> &SiPMspectra, map<int, SumSpectra> &bcalfADC)
{
   /// Loop over the CellSpectra objects and copy pointers to them into SumSpectra objects.
   ///
   /// For the BCAL, multiple SiPMs are summed together. This routine gathers individual
   /// SiPM spectra into single SumSpectra objects. Each SumSpectra represents a summed
   /// cell that is readout by an fADC channel. Since not every cell has signal in it, each
   /// SumSpectra object may not have as many input spectra as SiPMs that will actually be
   /// contributing. This is only an issue for dark hits in SiPMs with no signal. Those 
   /// non-signal SiPM dark hits are handled later in AddDarkHitsForNonHitSiPMs.
   ///
   /// To save space and, more importantly, time, the first histogram is promoted to be the
   /// histo for the SumSpectra and the pointer in the CellSpectra object set to NULL. (This 
   /// prevents it from being returned to the pool twice at the end of the event.)
   
   // Loop over SiPMspectra and copy a pointer to it to the correct SumSpectra
   // element in the bcalfADC container. The bcalfADC container is an STL map which is
   // convenient since it creates a new SumSpectra object for us if it doesn't exist,
   // but otherwise, returns a reference to the existing object.
   
   map<bcal_index, CellSpectra>::iterator iter = SiPMspectra.begin();
   for(; iter!=SiPMspectra.end(); iter++){
      
      // Get reference to SumSpectra object
      const bcal_index &idx = iter->first;
      int fADCId = DBCALGeometry::fADCId( idx.module, idx.layer, idx.sector);
      SumSpectra &sumspectra = bcalfADC[fADCId];
      
      // Add CellSpectra object to list in SumSpectra
      CellSpectra &cellspectra = iter->second;
      sumspectra.cellspectra.push_back(&cellspectra);
      
      // If this is the first cell added to the SumSpectra, re-assign its
      // histo to the SumSpectra. Otherwise, just add to it.

      // Upstream
      if(sumspectra.hup==NULL){
         sumspectra.hup = cellspectra.hup;
         cellspectra.hup = NULL;
      }else{
         if(cellspectra.hup)sumspectra.hup->Add(cellspectra.hup);
      }
      
      // Downstream
      if(sumspectra.hdn==NULL){
         sumspectra.hdn = cellspectra.hdn;
         cellspectra.hdn = NULL;
      }else{
         if(cellspectra.hdn)sumspectra.hdn->Add(cellspectra.hdn);
      }
   }
}

//-----------
// AddDarkHitsForNonHitSiPMs
//-----------
void AddDarkHitsForNonHitSiPMs(map<int, SumSpectra> &bcalfADC)
{
   /// Loop over the SumSpectra objects and add Dark Hits for channels
   /// that do not have any signal.
   ///
   /// Dark hits have already been added for cells that have signal in them
   /// (see AddDarkHits(...) ). Here, we add them in for the rest of the 
   /// channels. There are two categories of channels we need to add Dark Hits
   /// for:
   ///
   /// 1. Summed cells that already have at least one SiPM with signal, but
   ///    less than the number of summed SiPMs have signal
   ///
   /// 2. Summed cells with no SiPMs having signal. 
   ///
   /// Case 1 will add to the existing bcalfADC elements. Case 2 will
   /// need to create new elements. For dark-hit-only channels, a threshold
   /// is applied so only summed channels with at least BCAL_MIN_ENERGY_MEV
   /// equivalent deposited will have non-NULL histograms. All summed channels
   /// will have an entry in bcalfADC upon exit though, even if they have NULL
   /// histos.
   
   // In order to avoid potentially expensive mutex locks, we grab a couple
   // of histograms from the pool now and reuse them for the dark-hits-only
   // channels until we actually find one that has enough dark hits for us 
   // to keep it. Then we'll grab another from the pool. Since most channels
   // do NOT have signal and also do NOT have enough dark hits to make a
   // signal above threshold, it is worthwhile to reuse these temporary
   // histograms for multiple summed cells.

   if(NO_DARK_PULSES)return;

   DHistogram *hup_tmp = GetHistoFromPool();
   DHistogram *hdn_tmp = GetHistoFromPool();
   
   // Loop over all fADC readout cells
   for(int imodule=1; imodule<=DBCALGeometry::NBCALMODS; imodule++){

      int n_layers = DBCALGeometry::NBCALLAYSIN + DBCALGeometry::NBCALLAYSOUT;
      for(int fADC_lay=1; fADC_lay<=n_layers; fADC_lay++){

         int n_sectors = (fADC_lay <= DBCALGeometry::NBCALLAYSIN)? DBCALGeometry::NBCALSECSIN : DBCALGeometry::NBCALSECSOUT;
         for(int fADC_sec=1; fADC_sec<=n_sectors; fADC_sec++){

            // Use cellId(...) to convert fADC layer and sector into fADCId
            // (see DBCALGeometry::fADCId)
            int fADCId = DBCALGeometry::cellId(imodule, fADC_lay, fADC_sec);
            
            // Get SumSpectra object if it already exists or create new one 
            // if it doesn't.
            SumSpectra &sumspectra = bcalfADC[fADCId];
            
            // If we just created this, then grab some histos from pool.
            // n.b. we try first to use the pre-grabbed histos
            if(sumspectra.hup == NULL){
               if(hup_tmp == NULL){
                  sumspectra.hup = GetHistoFromPool();
               }else{
                  sumspectra.hup = hup_tmp;
                  hup_tmp = NULL;
               }
            }
            if(sumspectra.hdn == NULL){
               if(hdn_tmp == NULL){
                  sumspectra.hdn = GetHistoFromPool();
               }else{
                  sumspectra.hdn = hdn_tmp;
                  hdn_tmp = NULL;
               }
            }

            // Get the number of SiPMs in this readout channel.
            int NSiPMs = DBCALGeometry::NSiPMs(fADCId);
            
            // Add in dark hits for any additional SiPMs.
            int NSiPMs_no_signal =  NSiPMs - (int)sumspectra.cellspectra.size();
            if(NSiPMs_no_signal>0){
               double dEup = AddDarkHitsToOne(sumspectra.hup, NSiPMs_no_signal);
               double dEdn = AddDarkHitsToOne(sumspectra.hdn, NSiPMs_no_signal);
               
               // If this is a dark-hit-only channel, check if enough
               // equivalent energy was added to make it worthwhile to keep
               // these. Otherwise, recycle the histos for the next iteration
               // to avoid the expensive convolution and checking for thresholds
               // stages later.
               // BCAL_MIN_ENERGY_MEV (= 10 MeV) is set appropriately so that
               // it corresponds roughly with the later threshold
               // (BCAL_ADC_THRESHOLD = 4 mV), so that nothing eliminated here
               // would cross the later threshold.
               if(NSiPMs_no_signal == NSiPMs){
                  if(dEup < BCAL_MIN_ENERGY_MEV){
                     hup_tmp = sumspectra.hup;
                     hup_tmp->Reset();
                     sumspectra.hup = NULL;
                  }
                  if(dEdn < BCAL_MIN_ENERGY_MEV){
                     hdn_tmp = sumspectra.hdn;
                     hdn_tmp->Reset();
                     sumspectra.hdn = NULL;
                  }
               }
            }
         }
      }
   }
   
   // Return temporary histos to pool
   if(hup_tmp)ReturnHistoToPool(hup_tmp);
   if(hdn_tmp)ReturnHistoToPool(hdn_tmp);
   
}

//-----------
// ApplyElectronicPulseShape
//-----------
void ApplyElectronicPulseShape(map<int, SumSpectra> &bcalfADC)
{
   /// Loop over the SumSpectra objects and convolute any histograms
   /// with the electronic pulse shape.

   map<int, SumSpectra>::iterator iter = bcalfADC.begin();
   for(; iter!=bcalfADC.end(); iter++){
      
      SumSpectra &sumspectra = iter->second;
      
      if(sumspectra.hup)ApplyElectronicPulseShapeOneHisto(sumspectra.hup);
      if(sumspectra.hdn)ApplyElectronicPulseShapeOneHisto(sumspectra.hdn);
   }
}

//-----------
// ApplyElectronicPulseShapeOneHisto
//-----------
void ApplyElectronicPulseShapeOneHisto(DHistogram *h)
{
   /// Convolute the given histogram with the electronic pulse
   /// shape, overwriting the input histograms contents.
   
   // Make copy of input histogram and then reset it
   DHistogram tmp(*h);
   h->Reset();
   
   int Nbins = h->GetNbins();

   // Get pointers to the bin content arrays directly.
   // This will allow the two nested loops below to
   // work much more efficiently.
   // n.b. DHistogram content arrays allocate an extra
   // element at the beginning so it can use bin=1 as the
   // first bin. Thus, our loops run from 1 to Nbins rather than from
   // 0 to Nbins-1 (that first bin may eventually
   // be used for underflow).
   float* tmp_content = tmp.GetContentPointer();
   float* h_content = h->GetContentPointer();
   float* pulse_shape = BCAL_PULSE_SHAPE_HISTO->GetContentPointer();
   //BCAL_PULSE_SHAPE_HISTO will generally have fewer bins than Nbins
   int pulse_Nbins = BCAL_PULSE_SHAPE_HISTO->GetNbins();

   // Loop over bins of input histo
   for(int ibin=1; ibin<=Nbins; ibin++){
      float A = tmp_content[ibin];
      if(A==0.0)continue; // no need to continue for empty bins
      
      int end_bin = pulse_Nbins;
      if(end_bin > (Nbins-ibin+1))end_bin = (Nbins-ibin+1);
      for(int jbin=1; jbin<=end_bin; jbin++){
         float weight = A*(pulse_shape[jbin]);
         h_content[ibin+jbin-1] += weight;
      }      
   }
}

//-----------
// ApplyTimeSmearing
//-----------
void ApplyTimeSmearing(double sigma_ns, map<int, fADCHitList> &fADCHits)
{
   /// The fADC250 will extract a time from the samples by applying an algorithm
   /// to a few of the samples taken every 4ns. The time resolution of this is
   /// going to be worse than that obtained from our histograms which effectively
   /// sample every 0.1ns. This will apply additional smearing to the times in the
   /// given list of hits to account for this.

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
}

//-----------
// FindHits
//-----------
void FindHits(double thresh_mV, map<int, SumSpectra> &bcalfADC, map<int, fADCHitList> &fADCHits)
{
   /// Loop over SumSpectra objects and find places where it crosses threshold
   map<int, SumSpectra>::iterator iter = bcalfADC.begin();
   for(; iter!=bcalfADC.end(); iter++){
      
      int fADCId = iter->first;
      SumSpectra &sumspectra = iter->second;
      
      // Find threshold crossings for both upstream and downstream readout channels
      vector<fADCHit> uphits;
      vector<fADCHit> dnhits;

      if(sumspectra.hup)FindHitsOneHisto(thresh_mV, sumspectra.hup, uphits);
      if(sumspectra.hdn)FindHitsOneHisto(thresh_mV, sumspectra.hdn, dnhits);
      
      // If at least one readout channel has a hit, add the readout cell to fADCHits
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
   }
}

//-----------
// FindHitsOneHisto
//-----------
void FindHitsOneHisto(double thresh_mV, DHistogram *h, vector<fADCHit> &hits)
{
   /// Look through the given histogram and find places where the signal
   /// size crosses the given threshold. The signal is integrated before
   /// and after the crossing time by 20ns and 180ns respectively. The 
   /// times are approximate as they are measured in number of bins by
   /// dividing by the bin width of the histo.
   ///
   /// The amplitude of the signal is scaled down to account for the pre_amp
   /// gain difference between TDC and fADC signal splits. It is also converted
   /// into GeV-equivalent units.
   ///
   /// n.b. no quantization of the ADC counts is done at this stage since
   /// it is unclear if the conversions are correct. This can be done
   /// easily enough at a later stage after any additional scaling factors
   /// are applied.
   ///
   /// The time is taken from the center of the bin with the largest amplitude
   
   double bin_width = h->GetBinWidth();
   int Nbins_before = (int)(BCAL_FADC_INTEGRATION_WINDOW_PRE/bin_width);
   int Nbins_after = (int)(BCAL_FADC_INTEGRATION_WINDOW_POST/bin_width);

   // Loop 
   int start_bin = 1;
   int Nbins = h->GetNbins();
   while(start_bin<=Nbins){
      int ibin = h->FindFirstBinAbove(thresh_mV, start_bin);
      if(ibin<1 || ibin>Nbins)break;
      
      // Signal time. Start with the center of the bin as the time then
      // linearly interpolate (below) to the threshold crossing point.
      // Notice that we will have 2 times from the BCAL, one from the
      // fADC and one from the TDC. Here, we only specify one and make it
      // more closely match the TDC resolution.
      double t = h->GetBinCenter(ibin);
      
      // Linearly interpolate between this and previous bin
      if(ibin>1){
         double a1 = h->GetBinContent(ibin-1);
         double a2 = h->GetBinContent(ibin);
         double delta_t = bin_width*(a2-thresh_mV)/(a2-a1);
         t -= delta_t;
      }
      
      // Calculate integration limits for signal amplitude
      int istart = ibin - Nbins_before;
      int iend = ibin + Nbins_after;
      if(istart<1)istart=1;
      if(iend>Nbins)iend = Nbins;

      // Integrate signal
      // n.b. we use the start_bin variable so it is left
      // pointing to the end of the integration window which
      // is where we start looking for the next hit on the next iteration 
      double integral = 0.0;
      for(start_bin=istart; start_bin<=iend; start_bin++){
         integral += h->GetBinContent(start_bin);
      }
      
      //Previously we converted this integral into fADC units. This
      //conversion factor was simply bin_width/(4 ns)*4096/(2000 mV).
      //Now, however, for consistency with the rest of the software, we
      //output hit energy in GeV. Naively we would expect the conversion
      //factor between integrated mV and integrated MeV to be the same as the
      //mV_per_MeV factor calculated in bcal_init() above (numerically this
      //factor is about 0.88), but this is not
      //the case for some reason, probably because the integral of the
      //response function (spline) is not equal to one.
      //When performing the ApplyElectronicPulseShape() step,
      //one notices that the integral of
      //the histogram increases by a factor of about 303 (this can be
      //verified using the debug_hists option). This is effectively
      //the conversion factor between integrated mV and integrated MeV.
      //-WL 2014-07-25
      double adhoc_mV_per_MeV = 303.0;
      double E_GeV = integral/adhoc_mV_per_MeV;
      E_GeV /= 1000.0; //convert MeV to GeV

      // Store hit in container
      hits.push_back(fADCHit(E_GeV,t));
   }

}

//-----------
// FindTDCHits
//-----------
void FindTDCHits(double thresh_mV, map<int, SumSpectra> &bcalfADC, map<int, TDCHitList> &F1TDCHits)
{
   /// Loop over SumSpectra objects and find places where it crosses threshold
   map<int, SumSpectra>::iterator iter = bcalfADC.begin();
   for(; iter!=bcalfADC.end(); iter++){
      
      int fADCId = iter->first;

      //the outermost layer of the detector is not equipped with TDCs, so don't generate any TDC hits
      int layer = DBCALGeometry::layer(fADCId);
      if (layer == 4) continue;

      SumSpectra &sumspectra = iter->second;
      
      // Find threshold crossings for both upstream and downstream readout channels
      vector<double> uphits;
      vector<double> dnhits;

      if(sumspectra.hup)FindTDCHitsOneHisto(thresh_mV, sumspectra.hup, uphits);
      if(sumspectra.hdn)FindTDCHitsOneHisto(thresh_mV, sumspectra.hdn, dnhits);

      // If at least one readout channel has a hit, add the readout cell to fADCHits
      if(uphits.size()>0 || dnhits.size()>0){
         TDCHitList &hitlist = F1TDCHits[fADCId];
         
         // The module, fADC layer, and fADC sector are encoded in fADCId
         // (n.b. yes, these are the same methods used for extracting 
         // similar quantities from the cellId.)
         hitlist.module = DBCALGeometry::module(fADCId);
         hitlist.sumlayer = DBCALGeometry::layer(fADCId);
         hitlist.sumsector = DBCALGeometry::sector(fADCId);
         
         hitlist.uphits = uphits;
         hitlist.dnhits = dnhits;
      }
   }
}

//-----------
// FindTDCHitsOneHisto
//-----------
void FindTDCHitsOneHisto(double thresh_mV, DHistogram *h, vector<double> &t_hits)
{
   /// This is used for finding accurate leading edge times as would be
   /// reported by the F1TDCs.
   ///
   /// Looks through the given histogram and finds places where the signal
   /// size crosses the given threshold. The time is linearly interpolated
   /// between bins to more accurately determine it.
   ///
   /// Multiple hits may be found. Second pulses are a minimum of 20ns apart
   /// and require a rising edge crossing the threshold.

   double bin_width = h->GetBinWidth();
   double minimum_pulse_separation = 20.0; // ns
   int Nbins_to_skip = (int)(minimum_pulse_separation/bin_width);

   // The histogram should have the signal size for the ADC, but the
   // TADC leg will actually have a larger size since the pre-amp gain
   // will be set differently. Scale the threshold down here to accomodate
   // this.
   double preamp_gain_tdc = 5.0;
   thresh_mV /= preamp_gain_tdc;

   // Loop 
   int start_bin = 1;
   int Nbins = h->GetNbins();
   while(start_bin<=Nbins){
      int ibin = h->FindFirstBinAbove(thresh_mV, start_bin);
      if(ibin<1 || ibin>Nbins)break;
      
      // Signal time. Start with the center of the bin as the time then
      // linearly interpolate (below) to the threshold crossing point.
      // Notice that we will have 2 times from the BCAL, one from the
      // fADC and one from the TDC. Here, we only specify one and make it
      // more closely match the TDC resolution.
      double t = h->GetBinCenter(ibin);
      
      // Linearly interpolate between this and previous bin
      if(ibin>1){
         double a1 = h->GetBinContent(ibin-1);
         double a2 = h->GetBinContent(ibin);
         double delta_t = bin_width*(a2-thresh_mV)/(a2-a1);
         t -= delta_t;
      }

      // Store hit in container
      t_hits.push_back(t);
      
      // Skip the minimum number of bins before we can start looking
      // for a new edge.
      start_bin = ibin + Nbins_to_skip;
      
      // Advance until we are either below threshold, or hit the end of the histo
      for(; start_bin<=Nbins; start_bin++){
         if(h->GetBinContent(start_bin) < thresh_mV) break;
      }
   }
}

//-----------
// CopyBCALHitsToHDDM
//-----------
void CopyBCALHitsToHDDM(map<int, fADCHitList> &fADCHits,
                        map<int, TDCHitList> &F1TDCHits,
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
   if (bcals.size() == 0)
      bcals = record->getHitViews().begin()->addBarrelEMcals();
   hddm_s::BcalCellList cells = bcals().getBcalCells();
   hddm_s::BcalCellList::iterator iter;
   for (iter = cells.begin(); iter != cells.end(); ++iter) {

      // Delete any existing bcalfADCHit and bcalTDCHit structures
       iter->deleteBcalfADCHits();
       iter->deleteBcalTDCHits();
   }

   // If we have no cells over threshold, then bail now.
   if (fADCHits.size() == 0 && F1TDCHits.size() == 0)
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
      
      // Copy hits into HDDM
      for (unsigned int i = 0; i < hitlist.uphits.size(); i++) {
         hddm_s::BcalfADCHitList fadcs = iter->addBcalfADCHits();
         fadcs().setEnd(bcal_index::kUp);
         fadcs().setE(hitlist.uphits[i].E);
         fadcs().setT(hitlist.uphits[i].t);
      }
      for (unsigned int i = 0; i < hitlist.dnhits.size(); i++) {
         hddm_s::BcalfADCHitList fadcs = iter->addBcalfADCHits();
         fadcs().setEnd(bcal_index::kDown);
         fadcs().setE(hitlist.dnhits[i].E);
         fadcs().setT(hitlist.dnhits[i].t);
      }
   }

   // Create bcalTDCHit structures to hold our F1TDC hits
   map<int, TDCHitList>::iterator ittdc;
   for (ittdc = F1TDCHits.begin(); ittdc != F1TDCHits.end(); ittdc++) {
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
      
      // Add TDC hits to list
      for (unsigned int i = 0; i < hitlist.uphits.size(); ++i) {
         hddm_s::BcalTDCHitList tdcs = iter->addBcalTDCHits();
         tdcs().setEnd(bcal_index::kUp);
         tdcs().setT(hitlist.uphits[i]);
      }
      for (unsigned int i = 0; i < hitlist.dnhits.size(); i++) {
         hddm_s::BcalTDCHitList tdcs = iter->addBcalTDCHits();
         tdcs().setEnd(bcal_index::kDown);
         tdcs().setT(hitlist.dnhits[i]);
      }
   }
}

//-----------
// GetHistoFromPool
//-----------
DHistogram* GetHistoFromPool(void)
{
   DHistogram *h = NULL;

   pthread_mutex_lock(&histo_pool_mutex);
   if(histo_pool.size()>0){
      // Take histo out of pool
      h = histo_pool.front();
      histo_pool.pop();
   }
   pthread_mutex_unlock(&histo_pool_mutex);
   
   // Do this outside of mutex lock to speed things up
   if(h==NULL){
      // There wasn't a histo in the pool. Create a new one
      h = new DHistogram(4000, -100.0, 300.0);
   }else{
      // Histo came from pool. Reset it
      h->Reset();
   }

   return h;
}

//-----------
// ReturnHistoToPool
//-----------
void ReturnHistoToPool(DHistogram *h)
{
   pthread_mutex_lock(&histo_pool_mutex);
   if(h!=NULL){
      // Limit size of pool
      if(histo_pool.size() >= 200){
         delete h;
      }else{
         histo_pool.push(h);
      }
   }
   pthread_mutex_unlock(&histo_pool_mutex);
}

//---------------
// MakeTSpline
//---------------
TSpline* MakeTSpline(void)
{
   // These points are taken by mapping out figure 6
   // of GlueX-doc-1795 (the "5ns" rise time pulse).

   double t[] ={
      -7.00,
      -4.92,
      -0.75,
      4.24,
      9.37,      // 5
      
      9.71,
      10.05,
      10.39,
      10.73,
      
      11.07,
      12.03,
      12.92,
      14.08,
      16.05,   // 10
      18.56,
      21.22,
      24.31,
      28.44,
      32.92,   // 15
      38.64,
      42.07,
      //47.10,
      //52.59,
      //58.22,   // 20
      67.95,
      71.74,
      79.85,
      120.0,
      125.0,   // 25
      130.0,
      160.0,
      180.0,
      -1000.0 // flag end of data
   };
   
   double mV[]={
      0,
      0,
      0,
      0,
      0-1,               // 5

      -1.2,
      -1.4,
      -1.6,
      -1.8,

      0-2,
      -49.16039931,
      -159.0483507,
      -375.9324653,
      -711.3798959,   // 10
      -902.2379167,
      -988.9915625,
      -1020.801233,
      -960.0736806,
      -850.1857292,   // 15
      -694.0291667,
      -601.4919445,
      //-539.29,
      //-419.88,
      //-300.46,         // 20
      -142.53,
      -100.15,
      -61.63,
      -8.0,
      -4.0,            // 25
      -2.0,
      -1.0,
      0.0,
      -1000.0 // flag end of data
   };

   // Find number points
   int Npoints=0;
   while(true){
      if(t[Npoints]==-1000.0 && mV[Npoints]==-1000.0)break;
      Npoints++;
   }
   
   TSpline *s = new TSpline3("spline", t, mV, Npoints);
   
   return s;
}

//-----------
// SaveDebugHistos
//-----------
void SaveDebugHistos(int32_t eventNumber, map<bcal_index, CellSpectra> &SiPMspectra, const char *suffix, const char *yunits)
{
   /// Create ROOT histograms out of the CellSpectra objects. These would
   /// correspond to the individual SiPMs. These will be saved in the smear.root 
   /// file inside of a directory structure to allow debugging.
   
   // This should probably only be called when running single-threaded
   // (better safe than sorry)!
   pthread_mutex_lock(&root_mutex);

   // Save the current ROOT directory so we can restore it before returning
   TDirectory *savedir = gDirectory;
   
   // See if TDirectory for this event exists. If not, create one.
   char eventDirName[256];
   sprintf(eventDirName, "Event%04d", eventNumber);
   TDirectory *eventdir = (TDirectory*)gDirectory->FindObject(eventDirName);
   if(eventdir == NULL){
      eventdir = gDirectory->mkdir(eventDirName);
   }
   eventdir->cd();
   
   map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
   for(; iter!=SiPMspectra.end(); iter++){
      const bcal_index &idx = iter->first;
      CellSpectra &cellspectra = iter->second;
      
      // Skip empty spectra objects
      if(cellspectra.hup==NULL && cellspectra.hdn==NULL)continue;
      bool allempty=true;
      if(cellspectra.hup && cellspectra.hup->Integral()!=0.0)allempty = false;
      if(cellspectra.hdn && cellspectra.hdn->Integral()!=0.0)allempty = false;
      if(allempty)continue;
      
      // Extract location info
      int module = idx.module;
      int layer = idx.layer;
      int sector = idx.sector;
      int incident_id = idx.incident_id;
      
      // Create a directory to hold this readout cell's histos
      char dirname[256];
      sprintf(dirname, "SiPM_m%02dl%ds%dip%d%s", module, layer, sector, incident_id, suffix);
      
      // Directory may already exist since separate entries may be kept
      // for upstream and downstream. Check if it exists first and only
      // create it if necessary.
      eventdir->cd();
      TDirectory *rcdir = (TDirectory*)eventdir->FindObject(dirname);
      if(!rcdir)rcdir = eventdir->mkdir(dirname);
      rcdir->cd();
      
      // Make ROOT histograms out of upstream and downstream histos
      char hname[256];
      if(cellspectra.hup && cellspectra.hup->Integral()!=0.0){
         sprintf(hname, "Upstream summed cell mod=%d lay=%d sec=%d", module, layer, sector);
         TH1D *rh = cellspectra.hup->MakeTH1D("hup", hname);
         rh->SetXTitle("Time (ns)");
         rh->SetYTitle(yunits);
      }
      if(cellspectra.hdn && cellspectra.hdn->Integral()!=0.0){
         sprintf(hname, "Downstream summed cell mod=%d lay=%d sec=%d", module, layer, sector);
         TH1D *rh = cellspectra.hdn->MakeTH1D("hdn", hname);
         rh->SetXTitle("Time (ns)");
         rh->SetYTitle(yunits);
      }
      
   }
   
   // Restore ROOT working directory
   savedir->cd();

   pthread_mutex_unlock(&root_mutex);
   
}

//-----------
// SaveDebugHistos
//-----------
void SaveDebugHistos(int32_t eventNumber, map<int, SumSpectra> &bcalfADC, const char *suffix)
{
   /// Create ROOT histograms out of the histograms in the SumSpectra
   /// objects. These can be either electronics pulse shapes or attenuated
   /// energy spectra (depending upon what stage of the smearing has been
   /// done to the objects given us). These will be saved in the smear.root 
   /// file inside of a directory structure to allow debugging.
   ///
   /// WARNING: This will save several histograms for every event!
   /// Usually, you will only want to run mcsmear for one or two
   /// events when this is being called.
   
   // This should probably only be called when running single-threaded
   // (better safe than sorry)!
   pthread_mutex_lock(&root_mutex);

   // Save the current ROOT directory so we can restore it before returning
   TDirectory *savedir = gDirectory;

   // See if TDirectory for this event exists. If not, create one.
   char eventDirName[256];
   sprintf(eventDirName, "Event%04d", eventNumber);
   TDirectory *eventdir = (TDirectory*)gDirectory->FindObject(eventDirName);
   if(eventdir == NULL){
      eventdir = gDirectory->mkdir(eventDirName);
   }
   eventdir->cd();
   
   map<int, SumSpectra>::iterator iter=bcalfADC.begin();
   for(; iter!=bcalfADC.end(); iter++){
      int fADCId = iter->first;
      SumSpectra &sumspectra = iter->second;
      
      // Skip empty spectra objects
      if(sumspectra.hup==NULL && sumspectra.hdn==NULL)continue;
      bool allempty=true;
      if(sumspectra.hup && sumspectra.hup->Integral()!=0.0)allempty = false;
      if(sumspectra.hdn && sumspectra.hdn->Integral()!=0.0)allempty = false;
      if(allempty)continue;
      
      // Extract location info
      int module = DBCALGeometry::module(fADCId);
      int fADC_layer = DBCALGeometry::layer(fADCId);
      int fADC_sector = DBCALGeometry::sector(fADCId);
      
      // Create a directory to hold this readout cell's histos
      char dirname[256];
      sprintf(dirname, "Sum_m%02dl%ds%d%s", module, fADC_layer, fADC_sector, suffix);
      TDirectory *rcdir = eventdir->mkdir(dirname);
      rcdir->cd();
      
      // Make ROOT histograms out of upstream and downstream histos
      char hname[256];
      if(sumspectra.hup && sumspectra.hup->Integral()!=0.0){
         sprintf(hname, "Upstream summed cell mod=%d lay=%d sec=%d", module, fADC_layer, fADC_sector);
         TH1D *rh = sumspectra.hup->MakeTH1D("hup", hname);
         rh->SetXTitle("Time (ns)");
         rh->SetYTitle("Amplitude (mV)");
      }
      if(sumspectra.hdn && sumspectra.hdn->Integral()!=0.0){
         sprintf(hname, "Downstream summed cell mod=%d lay=%d sec=%d", module, fADC_layer, fADC_sector);
         TH1D *rh = sumspectra.hdn->MakeTH1D("hdn", hname);
         rh->SetXTitle("Time (ns)");
         rh->SetYTitle("Amplitude (mV)");
      }
      
   }
   
   // Restore ROOT working directory
   savedir->cd();

   pthread_mutex_unlock(&root_mutex);
}

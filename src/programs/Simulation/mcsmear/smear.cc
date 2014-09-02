// $Id$
//
// Created June 22, 2005  David Lawrence

#define JANA_ENABLED 1

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <FCAL/DFCALGeometry.h>
#include <CCAL/DCCALGeometry.h>
#include <BCAL/DBCALGeometry.h>

#include <math.h>
#include "units.h"
#include "HDDM/hddm_s.hpp"
#include <TF1.h>
#include <TH2.h>

#include "DRandom2.h"

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif


extern vector<vector<float> >fdc_smear_parms; 
extern TH2F *fdc_drift_time_smear_hist;
extern TH2F *fdc_drift_dist_smear_hist;
extern TH2F *fdc_drift_time;
extern TH1F *fdc_cathode_charge;
extern TH2F *cdc_drift_time;
extern TH1F *cdc_charge;
extern TH1F *fdc_anode_mult;
extern TH2F *cdc_drift_smear;

void GetAndSetSeeds(hddm_s::HDDM *record);
void SmearCDC(hddm_s::HDDM *record);
void AddNoiseHitsCDC(hddm_s::HDDM *record);
void SmearFDC(hddm_s::HDDM *record);
void AddNoiseHitsFDC(hddm_s::HDDM *record);
void SmearFCAL(hddm_s::HDDM *record);
void SmearCCAL(hddm_s::HDDM *record);
void SmearBCAL(hddm_s::HDDM *record);
void SmearTOF(hddm_s::HDDM *record);
void SmearSTC(hddm_s::HDDM *record);
void SmearCherenkov(hddm_s::HDDM *record);
void SmearTAGM(hddm_s::HDDM *record);
void SmearTAGH(hddm_s::HDDM *record);
void InitCDCGeometry(void);
void InitFDCGeometry(void);

pthread_mutex_t mutex_fdc_smear_function = PTHREAD_MUTEX_INITIALIZER;

bool CDC_GEOMETRY_INITIALIZED = false;
int CDC_MAX_RINGS=0;
extern vector<unsigned int> NCDC_STRAWS;
extern vector<double> CDC_RING_RADIUS;

DFCALGeometry *fcalGeom = NULL;
DCCALGeometry *ccalGeom = NULL;
bool FDC_GEOMETRY_INITIALIZED = false;
unsigned int NFDC_WIRES_PER_PLANE;
extern vector<double> FDC_LAYER_Z;
extern double FDC_RATE_COEFFICIENT;

double SampleGaussian(double sigma);
double SamplePoisson(double lambda);
double SampleRange(double x1, double x2);

// Do we or do we not add noise hits
extern bool ADD_NOISE;

// Do we or do we not smear real hits
extern bool SMEAR_HITS;

// Use or ignore random number seeds found in HDDM file
extern bool IGNORE_SEEDS;

// The error on the drift time in the CDC. The drift times
// for the actual CDC hits coming from the input file
// are smeared by a gaussian with this sigma.
extern double CDC_TDRIFT_SIGMA;

// The time window for which CDC hits are accumulated.
// This is used to determine the number of background
// hits in the CDC for a given event.
extern double CDC_TIME_WINDOW;

// The error in the energy deposition measurement in the CDC due to pedestal noise
extern double CDC_PEDESTAL_SIGMA;

// Number of sigmas above threshold for sparcification of CDC data
extern double CDC_THRESHOLD_FACTOR;
 
// If the following flag is true, then include the drift-distance
// dependency on the error in the FDC position. Otherwise, use a
// flat distribution given by the FDC_TDRIFT_SIGMA below.
extern bool FDC_USE_PARAMETERIZED_SIGMA;

// The error on the drift time in the FDC. The drift times
// for the actual FDC hits coming from the input file
// are smeared by a gaussian with this sigma.
extern double FDC_TDRIFT_SIGMA;

// The error in the distance along the wire as measured by
// the cathodes. This should NOT include the Lorentz
// effect which is already included in hdgeant. It
// should include any fluctuations due to ion trail density
// etc.
extern double FDC_CATHODE_SIGMA;

// The FDC pedestal noise is used to smear the cathode ADC
// values such that the position along the wire has the resolution
// specified by FDC_CATHODE_SIGMA.
extern double FDC_PED_NOISE; //pC (calculated from FDC_CATHODE_SIGMA in SmearFDC)

// Number of sigmas above threshold for sparcification of FDC data
extern double FDC_THRESHOLD_FACTOR;
 
// If energy loss was turned off in the FDC then the pedestal
// noise will not be scaled properly to give the nominal 200 micron
// resolution along the wires. This flag is used to indicated
// the magic scale factor should be applied to FDC_PED_NOISE
// when it is calculated below to get the correct resolution.
extern bool FDC_ELOSS_OFF;

// Time window for acceptance of FDC hits
extern double FDC_TIME_WINDOW;

// Fraction of FDC hits to randomly drop (0=drop nothing 1=drop everything)
extern double FDC_HIT_DROP_FRACTION;

// FDC discriminator threshold in keV
extern double FDC_THRESH_KEV;

// Single FTOF bar energy threshold (applied after smearing)
extern double FTOF_BAR_THRESHOLD;

// Single STC paddle energy threshold (applied after smearing)
extern double STC_PADDLE_THRESHOLD;

// Tagging counter time resolutions (s)
extern double TAGM_TSIGMA;
extern double TAGH_TSIGMA;

// Tagging counter fADC responses
extern double TAGM_FADC_TSIGMA;
extern double TAGH_FADC_TSIGMA;
extern double TAGM_NPIX_PER_GEV;
extern double TAGH_NPE_PER_GEV;

// Photon-statistics factor for smearing hit energy (from Criss's MC)
// (copied from DFCALMCResponse_factory.cc 7/2/2009 DL)
extern double FCAL_PHOT_STAT_COEF;

// Single block energy threshold (applied after smearing)
extern double FCAL_BLOCK_THRESHOLD;

// Photon-statistics factor for smearing hit energy for CompCal
// (This is just a rough estimate 11/30/2010 DL)
double CCAL_PHOT_STAT_COEF = 0.035/2.0;

// Single block energy threshold (applied after smearing)
// (This is just a rough estimate 11/30/2010 DL)
double CCAL_BLOCK_THRESHOLD = 20.0*k_MeV;


// Forward TOF resolution
extern double TOF_SIGMA;
extern double TOF_PHOTONS_PERMEV;

// Start counter resolution
extern double START_SIGMA ;
extern double START_PHOTONS_PERMEV;

extern bool SMEAR_BCAL;

// Beginning of readout window
extern double TRIGGER_LOOKBACK_TIME;

extern bool DROP_TRUTH_HITS;

// Polynomial interpolation on a grid.
// Adapted from Numerical Recipes in C (2nd Edition), pp. 121-122.
void polint(float *xa, float *ya,int n,float x, float *y,float *dy){
  int i,m,ns=0;
  float den,dif,dift,ho,hp,w;

  float *c=(float *)calloc(n,sizeof(float));
  float *d=(float *)calloc(n,sizeof(float));

  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++){
    if ((dift=fabs(x-xa[i]))<dif){
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];

  for (m=1;m<n;m++){
    for (i=1;i<=n-m;i++){
      ho=xa[i-1]-x;
      hp=xa[i+m-1]-x;
      w=c[i+1-1]-d[i-1];
      if ((den=ho-hp)==0.0){
   free(c);
   free(d);
   return;
      }
  
      den=w/den;
      d[i-1]=hp*den;
      c[i-1]=ho*den;
      
    }
    
    *y+=(*dy=(2*ns<(n-m) ?c[ns+1]:d[ns--]));
  }
  free(c);
  free(d);
}

//-----------
// Smear
//-----------
void Smear(hddm_s::HDDM *record)
{
   GetAndSetSeeds(record);

   if (SMEAR_HITS) {
      SmearCDC(record);
      SmearFDC(record);
      SmearFCAL(record);
      SmearCCAL(record);
      SmearTOF(record);
      SmearSTC(record);
      SmearCherenkov(record);
      SmearTAGM(record);
      SmearTAGH(record);
   }
   if (ADD_NOISE) {
      AddNoiseHitsCDC(record);
      AddNoiseHitsFDC(record);
   }
   if (SMEAR_BCAL) 
      SmearBCAL(record);
}

//-----------
// SetSeeds
//-----------
void SetSeeds(const char *vals)
{
   /// This is called from the command line parser to
   /// set the initial seeds based on user input from
   /// the command line.
   //
   //
   stringstream ss(vals);
   Int_t seed1, seed2, seed3;
   ss >> seed1 >> seed2 >> seed3;
   UInt_t *useed1 = reinterpret_cast<UInt_t*>(&seed1);
   UInt_t *useed2 = reinterpret_cast<UInt_t*>(&seed2);
   UInt_t *useed3 = reinterpret_cast<UInt_t*>(&seed3);
   gDRandom.SetSeeds(*useed1, *useed2, *useed3);

   cout << "Seeds set from command line. Any random number" << endl;
   cout << "seeds found in the input file will be ignored!" << endl;
   IGNORE_SEEDS = true;
}

//-----------
// GetAndSetSeeds
//-----------
void GetAndSetSeeds(hddm_s::HDDM *record)
{
   // Check if non-zero seed values exist in the input HDDM file.
   // If so, use them to set the seeds for the random number
   // generator. Otherwise, make sure the seeds that are used
   // are stored in the output event.
   
   if (record == 0)
      return;
   else if (record->getReactions().size() == 0)
      return;

   hddm_s::ReactionList::iterator reiter = record->getReactions().begin();
   if (reiter->getRandoms().size() == 0) {
      // No seeds stored in event. Add them
      hddm_s::RandomList blank_rand = reiter->addRandoms();
      blank_rand().setSeed1(0);
      blank_rand().setSeed2(0);
      blank_rand().setSeed3(0);
      blank_rand().setSeed4(0);
   }

   UInt_t seed1, seed2, seed3;
   hddm_s::Random my_rand = reiter->getRandom();

   if (!IGNORE_SEEDS) {
      // Copy seeds from event record to local variables
      seed1 = my_rand.getSeed1();
      seed2 = my_rand.getSeed2();
      seed3 = my_rand.getSeed3();
      
      // If the seeds in the event are all zeros it means they
      // were not set. In this case, initialize seeds to constants
      // to guarantee the seeds are used if this input file were
      // smeared again with the same command a second time. These
      // are set here to the fractional part of the cube roots of
      // the first three primes, truncated to 9 digits.
      if ((seed1 == 0) || (seed2 == 0) || (seed3 == 0)){
         int eventNo = record->getPhysicsEvent().getEventNo();
         seed1 = 259921049 + eventNo;
         seed2 = 442249570 + eventNo;
         seed3 = 709975946 + eventNo;
      }
      
      // Set the seeds in the random generator.
      gDRandom.SetSeeds(seed1, seed2, seed3);
   }

   // Copy seeds from generator to local variables
   gDRandom.GetSeeds(seed1, seed2, seed3);

   // Copy seeds from local variables to event record
   my_rand.setSeed1(seed1);
   my_rand.setSeed2(seed2);
   my_rand.setSeed3(seed3);
}

//-----------
// SmearCDC
//-----------
void SmearCDC(hddm_s::HDDM *record)
{
   /// Smear the drift times of all CDC hits.
   /// This will add cdcStrawHit objects generated by smearing values in the
   /// cdcStrawTruthHit objects that hdgeant outputs. Any existing cdcStrawHit
   /// objects will be replaced.

   double t_max = TRIGGER_LOOKBACK_TIME + CDC_TIME_WINDOW;
   double threshold = CDC_THRESHOLD_FACTOR * CDC_PEDESTAL_SIGMA; // for sparcification

   // Loop over all cdcStraw tags
   hddm_s::CdcStrawList straws = record->getCdcStraws();
   hddm_s::CdcStrawList::iterator iter;
   for (iter = straws.begin(); iter != straws.end(); ++iter) {
 
      // If the element already contains a cdcStrawHit list then delete it.
      hddm_s::CdcStrawHitList hits = iter->getCdcStrawHits();
      if (hits.size() > 0) {
         static bool warned = false;
         iter->deleteCdcStrawHits();
         if (!warned) {
            warned = true;
            cerr << endl;
            cerr << "WARNING: CDC hits already exist in input file! Overwriting!"
                 << endl << endl;
         }
      }

      // Create new cdcStrawHit from cdcStrawTruthHit information
      hddm_s::CdcStrawTruthHitList thits = iter->getCdcStrawTruthHits();
      hddm_s::CdcStrawTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++ titer) {
         // Pedestal-smeared charge
         double q = titer->getQ() + SampleGaussian(CDC_PEDESTAL_SIGMA);
         // Smear out the CDC drift time using the specified sigma.
         // This is for timing resolution from the electronics;
         // diffusion is handled in hdgeant.
         double t = titer->getT() + SampleGaussian(CDC_TDRIFT_SIGMA)*1.0e9;
         if (t > TRIGGER_LOOKBACK_TIME && t < t_max && q > threshold) {
            if (iter->getRing() == 1) {
               double t_corr = t-0.33;
               cdc_drift_time->Fill(t_corr, titer->getD());
               cdc_drift_smear->Fill(t_corr, titer->getD()-
                                             (0.0285*sqrt(t_corr)+0.014));
               cdc_charge->Fill(q);
            }
            hits = iter->addCdcStrawHits();
            hits().setT(t);
            hits().setQ(q);
         }

         if (DROP_TRUTH_HITS) {
            iter->deleteCdcStrawTruthHits();
         }
      }
   }
}

//-----------
// AddNoiseHitsCDC
//-----------
void AddNoiseHitsCDC(hddm_s::HDDM *record)
{
   // Calculate the number of noise hits for each straw and store
   // them in a sparse map, then copy into the output hddm record.
   //
   // The straw rates are obtained using a parameterization done
   // to calculate the event size for the August 29, 2007 online
   // meeting. This parameterization is almost already obsolete.
   // 10/12/2007 D. L.

#if ! JANA_ENABLED
   if(!CDC_GEOMETRY_INITIALIZED)InitCDCGeometry();
#endif   

   vector<int> Nstraw_hits;
   vector<int> straw_number;
   vector<int> ring_number;
   int Nnoise_straws = 0;
   int Nnoise_hits = 0;
   for(unsigned int ring=1; ring <= NCDC_STRAWS.size(); ring++){
      double p[2] = {10.4705, -0.103046};
      double r_prime = (double)(ring+3);
      double N = exp(p[0] + r_prime*p[1]);
      N *= CDC_TIME_WINDOW;
      for (unsigned int straw=1; straw<=NCDC_STRAWS[ring-1]; straw++) {
         // Indivdual straw rates should be way less than 1/event so
         // we just use the rate as a probablity.
         double Nhits = SampleRange(0.0, 1.0)<N ? 1.0:0.0;
         if(Nhits<1.0)continue;
         int iNhits = (int)floor(Nhits);
         Nstraw_hits.push_back(iNhits);
         straw_number.push_back(straw);
         ring_number.push_back(ring);
         Nnoise_straws++;
         Nnoise_hits+=iNhits;
      }
   }

   double t_max = TRIGGER_LOOKBACK_TIME + CDC_TIME_WINDOW;
   double threshold = CDC_THRESHOLD_FACTOR * CDC_PEDESTAL_SIGMA; // for sparcification
   
   // Loop over straws with noise hits
   hddm_s::CentralDCList cdc = record->getCentralDCs();
   hddm_s::CdcStrawList straws = record->getCdcStraws();
   for (unsigned int j=0; j < Nstraw_hits.size(); j++) {
      hddm_s::CdcStrawList::iterator iter;
      for (iter = straws.begin(); iter != straws.end(); ++iter) {
         if (iter->getRing() == ring_number[j] &&
             iter->getStraw() == straw_number[j])
            break;
      }
      if (iter == straws.end()) {
         if (cdc.size() == 0)
            cdc = record->getHitViews().begin()->addCentralDCs();
         iter = cdc().addCdcStraws().begin();
         iter->setRing(ring_number[j]);
         iter->setStraw(straw_number[j]);
      }
      for (int k=0; k < Nstraw_hits[j]; k++) {
         double q = SampleGaussian(CDC_PEDESTAL_SIGMA);
         if (q > threshold) {
            hddm_s::CdcStrawHitList hits = iter->addCdcStrawHits();
            hits().setQ(q);
            hits().setT(SampleRange(TRIGGER_LOOKBACK_TIME,t_max));
            cdc_charge->Fill(hits().getQ());   
            cdc_drift_time->Fill(hits().getT(), 0.);
         }
      }
   }
}

//-----------
// SmearFDC
//-----------
void SmearFDC(hddm_s::HDDM *record)
{
   // Calculate ped noise level based on position resolution
   //   FDC_PED_NOISE = -0.004594 + 0.008711*FDC_CATHODE_SIGMA +
   //                    0.000010*FDC_CATHODE_SIGMA*FDC_CATHODE_SIGMA; //pC
   FDC_PED_NOISE = -0.0938 + 0.0485*FDC_CATHODE_SIGMA;
   if (FDC_ELOSS_OFF)
      FDC_PED_NOISE *= 7.0; // empirical  4/29/2009 DL

   double t_max = TRIGGER_LOOKBACK_TIME + FDC_TIME_WINDOW;
   double threshold = FDC_THRESHOLD_FACTOR * FDC_PED_NOISE; // for sparcification

   hddm_s::FdcChamberList chambers = record->getFdcChambers();
   hddm_s::FdcChamberList::iterator iter;
   for (iter = chambers.begin(); iter != chambers.end(); ++iter) {

      // Add pedestal noise to strip charge data
      hddm_s::FdcCathodeStripList strips = iter->getFdcCathodeStrips();
      hddm_s::FdcCathodeStripList::iterator siter;
      for (siter = strips.begin(); siter != strips.end(); ++siter) {
          // If a fdcCathodeHit already exists delete it
          siter->deleteFdcCathodeHits();
          hddm_s::FdcCathodeTruthHitList thits = 
                                         siter->getFdcCathodeTruthHits();
          hddm_s::FdcCathodeTruthHitList::iterator titer;
          for (titer = thits.begin(); titer != thits.end(); ++titer) {
            //if (SampleRange(0.0, 1.0) <= FDC_HIT_DROP_FRACTION)
            //   continue;
            double q = titer->getQ() + SampleGaussian(FDC_PED_NOISE);
            double t = titer->getT() +
                       SampleGaussian(FDC_TDRIFT_SIGMA)*1.0e9;
            if (q > threshold && t > TRIGGER_LOOKBACK_TIME && t < t_max) {
               hddm_s::FdcCathodeHitList hits = siter->addFdcCathodeHits();
               hits().setQ(q);
               hits().setT(t);
            }
            fdc_cathode_charge->Fill(q);
         }

         if (DROP_TRUTH_HITS)
            siter->deleteFdcCathodeTruthHits();
      }

      // Add drift time varation to the anode data 
      hddm_s::FdcAnodeWireList wires = iter->getFdcAnodeWires();
      hddm_s::FdcAnodeWireList::iterator witer;
      for (witer = wires.begin(); witer != wires.end(); ++witer) {
         // If a fdcAnodeHit exists already delete it
         witer->deleteFdcAnodeHits();
         hddm_s::FdcAnodeTruthHitList thits = witer->getFdcAnodeTruthHits();
         hddm_s::FdcAnodeTruthHitList::iterator titer;
         for (titer = thits.begin(); titer != thits.end(); ++titer) {
            double t = titer->getT() + SampleGaussian(FDC_TDRIFT_SIGMA)*1.0e9;
            if (t > TRIGGER_LOOKBACK_TIME && t < t_max) {
               hddm_s::FdcAnodeHitList hits = witer->addFdcAnodeHits();
               hits().setT(t);
               hits().setDE(titer->getDE());
            }
            if (witer->getLayer() == 1 && witer->getModule() == 1) {
               // 3.7 ns flight time for c=1 to first fdc plane
               double dt = t - titer->getT() - 3.7;
               fdc_drift_time_smear_hist->Fill(0., dt);
               fdc_drift_dist_smear_hist->Fill(0.,
                         dt*(0.5*0.02421/sqrt(titer->getT())+5.09e-4));
               fdc_drift_time->Fill(t - 3.7, 0.);
            }
         }
         fdc_anode_mult->Fill(witer->getFdcAnodeHits().size());

         if (DROP_TRUTH_HITS)
            witer->deleteFdcAnodeTruthHits();
      }
   }
}

//-----------
// AddNoiseHitsFDC
//-----------
void AddNoiseHitsFDC(hddm_s::HDDM *record)
{
   // Calculate the number of noise hits for each FDC wire and store
   // them in a sparse map.
   //
   // We do this using the individual wire rates to calculate the probability
   // of the wire firing for a single event. For the FDC, we calculate the
   // wire rates as a function of both wire number (distance from beam line)
   // and layer (position in z). We want a roughly 1/r distribution in the
   // radial direction and a roughly exponential rise in rate in the
   // +z direction.
   //
   // The wire rates are obtained using a parameterization done
   // to calculate the event size for the August 29, 2007 online
   // meeting. This parameterization is almost already obsolete.
   // In rough terms, the layer rate (integrated over all wires)
   // is about 1 MHz. For a 24 layer chamber with a 1us time window,
   // we should have approximately 24 background hits per event.
   // 11/9/2007 D. L.
   //
   // Note by RTJ: Where are the noise hits in the FDC cathode strips??
   // From the Ncathode_hits vector, it seems that this was originally
   // intended to be added, once the wire noise model is validated.
   // The FDC noise algorithm represented below is incomplete.
  
#if ! JANA_ENABLED
   if(!FDC_GEOMETRY_INITIALIZED)InitFDCGeometry();
#endif   

   vector<int> Nwire_hits;
   vector<int> Ncathode_hits;
   vector<int> wire_number;
   vector<int> layer_number;
   int Nnoise_wires = 0;
   int Nnoise_hits = 0;
   for (unsigned int layer=1; layer <= FDC_LAYER_Z.size(); layer++) {
      double No = FDC_RATE_COEFFICIENT*exp((double)layer*log(4.0)/24.0);
      for (unsigned int wire=1; wire <= 96; wire++) {
         double rwire = fabs(96.0/2.0 - (double)wire);
         double N = No*log((rwire+0.5)/(rwire-0.5));

         // Indivdual wire rates should be way less than 1/event so
         // we just use the rate as a probablity.
         double Nhits = (SampleRange(0.0, 1.0) < N)? 1.0 : 0.0;
         if (Nhits < 1.0)
            continue;
         int iNhits = (int)floor(Nhits);
         Nwire_hits.push_back(iNhits);
         wire_number.push_back(wire);
         layer_number.push_back(layer);
         Nnoise_wires++;
         Nnoise_hits+=iNhits;
      }
   }

   double t_max = TRIGGER_LOOKBACK_TIME + FDC_TIME_WINDOW;
   //double threshold = FDC_THRESHOLD_FACTOR * FDC_PED_NOISE; // for sparcification

   hddm_s::ForwardDCList fdc = record->getForwardDCs();
   hddm_s::FdcCathodeStripList strips = record->getFdcCathodeStrips();

   // Loop over wires with noise hits
   for (unsigned int j=0; j < Nwire_hits.size(); j++) {
      hddm_s::FdcAnodeWireList wires = record->getFdcAnodeWires();
      hddm_s::FdcAnodeWireList::iterator witer;
      for (witer = wires.begin(); witer != wires.end(); ++witer) {
         int layerId = witer->getLayer() + (witer->getModule() - 1)*8;
         if (layer_number[j] == layerId && wire_number[j] == witer->getWire())
            break;
      }
      if (witer == wires.end()) {
         if (fdc.size() == 0)
            fdc = record->getHitViews().begin()->addForwardDCs();
         hddm_s::FdcChamberList chambers = record->getFdcChambers();
         hddm_s::FdcChamberList::iterator citer;
         for (citer = chambers.begin(); citer != chambers.end(); ++citer) {
            int layerId = witer->getLayer() + (witer->getModule() - 1)*8;
            if (layer_number[j] == layerId)
               break;
         }
         if (citer == chambers.end()) {
            citer = fdc().addFdcChambers().begin();
            citer->setModule((layer_number[j] - 1)/3 + 1);
            citer->setLayer((layer_number[j] - 1)%3 + 1);
         }
         witer = citer->addFdcAnodeWires().begin();
         witer->setWire(wire_number[j]);
      }

      double dEsigma=FDC_THRESH_KEV/FDC_THRESHOLD_FACTOR;
      for (int k=0; k < Nwire_hits[j]; k++) {
         // Simulated random hit as pedestal noise 
         double dE = SampleGaussian(dEsigma);
         if (dE > FDC_THRESH_KEV) {
            hddm_s::FdcAnodeHitList hits = witer->addFdcAnodeHits();
            hits().setDE(dE); // what should this be?
            hits().setT(SampleRange(TRIGGER_LOOKBACK_TIME, t_max));
            fdc_drift_time->Fill(hits().getT(), 0.);
         }
      }
   }
}

//-----------
// SmearFCAL
//-----------
void SmearFCAL(hddm_s::HDDM *record)
{
   /// Smear the FCAL hits using the nominal resolution of the individual blocks.
   /// The way this works is a little funny and warrants a little explanation.
   /// The information coming from hdgeant is truth information indexed by 
   /// row and column, but containing energy deposited and time. The mcsmear
   /// program will copy the truth information from the fcalTruthHit element
   /// to a new fcalHit element, smearing the values with the appropriate detector
   /// resolution.
   ///
   /// To access the "truth" values in DANA, get the DFCALHit objects using the
   /// "TRUTH" tag.
   
   if (!fcalGeom)
      fcalGeom = new DFCALGeometry();

   hddm_s::FcalBlockList blocks = record->getFcalBlocks();
   hddm_s::FcalBlockList::iterator iter;
   for (iter = blocks.begin(); iter != blocks.end(); ++iter) {
      iter->deleteFcalHits();
      hddm_s::FcalTruthHitList thits = iter->getFcalTruthHits();
      hddm_s::FcalTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // Simulation simulates a grid of blocks for simplicity. 
         // Do not bother smearing inactive blocks. They will be
         // discarded in DEventSourceHDDM.cc while being read in
         // anyway.
         if (!fcalGeom->isBlockActive(iter->getRow(), iter->getColumn()))
            continue;
         // Smear the energy and timing of the hit
         double sigma = FCAL_PHOT_STAT_COEF/sqrt(titer->getE());
         double E = titer->getE() * (1.0 + SampleGaussian(sigma));
         // Smear the time by 200 ps (fixed for now) 7/2/2009 DL
         double t = titer->getT() + SampleGaussian(200.0e-3);
         // Apply a single block threshold. 
         if (E >= FCAL_BLOCK_THRESHOLD) {
            hddm_s::FcalHitList hits = iter->addFcalHits();
            hits().setE(E);
            hits().setT(t);
         }
      }

      if (DROP_TRUTH_HITS)
         iter->deleteFcalTruthHits();
   }
}

//-----------
// SmearCCAL
//-----------
void SmearCCAL(hddm_s::HDDM *record)
{
   /// Smear the CCAL hits using the same procedure as the FCAL above.
   /// See those comments for details.

   if (!ccalGeom)
      ccalGeom = new DCCALGeometry();

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
         double sigma = CCAL_PHOT_STAT_COEF/sqrt(titer->getE()) ;
         double E = titer->getE() * (1.0 + SampleGaussian(sigma));
         // Smear the time by 200 ps (fixed for now) 7/2/2009 DL
         double t = titer->getT() + SampleGaussian(200.0e-3);
         // Apply a single block threshold. If the (smeared) energy is below this,
         // then set the energy and time to zero. 
         if (E > CCAL_BLOCK_THRESHOLD) {
            hddm_s::CcalHitList hits = iter->addCcalHits();
            hits().setE(E);
            hits().setT(t);
         }
      }

      if (DROP_TRUTH_HITS)
         iter->deleteCcalTruthHits();
   }
}

//-----------
// SmearTOF
//-----------
void SmearTOF(hddm_s::HDDM *record)
{
   hddm_s::FtofCounterList tofs = record->getFtofCounters();
   hddm_s::FtofCounterList::iterator iter;
   for (iter = tofs.begin(); iter != tofs.end(); ++iter) {
      // take care of hits
      iter->deleteFtofHits();
      hddm_s::FtofTruthHitList thits = iter->getFtofTruthHits();
      hddm_s::FtofTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // Smear the time
         double t = titer->getT() + SampleGaussian(TOF_SIGMA);
         // Smear the energy
         double npe = titer->getDE() * 1000. * TOF_PHOTONS_PERMEV;
         npe = npe +  SampleGaussian(sqrt(npe));
         float NewE = npe/TOF_PHOTONS_PERMEV/1000.;
         if (NewE > FTOF_BAR_THRESHOLD) {
            hddm_s::FtofHitList hits = iter->addFtofHits();
            hits().setEnd(titer->getEnd());
            hits().setT(t);
            hits().setDE(NewE);
         }
      }
    
      if (DROP_TRUTH_HITS) {
         iter->deleteFtofTruthHits();
      }
   }
}

//-----------
// SmearSTC - smear hits in the start counter
//-----------
void SmearSTC(hddm_s::HDDM *record)
{
   hddm_s::StcPaddleList pads = record->getStcPaddles();
   hddm_s::StcPaddleList::iterator iter;
   for (iter = pads.begin(); iter != pads.end(); ++iter) {
      iter->deleteStcHits();
      hddm_s::StcTruthHitList thits = iter->getStcTruthHits();
      hddm_s::StcTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time
         double t = titer->getT() + SampleGaussian(START_SIGMA);
         // smear the energy
         double npe = titer->getDE() * 1000. *  START_PHOTONS_PERMEV;
         npe = npe +  SampleGaussian(sqrt(npe));
         double NewE = npe/START_PHOTONS_PERMEV/1000.;
         if (NewE > STC_PADDLE_THRESHOLD) {
            hddm_s::StcHitList hits = iter->addStcHits();
            hits().setT(t);
            hits().setDE(NewE);
         }
      }

      if (DROP_TRUTH_HITS)
         iter->deleteStcTruthHits();
   }
}

//-----------
// SmearCherenkov
//-----------
void SmearCherenkov(hddm_s::HDDM *record)
{
}

//-----------
// SmearTAGM
//-----------
void SmearTAGM(hddm_s::HDDM *record)
{
   hddm_s::MicroChannelList tagms = record->getMicroChannels();
   hddm_s::MicroChannelList::iterator iter;
   for (iter = tagms.begin(); iter != tagms.end(); ++iter) {
      iter->deleteTaggerHits();
      hddm_s::TaggerTruthHitList thits = iter->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time
         double t = titer->getT() + SampleGaussian(TAGM_TSIGMA);
         double tADC = titer->getT() + SampleGaussian(TAGM_FADC_TSIGMA);
         double npe = SamplePoisson(titer->getDE() * TAGM_NPIX_PER_GEV);
         hddm_s::TaggerHitList hits = iter->addTaggerHits();
         hits().setT(t);
         hits().setTADC(tADC);
         hits().setNpe(npe);
      }

      if (DROP_TRUTH_HITS)
         iter->deleteTaggerTruthHits();
   }
}

//-----------
// SmearTAGH
//-----------
void SmearTAGH(hddm_s::HDDM *record)
{
   hddm_s::HodoChannelList taghs = record->getHodoChannels();
   hddm_s::HodoChannelList::iterator iter;
   for (iter = taghs.begin(); iter != taghs.end(); ++iter) {
      iter->deleteTaggerHits();
      hddm_s::TaggerTruthHitList thits = iter->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // smear the time
         double t = titer->getT() + SampleGaussian(TAGH_TSIGMA);
         double tADC = titer->getT() + SampleGaussian(TAGH_FADC_TSIGMA);
         double npe = SamplePoisson(titer->getDE() * TAGH_NPE_PER_GEV);
         hddm_s::TaggerHitList hits = iter->addTaggerHits();
         hits().setT(t);
         hits().setTADC(tADC);
         hits().setNpe(npe);
      }

      if (DROP_TRUTH_HITS)
         iter->deleteTaggerTruthHits();
   }
}

//-----------
// InitCDCGeometry
//-----------
void InitCDCGeometry(void)
{
  CDC_GEOMETRY_INITIALIZED = true;
  
  CDC_MAX_RINGS = 28;
  
  //-- This was cut and pasted from DCDCTrackHit_factory.cc on 10/11/2007 --
  
  float degrees0 = 0.0;
  float degrees6 = 6.0*M_PI/180.0;
  
  for(int ring=1; ring<=CDC_MAX_RINGS; ring++){
    int myNstraws=0;
    float radius = 0.0;
    float stereo=0.0;
    float phi_shift=0.0;
    float deltaX=0.0, deltaY=0.0;
    float rotX=0.0, rotY=0.0;
    switch(ring){
      // axial
    case  1:   myNstraws=  42;   radius= 10.7219;   stereo=  degrees0; phi_shift= 0.00000;   break;
    case  2:   myNstraws=  42;   radius= 12.097;           stereo=  degrees0; phi_shift= 4.285714;   break;
    case  3:   myNstraws=  54;   radius= 13.7803;   stereo=  degrees0; phi_shift= 2.00000;   break;
    case  4:   myNstraws=  54;   radius= 15.1621;   stereo=  degrees0; phi_shift= 5.3333333;   break;
      
      // -stereo
    case  5:   myNstraws=  66;   radius= 16.9321;   stereo= -degrees6; phi_shift= 0.33333;   break;
    case  6:   myNstraws=  66;   phi_shift= 0.33333;   deltaX= 18.2948 ;   deltaY= 0.871486;   rotX=-6.47674;   rotY=-0.302853;   break;
    case  7:   myNstraws=  80;   radius= 20.5213;   stereo= -degrees6; phi_shift= -0.5000;   break;
    case  8:   myNstraws=  80;   phi_shift= -0.5000;   deltaX= 21.8912;   deltaY= 0.860106;   rotX=-6.39548;   rotY=-0.245615;   break;
      
      // +stereo
    case  9:   myNstraws=  93;   radius= 23.8544;   stereo= +degrees6; phi_shift= 1.1000;   break;
    case 10:   myNstraws=  93;   phi_shift= 1.1000;   deltaX= 25.229;   deltaY= 0.852573;   rotX=+6.34142;   rotY=+0.208647;   break;
    case 11:   myNstraws= 106;   radius= 27.1877;   stereo= +degrees6; phi_shift= -1.40;   break;
    case 12:   myNstraws= 106;   phi_shift= -1.400;   deltaX= 28.5658;   deltaY= 0.846871;   rotX=+6.30035;   rotY=+0.181146;   break;
      
      // axial
    case 13:   myNstraws= 123;   radius= 31.3799;   stereo=  degrees0; phi_shift= 0.5000000;   break;
    case 14:   myNstraws= 123;   radius= 32.7747;   stereo=  degrees0; phi_shift= 1.9634146;   break;
    case 15:   myNstraws= 135;   radius= 34.4343;   stereo=  degrees0; phi_shift= 1.0000000;   break;
    case 16:   myNstraws= 135;   radius= 35.8301;   stereo=  degrees0; phi_shift= 2.3333333;   break;
      
      // -stereo
    case 17:   myNstraws= 146;   radius= 37.4446;   stereo= -degrees6; phi_shift= 0.2;   break;
    case 18:   myNstraws= 146;   phi_shift= 0.2;   deltaX= 38.8295;   deltaY= 0.835653;   rotX=-6.21919 ;   rotY=-0.128247;   break;
    case 19:   myNstraws= 158;   radius= 40.5364;   stereo= -degrees6; phi_shift= 0.7;   break;
    case 20:   myNstraws= 158;   phi_shift= 0.7;   deltaX=41.9225 ;   deltaY= 0.833676;   rotX=-6.20274;   rotY=-0.118271;   break;
      
      // +stereo
    case 21:   myNstraws= 170;   radius= 43.6152;   stereo= +degrees6; phi_shift= 1.1000;   break;
    case 22:   myNstraws= 170;   phi_shift= 1.1000;   deltaX=45.0025  ;   deltaY= 0.831738;   rotX=+6.18859;   rotY=+0.109325;   break;
    case 23:   myNstraws= 182;   radius= 46.6849;   stereo= +degrees6; phi_shift= 1.40;   break;
    case 24:   myNstraws= 182;   phi_shift= 1.400;   deltaX= 48.0733;   deltaY= 0.829899;   rotX=+6.1763;   rotY=+0.101315;   break;
      
      // axial
    case 25:   myNstraws= 197;   radius= 50.37;   stereo=  degrees0; phi_shift= 0.200000000;   break;
    case 26:   myNstraws= 197;   radius= 51.77;   stereo=  degrees0; phi_shift= 1.113705000;   break;
    case 27:   myNstraws= 209;   radius= 53.363;   stereo=  degrees0; phi_shift= 0.800000000;   break;
    case 28:   myNstraws= 209;   radius= 54.76;   stereo=  degrees0; phi_shift= 1.661244;   break;
    default:
      cerr<<__FILE__<<":"<<__LINE__<<" Invalid value for CDC ring ("<<ring<<") should be 1-28 inclusive!"<<endl;
    }
    NCDC_STRAWS.push_back(myNstraws);
    CDC_RING_RADIUS.push_back(radius);
  }
  
  double Nstraws = 0;
  double alpha = 0.0;
  for(unsigned int i=0; i<NCDC_STRAWS.size(); i++){
    Nstraws += (double)NCDC_STRAWS[i];
    alpha += (double)NCDC_STRAWS[i]/CDC_RING_RADIUS[i];
  }
}


//-----------
// InitFDCGeometry
//-----------
void InitFDCGeometry(void)
{
  FDC_GEOMETRY_INITIALIZED = true;
  
  int FDC_NUM_LAYERS = 24;
  //int WIRES_PER_PLANE = 96;
  //int WIRE_SPACING = 1.116;
  
  for (int layer=1; layer <= FDC_NUM_LAYERS; layer++) {
    
    float degrees00 = 0.0;
    float degrees60 = M_PI*60.0/180.0;
    
    float angle = 0.0;
    float z_anode = 212.0 + 95.5;
    switch (layer) {
     case  1:
       z_anode += -92.5 - 2.0;
       angle =  degrees00;
       break;
     case  2:
       z_anode += -92.5 + 0.0;
       angle = +degrees60;
       break;
     case  3:
       z_anode += -92.5 + 2.0;
       angle = -degrees60;
       break;
     case  4:
       z_anode += -86.5 - 2.0;
       angle =  degrees00;
       break;
     case  5:
       z_anode += -86.5 + 0.0;
       angle = +degrees60;
       break;
     case  6:
       z_anode += -86.5 + 2.0;
       angle = -degrees60;
       break;
     case  7:
       z_anode += -32.5 - 2.0;
       angle =  degrees00;
       break;
     case  8:
       z_anode += -32.5 + 0.0;
       angle = +degrees60;
       break;
     case  9:
       z_anode += -32.5 + 2.0;
       angle = -degrees60;
       break;
     case 10:
       z_anode += -26.5 - 2.0;
       angle =  degrees00;
       break;
     case 11:
       z_anode += -26.5 + 0.0;
       angle = +degrees60;
       break;
     case 12:
       z_anode += -26.5 + 2.0;
       angle = -degrees60;
       break;
     case 13:
       z_anode += +26.5 - 2.0;
       angle =  degrees00;
       break;
     case 14:
       z_anode += +26.5 + 0.0;
       angle = +degrees60;
       break;
     case 15:
       z_anode += +26.5 + 2.0;
       angle = -degrees60;
       break;
     case 16:
       z_anode += +32.5 - 2.0;
       angle =  degrees00;
       break;
     case 17:
       z_anode += +32.5 + 0.0;
       angle = +degrees60;
       break;
     case 18:
       z_anode += +32.5 + 2.0;
       angle = -degrees60;
       break;
     case 19:
       z_anode += +86.5 - 2.0;
       angle =  degrees00;
       break;
     case 20:
       z_anode += +86.5 + 0.0;
       angle = +degrees60;
       break;
     case 21:
       z_anode += +86.5 + 2.0;
       angle = -degrees60;
       break;
     case 22:
       z_anode += +92.5 - 2.0;
       angle =  degrees00;
       break;
     case 23:
       z_anode += +92.5 + 0.0;
       angle = +degrees60;
       break;
     case 24:
       z_anode += +92.5 + 2.0;
       angle = -degrees60;
       break;
    }
    FDC_LAYER_Z.push_back(z_anode);
  }
  
  // Coefficient used to calculate FDCsingle wire rate. We calculate
  // it once here just to save calculating it for every wire in every event
  FDC_RATE_COEFFICIENT = exp(-log(4.0)/23.0) / 2.0 /
                         log(24.0) * FDC_TIME_WINDOW/1000.0e-9;
  
  // Something is a little off in my calculation above so I scale it down via
  // an emprical factor:
  FDC_RATE_COEFFICIENT *= 0.353;
}

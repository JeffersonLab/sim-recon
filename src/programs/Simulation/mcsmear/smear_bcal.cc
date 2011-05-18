// $Id: smear.cc 7650 2011-03-29 22:52:30Z shepherd $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
using namespace std;

#include <BCAL/DBCALGeometry.h>

#include <math.h>
#include "units.h"
#include "HDDM/hddm_s.h"
#include <TF1.h>
#include <TH2.h>
#include <TRandom3.h>

#include "DBCALReadoutChannel.h"

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif


void bcalInit(DBCALGeometry &bcalGeom);
float bcalSamplingSmear( float E );
float bcalTimeSmear(float t, float e);
int getDarkHits();

// The following are set in bcalInit below
bool BCAL_INITIALIZED = false;
double BCAL_mevPerPE=0.0;
float BCAL_UNATTENUATE_TO_CENTER = 0.0;
int BCAL_Nsum_inner = 0;
int BCAL_Nsum_outer = 0;
double BCAL_sigma_singlePE_sq = 0.0;
double BCAL_sigma_ped_sq = 0.0;

map<int, DBCALReadoutChannel> bcal_fADCs; // key is DBCALGeometry::fADCId()

double SampleGaussian(double sigma);
int SamplePoisson(float lambda);
double SampleRange(double x1, double x2);


// Flag used specifically for BCAL
extern bool SMEAR_BCAL;

extern double BCAL_TIME_WINDOW;

// setup response parameters
extern float BCAL_DARKRATE_GHZ        ;// 0.041;
extern float BCAL_XTALK_FRACT         ;//0.03;
extern float BCAL_INTWINDOW_NS        ;//100;
extern float BCAL_DEVICEPDE           ;//0.12;
extern float BCAL_SAMPLING_FRACT      ;//0.15;
extern float BCAL_MAXOCCUPANCY_FRACT  ;//0.05;

// GX-doc 1069, Table 1 -- try to extract back to
// photons per side per MeV in fiber
// 4.6 / PDE / attentuation  (meaurements performed in center)
// 75 = 4.6  / 0.12 / exp( -200 / 300 )
extern float BCAL_PHOTONSPERSIDEPERMEV_INFIBER; // 75

// set the sampling smearing coefficients:
// (from GlueX-doc 827 v3 Figure 13 )
extern float BCAL_SAMPLINGCOEFA ; //= 0.042;
extern float BCAL_SAMPLINGCOEFB ; //= 0.013;

// time smearing comes from beam test time difference resolution with
// beam incident on the center of the module
extern float BCAL_TIMEDIFFCOEFA ; //= 0.07 * sqrt( 2 );

// no floor term, but leave the option here:
extern float BCAL_TIMEDIFFCOEFB ; //= 0.0 * sqrt( 2 );

// calculated later in based on number of photons and ph threshold
extern float Bcal_CellInnerThreshold ;

// The following are all false by default, but can be
// set to true via command line parameters. Setting
// one of these to true will turn off the feature.
extern bool NO_E_SMEAR;
extern bool NO_T_SMEAR;
extern bool NO_DARK_PULSES;
extern bool NO_THRESHOLD_CUT;


//-----------
// SmearBCAL
//-----------
void SmearBCAL(s_HDDM_t *hddm_s)
{
	DBCALGeometry bcalGeom;

	// Initialize BCAL globals on first call
	if(!BCAL_INITIALIZED)bcalInit(bcalGeom);
	
	// Container to keep track of the elements in bcal_fADCs that
	// have hits so we can sparsely clear them at the end of this
	set<int> fADC_ids_with_hits;
	
	// Loop over PhysicsEvents
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	for(unsigned int iphysics_event=0; iphysics_event<PE->mult; iphysics_event++){
	
		// Verify BCAL hits exist
		s_HitView_t *hits = PE->in[iphysics_event].hitView;
		if (hits == HDDM_NULL || hits->barrelEMcal == HDDM_NULL || hits->barrelEMcal->bcalCells == HDDM_NULL)continue;
    
		// Loop over GEANT hits in BCAL
		s_BcalCells_t *cells = hits->barrelEMcal->bcalCells;
		for(unsigned int j=0; j<cells->mult; j++){
			s_BcalCell_t *cell = &cells->in[j];

			// If data exists in HDDM tree for SiPM hits, then delete it
			// so we can write in new hits below.
			if(cell->bcalSiPMUpHits!=HDDM_NULL)free(cell->bcalSiPMUpHits);
			if(cell->bcalSiPMDownHits!=HDDM_NULL)free(cell->bcalSiPMDownHits);

			// We make one SiPM hit for every GEANT hit in both the up and
			// down stream ends
			unsigned int Nbcalhits = cell->bcalHits->mult;
			cell->bcalSiPMUpHits = make_s_BcalSiPMUpHits(Nbcalhits);
			cell->bcalSiPMDownHits = make_s_BcalSiPMDownHits(Nbcalhits);

			cell->bcalSiPMUpHits->mult = Nbcalhits;
			cell->bcalSiPMDownHits->mult = Nbcalhits;			

			// Loop over GEANT hits in BCAL for this cell
			for(unsigned int k=0; k<Nbcalhits; k++){
				s_BcalHit_t *bcalhit = &cell->bcalHits->in[k];
				s_BcalSiPMUpHit_t *bcaluphit = &cell->bcalSiPMUpHits->in[k];
				s_BcalSiPMDownHit_t *bcaldownhit = &cell->bcalSiPMDownHits->in[k];
	
				// Distance from location of energy deposition to each end of BCAL
				float upDist   = ( bcalGeom.BCALFIBERLENGTH / 2 ) + bcalhit->zLocal;
				float downDist = ( bcalGeom.BCALFIBERLENGTH / 2 ) - bcalhit->zLocal;
	
				// sampling fluctuations are correlated between ends
				float smearedE = bcalSamplingSmear( bcalhit->E );
				float upEnergy   = smearedE * exp( -upDist / bcalGeom.ATTEN_LENGTH );
				float downEnergy = smearedE * exp( -downDist / bcalGeom.ATTEN_LENGTH );

				// independently smear time for both ends -- time smearing 
				// parameters come from data taken with beam at the center of 
				// the module so there is an implicit exp( ( -L / 2 ) / lambda ) 
				// that needs to be canceled out since we are working
				// at this stage with attenuated energies
				float smearedtUp   = bcalTimeSmear( bcalhit->t,   upEnergy*BCAL_UNATTENUATE_TO_CENTER );
				float smearedtDown = bcalTimeSmear( bcalhit->t, downEnergy*BCAL_UNATTENUATE_TO_CENTER );

				// To get a realistic distribution we have to convert the smeared energy into 
				// an integer number of photoelectrons. We'll add the dark hits to that and
				// then calculate a sigma in energy based on the total number of pixels that
				// fired. The number of photoelectrons will be converted back into energy
				// and then smeared using this sigma.
				double upPE = floor(0.5 + upEnergy/(BCAL_mevPerPE * k_MeV)) + getDarkHits();
				double downPE = floor(0.5 + downEnergy/(BCAL_mevPerPE * k_MeV)) + getDarkHits();

				double sigma_up = sqrt(BCAL_sigma_ped_sq + upPE*BCAL_sigma_singlePE_sq);
				double sigma_down = sqrt(BCAL_sigma_ped_sq + downPE*BCAL_sigma_singlePE_sq);
				
				// Convert integer # of PE back into energy
				upEnergy = upPE * BCAL_mevPerPE * k_MeV;
				downEnergy = downPE * BCAL_mevPerPE * k_MeV;

				// Add in pedestal widths due to SiPMs
				upEnergy += SampleGaussian(sigma_up);
				downEnergy += SampleGaussian(sigma_down);
	
				// now offset times for propagation distance
				float upTime = smearedtUp + upDist / bcalGeom.C_EFFECTIVE;
				float downTime = smearedtDown + downDist / bcalGeom.C_EFFECTIVE;
	
				// If energy is smeared to negative or time is nan, set to 0.
				if(upEnergy <= 0 || !finite(upTime)){
					upEnergy = 0;
					upTime = 0;
				}
				if(downEnergy <= 0 || !finite(downTime)){
					downEnergy = 0;
					downTime = 0;
				}

				// Record upstream SiPM values in HDDM
				bcaluphit->E = upEnergy;
				bcaluphit->t = upTime;

				// Record downstream SiPM values in HDDM
				bcaldownhit->E = downEnergy;
				bcaldownhit->t = downTime;
			
				// Remember these hits in appropriate readout channel
				int fADC_Id = bcalGeom.fADCId( cell->module, cell->layer, cell->sector);
				bcal_fADCs[fADC_Id].uphits.push_back(bcaluphit);
				bcal_fADCs[fADC_Id].downhits.push_back(bcaldownhit);
				fADC_ids_with_hits.insert(fADC_Id);
			} //k (bcal hits)
		} //j (cells)

		// Above, we looped over GEANT hits to create SiPM hits. Below we loop
		// over fADC (readout) channels which may contain multiple SiPMs.
		// For each channel, we determine how many hits occurred in
		// BCAL_TIME_WINDOW, including cumulative dark hits. The number of 
		// cells over threshold are counted and the values recorded so that
		// we can create the correct number of bcalfADCCell stuctures in HDDM
		// later.

		// Loop over readout channels
		set<int> fADC_ids_over_thresh;
		map<int, DBCALReadoutChannel>::iterator iter;
		for(iter = bcal_fADCs.begin(); iter!=bcal_fADCs.end(); iter++){

			int fADC_id = iter->first;
			DBCALReadoutChannel &bcalfADC = iter->second;
			vector<s_BcalSiPMUpHit_t*> &uphits = bcalfADC.uphits;
			vector<s_BcalSiPMDownHit_t*> &downhits = bcalfADC.downhits;

			//--------- Upstream -------------
			// Initialize
			bcalfADC.Eup = 0.0;
			bcalfADC.tup = 0.0;

			// Sum upstream hits
			for(unsigned int j=0; j<uphits.size(); j++){
				bcalfADC.Eup += uphits[j]->E;
				bcalfADC.tup += uphits[j]->t * uphits[j]->E; // energy weighted time for sum
			}
			
			// Add in dark pulses for upstream SiPMs not hit, but that are included in sum
			unsigned int Ndark_channels_up = bcalfADC.NSiPM - uphits.size();
			if(uphits.size() > bcalfADC.NSiPM)Ndark_channels_up = 0;
			for(unsigned int j=0; j<Ndark_channels_up; j++){
				double darkPE = getDarkHits();
				double sigma_up = sqrt(BCAL_sigma_ped_sq + darkPE*BCAL_sigma_singlePE_sq);				
				double Edark = (darkPE * BCAL_mevPerPE * k_MeV) + SampleGaussian(sigma_up);
				double tdark = SampleRange( -0.25 * BCAL_INTWINDOW_NS, 0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
				bcalfADC.Eup += Edark;
				bcalfADC.tup += tdark * Edark;
			}
			
			// normalize time
			bcalfADC.tup /= bcalfADC.Eup;

			//--------- Downstream -------------
			// Initialize
			bcalfADC.Edown = 0.0;
			bcalfADC.tdown = 0.0;

			// Sum downstream hits
			for(unsigned int j=0; j<downhits.size(); j++){
				bcalfADC.Edown += downhits[j]->E;
				bcalfADC.tdown += downhits[j]->t * downhits[j]->E; // energy weighted time for sum
			}
			
			// Add in dark pulses for downstream SiPMs not hit, but that are included in sum
			unsigned int Ndark_channels_down = bcalfADC.NSiPM - downhits.size();
			if(downhits.size() > bcalfADC.NSiPM)Ndark_channels_down = 0;
			for(unsigned int j=0; j<Ndark_channels_down; j++){
				double darkPE = getDarkHits();
				double sigma_down = sqrt(BCAL_sigma_ped_sq + darkPE*BCAL_sigma_singlePE_sq);				
				double Edark = (darkPE * BCAL_mevPerPE * k_MeV) + SampleGaussian(sigma_down);
				double tdark = SampleRange( -0.25 * BCAL_INTWINDOW_NS, 0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
				bcalfADC.Edown += Edark;
				bcalfADC.tdown += tdark * Edark;
			}
			
			// normalize time
			bcalfADC.tdown /= bcalfADC.Edown;


			// If either the upstream or downstream hit is over threshold, then
			// a bcalfADCCell will need to be created so remember this id if needed.
			if((bcalfADC.Eup > bcalfADC.threshold) || (bcalfADC.Edown > bcalfADC.threshold)){
				fADC_ids_over_thresh.insert(fADC_id);
			}
		} // readout channel
		
		
		// At this point we have summed all SiPMs with dark hits included. Any cell with
		// either the upstream or downstream fADC over threshold has its id stored in
		// the fADC_ids_over_thresh container. We now need to create bcalfADCCell,
		// bcalfADCUpHit, and bcalfADCDownHit structures within the HDDM tree to hold
		// the fADC values.
		
		
		// Delete any existing bcalfADCCell structures (along with their bcalfADCUpHit
		// and bcalfADCDownHit structures).
		if(hits->barrelEMcal->bcalfADCCells!=HDDM_NULL){
			for(unsigned int i=0; i<hits->barrelEMcal->bcalfADCCells->mult; i++){
				if(hits->barrelEMcal->bcalfADCCells->in[i].bcalfADCUpHits)
					free(hits->barrelEMcal->bcalfADCCells->in[i].bcalfADCUpHits);
				if(hits->barrelEMcal->bcalfADCCells->in[i].bcalfADCDownHits)
					free(hits->barrelEMcal->bcalfADCCells->in[i].bcalfADCDownHits);
			}
			free(hits->barrelEMcal->bcalfADCCells);
		}
		
		// Make sure bcalfADCCells pointer is empty in case we return early below
		hits->barrelEMcal->bcalfADCCells = (s_BcalfADCCells_t*)HDDM_NULL;

		// If we have no cells over threshold, then bail now.
		if(fADC_ids_over_thresh.size()==0) continue; // next iphysics_event
		
		// Create bcalfADCCell structures to hold all of our hits
		hits->barrelEMcal->bcalfADCCells = make_s_BcalfADCCells(fADC_ids_over_thresh.size());
		unsigned int &mult = hits->barrelEMcal->bcalfADCCells->mult;
		set<int>::iterator it = fADC_ids_over_thresh.begin();
		for(mult=0; mult<fADC_ids_over_thresh.size(); mult++, it++){
		
			// Get pointer to this fADC cell in the HDDM tree
			s_BcalfADCCell_t *fADCcell = &hits->barrelEMcal->bcalfADCCells->in[mult];
			
			// Get pointer to our fADC cell information that needs to be copied to HDDM
			DBCALReadoutChannel &bcalfADC = bcal_fADCs[*it];
			
			fADCcell->module = bcalfADC.module;
			fADCcell->layer  = bcalfADC.layer;
			fADCcell->sector = bcalfADC.sector;

			// Upstream hit
			if(bcalfADC.Eup > bcalfADC.threshold){
				fADCcell->bcalfADCUpHits = make_s_BcalfADCUpHits(1);
				fADCcell->bcalfADCUpHits->mult = 1;
				s_BcalfADCUpHit_t *fadcuphit = &fADCcell->bcalfADCUpHits->in[0];
				
				fadcuphit->E = bcalfADC.Eup;
				fadcuphit->t = bcalfADC.tup;
			}else{
				fADCcell->bcalfADCUpHits = (s_BcalfADCUpHits_t*)HDDM_NULL;
			}

			// Downstream hit
			if(bcalfADC.Edown > bcalfADC.threshold){
				fADCcell->bcalfADCDownHits = make_s_BcalfADCDownHits(1);
				fADCcell->bcalfADCDownHits->mult = 1;
				s_BcalfADCDownHit_t *fadcdownhit = &fADCcell->bcalfADCDownHits->in[0];
				
				fadcdownhit->E = bcalfADC.Edown;
				fadcdownhit->t = bcalfADC.tdown;
			}else{
				fADCcell->bcalfADCDownHits = (s_BcalfADCDownHits_t*)HDDM_NULL;
			}
		} // fADC_ids_over_thresh
	} // iphysics_event

	
	// Clear the bcal_fADCs elements that have hits.
	// We keep track above so we can do a sparse clear here to
	// save time clearing every fADC array for every event.
	set<int>::iterator iter=fADC_ids_with_hits.begin();
	for(; iter!=fADC_ids_with_hits.end(); iter++){
		bcal_fADCs[*iter].Clear();
	}
}


void bcalInit(DBCALGeometry &bcalGeom)
{

// this the average number of (assumed single PE) 
  // dark pulses in a window

  double nAvg = BCAL_DARKRATE_GHZ * BCAL_INTWINDOW_NS;

  // now we need to find n such that 
  // sum_i=1^n P(n) > ( 1 - maxOccupancy )
  // P(n) is the probability to have n pulses
  //
  // P(n) is really a convolution since we have:
  // P(n) = P( n_d + n_x ) = 
  //    P( n_d; nAvg ) * P( n_x; n_d * x )
  // 
  // n_d number of dark pulses
  // x is the cross talk rate

  // numerically build int_0^n P(n)

  // in practice the cutoff is going to be < 100 PE
  // we can build an accurate pdf up to that
  double pdf[100];
  for( int i = 0; i < 100; ++i ){ pdf[i] = 0; }

  // probability for zero is easy:
  double darkTerm = exp( -nAvg );
  pdf[0] = darkTerm;

  for( int n_d = 1; n_d < 100; ++n_d ){

    darkTerm *= ( nAvg / n_d );

    double xTalkAvg = n_d * BCAL_XTALK_FRACT;

    // probability for zero x-talk pulses
    double xTerm = exp( -xTalkAvg );
    pdf[n_d] += ( xTerm * darkTerm );

    // now include probability for additional
    // cross talk pulses
    for( int n_x = 1; n_x + n_d < 100; ++n_x ){

      xTerm *= ( xTalkAvg / n_x );

      pdf[n_d+n_x] += ( xTerm * darkTerm );
    }
  }

  double integral = 0;
  int nPEThresh = 0;
  while( integral < ( 1 - BCAL_MAXOCCUPANCY_FRACT ) ){

    // post increment includes zero and requires
    // one more PE than what breaks the loop
    integral += pdf[nPEThresh];
    ++nPEThresh;
  }

  // now get the photon theshold
  double photonThresh = nPEThresh / BCAL_DEVICEPDE;

  // now convert this into a energy threshold
  Bcal_CellInnerThreshold = 
    ( photonThresh / BCAL_PHOTONSPERSIDEPERMEV_INFIBER ) / 
    BCAL_SAMPLING_FRACT * k_MeV; 

	// When comparing the new to old verison of the BCAL smearing
	// code it was noticed that a number of hits right at threshold
	// were being accepted by the new and rejected by the old. It
	// is believed this is due to round-off errors arising from the
	// former using doubles in several places. That coupled with the
	// threshold being determined by the number of photo-electrons in
	// the dark hits. Therefore, for channels with only dark hits and
	// that number happening to be exactly the same as the threshold,
	// the round-off determined whether the hit was accepted or not.
	// This turned out to be a significant number. Thus,we add 10keV
	// to the threshold to nudge it out of range of the round-off,
	// replicating the result from the old code.
	Bcal_CellInnerThreshold += 0.00001;

	if(NO_THRESHOLD_CUT){
		Bcal_CellInnerThreshold = 0.0; // n.b. the outer threshold is derived from this too
	}

	// MeV per PhotoElectron in SiPM
	BCAL_mevPerPE = 1 /
		( BCAL_PHOTONSPERSIDEPERMEV_INFIBER * BCAL_DEVICEPDE * BCAL_SAMPLING_FRACT );
	
	// Factor to unattenuate the energy at the end to what it *would* be if it came
	// from the center of the module. This is so the proper time smearing can be
	// applied based on the 2006 beam test results which shot photons into the center
	// of the module.
	BCAL_UNATTENUATE_TO_CENTER = exp(bcalGeom.BCALFIBERLENGTH/2/ bcalGeom.ATTEN_LENGTH);
	
	// Number SiPMs added to give one fADC signal in inner and outer regions
	BCAL_Nsum_inner = DBCALGeometry::NSUMLAYS1 * DBCALGeometry::NSUMSECS1;
	BCAL_Nsum_outer = DBCALGeometry::NSUMLAYS2 * DBCALGeometry::NSUMSECS2;

	// The values of 0.108 and 0.143 came from the plot on slide 10 of GlueX-doc-1754-v3
	// 6.489/602 = 0.108
	// 8.585/602 = 0.143
	double sigma_singlePE = 0.108*BCAL_mevPerPE * k_MeV;
	double sigma_ped = sqrt(16.0)*0.143*BCAL_mevPerPE * k_MeV; // 0.143 is for single tile. sqrt(16) scales up to 16 tile SiPM array
	BCAL_sigma_singlePE_sq = sigma_singlePE*sigma_singlePE;
	BCAL_sigma_ped_sq = sigma_ped*sigma_ped;

	//
	// NOTE!!! This scheme needs to be changed if we want to get any
	// advantage from multi-threading! Otherwise, a mutex lock needs
	// to be added that encloses the entire BCAL smearing code.
	//
	// Create DBCALReadoutChannel objects for all readout channels.
	// It's easiest to do this for every cell regardless if summing is
	// turned on. Since the values are stored in a map, only one
	// DBCALReadoutChannel object will exist for every readout channel
	// in the end. That leaves the map such that it is fast to access
	// during event time since it is indexed by DBCALGeometry::fADCId(). 
	for(int i = 1; i<=48;i++){

		// inner
		unsigned int NSiPM_up = DBCALGeometry::NSUMLAYS1 * DBCALGeometry::NSUMSECS1;
		for(int j = 1; j<=6;j++){
			for(int k = 1; k<=4; k++){
				int cellId = DBCALGeometry::cellId( i, j, k);
				int fADCId = DBCALGeometry::fADCId( i, j, k);
				double thresh = Bcal_CellInnerThreshold * (double)BCAL_Nsum_inner;
				bcal_fADCs[fADCId] = DBCALReadoutChannel(NSiPM_up, thresh, i, DBCALGeometry::fADC_layer(cellId), DBCALGeometry::fADC_sector(cellId));
			}
		}

		// outer
		unsigned int NSiPM_down = DBCALGeometry::NSUMLAYS2 * DBCALGeometry::NSUMSECS2;
		for(int j = 1; j<=4;j++){
			for(int k = 1; k<=4; k++){
				int cellId = DBCALGeometry::cellId( i, j+6, k);
				int fADCId = DBCALGeometry::fADCId( i, j+6, k);
				double thresh = Bcal_CellInnerThreshold * (double)BCAL_Nsum_outer;
				bcal_fADCs[fADCId] = DBCALReadoutChannel(NSiPM_down, thresh, i, DBCALGeometry::fADC_layer(cellId), DBCALGeometry::fADC_sector(cellId));
			}
		}
	}

	// Flag that we have been called
	BCAL_INITIALIZED = true;
}


//-----------
// bcalSamplingSmear
//-----------
float bcalSamplingSmear( float E )
{
	if(NO_E_SMEAR)return E;

    double sigmaSamp = BCAL_SAMPLINGCOEFA / sqrt( E ) + BCAL_SAMPLINGCOEFB;
    
    return( E * (1.0 + SampleGaussian(sigmaSamp)) );
}

//-----------
// bcalTimeSmear
//-----------
float bcalTimeSmear( float t, float E )
{
	if(NO_T_SMEAR)return t;

  double sigmaT = BCAL_TIMEDIFFCOEFA / sqrt( E ) + BCAL_TIMEDIFFCOEFB;

  return( t + SampleGaussian(sigmaT) );
}

//-----------
// getDarkHits
//-----------
int getDarkHits()
{
	if(NO_DARK_PULSES)return 0;

  int darkPulse = SamplePoisson( BCAL_DARKRATE_GHZ* BCAL_INTWINDOW_NS );

  int xTalk = SamplePoisson( darkPulse * BCAL_XTALK_FRACT );

  return( xTalk + darkPulse );
}


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
#include <TMath.h>
#include <TMatrix.h>

#include "DBCALReadoutChannel.h"

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

// From SampleGaussian.cc
double SampleGaussian(double sigma);
double SamplePoisson(double lambda);
double SampleRange(double x1, double x2);

// Defined in this file
double bcalSamplingSmear( double E );
double bcalTimeSmear(double t, double e);
int getDarkHits();
void bcalInit(void);
void CalculateThresholds(void);
TH1D* CombineSiPM_PDFs(const char *hname, TH1D *single, int NSiPM);


// The following are set in bcalInit below
bool BCAL_INITIALIZED = false;
double BCAL_mevPerPE=0.0;
double BCAL_UNATTENUATE_TO_CENTER = 0.0;
int BCAL_Nsum_inner = 0;
int BCAL_Nsum_outer = 0;
double BCAL_sigma_singlePE_sq = 0.0;
double BCAL_sigma_ped_sq = 0.0;
double BCAL_inner_thresh = 0.0;
double BCAL_outer_thresh = 0.0;
double BCAL_inner_mean = 0.0; // Mean of dark hit distribution in GeV
double BCAL_outer_mean = 0.0; // Mean of dark hit distribution in GeV
int Nsum_inner;
int Nsum_outer;
map<int, DBCALReadoutChannel> bcal_fADCs; // key is DBCALGeometry::fADCId()



// Flag used specifically for BCAL
extern bool SMEAR_BCAL;

// The following are all false by default, but can be
// set to true via command line parameters. Setting
// one of these to true will turn OFF the feature.
extern bool NO_E_SMEAR;
extern bool NO_T_SMEAR;
extern bool NO_DARK_PULSES;
extern bool NO_THRESHOLD_CUT;


// setup response parameters
extern double BCAL_DARKRATE_GHZ;                // 0.0176 (from calibDB BCAL/bcal_parms) for 4x4 array
extern double BCAL_SIGMA_SIG_RELATIVE;          // 0.105  (from calibDB BCAL/bcal_parms)
extern double BCAL_SIGMA_PED_RELATIVE;          // 0.139  (from calibDB BCAL/bcal_parms)
extern double BCAL_SIPM_GAIN_VARIATION;         // 0.04   (from calibDB BCAL/bcal_parms)
extern double BCAL_XTALK_FRACT;                 // 0.157  (from calibDB BCAL/bcal_parms)
extern double BCAL_INTWINDOW_NS;                // 100    (from calibDB BCAL/bcal_parms)
extern double BCAL_DEVICEPDE;                   // 0.21   (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLING_FRACT;              // 0.095  (from calibDB BCAL/bcal_parms)
extern double BCAL_PHOTONSPERSIDEPERMEV_INFIBER;// 75 (from calibDB BCAL/bcal_parms) 
extern double BCAL_AVG_DARK_DIGI_VALS_PER_EVENT;// 240 used to set thresholds
extern double BCAL_SAMPLINGCOEFA;               // 0.042 (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLINGCOEFB;               // 0.013 (from calibDB BCAL/bcal_parms)
extern double BCAL_TIMEDIFFCOEFA;               // 0.07 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
extern double BCAL_TIMEDIFFCOEFB;               // 0.00 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)



//-----------
// SmearBCAL
//-----------
void SmearBCAL(s_HDDM_t *hddm_s)
{
	DBCALGeometry bcalGeom;

	// Initialize BCAL globals on first call
	if(!BCAL_INITIALIZED)bcalInit();
	
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
				double upDist   = ( bcalGeom.BCALFIBERLENGTH / 2 ) + bcalhit->zLocal;
				double downDist = ( bcalGeom.BCALFIBERLENGTH / 2 ) - bcalhit->zLocal;
	
				// sampling fluctuations are correlated between ends
				double smearedE = bcalSamplingSmear( bcalhit->E );
				double upEnergy   = smearedE * exp( -upDist / bcalGeom.ATTEN_LENGTH );
				double downEnergy = smearedE * exp( -downDist / bcalGeom.ATTEN_LENGTH );

				// independently smear time for both ends -- time smearing 
				// parameters come from data taken with beam at the center of 
				// the module so there is an implicit exp( ( -L / 2 ) / lambda ) 
				// that needs to be canceled out since we are working
				// at this stage with attenuated energies
				double smearedtUp   = bcalTimeSmear( bcalhit->t,   upEnergy*BCAL_UNATTENUATE_TO_CENTER );
				double smearedtDown = bcalTimeSmear( bcalhit->t, downEnergy*BCAL_UNATTENUATE_TO_CENTER );

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
				double upTime = smearedtUp + upDist / bcalGeom.C_EFFECTIVE;
				double downTime = smearedtDown + downDist / bcalGeom.C_EFFECTIVE;
	
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

//-----------
// bcalSamplingSmear
//-----------
double bcalSamplingSmear( double E )
{
	if(NO_E_SMEAR)return E;

    double sigmaSamp = BCAL_SAMPLINGCOEFA / sqrt( E ) + BCAL_SAMPLINGCOEFB;
    
    return( E * (1.0 + SampleGaussian(sigmaSamp)) );
}

//-----------
// bcalTimeSmear
//-----------
double bcalTimeSmear( double t, double E )
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

  int darkPulse = (int)floor(0.5+SamplePoisson( BCAL_DARKRATE_GHZ* BCAL_INTWINDOW_NS ));

  int xTalk = (int)floor(0.5+SamplePoisson( (double)darkPulse * BCAL_XTALK_FRACT ));

  return( xTalk + darkPulse );
}

//-----------
// bcalInit
//-----------
void bcalInit(void)
{
	// Calculate the inner and outer thresholds based on
	// the parameters from the calibDB. The threshold values
	// are left in the global variables BCAL_inner_thresh
	// and BCAL_outer_thresh.
	CalculateThresholds();

	if(NO_THRESHOLD_CUT){
		BCAL_inner_thresh = 0.0;
		BCAL_outer_thresh = 0.0;
	}

	cout<<"BCAL inner threshold: "<<BCAL_inner_thresh*1000.0<<" MeV"<<" ("<<Nsum_inner<<" SiPMs summed)"<<endl;
	cout<<"BCAL outer threshold: "<<BCAL_outer_thresh*1000.0<<" MeV"<<" ("<<Nsum_outer<<" SiPMs summed)"<<endl;
	
	// Factor to unattenuate the energy at the end to what it *would* be if it came
	// from the center of the module. This is so the proper time smearing can be
	// applied based on the 2006 beam test results which shot photons into the center
	// of the module.
	BCAL_UNATTENUATE_TO_CENTER = exp((DBCALGeometry::BCALFIBERLENGTH/2.0)/DBCALGeometry::ATTEN_LENGTH);
	
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
		for(int j = 1; j<=6;j++){
			for(int k = 1; k<=4; k++){
				int cellId = DBCALGeometry::cellId( i, j, k);
				int fADCId = DBCALGeometry::fADCId( i, j, k);
				double thresh = BCAL_inner_thresh;
				bcal_fADCs[fADCId] = DBCALReadoutChannel(Nsum_inner, thresh, i, DBCALGeometry::fADC_layer(cellId), DBCALGeometry::fADC_sector(cellId));
			}
		}

		// outer
		for(int j = 1; j<=4;j++){
			for(int k = 1; k<=4; k++){
				int cellId = DBCALGeometry::cellId( i, j+6, k);
				int fADCId = DBCALGeometry::fADCId( i, j+6, k);
				double thresh = BCAL_outer_thresh;
				bcal_fADCs[fADCId] = DBCALReadoutChannel(Nsum_outer, thresh, i, DBCALGeometry::fADC_layer(cellId), DBCALGeometry::fADC_sector(cellId));
			}
		}
	}

	// Flag that we have been called
	BCAL_INITIALIZED = true;
}

//-----------
// CalculateThreshold
//-----------
void CalculateThresholds(void)
{
	/// Calculate the BCAL thresholds based on the parameters in 
	/// the calibDB in BCAL/bcal_parms.
	///
	/// To do this, the electronic signal distribution (in MeV) is
	/// calculated for a single SiPM. This is used to determine a 
	/// cumulative distribution for multiple SiPMs in the case of
	/// summing.
	///
	/// For the course segmentation scheme, 3 inner SiPMs are summed
	/// while 4 outer SiPMs are summed. Because of this, different
	/// thresholds are calculated for the inner and outer regions.
	///
	/// The average number of "hits" per event due to dark pulses
	/// is given by the calibration parameter BCAL_AVG_DARK_DIGI_VALS_PER_EVENT.
	/// This sets a limit on the BCAL contribution to the average
	/// event size due to dark pulses and the thresholds are set 
	/// based on that.
	///
	/// A choice must be made on how to distribute the allowed,
	/// average number of hits between the inner and outer regions
	/// because they sum different numbers of SiPMs. This means
	/// they have different contributions to the signal spectrum
	/// due to the dark hits. I.e. for a given energy, there will
	/// be more dark pulses when summing 4 SiPMs together than
	/// when summing 3.
	/// 
	/// The method implemented here is to set the thresholds such that
	/// the fraction of SimPMs that are over threshold in a given event
	/// is the same for the inner SiPMs as for the outer.
	
	// Number of readout channels total
	DBCALGeometry g;
	int Nchan_inner = 2*g.NBCALMODS * g.NBCALLAYS1 * g.NBCALSECS1;
	int Nchan_outer = 2*g.NBCALMODS * g.NBCALLAYS2 * g.NBCALSECS2;
	
	// The inner BCAL region will have both an fADC and a TDC with
	// the fADC supplying both an amplitude and a time. The outer
	// region will have only an fADC so only 2 values will come
	// from those hits.
	//
	// To find the fraction of channels above threshold for a given
	// BCAL_AVG_DARK_DIGI_VALS_PER_EVENT, we use the fact that the
	// fraction of inner channels firing is equal to the fraction
	// of outer channels firing and that the sum of digitized values
	// is given by BCAL_AVG_DARK_DIGI_VALS_PER_EVENT:
	//
	// Ninner/Nchan_inner = Nouter/Nchan_outer
	//
	// BCAL_AVG_DARK_DIGI_VALS_PER_EVENT = 3*Ninner + 2*Nouter
	//
	// Solving for Nouter...
	//
	// Nouter = BCAL_AVG_DARK_DIGI_VALS_PER_EVENT/(3*Nchan_inner/Nchan_outer + 2)
	//
	// The desired fraction of readout channels firing is then:
	//
	// fraction_above_threshold = Nouter/Nchan_outer
	double Nouter = BCAL_AVG_DARK_DIGI_VALS_PER_EVENT/(3.0*(double)Nchan_inner/(double)Nchan_outer + 2.0);
	double fraction_above_threshold = Nouter/(double)Nchan_outer; // also equals Ninner/Nchan_inner

	// Some quantities are for single SiPM tiles, and some for 4x4 arrays
	double Ntiles = 16.0;

	// MeV per PhotoElectron in SiPM
	BCAL_mevPerPE = 1.0/( BCAL_PHOTONSPERSIDEPERMEV_INFIBER * BCAL_DEVICEPDE * BCAL_SAMPLING_FRACT );

	// Single PE signal and pedestal widths (values in calibDB are for single tile)
	double sigma_singlePE = BCAL_SIGMA_SIG_RELATIVE * BCAL_mevPerPE * k_MeV;
	double sigma_ped = sqrt(Ntiles) * BCAL_SIGMA_PED_RELATIVE * BCAL_mevPerPE * k_MeV;
	double sigma_gain_variation = BCAL_SIPM_GAIN_VARIATION * BCAL_mevPerPE * k_MeV;
	
	BCAL_sigma_singlePE_sq = sigma_singlePE*sigma_singlePE + sigma_gain_variation*sigma_gain_variation; // gain variation contributes in same way as single PE width
	BCAL_sigma_ped_sq = sigma_ped*sigma_ped;

	// Mean number of dark pulses from one SiPM inside integration window
	double mean_dark_pulses_array = BCAL_DARKRATE_GHZ * BCAL_INTWINDOW_NS; // BCAL_DARKRATE_GHZ is for whole 4x4 array


	// The following is code originally written by Dan Bennet at IU (I think).
	// The gist of it is that it uses the recursion relation for Poisson
	// distributions to build up a probability distribution for the number of
	// dark pulse photoelectrons per event (i.e. BCAL_INTWINDOW_NS time window).
	// This includes the cross-talk contribution, but not the cross-talk due
	// to cross talk, etc...
	// The probability of getting N photoelectrons from a SiPM array in a 
	// random BCAL_INTWINDOW_NS time window will be given by pdf[N].

	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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
	double darkTerm = exp( -mean_dark_pulses_array );
	pdf[0] = darkTerm;

	for( int n_d = 1; n_d < 100; ++n_d ){

		darkTerm *= ( mean_dark_pulses_array / n_d );

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
	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><>

	// At this point we add in electronic noise using the sigma_singlePE
	// and sigma_ped calculated above. This is done by creating a histogram
	// and filling it with contributions from Npe=0 to Npe=99
	TH1D *bcal_dark_signal = new TH1D("bcal_dark_signal", "fADC signal", 1001, -0.05, 100.05);
	bcal_dark_signal->SetXTitle("fADC (MeV)");
	for(int bin=1; bin<=bcal_dark_signal->GetNbinsX(); bin++){

		double x = bcal_dark_signal->GetBinCenter(bin);
		
		
		// Loop over Npe
		double prob = 0.0;
		for(int Npe=0; Npe<100; Npe++){
			double sigma_tot = sqrt(BCAL_sigma_ped_sq + (double)Npe*BCAL_sigma_singlePE_sq)*1000.0; // in MeV
			double E = (double)Npe * BCAL_mevPerPE; // in MeV
			
			prob += TMath::Gaus(x, E , sigma_tot)*pdf[Npe];
		}

		bcal_dark_signal->SetBinContent(bin, prob);
	}

	// The histogram bcal_dark_signal now contains the probability 
	// distribution for the electronic response of a single SiPM array
	// in units of MeV. We wish to convolute that with itself
	// multiple times to get the PDF for multiple SiPMs added together
	// (in the case of summing).
	Nsum_inner = g.NSUMLAYS1 * g.NSUMSECS1;
	Nsum_outer = g.NSUMLAYS2 * g.NSUMSECS2;
	
	TH1D *bcal_dark_signal_sum_inner = CombineSiPM_PDFs("bcal_dark_signal_sum_inner", bcal_dark_signal, Nsum_inner);
	TH1D *bcal_dark_signal_sum_outer = CombineSiPM_PDFs("bcal_dark_signal_sum_outer", bcal_dark_signal, Nsum_outer);
	
	// Integrate the combined PDFs
	TH1D *ibcal_dark_signal_sum_inner = (TH1D*)bcal_dark_signal_sum_inner->Clone("ibcal_dark_signal_sum_inner");
	TH1D *ibcal_dark_signal_sum_outer = (TH1D*)bcal_dark_signal_sum_outer->Clone("ibcal_dark_signal_sum_outer");
	for(int bin=2; bin<=ibcal_dark_signal_sum_inner->GetNbinsX(); bin++){
		double sum1 = ibcal_dark_signal_sum_inner->GetBinContent(bin-1) + ibcal_dark_signal_sum_inner->GetBinContent(bin);
		ibcal_dark_signal_sum_inner->SetBinContent(bin, sum1);

		double sum2 = ibcal_dark_signal_sum_outer->GetBinContent(bin-1) + ibcal_dark_signal_sum_outer->GetBinContent(bin);
		ibcal_dark_signal_sum_outer->SetBinContent(bin, sum2);
	}

	
	// Normalize
	int last_bin = ibcal_dark_signal_sum_inner->GetNbinsX();
	ibcal_dark_signal_sum_inner->Scale(1.0/ibcal_dark_signal_sum_inner->GetBinContent(last_bin));
	ibcal_dark_signal_sum_outer->Scale(1.0/ibcal_dark_signal_sum_outer->GetBinContent(last_bin));

	// Finally, search the integrated, normalized histograms for the bin
	// where the integral fraction above it is less than fraction_above_threshold
	for(int bin=1 ; bin<=last_bin; bin++){
		if(ibcal_dark_signal_sum_inner->GetBinContent(bin) < (1.0 - fraction_above_threshold)){
			BCAL_inner_thresh = ibcal_dark_signal_sum_inner->GetBinCenter(bin) * k_MeV;
		}
		if(ibcal_dark_signal_sum_outer->GetBinContent(bin) < (1.0 - fraction_above_threshold)){
			BCAL_outer_thresh = ibcal_dark_signal_sum_inner->GetBinCenter(bin) * k_MeV;
		}
	}
}


//-----------
// CombineSiPM_PDFs
//-----------
TH1D* CombineSiPM_PDFs(const char *hname, TH1D *single, int NSiPM)
{
	// It is quicker (and easier) to do this in stages by calculating the
	// probability of signals adding to a given energy and then 
	// subsequently calculating the probability of adding one more.
	
	TH1D *summed = (TH1D*)single->Clone(hname);
	
	for(int i=2; i<=NSiPM; i++){
	
		// Make copy of summed histogram since we'll be overwritting it
		TH1D *tmp = (TH1D*)summed->Clone("tmp");
		
		// Reset summed histo (for safety)
		summed->Reset();
		
		// Loop over all combinations of bins in the previous iteration's
		// summed probability histo and the single SiPM probability
		// histo. Fill in the new summed probability histo.
		for(int ibin=1; ibin<=tmp->GetNbinsX(); ibin++){

			double NPE1 = tmp->GetBinCenter(ibin);
			double prob1 = tmp->GetBinContent(ibin);

			for(int jbin=1; jbin<=single->GetNbinsX(); jbin++){

				double NPE2 = single->GetBinCenter(jbin);
				double prob2 = single->GetBinContent(jbin);
				
				// Find the bin in the summed histogram corresponding to
				// the sum of the energies from the bin in the previous
				// iteration's summed prob. histo and the single SiPM
				// prob. histo. If the energy sum is in range, then set
				// the bin content of the current iteration's summed prob.
				// histo.
				int kbin = summed->FindBin(NPE1 + NPE2);
				if(kbin>=1 && kbin<=summed->GetNbinsX()){
					summed->Fill(NPE1 + NPE2, prob1*prob2);
				}
			}
		}
		
		delete tmp;
	}
	
	return summed;
}



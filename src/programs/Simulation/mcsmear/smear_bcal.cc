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


// Uncomment any of the following to turn off the specific feature
//#define NO_E_SMEAR
//#define NO_T_SMEAR
//#define NO_DARK_PULSES
//#define NO_THRESHOLD_CUT


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
double BCAL_mean_darkpulse_energy_SiPM = 0.0;
double BCAL_mean_darkpulse_energy_inner = 0.0;
double BCAL_mean_darkpulse_energy_outer = 0.0;

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


// Below are two versions of the BCAL smearing code. The first
// is Beni's version copied from smear.cc and the second is
// Dave's (incomplete) version. The following #define can be
// used to switch between the two with 1=Beni and 0=Dave.
//
// You do not want to use the Dave version. It is being included
// here since some non-trivial amount of work has been put into it
// and I want to back it up in the repository, but without
// making it the default yet. Once it is completed and tested,
// The Dave version will be made the default and the Beni version
// removed.
//
// The "Dave" version is based on the "Beni" version, but tries
// to address a few things:
//
// 1. Multiple hits in the same readout channel that are out 
//    of time are made into separate fADC hits
//
// 2. Multiple dark pulses whose sum exceeds threshold to create
//    purely dark noise hits
//
// 3. Speed optimization.
//
#if 0

//-----------
// SmearBCAL
//-----------
void SmearBCAL(s_HDDM_t *hddm_s)
{
  DBCALGeometry *bcalGeom = new DBCALGeometry();
  
  map< int, pair< int, int > > darkHits;
  map< int, pair< int, int > >::iterator darkHitItr;
  
  
  for( int m = 1; m <= 48; ++m ){
    for( int l = 1; l <= 10; ++l ){
      for( int s = 1; s <= 4; ++s ){
	
	pair< int, int > nHits( getDarkHits(), getDarkHits() );
	
	darkHits[DBCALGeometry::cellId( m, l, s )] = nHits;
      }
    }
  }
  
  double mevPerPE = 1 / 
    ( BCAL_PHOTONSPERSIDEPERMEV_INFIBER * BCAL_DEVICEPDE *
      BCAL_SAMPLING_FRACT );
  
  bcalInit(*bcalGeom); //find inner cell threshold E
  
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return;
  
  for(unsigned int i=0; i<PE->mult; i++){ 
    
    float eUpStore[49][11][5] = {{{0}}}; //store energies [mod][lay][sect]
    float eDownStore[49][11][5]= {{{0}}};	  
    
    float eUpfADC[49][11][5]={{{0}}};
    float eDownfADC[49][11][5]={{{0}}};
    
    float tUpStore[49][11][5] = {{{0}}}; //store times [mod][lay][sect]
    float tDownStore[49][11][5]= {{{0}}};
    
    float tUpfADC[49][11][5]={{{0}}};
    float tDownfADC[49][11][5]={{{0}}};
    
    int fADCCellCount = 0;
    
    s_HitView_t *hits = PE->in[i].hitView;
    if (hits == HDDM_NULL ||
	hits->barrelEMcal == HDDM_NULL ||
	hits->barrelEMcal->bcalCells == HDDM_NULL)continue;
    
    
    s_BcalCells_t *cells = hits->barrelEMcal->bcalCells;
    
    
    
    for(unsigned int j=0; j<cells->mult; j++){
      s_BcalCell_t *cell = &cells->in[j];
      
      // dcell is needed for dark hits cellID
      int dcell = DBCALGeometry::cellId( cell->module, cell->layer, cell->sector);
      darkHitItr = darkHits.find( dcell );
      
      
      //Create BCAL hits structure to put smeared data into
      
      if(cell->bcalSiPMUpHits!=HDDM_NULL)free(cell->bcalSiPMUpHits);
      cell->bcalSiPMUpHits = make_s_BcalSiPMUpHits(cell->bcalHits->mult);
      cell->bcalSiPMUpHits->mult = cell->bcalHits->mult;
      
      if(cell->bcalSiPMDownHits!=HDDM_NULL)free(cell->bcalSiPMDownHits);
      cell->bcalSiPMDownHits = make_s_BcalSiPMDownHits(cell->bcalHits->mult);
      cell->bcalSiPMDownHits->mult = cell->bcalHits->mult;			
      
      for(unsigned int k=0; k<cell->bcalHits->mult; k++){
	s_BcalHit_t *bcalhit = &cell->bcalHits->in[k];
	s_BcalSiPMUpHit_t *bcaluphit = &cell->bcalSiPMUpHits->in[k];
	s_BcalSiPMDownHit_t *bcaldownhit = &cell->bcalSiPMDownHits->in[k];
	
	float smearedE = bcalSamplingSmear( bcalhit->E );
	
	float upDist = ( bcalGeom->BCALFIBERLENGTH / 2 ) + bcalhit->zLocal;
	float downDist = ( bcalGeom->BCALFIBERLENGTH / 2 ) - bcalhit->zLocal;
	
	// sampling fluctuations are correlated between ends
	float upEnergy = smearedE * exp( -upDist / bcalGeom->ATTEN_LENGTH );
	float downEnergy = smearedE * exp( -downDist / bcalGeom->ATTEN_LENGTH );
	
	// independently smear time for both ends -- time smearing 
	// parameters come from data taken with beam at the center of 
	// the module so there is an implicit exp( ( -L / 2 ) / lambda ) 
	// that needs to be canceled out since we are working
	// at this stage with attenuated energies 
	float smearedtUp = 
	  bcalTimeSmear( bcalhit->t, 
			 upEnergy * exp( ( bcalGeom->BCALFIBERLENGTH / 2 ) / 
					 bcalGeom->ATTEN_LENGTH ) );
	float smearedtDown = 
	  bcalTimeSmear( bcalhit->t, 
			 downEnergy * exp( ( bcalGeom->BCALFIBERLENGTH / 2 ) / 
					   bcalGeom->ATTEN_LENGTH ) );
	
	darkHitItr = darkHits.find( dcell );
	if( darkHitItr != darkHits.end() ){
	  
	  upEnergy += ( darkHitItr->second.first * mevPerPE * k_MeV );
	  downEnergy += ( darkHitItr->second.second * mevPerPE * k_MeV );
	  
	  // now delete this from the map so we don't create
	  // additional hits later
	  darkHits.erase( darkHitItr );
	}
	
	// now offset times for propagation distance
	float upTime = smearedtUp + upDist / bcalGeom->C_EFFECTIVE;
	float downTime = smearedtDown + downDist / bcalGeom->C_EFFECTIVE;
	
	//If energy is smeared to negative, set to 0.
	if(upEnergy <= 0 || upTime!=upTime)
	  {
	    upEnergy = 0;
	    upTime = 0;
	  }
	if(downEnergy <= 0|| downTime!=downTime)
	  {
	    downEnergy = 0;
	    downTime = 0;
	  }
	
	bcaluphit->E = upEnergy;
	eUpStore[cell->module][cell->layer][cell->sector]+= upEnergy;
	bcaluphit->t = upTime;
	tUpStore[cell->module][cell->layer][cell->sector]= upTime;
	
	
	bcaldownhit->E = downEnergy;
	eDownStore[cell->module][cell->layer][cell->sector]+= downEnergy;
	bcaldownhit->t = downTime;
	tDownStore[cell->module][cell->layer][cell->sector]= downTime;
	
	
      } //k (bcal hits)
      
    } //j (cells)
    
    //Add in dark hits to empty cells
    for(int i = 1; i<=48;i++)
      {
	for(int j = 1; j<=10;j++)
	  {
	    for(int k = 1; k<=4; k++)
	      {
		int dcell = DBCALGeometry::cellId( i, j, k);
		darkHitItr = darkHits.find( dcell );
		if( darkHitItr != darkHits.end() ){
		  eUpStore[i][j][k] += ( darkHitItr->second.first * mevPerPE * k_MeV );
		  eDownStore[i][j][k] += ( darkHitItr->second.second * mevPerPE * k_MeV );
		  tUpStore[i][j][k] = SampleRange( -0.25 * BCAL_INTWINDOW_NS,
						   0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
		  tDownStore[i][j][k] = SampleRange( -0.25 * BCAL_INTWINDOW_NS,
						     0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
		}
	      }
	  }
      }		     
    
    //Inner Cells, Summing
    for(int i=1;i<=48;i++)
      {
	for(int j=1;j<=bcalGeom->NBCALLAYS1;j++)
	  {
	    for(int k=1;k<=bcalGeom->NBCALSECS1;k++)
	      {
		for(int l=1;l<=(bcalGeom->BCALMID-1)/bcalGeom->NBCALLAYS1;l++)
		  {
		    for(int m=1;m<=4/bcalGeom->NBCALSECS1;m++)
		      {			 
			eUpfADC[i][j][k]+= eUpStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m];//Sum energies
			eDownfADC[i][j][k]+= eDownStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m];
			
			tUpfADC[i][j][k]+= eUpStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m] //Sum times weighted by Energy
			  *tUpStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m];
			tDownfADC[i][j][k]+= eDownStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m]
			  *tDownStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m];
			
		      }
		  }
		
		if(eUpfADC[i][j][k]!=0) //Divide time by Energy, to average weighted by energy, but make sure nonzero energy present.
		  tUpfADC[i][j][k] = tUpfADC[i][j][k]/eUpfADC[i][j][k];
		else
		  tUpfADC[i][j][k] = 0;
		
		if(eDownfADC[i][j][k]!=0)
		  tDownfADC[i][j][k] = tDownfADC[i][j][k]/eDownfADC[i][j][k];
		else
		  tDownfADC[i][j][k] = 0;
	      }
	  }
      }
    //Outer Cells, Summing
    for(int i=1;i<=48;i++)
      {
	for(int j=bcalGeom->NBCALLAYS1+1;j<=bcalGeom->NBCALLAYS2+bcalGeom->NBCALLAYS1;j++)
	  {
	    for(int k=1;k<=bcalGeom->NBCALSECS2;k++)
	      {
		for(int l=1;l<=(11-bcalGeom->BCALMID)/bcalGeom->NBCALLAYS2;l++)
		  {
		    for(int m=1;m<=4/bcalGeom->NBCALSECS2;m++)
		      {			 
			eUpfADC[i][j][k]+= eUpStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m];//Sum energies
			eDownfADC[i][j][k]+= eDownStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m];
			
			tUpfADC[i][j][k]+= eUpStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m] //Sum times weighted by energy
			  *tUpStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m];
			tDownfADC[i][j][k]+= eDownStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m]
			  *tDownStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m];			 
		      }
		  }
		
		if(eUpfADC[i][j][k]!=0) //Divide time by Energy, to average weighted by energy, but make sure nonzero energy present.
		  tUpfADC[i][j][k] = tUpfADC[i][j][k]/eUpfADC[i][j][k];
		else
		  tUpfADC[i][j][k] = 0;
		
		if(eDownfADC[i][j][k]!=0)
		  tDownfADC[i][j][k] = tDownfADC[i][j][k]/eDownfADC[i][j][k];
		else
		  tDownfADC[i][j][k] = 0;
	      }
	  }
      }
    
    //Naively multiply threshold by the square root of the number summed cells
    //in order to apply an effective 95% cut on fADC dark counts
    //passing.  Needs more thought.	
    float InnerThreshold = Bcal_CellInnerThreshold * 
      ((bcalGeom->BCALMID-1)/bcalGeom->NBCALLAYS1) * (4/bcalGeom->NBCALSECS1);
    float OuterThreshold = Bcal_CellInnerThreshold * 
      ((11-bcalGeom->BCALMID)/bcalGeom->NBCALLAYS2) * (4/bcalGeom->NBCALSECS2);	
    
    for(int i = 1;i<=48;i++)
      {
	for(int j = 1;j<=bcalGeom->NBCALLAYS1 ;j++)
	  {
	    for(int k = 1;k<=bcalGeom->NBCALSECS1;k++)
	      {
		if(eUpfADC[i][j][k]>InnerThreshold||eDownfADC[i][j][k]>InnerThreshold)
		  {		
		    fADCCellCount++;
		  }
	      }
	  }
      }
    for(int i = 1;i<=48;i++)
      {
	for(int j = bcalGeom->NBCALLAYS1+1;j<=bcalGeom->NBCALLAYS2+bcalGeom->NBCALLAYS1 ;j++)
	  {
	    for(int k = 1;k<=bcalGeom->NBCALSECS2 ;k++)
	      {
		if(eUpfADC[i][j][k]> OuterThreshold ||eDownfADC[i][j][k]>OuterThreshold)
		  {	
		    fADCCellCount++;
		  }
	      }
	  }
      }
    
    if(hits->barrelEMcal->bcalfADCCells!=HDDM_NULL)free(hits->barrelEMcal->bcalfADCCells);
    hits->barrelEMcal->bcalfADCCells = make_s_BcalfADCCells(fADCCellCount);
    hits->barrelEMcal->bcalfADCCells->mult = fADCCellCount;
    
    for(unsigned int j=0; j<hits->barrelEMcal->bcalfADCCells->mult; j++)
      {
	s_BcalfADCCell_t *fADCcell = &hits->barrelEMcal->bcalfADCCells->in[j];
	
	bool Etrigg = FALSE;		 
	
	//Inner Cells
	for(int i = 1;i<=48 && !Etrigg ;i++)
	  {
	    for(int j = 1;j<=bcalGeom->NBCALLAYS1 && !Etrigg ;j++)
	      {
		for(int k = 1;k<=bcalGeom->NBCALSECS1 && !Etrigg ;k++)
		  {
		    if(eUpfADC[i][j][k]>InnerThreshold||eDownfADC[i][j][k]>InnerThreshold)
		      {			       
			fADCcell->module = i;
			fADCcell->layer = j;
			fADCcell->sector = k;
			
			if(fADCcell->bcalfADCUpHits!=HDDM_NULL)free(fADCcell->bcalfADCUpHits);
			fADCcell->bcalfADCUpHits = make_s_BcalfADCUpHits(1);
			fADCcell->bcalfADCUpHits->mult = 1;
			
			if(fADCcell->bcalfADCDownHits!=HDDM_NULL)free(fADCcell->bcalfADCDownHits);
			fADCcell->bcalfADCDownHits = make_s_BcalfADCDownHits(1);
			fADCcell->bcalfADCDownHits->mult = 1;
			
			s_BcalfADCUpHit_t *fadcuphit = &fADCcell->bcalfADCUpHits->in[0];
			s_BcalfADCDownHit_t *fadcdownhit = &fADCcell->bcalfADCDownHits->in[0];				   
			
			if( eUpfADC[i][j][k] > InnerThreshold){
			  fadcuphit->E = eUpfADC[i][j][k];
			  fadcuphit->t = tUpfADC[i][j][k];
			}
			else {fadcuphit->E = 0; fadcuphit->t = 0;}
			
			if( eDownfADC[i][j][k] > InnerThreshold){
			  fadcdownhit->E = eDownfADC[i][j][k];
			  fadcdownhit->t = tDownfADC[i][j][k];
			}
			else {fadcdownhit->E = 0; fadcdownhit->t = 0;}			       
			
			Etrigg = TRUE;
			eUpfADC[i][j][k]=eDownfADC[i][j][k]=0;
		      }			        	   
		  }
	      }
	  }
	//Outer Cells
	for(int i = 1;i<=48 && !Etrigg ;i++)
	  {
	    for(int j = bcalGeom->NBCALLAYS1+1;j<=bcalGeom->NBCALLAYS2+bcalGeom->NBCALLAYS1 && !Etrigg ;j++)
	      {
		for(int k = 1;k<=bcalGeom->NBCALSECS2 && !Etrigg ;k++)
		  {
		    if(eUpfADC[i][j][k]> OuterThreshold ||eDownfADC[i][j][k]>OuterThreshold)
		      {			  
			fADCcell->module = i;
			fADCcell->layer = j;
			fADCcell->sector = k;
			
			if(fADCcell->bcalfADCUpHits!=HDDM_NULL)free(fADCcell->bcalfADCUpHits);
			fADCcell->bcalfADCUpHits = make_s_BcalfADCUpHits(1);
			fADCcell->bcalfADCUpHits->mult = 1;
			
			if(fADCcell->bcalfADCDownHits!=HDDM_NULL)free(fADCcell->bcalfADCDownHits);
			fADCcell->bcalfADCDownHits = make_s_BcalfADCDownHits(1);
			fADCcell->bcalfADCDownHits->mult = 1;
			
			s_BcalfADCUpHit_t *fadcuphit = &fADCcell->bcalfADCUpHits->in[0];
			s_BcalfADCDownHit_t *fadcdownhit = &fADCcell->bcalfADCDownHits->in[0];
			
			if( eUpfADC[i][j][k] > OuterThreshold){
			  fadcuphit->E = eUpfADC[i][j][k];
			  fadcuphit->t = tUpfADC[i][j][k];
			}
			else {fadcuphit->E = 0; fadcuphit->t = 0;}
			
			if( eDownfADC[i][j][k] > OuterThreshold){
			  fadcdownhit->E = eDownfADC[i][j][k];
			  fadcdownhit->t = tDownfADC[i][j][k];
			}
			else {fadcdownhit->E = 0; fadcdownhit->t = 0;}			       
			
			Etrigg = TRUE;
			eUpfADC[i][j][k]=eDownfADC[i][j][k]=0;
		      }			        	   
		  }
	      }
	  }
      }		                   		       
  } //i (PhysicsEvents)
  
}

#else // Beni's version above. Dave's version below

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
	
	// Clear DBCALReadoutChannel objects
	//map<int, DBCALReadoutChannel>::iterator iter=bcal_fADCs.begin();
	//for(; iter!= bcal_fADCs.end(); iter++)iter->second.Clear();

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

				// Add dark noise hits
				upEnergy   += ( getDarkHits() * BCAL_mevPerPE * k_MeV );
				downEnergy += ( getDarkHits() * BCAL_mevPerPE * k_MeV );
	
				// now offset times for propagation distance
				float upTime = smearedtUp + upDist / bcalGeom.C_EFFECTIVE;
				float downTime = smearedtDown + downDist / bcalGeom.C_EFFECTIVE;
	
				// If energy is smeared to negative, set to 0.
				if(upEnergy <= 0){
					upEnergy = 0;
					upTime = 0;
				}
				if(downEnergy <= 0){
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
				double Edark = ( getDarkHits() * BCAL_mevPerPE * k_MeV );
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
				double Edark = ( getDarkHits() * BCAL_mevPerPE * k_MeV );
				double tdark = SampleRange( -0.25 * BCAL_INTWINDOW_NS, 0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
				bcalfADC.Edown += Edark;
				bcalfADC.tdown += tdark * Edark;
			}
			
			// normalize time
			bcalfADC.tdown /= bcalfADC.Edown;


			// If either the upstream or downstream hit is over threshold, then
			// a bcalfADCCell will need to be created so remember this id if needed.
			if((bcalfADC.Eup >= bcalfADC.threshold) || (bcalfADC.Edown >= bcalfADC.threshold)){
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
			if(bcalfADC.Eup >= bcalfADC.threshold){
				fADCcell->bcalfADCUpHits = make_s_BcalfADCUpHits(1);
				fADCcell->bcalfADCUpHits->mult = 1;
				s_BcalfADCUpHit_t *fadcuphit = &fADCcell->bcalfADCUpHits->in[0];
				
				fadcuphit->E = bcalfADC.Eup;
				fadcuphit->t = bcalfADC.tup;
			}else{
				fADCcell->bcalfADCUpHits = (s_BcalfADCUpHits_t*)HDDM_NULL;
			}

			// Downstream hit
			if(bcalfADC.Edown >= bcalfADC.threshold){
				fADCcell->bcalfADCDownHits = make_s_BcalfADCDownHits(1);
				fADCcell->bcalfADCDownHits->mult = 1;
				s_BcalfADCDownHit_t *fadcdownhit = &fADCcell->bcalfADCDownHits->in[0];
				
				fadcdownhit->E = bcalfADC.Edown;
				fadcdownhit->t = bcalfADC.Edown;
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

#endif // end of Dave's verison of SmearBCAL


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

#ifdef NO_THRESHOLD_CUT
	Bcal_CellInnerThreshold = 0.0; // n.b. the outer threshold is derived from this too
#endif

	// MeV per PhotoElectron in SiPM
	BCAL_mevPerPE = 1 /
		( BCAL_PHOTONSPERSIDEPERMEV_INFIBER * BCAL_DEVICEPDE * BCAL_SAMPLING_FRACT );
	
	// Factor to unattenuate the energy at the end to what it *would* be if it came
	// from the center of the module. This is so the proper time smearing can be
	// applied based on the 2006 beam test results which shot photons into the center
	// of the module.
	BCAL_UNATTENUATE_TO_CENTER = exp(bcalGeom.BCALFIBERLENGTH/2/ bcalGeom.ATTEN_LENGTH);
	
	// Calculate the mean energy (effective) due to dark pulses. The dark pulses
	// effectively add to the pedestal and so this value should be subtracted before
	// the threshold is applied.
	BCAL_Nsum_inner = DBCALGeometry::NSUMLAYS1 * DBCALGeometry::NSUMSECS1;
	BCAL_Nsum_outer = DBCALGeometry::NSUMLAYS2 * DBCALGeometry::NSUMSECS2;
	BCAL_mean_darkpulse_energy_SiPM = BCAL_DARKRATE_GHZ * BCAL_INTWINDOW_NS * (1.0 + BCAL_XTALK_FRACT) * BCAL_mevPerPE;
	BCAL_mean_darkpulse_energy_inner = BCAL_mean_darkpulse_energy_SiPM * (double)BCAL_Nsum_inner;
	BCAL_mean_darkpulse_energy_outer = BCAL_mean_darkpulse_energy_SiPM * (double)BCAL_Nsum_outer;
	
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
				double thresh = Bcal_CellInnerThreshold; /// FIXME!!!!
				bcal_fADCs[fADCId] = DBCALReadoutChannel(NSiPM_up, thresh, i, DBCALGeometry::fADC_layer(cellId), DBCALGeometry::fADC_sector(cellId));
			}
		}

		// outer
		unsigned int NSiPM_down = DBCALGeometry::NSUMLAYS2 * DBCALGeometry::NSUMSECS2;
		for(int j = 1; j<=4;j++){
			for(int k = 1; k<=4; k++){
				int cellId = DBCALGeometry::cellId( i, j+6, k);
				int fADCId = DBCALGeometry::fADCId( i, j+6, k);
				double thresh = Bcal_CellInnerThreshold; /// FIXME!!!!
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
#ifdef NO_E_SMEAR
return E;
#endif

    double sigmaSamp = BCAL_SAMPLINGCOEFA / sqrt( E ) + BCAL_SAMPLINGCOEFB;
    
    return( E * (1.0 + SampleGaussian(sigmaSamp)) );
}

//-----------
// bcalTimeSmear
//-----------
float bcalTimeSmear( float t, float E )
{
#ifdef NO_T_SMEAR
return t;
#endif
  double sigmaT = BCAL_TIMEDIFFCOEFA / sqrt( E ) + BCAL_TIMEDIFFCOEFB;

  return( t + SampleGaussian(sigmaT) );
}

//-----------
// getDarkHits
//-----------
int getDarkHits()
{
#ifdef NO_DARK_PULSES
	return 0;
#endif

  int darkPulse = SamplePoisson( BCAL_DARKRATE_GHZ* BCAL_INTWINDOW_NS );

  int xTalk = SamplePoisson( darkPulse * BCAL_XTALK_FRACT );

  return( xTalk + darkPulse );
}


// $Id$
//
//    File: DCDCTrackHit_factory.cc
// Created: Mon Oct 16 10:20:07 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#include <cmath>
#include <pthread.h>
using namespace std;

#include "DCDCTrackHit_factory.h"
#include "DCDCHit.h"
#include <HDGEOMETRY/DGeometry.h>
#include <TRACKING/DTrackHitSelectorTHROWN.h>
#include <TRACKING/DMCTrackHit.h>
 
DCDCTrackHit_factory::~DCDCTrackHit_factory(){
  if (cdcwires.size()){
    for (unsigned int i=0;i<cdcwires.size();i++){
      for (unsigned int j=0;j<cdcwires[i].size();j++){
	delete cdcwires[i][j];
      }
    }    
  }
}

//------------------
// init
//------------------
jerror_t DCDCTrackHit_factory::init(void)
{
   MATCH_TRUTH_HITS=false;
  
   gPARMS->SetDefaultParameter("CDC:MATCH_TRUTH_HITS",MATCH_TRUTH_HITS,"Match truth hits to CDC hits (DEF=false)");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCDCTrackHit_factory::brun(JEventLoop *loop, int runnumber)
{
  // Get pointer to DGeometry object
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  dgeom  = dapp->GetDGeometry(runnumber);
  
  // Get the CDC wire table from the XML
  //jout<< "Getting map of cdc wires from the XML" <<endl;
  dgeom->GetCDCWires(cdcwires);
  
  // Fill array with the number of straws for each layer
  for (unsigned int i=0;i<cdcwires.size();i++){
    Nstraws[i]=cdcwires[i].size();
  }

	// Get drift time parameters
	JCalibration *jcalib = dapp->GetJCalibration((loop->GetJEvent()).GetRunNumber());
	map<string, double> cdc_drift_parms;
	jcalib->Get("CDC/cdc_drift_parms", cdc_drift_parms);
	CDC_DRIFT_BSCALE_PAR1 = cdc_drift_parms["bscale_par1"];
	CDC_DRIFT_BSCALE_PAR2 = cdc_drift_parms["bscale_par2"];

	typedef map<string,double>::iterator iter_double;
  	vector< map<string, double> > tvals;
	if (jcalib->Get("CDC/cdc_drift_table", tvals)==false){    
		for(unsigned int i=0; i<tvals.size(); i++){
			map<string, double> &row = tvals[i];
			iter_double iter = row.find("t");
			cdc_drift_table.push_back(1000.*iter->second);
		}
	}
	if(cdc_drift_table.empty()){
		jerr << endl;
		jerr << " No values found for \"CDC/cdc_drift_table\"!" <<endl;
		jerr << endl;
		jerr << " This probably means you'r using an old calibration DB." << endl;
		jerr << " Check your JANA_CALIB_URL environment variable." << endl;
		jerr << " (This message printed from DCDCTrackHit_factory::brun())" << endl;
		exit(-1);
	}
	cdc_drift_table_min = cdc_drift_table[0];
	cdc_drift_table_max = cdc_drift_table[cdc_drift_table.size()-1];

  return NOERROR;
}

jerror_t DCDCTrackHit_factory::erun(void){
  if (cdcwires.size()){
    for (unsigned int i=0;i<cdcwires.size();i++){
      for (unsigned int j=0;j<cdcwires[i].size();j++){
	delete cdcwires[i][j];
      }
    }    
  }
  cdcwires.clear();
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DCDCTrackHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Convert from ring/straw indexing to x/y position
	/// of wire center and stereo angle.
	vector<const DCDCHit*> cdchits;
	loop->Get(cdchits);
	if (cdchits.size()==0) return NOERROR;
	
	// If this is simulated data then we want to match up the truth hit
	// with this "real" hit. Ideally, this would be done at the
	// DCDCHit object level, but the organization of the data in HDDM
	// makes that difficult. Here we have the full wire definition so
	// we make the connection here.
	vector<const DMCTrackHit*> mctrackhits;
	if (MATCH_TRUTH_HITS)loop->Get(mctrackhits);
	
	for(unsigned int i=0; i<cdchits.size(); i++){
		const DCDCHit* cdchit = cdchits[i];

		if(cdchit->ring>CDC_MAX_RINGS || cdchit->ring<1
			|| cdchit->straw>Nstraws[cdchit->ring-1] || cdchit->straw<1){
			cerr<<__FILE__<<":"<<__LINE__<<" Ring or straw out of range! ring="
				<<cdchit->ring<<" (should be 1-"<<CDC_MAX_RINGS<<")  straw="
				<<cdchit->straw<<" (should be 1-"<<Nstraws[cdchit->ring-1]<<")"<<endl;
			continue;
		}

		DCDCTrackHit *hit = new DCDCTrackHit;
		hit->wire = cdcwires[cdchit->ring-1][cdchit->straw-1];
		hit->is_stereo=((cdchit->ring>4&&cdchit->ring<13)
				||(cdchit->ring>16&&cdchit->ring<25))
				?true:false;
		hit->tdrift = cdchit->t;

		double w_eff=29.5e-9;
		double gas_gain=1e5;
		double electron_charge=1.6022e-4; /* fC */
		hit->dE=cdchit->q*w_eff/(gas_gain*electron_charge);
		
		// Calculate drift distance assuming no tof to wire for now.
		// The tracking package will replace this once an estimate
		// of the tof to the wire is known. This algorithm was copied
		// from the method DTrackFitterKalmanSIMD::ComputeCDCDrift
		// so that a distance could be used for drawing the CDC 
		// drift times in hdview2, even when track fitting isn't done.
		double B = 2.0; // just assume 2T B-field everywhere
		double dtc =(CDC_DRIFT_BSCALE_PAR1 + CDC_DRIFT_BSCALE_PAR2 * B)* hit->tdrift;
    	double tcorr = hit->tdrift - dtc;
		unsigned int index=0;
		if(tcorr < cdc_drift_table_min){
			index = 0;
		}else if (tcorr >= cdc_drift_table_max){
			index = cdc_drift_table.size()-1;
		}else{
			index=locate(cdc_drift_table,tcorr);
		}
		double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
		double frac=(tcorr-cdc_drift_table[index])/dt;
		hit->dist = 0.01*(double(index)+frac);  // the actual drift distance is calculated later, use a placeholder value here

		hit->AddAssociatedObject(cdchit);

		if (MATCH_TRUTH_HITS==true&&mctrackhits.size()>0){
		  // Estimate for drift distance, ignoring flight time
		  double d=0.;
		  if (cdchit->t>0) d=0.0279*sqrt(cdchit->t);

		  // Try matching truth hit with this "real" hit.
		  const DMCTrackHit *mctrackhit = DTrackHitSelectorTHROWN::GetMCTrackHit(hit->wire,d, mctrackhits);
		

		  if(mctrackhit)hit->AddAssociatedObject(mctrackhit);
		}

		_data.push_back(hit);
	}
	
	return NOERROR;
}

//------------------
// locate
// Locate a position in vector xx given x
// (this copied from DTrackFitterKalmanSIMD.cc)
//------------------
unsigned int DCDCTrackHit_factory::locate(vector<double>&xx,double x){
  int ju,jm,jl;
  int ascnd;

  int n=xx.size();

  jl=-1;
  ju=n;
  ascnd=(xx[n-1]>=xx[0]);
  while(ju-jl>1){
    jm=(ju+jl)>>1;
    if ( (x>=xx[jm])==ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x==xx[0]) return 0;
  else if (x==xx[n-1]) return n-2;
  return jl;
}


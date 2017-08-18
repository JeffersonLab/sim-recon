// Test version V0 of the L1 trigger simulation, Aug 18, 2017 (asomov)

#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JApplication.h>
#include <DAQ/DCODAROCInfo.h>
#include <DAQ/DL1Info.h>

using namespace jana;

#include "DL1MCTrigger_factory.h"


//------------------
// init
//------------------
jerror_t DL1MCTrigger_factory::init(void)
{

  // Default parameters for the main production trigger are taken from the 
  // spring run of 2017: 25 F + B > 45000

  FCAL_ADC_PER_MEV =  3.73;
  FCAL_CELL_THR    =  65;
  FCAL_EN_SC       =  25;
  FCAL_NSA         =  10;
  FCAL_NSB         =  3;
  FCAL_WINDOW      =  10;

  BCAL_ADC_PER_MEV =  34.48276;
  BCAL_CELL_THR    =  20;
  BCAL_EN_SC       =  1;
  BCAL_NSA         =  19;
  BCAL_NSB         =  3;
  BCAL_WINDOW      =  20;

  FCAL_BCAL_EN     =  45000; 

  ST_ADC_PER_MEV   =  1.;
  ST_CELL_THR      =  40;
  ST_NSA           =  10;
  ST_NSB           =  3;
  ST_WINDOW        =  10;
  ST_NHIT          =  1;

  BCAL_OFFSET      =  2;


  gPARMS->SetDefaultParameter("TRIG:FCAL_ADC_PER_MEV", FCAL_ADC_PER_MEV,
			      "FCAL energy calibration for the Trigger");
  gPARMS->SetDefaultParameter("TRIG:FCAL_CELL_THR", FCAL_CELL_THR,
			      "FCAL energy threshold per cell");
  gPARMS->SetDefaultParameter("TRIG:FCAL_EN_SC", FCAL_EN_SC,
			      "FCAL energy threshold");
  gPARMS->SetDefaultParameter("TRIG:FCAL_NSA", FCAL_NSA,
			      "FCAL NSA");
  gPARMS->SetDefaultParameter("TRIG:FCAL_NSB", FCAL_NSB,
			      "FCAL NSB");
  gPARMS->SetDefaultParameter("TRIG:FCAL_WINDOW", FCAL_WINDOW,
			      "FCAL GTP integration window");

  gPARMS->SetDefaultParameter("TRIG:BCAL_ADC_PER_MEV", BCAL_ADC_PER_MEV,
			      "BCAL energy calibration for the Trigger");
  gPARMS->SetDefaultParameter("TRIG:BCAL_CELL_THR", BCAL_CELL_THR,
			      "BCAL energy threshold per cell");
  gPARMS->SetDefaultParameter("TRIG:BCAL_EN_SC", BCAL_EN_SC,
			      "BCAL energy threshold");
  gPARMS->SetDefaultParameter("TRIG:BCAL_NSA", BCAL_NSA,
			      "BCAL NSA");
  gPARMS->SetDefaultParameter("TRIG:BCAL_NSB", BCAL_NSB,
			      "BCAL NSB");
  gPARMS->SetDefaultParameter("TRIG:BCAL_WINDOW", BCAL_WINDOW,
			      "BCAL GTP integration window");

  gPARMS->SetDefaultParameter("TRIG:ST_ADC_PER_MEV", ST_ADC_PER_MEV,
			      "ST energy calibration for the Trigger");
  gPARMS->SetDefaultParameter("TRIG:ST_CELL_THR", ST_CELL_THR,
			      "ST energy threshold per cell");
  gPARMS->SetDefaultParameter("TRIG:ST_NSA", ST_NSA,
			      "ST NSA");
  gPARMS->SetDefaultParameter("TRIG:ST_NSB", ST_NSB,
			      "ST NSB");  
  gPARMS->SetDefaultParameter("TRIG:ST_WINDOW", ST_WINDOW,
			      "ST window for merging hits (GTP)");
  gPARMS->SetDefaultParameter("TRIG:ST_NHIT", ST_NHIT,
			      "Number of hits in ST");
  
  gPARMS->SetDefaultParameter("TRIG:FCAL_BCAL_EN", FCAL_BCAL_EN,
			      "Energy threshold for the FCAL & BCAL trigger");

  gPARMS->SetDefaultParameter("TRIG:BCAL_OFFSET", BCAL_OFFSET,
			      "Timing offset between BCAL and FCAL energies at GTP (sampels)");
  

  time_shift = 100;

  time_min  =  0;
  time_max  =  (sample - 1)*max_adc_bins;


  return NOERROR;
}


//------------------
// brun
//------------------
jerror_t DL1MCTrigger_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{

  int use_rcdb = 1;

  int status   = 0;

  fcal_trig_mask.clear();
  bcal_trig_mask.clear();

  triggers_enabled.clear();

  // Only print messages for one thread whenever run number change

  static pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;
  static set<int> runs_announced;
  pthread_mutex_lock(&print_mutex);
  bool print_messages = false;

  if(runs_announced.find(runnumber) == runs_announced.end()){
    print_messages = true;
    runs_announced.insert(runnumber);
  }
  
  pthread_mutex_unlock(&print_mutex);
  
  // Don't use RCDB for mc_generic: ideal MC simulation

  string JANA_CALIB_CONTEXT = getenv("JANA_CALIB_CONTEXT");

  if(getenv("JANA_CALIB_CONTEXT") != NULL ){ 
    if(JANA_CALIB_CONTEXT.find("mc_generic") != string::npos){
      use_rcdb = 0;
    }
  }

  //  runnumber = 30942;

  if(use_rcdb == 1){
    status = Read_RCDB(runnumber);
    PrintTriggers();
  }

  
  if( (use_rcdb == 0) || (status > 0) || (triggers_enabled.size() == 0)){

    // Simulate FCAL & BCAL main production trigger only
  
    trigger_conf trig_tmp;
    trig_tmp.bit = 0;
    trig_tmp.gtp.fcal   = FCAL_EN_SC;
    trig_tmp.gtp.bcal   = BCAL_EN_SC;
    trig_tmp.gtp.en_thr = FCAL_BCAL_EN;
    triggers_enabled.push_back(trig_tmp);

    cout << " Do not use RCDB for the trigger simulation. Default (spring 2017) trigger settings are used " << endl;
  }

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DL1MCTrigger_factory::evnt(JEventLoop *loop, uint64_t eventnumber){

        int l1_found = 1;  

	int status = 0;

	fcal_signal_hits.clear();
	bcal_signal_hits.clear();

	fcal_merged_hits.clear();
	bcal_merged_hits.clear();

	memset(fcal_ssp,0,sizeof(fcal_ssp));
	memset(fcal_gtp,0,sizeof(fcal_gtp));

	memset(bcal_ssp,0,sizeof(bcal_ssp));
	memset(bcal_gtp,0,sizeof(bcal_gtp));


        vector<const DFCALHit*>  fcal_hits;
	vector<const DBCALHit*>  bcal_hits;

	loop->Get(fcal_hits);
	loop->Get(bcal_hits);


	DL1MCTrigger *trigger = new DL1MCTrigger;

	//  FCAL energy sum	
	double fcal_hit_en = 0;
	
	for (unsigned int ii = 0; ii < fcal_hits.size(); ii++){
	  int row  = fcal_hits[ii]->row;
	  int col  = fcal_hits[ii]->column;

	  // Shift time to simulate pile up hits
	  double time = fcal_hits[ii]->t + time_shift;
	  if((time < time_min) || (time > time_max)){
	    continue;
	  }

	  // Check channels masked for trigger
	  int ch_masked = 0;
	  for(unsigned int jj = 0; jj < fcal_trig_mask.size(); jj++){
	    if( (row == fcal_trig_mask[jj].row) && (col == fcal_trig_mask[jj].col)){	      
	      ch_masked = 1;
	      break;
	    } 
	  }
	  
	  if(ch_masked == 0){
	    fcal_hit_en += fcal_hits[ii]->E;
	    
	    fcal_signal fcal_tmp;
	    
	    fcal_tmp.row     = row;
	    fcal_tmp.column  = col;
	    
	    fcal_tmp.energy  = fcal_hits[ii]->E;
	    fcal_tmp.time    = time;
	    memset(fcal_tmp.adc_amp,0,sizeof(fcal_tmp.adc_amp));
	    
	    double fcal_adc_en  = fcal_tmp.energy*FCAL_ADC_PER_MEV*1000;

	    status = Digitize(fcal_adc_en, fcal_tmp.time, fcal_tmp.adc_amp, 1);
	    status = 0;

	    fcal_signal_hits.push_back(fcal_tmp);
	  }
	  
	}       
	

	// Merge FCAL hits
	for(unsigned int ii = 0; ii < fcal_signal_hits.size(); ii++){	  

	  if(fcal_signal_hits[ii].merged == 1) continue;
	  
	  fcal_signal fcal_tmp;
	  fcal_tmp.row     = fcal_signal_hits[ii].row;
	  fcal_tmp.column  = fcal_signal_hits[ii].column;
	  fcal_tmp.energy  = 0.;
	  fcal_tmp.time    = 0.;
	  for(int kk = 0; kk < sample; kk++)	      
	    fcal_tmp.adc_amp[kk] = fcal_signal_hits[ii].adc_amp[kk];
	  
	  for(unsigned int jj = ii + 1; jj < fcal_signal_hits.size(); jj++){
	    if((fcal_signal_hits[ii].row  == fcal_signal_hits[jj].row) &&
	       (fcal_signal_hits[ii].column == fcal_signal_hits[jj].column)){

	      fcal_signal_hits[jj].merged = 1;

	      for(int kk = 0; kk < sample; kk++)	      
		fcal_tmp.adc_amp[kk] += fcal_signal_hits[jj].adc_amp[kk];	
	    }
	  }
	  
	  fcal_merged_hits.push_back(fcal_tmp);
	}	
	
	
	int fcal_hit_adc_en = 0;

	for(unsigned int ii = 0; ii < fcal_signal_hits.size(); ii++)
	  for(int jj = 0; jj < sample; jj++)
	    fcal_hit_adc_en += fcal_signal_hits[ii].adc_amp[jj];	  
	
	
	status += FADC_SSP(fcal_merged_hits, 1);

	status += GTP(1);


	// BCAL	energy sum
	double bcal_hit_en = 0;
	
	for (unsigned int ii = 0; ii < bcal_hits.size(); ii++){

	  // Shift time to simulate pile up hits
	  double time = bcal_hits[ii]->t + time_shift;	 
	  if((time < time_min) || (time > time_max)){
	    continue;
	  }
	  
	  int module = bcal_hits[ii]->module;
	  int layer  = bcal_hits[ii]->layer;
	  int sector = bcal_hits[ii]->sector;
	  int end    = bcal_hits[ii]->end;
	  
	  // Check trigger masks
	  int ch_masked = 0;

	  for(unsigned int jj = 0; jj < bcal_trig_mask.size(); jj++){
	    if( (module == bcal_trig_mask[jj].module) && (layer == bcal_trig_mask[jj].layer) && 
		(sector == bcal_trig_mask[jj].sector) && (end == bcal_trig_mask[jj].end)  ){ 
	      ch_masked = 1;
	      break;
	    } 
	  }
	  

	  if(ch_masked == 0){
	    bcal_hit_en += bcal_hits[ii]->E;
	    
	    bcal_signal bcal_tmp;	    
	    bcal_tmp.module  = module;
	    bcal_tmp.layer   = layer;
  	    bcal_tmp.sector  = sector;
	    bcal_tmp.end     = end;

	    bcal_tmp.energy  = bcal_hits[ii]->E;
	    bcal_tmp.time    = time;
	    memset(bcal_tmp.adc_amp,0,sizeof(bcal_tmp.adc_amp));
	    
	    double bcal_adc_en  = bcal_tmp.energy*BCAL_ADC_PER_MEV*1000;

	    status = Digitize(bcal_adc_en, bcal_tmp.time, bcal_tmp.adc_amp, 2);
	    status = 0;
 
	    bcal_signal_hits.push_back(bcal_tmp);
	  }	  
	}       


	// Merge BCAL hits	
	for(unsigned int ii = 0; ii < bcal_signal_hits.size(); ii++){	  
	  
	  if(bcal_signal_hits[ii].merged == 1) continue;
	  
	  bcal_signal bcal_tmp;
	  bcal_tmp.module  =  bcal_signal_hits[ii].module;
	  bcal_tmp.layer   =  bcal_signal_hits[ii].layer;
	  bcal_tmp.sector  =  bcal_signal_hits[ii].sector;
	  bcal_tmp.end    =  bcal_signal_hits[ii].end;
	  
	  bcal_tmp.energy  = 0.;
	  bcal_tmp.time    = 0.;
	  
	  for(int kk = 0; kk < sample; kk++)	      
	    bcal_tmp.adc_amp[kk] = bcal_signal_hits[ii].adc_amp[kk];
	  
	  for(unsigned int jj = ii + 1; jj < bcal_signal_hits.size(); jj++){
	    if((bcal_signal_hits[ii].module == bcal_signal_hits[jj].module) &&
	       (bcal_signal_hits[ii].layer  == bcal_signal_hits[jj].layer) &&
	       (bcal_signal_hits[ii].sector == bcal_signal_hits[jj].sector) &&    
	       (bcal_signal_hits[ii].end   == bcal_signal_hits[jj].end)){

	      bcal_signal_hits[jj].merged = 1;

	      for(int kk = 0; kk < sample; kk++)	      
		bcal_tmp.adc_amp[kk] += bcal_signal_hits[jj].adc_amp[kk];	
	    }
	  }
	  
	  bcal_merged_hits.push_back(bcal_tmp);
	}
	

	int bcal_hit_adc_en = 0;
	for(unsigned int ii = 0; ii < bcal_signal_hits.size(); ii++)	  
	  for(int jj = 0; jj < sample; jj++)
	    bcal_hit_adc_en += bcal_signal_hits[ii].adc_amp[jj];	  
	
	
	status = FADC_SSP(bcal_merged_hits, 2);

	status = GTP(2);

	// Search for triggers

	l1_found = FindTriggers(trigger);
	
	if(l1_found){
	  
	  int fcal_gtp_max = 0;
	  int bcal_gtp_max = 0;
	  
	  for(unsigned int ii = 0; ii < sample; ii++){
	    if(fcal_gtp[ii] > fcal_gtp_max) fcal_gtp_max = fcal_gtp[ii];
	    if(bcal_gtp[ii] > bcal_gtp_max) bcal_gtp_max = bcal_gtp[ii];
	  }	
	  	  
	  trigger->fcal_en      =  fcal_hit_en;
	  trigger->fcal_adc     =  fcal_hit_adc_en;
	  trigger->fcal_adc_en  =  fcal_hit_adc_en/FCAL_ADC_PER_MEV/1000.;
	  trigger->fcal_gtp     =  fcal_gtp_max;
	  trigger->fcal_gtp_en  =  fcal_gtp_max/FCAL_ADC_PER_MEV/1000.;
	  
	  trigger->bcal_en      =  bcal_hit_en;
	  trigger->bcal_adc     =  bcal_hit_adc_en;
	  trigger->bcal_adc_en  =  bcal_hit_adc_en/BCAL_ADC_PER_MEV/2./1000.;
	  trigger->bcal_gtp     =  bcal_gtp_max;
	  trigger->bcal_gtp_en  =  bcal_gtp_max/BCAL_ADC_PER_MEV/2./1000.;	  	  
	  
	  _data.push_back(trigger);	 
	  
	} else{
	  delete trigger;
	}
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DL1MCTrigger_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DL1MCTrigger_factory::fini(void)
{
	return NOERROR;
}

//*********************
// Read RCDB 
//*********************

int  DL1MCTrigger_factory::Read_RCDB(int32_t runnumber){

  vector<const DTranslationTable*> ttab;  
  eventLoop->Get(ttab);
  
  vector<string> SectionNames = {"TRIGGER", "GLOBAL", "FCAL", "BCAL", "TOF", "ST", "TAGH",
				 "TAGM", "PS", "PSC", "TPOL", "CDC", "FDC"};
  
  string RCDB_CONNECTION;

  if( getenv("RCDB_CONNECTION")!= NULL )
    RCDB_CONNECTION = getenv("RCDB_CONNECTION");
  else
    RCDB_CONNECTION = "mysql://rcdb@hallddb.jlab.org/rcdb";   // default to outward-facing MySQL DB


  rcdb::Connection connection(RCDB_CONNECTION);
  
  auto  rtvsCnd  =  connection.GetCondition(runnumber, "rtvs");
  
  if( !rtvsCnd ) {
    cout<<"'rtvs' condition is not set for run " << runnumber << endl;
    return 2;
  }
  
  
  auto json = rtvsCnd->ToJsonDocument();               // The CODA rtvs is serialized as JSon dictionary.  
  string fileName(json["%(config)"].GetString());      // We need item with name '%(config)'
  
  auto file = connection.GetFile(runnumber, fileName);
  
  if(!file) {                                                         // If there is no such file, null is returned
    cout<<"File with name: "<< fileName
	<<" doesn't exist (not associated) with run: "<< runnumber << endl;
    return 3;
  }
  
  string fileContent = file->GetContent();                               // Get file content
  auto result = rcdb::ConfigParser::Parse(fileContent, SectionNames);    // Parse it!
  
  
  // Pulse parameters for trigger

  auto trig_thr = result.Sections["FCAL"].NameValues["FADC250_TRIG_THR"];
  auto trig_nsb = result.Sections["FCAL"].NameValues["FADC250_TRIG_NSB"];
  auto trig_nsa = result.Sections["FCAL"].NameValues["FADC250_TRIG_NSA"];
  
  if(trig_thr.size() > 0){
    FCAL_CELL_THR  =  stoi(trig_thr) - 100;
    if(FCAL_CELL_THR < 0) FCAL_CELL_THR = 0;
  }

  if(trig_nsb.size() > 0)
    FCAL_NSB  =  stoi(trig_nsb);

  if(trig_thr.size() > 0)
    FCAL_NSA   =  stoi(trig_nsa);

  trig_thr = result.Sections["BCAL"].NameValues["FADC250_TRIG_THR"];
  trig_nsb = result.Sections["BCAL"].NameValues["FADC250_TRIG_NSB"];
  trig_nsa = result.Sections["BCAL"].NameValues["FADC250_TRIG_NSA"];
  
  if(trig_thr.size() > 0){
    BCAL_CELL_THR   =  stoi(trig_thr) - 100;
    if(BCAL_CELL_THR < 0) BCAL_CELL_THR = 0;
  }

  if(trig_nsb.size() > 0)
    BCAL_NSB  =  stoi(trig_nsb);

  if(trig_thr.size() > 0)
    BCAL_NSA  =  stoi(trig_nsa);


  // List of enabled GTP equations and triggers
  vector<vector<string>> triggerTypes;

  for(auto row : result.Sections["TRIGGER"].Rows) {    

    if(row[0] == "TRIG_TYPE") {
      
      if(row.size() >= 9){
	triggerTypes.push_back(row);       // The line in config file starts with TRIG_TYPE
      } else {
	cout <<  " Cannot parse TRIG_TYPE. Insufficient number of parameters " <<  row.size() << endl; 
      }      
    }    
   

    // Integration windows for BCAL and FCAL 
    if(row[0] == "TRIG_EQ") {
      if(row.size() >= 5){
	if(stoi(row[4]) == 1){              // Trigger enabled 
	  if(row[1] == "FCAL")
	    FCAL_WINDOW = stoi(row[3]);  
	  if(row[1] == "BCAL_E")
	    BCAL_WINDOW = stoi(row[3]);  
	}
      }            
    }
  }
  
  

  
  for(int ii = 0; ii < 32; ii++){
    
    int trig_found = 0;
    
    trigger_conf trigger_tmp;

    memset(&trigger_tmp,0,sizeof(trigger_tmp));
    
    trigger_tmp.bit = ii;   
    
    for(unsigned int jj = 0; jj < triggerTypes.size(); jj++){
      
      if(triggerTypes[jj].size() < 9) continue;   // Minimum  9 parameters are required in the config file 

      int bit = stoi(triggerTypes[jj][8]);  // Trigger bit
      
      if( bit == ii){      
		
	// FCAL, BCAL triggers 
	if(triggerTypes[jj][1] == "BFCAL"){
	  
	  int fcal = 0; 
	  int bcal = 0; 

	  if(triggerTypes[jj][4].size() > 0) fcal = stoi(triggerTypes[jj][4]); 
	  if(triggerTypes[jj][5].size() > 0) bcal = stoi(triggerTypes[jj][5]); 


	  if( (fcal > 0) && (bcal > 0)){             // FCAL & BCAL 
	    trigger_tmp.type  =  0x3;
	    trigger_tmp.gtp.fcal      =   fcal;
	    trigger_tmp.gtp.bcal      =   bcal;
	    if(triggerTypes[jj][6].size() > 0) trigger_tmp.gtp.en_thr = stoi(triggerTypes[jj][6]);
	    
	    if(triggerTypes[jj].size() > 9){
	      if((triggerTypes[jj].size() >= 10) && (triggerTypes[jj][9].size() > 0))  trigger_tmp.gtp.fcal_min  =   stoi(triggerTypes[jj][9]);
	      if((triggerTypes[jj].size() >= 11) && (triggerTypes[jj][10].size() > 0)) trigger_tmp.gtp.fcal_max  =   stoi(triggerTypes[jj][10]);
	      if((triggerTypes[jj].size() >= 12) && (triggerTypes[jj][11].size() > 0)) trigger_tmp.gtp.bcal_min  =   stoi(triggerTypes[jj][11]);
	      if((triggerTypes[jj].size() >= 13) && (triggerTypes[jj][12].size() > 0)) trigger_tmp.gtp.bcal_max  =   stoi(triggerTypes[jj][12]);
	    }
	    
	    trig_found++;
	  } else if ((fcal > 0) && (bcal == 0)){      // FCAL only
	    trigger_tmp.type  =  0x1;
	    trigger_tmp.gtp.fcal      =   fcal;
	    if(triggerTypes[jj][6].size() > 0)  trigger_tmp.gtp.en_thr      =   stoi(triggerTypes[jj][6]);
	    if(triggerTypes[jj].size() > 9){
	      if(triggerTypes[jj][9].size() > 0) trigger_tmp.gtp.fcal_min   =   stoi(triggerTypes[jj][9]);
	      if(triggerTypes[jj][10].size() > 0) trigger_tmp.gtp.fcal_max  =   stoi(triggerTypes[jj][10]);
	    }
	    trig_found++;
	  } else if ((bcal > 0) && (fcal == 0)){      // BCAL only
	    trigger_tmp.type  =  0x2;
	    trigger_tmp.gtp.bcal      =   bcal;
	    if(triggerTypes[jj][6].size() > 0)  trigger_tmp.gtp.en_thr    =   stoi(triggerTypes[jj][6]);
	    if((triggerTypes[jj].size() >= 12) && (triggerTypes[jj][11].size() > 0))  trigger_tmp.gtp.bcal_min  =   stoi(triggerTypes[jj][11]);
	    if((triggerTypes[jj].size() >= 13) && (triggerTypes[jj][12].size() > 0))  trigger_tmp.gtp.bcal_max  =   stoi(triggerTypes[jj][12]);

	    
	    trig_found++;
	  } else {
	    cout << " Incorrect parameters for BFCAL trigger " << endl;
	  }
	  
	} else if(triggerTypes[jj][1] == "ST"){          //  ST 
	  trigger_tmp.type    =  0x4;
	  if(triggerTypes[jj][7].size() > 0)   trigger_tmp.gtp.st_nhit     =  stoi(triggerTypes[jj][7]);
	  if((triggerTypes[jj].size() >= 14) && (triggerTypes[jj][13].size() > 0))  trigger_tmp.gtp.st_pattern  =  stoi(triggerTypes[jj][13]);

	  trig_found++;
	}   else if(triggerTypes[jj][1] == "PS"){        //  PS
	  trigger_tmp.type    =  0x8;
	  if(triggerTypes[jj][7].size() > 0)  trigger_tmp.gtp.ps_nhit      =  stoi(triggerTypes[jj][7]);
	  if((triggerTypes[jj].size() >= 14) && (triggerTypes[jj][13].size() > 0)) trigger_tmp.gtp.ps_pattern   =  stoi(triggerTypes[jj][13]);
	  trig_found++;
	}  else if(triggerTypes[jj][1] == "TAGH"){       //  TAGH
	  trigger_tmp.type     =  0x10;
	  if(triggerTypes[jj][7].size() > 0)  trigger_tmp.gtp.tof_nhit      =  stoi(triggerTypes[jj][7]);
	  if((triggerTypes[jj].size() >= 14) && (triggerTypes[jj][13].size() > 0)){
	    trigger_tmp.gtp.tagh_pattern   =  stoi(triggerTypes[jj][13],nullptr,0);
	  }

	  trig_found++;
	}  else if(triggerTypes[jj][1] == "TOF"){        //  TOF
	  trigger_tmp.type     =  0x20;
	  if(triggerTypes[jj][7].size() > 0)  trigger_tmp.gtp.ps_nhit      =  stoi(triggerTypes[jj][7]);
	  if((triggerTypes[jj].size() >= 14) && (triggerTypes[jj][13].size() > 0)) trigger_tmp.gtp.tof_pattern   =  stoi(triggerTypes[jj][13]);
	  trig_found++;
	} else {
	  cout << " Incorrect Trigger type " << triggerTypes[jj][1] << endl;
	}
	
	if(trig_found > 0){
	  triggers_enabled.push_back(trigger_tmp);
	}

      }  // Trigger lane is enabled      
    }    // Trigger types
      
  }
  


  // Load FCAL Trigger Masks

  string comDir  = result.Sections["FCAL"].NameValues["FADC250_COM_DIR"];
  string comVer  = result.Sections["FCAL"].NameValues["FADC250_COM_VER"];
  string userDir = result.Sections["FCAL"].NameValues["FADC250_USER_DIR"];
  string userVer = result.Sections["FCAL"].NameValues["FADC250_USER_VER"];
  
  
  for(int crate = 1; crate <= 12; crate++){
    
    std::string s = std::to_string(crate);
    
    string comFileName   =  comDir  +  "/rocfcal" +s + "_fadc250_" + comVer + ".cnf";
    string userFileName  =  userDir + "/rocfcal"  +s + "_" + userVer + ".cnf";
        
    auto userFile = connection.GetFile(runnumber, userFileName);
    
    if(!userFile) {                                                         
      //      cout<<" USER File with name: "<< userFileName
      //	  <<" doesn't exist (not associated) with run: "<< runnumber <<endl;
      continue;      
    }
    
    
    auto userParseResult = rcdb::ConfigParser::ParseWithSlots(userFile->GetContent(), "FADC250_SLOTS");
    
    
    for(unsigned int slot = 3; slot <= 21; slot++){
      
      auto userValues = userParseResult.SectionsBySlotNumber[slot].NameVectors["FADC250_TRG_MASK"];  // Parse it and return
      
      if(userValues.size() > 0){
	for (unsigned int ch = 0; ch < userValues.size(); ++ch) {

	  if(userValues[ch].size() == 0) continue;
	  
	  if(stoi(userValues[ch]) > 0){

	    uint32_t roc_id = 10 + crate;

	    DTranslationTable::csc_t daq_index = {roc_id, slot, ch };
	    
	    DTranslationTable::DChannelInfo channel_info;
	    
	    try {		
	      channel_info = ttab[0]->GetDetectorIndex(daq_index);	      
	    }
	    
	    catch(...){
	      cout << "Exception: FCAL channel is not in the translation table  " <<  " Crate = " << 10 + crate << "  Slot = " << slot << 
		" Channel = " << ch << endl;
	      continue;
	    }
	    	    
	    fcal_mod tmp;

	    tmp.roc     =  crate;	    
	    tmp.slot    =  slot;
	    tmp.ch      =  ch;
	    
	    tmp.row = channel_info.fcal.row;
	    tmp.col = channel_info.fcal.col;
	    
	    fcal_trig_mask.push_back(tmp);
	  }
	  
	}
	
      }	
      
    }  // Loop over slots
  }    // Loop over crates       
  
  
  

  // Load BCAL Trigger Masks

  comDir  = result.Sections["BCAL"].NameValues["FADC250_COM_DIR"];
  comVer  = result.Sections["BCAL"].NameValues["FADC250_COM_VER"];
  userDir = result.Sections["BCAL"].NameValues["FADC250_USER_DIR"];
  userVer = result.Sections["BCAL"].NameValues["FADC250_USER_VER"];
  

  for(int crate = 1; crate <= 12; crate++){
    
    if( (crate == 3) || (crate == 6) || (crate == 9) || (crate == 12)) continue;

    std::string s = std::to_string(crate);
    
    string comFileName   =  comDir  +  "/rocbcal" +s + "_fadc250_" + comVer + ".cnf";
    string userFileName  =  userDir + "/rocbcal"  +s + "_" + userVer + ".cnf";
    
    auto userFile = connection.GetFile(runnumber, userFileName);
    
    if(!userFile){
      //      cout<<" USER File with name: "<< userFileName
      //	  <<" doesn't exist (not associated) with run: "<< runnumber <<endl;                                               
      continue;
    }
    
    
    auto userParseResult = rcdb::ConfigParser::ParseWithSlots(userFile->GetContent(), "FADC250_SLOTS");
      
    for(unsigned int slot = 3; slot <= 21; slot++){
      
      auto userValues = userParseResult.SectionsBySlotNumber[slot].NameVectors["FADC250_TRG_MASK"];  // Parse it and return
      
      if(userValues.size() > 0){

	for (unsigned int ch = 0; ch < userValues.size(); ++ch) {

	  if(userValues[ch].size() == 0) continue;
	  
	  if(stoi(userValues[ch]) > 0){

	    uint32_t roc_id = 30 + crate;

	    DTranslationTable::csc_t daq_index = {roc_id, slot, ch };
	    
	    DTranslationTable::DChannelInfo channel_info;
	    
	    try {		
	      channel_info = ttab[0]->GetDetectorIndex(daq_index);	      
	    }
	    
	    catch(...){
	      cout << "Exception: BCAL channel is not in the translation table  " <<  " Crate = " << 30 + crate << "  Slot = " << slot << 
		" Channel = " << ch << endl;
	      continue;
	    }
	    
	    
	    bcal_mod tmp;

	    tmp.roc     =  crate;	    
	    tmp.slot    =  slot;
	    tmp.ch      =  ch;

	    tmp.module  =  channel_info.bcal.module;
	    tmp.layer   =  channel_info.bcal.layer;
	    tmp.sector  =  channel_info.bcal.sector;	      
	    tmp.end     =  channel_info.bcal.end;

	    bcal_trig_mask.push_back(tmp);
	  }
	  
	}
	
      }	
      
    }  // Loop over slots
  }    // Loop over crates       
  
  return 0;
}


int  DL1MCTrigger_factory::Digitize(double en, double time, int amp_array[sample], int type){
   

  // Parameterize and digitize pulse shapes. Sum up amplitudes
  // type = 1 - FCAL
  // type = 2 - BCAL

  float exp_par = 0.358;

  int pulse_length = 20;
  
  if(type == 2) exp_par = 0.18;
  
  int sample_first = (int)floor(time/time_stamp);
  
  // digitization range
  int ind_min = sample_first + 1;
  int ind_max = ind_min + pulse_length + 1;
    
  if( (ind_min > sample) || (ind_min < 0)){
    //    cout << " Digitize() FATAL error: time out of range   "  <<  time <<  "     " << ind_min << endl;
    return 1;
  }
  
  if(ind_max > sample){
    //    cout << " Digitize: ind_max set to maximum" << time <<  "  "  << ind_max << endl;
    ind_max = sample - 1;
  }
  
  for(int i = ind_min; i < ind_max; i++ ){
    double adc_t  =  time_stamp*i - time;
    double amp    =  exp_par*exp_par*exp(-adc_t*exp_par)*adc_t;
    
    amp_array[i] += (int)(amp*time_stamp*en + 0.5);

    if(amp_array[i] > max_adc_bins){
      amp_array[i] = max_adc_bins;
    }    
  }
  
  return 0;  
}

int DL1MCTrigger_factory::GTP(int detector){

  // 1  - FCAL
  // 2  - BCAL

  int INT_WINDOW = 20; 

  switch(detector){
  case 1:
    INT_WINDOW =  FCAL_WINDOW;
    break;
  case 2:
    INT_WINDOW =  BCAL_WINDOW;
    break;
  default:
    break;
  }

  int index_min = 0;
  int index_max = 0;

  for(unsigned int samp = 0; samp < sample; samp++){

    index_max = samp;
    index_min = samp - INT_WINDOW;
    
    if(index_min < 0) index_min = 0;

    int energy_sum = 0;

    for(int ii = index_min; ii <= index_max; ii++){
      if(detector == 1)
	energy_sum += fcal_ssp[ii];
      else 
	energy_sum += bcal_ssp[ii];
    }

    if(detector == 1)
      fcal_gtp[samp] = energy_sum;
    else 
      bcal_gtp[samp] = energy_sum;
  }  

  return 0;

}


template <typename T>  int DL1MCTrigger_factory::FADC_SSP(vector<T> merged_hits, int detector){
  
  //  1  - FCAL
  //  2  - BCAL

  int EN_THR =  4096; 
  int NSA    =  10;
  int NSB    =  3;;
  
  switch(detector){
  case 1:
    EN_THR =  FCAL_CELL_THR;
    NSA    =  FCAL_NSA;
    NSB    =  FCAL_NSB;
    break;
  case 2:
    EN_THR =  BCAL_CELL_THR;
    NSA    =  BCAL_NSA;
    NSB    =  BCAL_NSB;
    break;
  default:
    break;
  }
  
  for(unsigned int hit = 0; hit < merged_hits.size(); hit++){

    int index_min  =  -10;
    int index_max  =  -10;
    
    for(int ii = 0; ii < sample; ii++){

      int pulse_found = 0;

      if(merged_hits[hit].adc_amp[ii] >= EN_THR){
	pulse_found = 1;
      } else {
	continue;
      }
      
      index_min = ii - NSB;

      if(index_max > index_min) index_min = index_max + 1;

      index_max = ii + NSA - 1;

      if(index_min < 0) index_min = 0;

      if(index_max >= sample){
	index_max = sample - 1;
      }

      int extend_nsa = 1;

      // Extend FADC range if needed

      while(extend_nsa){
	
	int index_tmp = index_max + 1;
	
	if(index_tmp < sample){
	  if(merged_hits[hit].adc_amp[index_tmp] >= EN_THR){
	    index_max += NSA;
	  } else extend_nsa = 0;
	} else extend_nsa = 0;	
      }
            
      if(index_max >= sample)
	index_max = sample - 1;
          
      for(int kk = index_min; kk <= index_max; kk++){
	if(detector == 1)
	  fcal_ssp[kk] += merged_hits[hit].adc_amp[kk];
	else if(detector == 2)
	  bcal_ssp[kk] += merged_hits[hit].adc_amp[kk];
      }

      if(pulse_found == 1){
	ii = index_max + 1;
	pulse_found = 0;
      }
      
    }        
    
  }

  return 0;
}

void DL1MCTrigger_factory::PrintTriggers(){
  
  string detector;
  int nhit = 0;
  unsigned int pattern = 0;

  cout << endl << endl;
  cout << " ------------  Trigger Settings --------------- " << endl;
  cout << endl << endl;

  cout << "----------- FCAL ----------- " << endl << endl;

  cout << "FCAL_CELL_THR  = " <<  setw(10) <<  FCAL_CELL_THR <<  endl;
  cout << "FCAL_NSA       = " <<  setw(10) <<  FCAL_NSA      <<  endl;
  cout << "FCAL_NSB       = " <<  setw(10) <<  FCAL_NSB      <<  endl;
  cout << "FCAL_WINDOW    = " <<  setw(10) <<  FCAL_WINDOW   <<  endl;

  cout << endl;
  
  cout << "----------- BCAL ----------- " << endl << endl;

  cout << "BCAL_CELL_THR  = " <<  setw(10) <<  BCAL_CELL_THR <<  endl;
  cout << "BCAL_NSA       = " <<  setw(10) <<  BCAL_NSA      <<  endl;
  cout << "BCAL_NSB       = " <<  setw(10) <<  BCAL_NSB      <<  endl;
  cout << "BCAL_WINDOW    = " <<  setw(10) <<  BCAL_WINDOW   <<  endl;

  cout << endl << endl;


  if(triggers_enabled.size() > 0){
    cout  << "TYPE     "  << "FCAL_E   " << "BCAL_E   " << 
      "EN_THR   " << "NHIT     "  << "LANE     " << "FCAL_EMIN    " << "FCAL_EMAX   "  <<
      "BCAL_EMIN   " << "BCAL_EMAX   " << "PATTERN   " << endl;
  }


  for(unsigned int ii = 0; ii < triggers_enabled.size(); ii++){
    
    switch(triggers_enabled[ii].type){
    case 1: detector =  "BFCAL  ";
      break;
    case 2: detector =  "BFCAL  ";
      break;
    case 3: detector =  "BFCAL  ";
      break;
    case 4: detector =  "ST     ";
      nhit    = triggers_enabled[ii].gtp.st_nhit;
      pattern = triggers_enabled[ii].gtp.st_pattern;
      break;
    case 8: detector = "PS     ";
      nhit    = triggers_enabled[ii].gtp.ps_nhit;
      pattern = triggers_enabled[ii].gtp.ps_pattern;
      break;
    case 16: detector = "TAGH   ";
      nhit    = triggers_enabled[ii].gtp.tagh_nhit;
      pattern = triggers_enabled[ii].gtp.tagh_pattern;
      break;
    case 32: detector = "TOF   ";
      nhit    = triggers_enabled[ii].gtp.tof_nhit;
      pattern = triggers_enabled[ii].gtp.tof_pattern;
      break;
    default: detector = "NONE ";
      cout << " Unknown detector ===== " << triggers_enabled[ii].type << endl;
      break;
    }
    
    cout  << detector  <<  setw(6) <<  
      triggers_enabled[ii].gtp.fcal << setw(9) <<  
      triggers_enabled[ii].gtp.bcal << setw(11) <<  
      triggers_enabled[ii].gtp.en_thr << setw(6) <<  
      nhit << setw(9) << 
      triggers_enabled[ii].bit << setw(12) <<  
      triggers_enabled[ii].gtp.fcal_min << setw(14) <<
      triggers_enabled[ii].gtp.fcal_max << setw(10) <<
      triggers_enabled[ii].gtp.bcal_min << setw(14) <<
      triggers_enabled[ii].gtp.bcal_max << setw(8) <<
      hex << uppercase << "0x" << pattern << nouppercase << dec <<  endl;

  }

  cout << endl << endl;

}



int DL1MCTrigger_factory::FindTriggers(DL1MCTrigger *trigger){
  
  int trig_found = 0;
  
  // Main production trigger  
  for(unsigned int ii = 0; ii < triggers_enabled.size(); ii++){
    
    if(triggers_enabled[ii].bit == 0){    // Main production trigger found
      
      int gtp_energy  = 0;
      int bcal_energy = 0;
      
      for(unsigned int samp = 0; samp < sample; samp++){
	
	int bcal_samp = samp - BCAL_OFFSET;
	
	if(bcal_samp < 0){
	  bcal_energy = 0;
	} else if(bcal_samp >= sample){ 
	  bcal_energy = 0;
	} else{ 
	  bcal_energy = bcal_gtp[bcal_samp];
	}


	gtp_energy = triggers_enabled[ii].gtp.fcal*fcal_gtp[samp] + 
	  triggers_enabled[ii].gtp.bcal*bcal_energy;
	
	if(gtp_energy >= triggers_enabled[ii].gtp.en_thr){

	  if(triggers_enabled[ii].gtp.fcal_min > 0) {     // FCAL > 0
	    if(fcal_gtp[samp] > triggers_enabled[ii].gtp.fcal_min){
	      trigger->trig_mask   = (trigger->trig_mask | 0x1);
	      trigger->trig_time[0]   = samp - 25;
	      trig_found++;
	    }
	  } else {
	    
	    trigger->trig_mask   = (trigger->trig_mask | 0x1);
	    trigger->trig_time[0]   = samp - 25; 
	    trig_found++;
	  }
	  
	  break;

	}  // Check energy threshold	
      }      
    }     // Trigger Bit 0
  

    if(triggers_enabled[ii].bit == 2){    //  BCAL trigger found
            
      int gtp_energy  = 0;
      int bcal_energy = 0;
      
      for(unsigned int samp = 0; samp < sample; samp++){
	
	int bcal_samp = samp - BCAL_OFFSET;
	
	if(bcal_samp < 0){
	  bcal_energy = 0;
	} else if(bcal_samp >= sample){ 
	  bcal_energy = 0;
	} else{ 
	  bcal_energy = bcal_gtp[bcal_samp];
	}
	      
	gtp_energy = triggers_enabled[ii].gtp.bcal*bcal_energy;
      
	if(gtp_energy >= triggers_enabled[ii].gtp.en_thr){
	  
	  trigger->trig_mask   = (trigger->trig_mask | 0x4);
	  trigger->trig_time[2]   = samp - 25; 
	  trig_found++;

	  break;
	  
	}  // Energy threshold			
      }
    }      // Trigger Bit  3

  }    // Loop over triggers found in the config file
  
  return trig_found;
  
}




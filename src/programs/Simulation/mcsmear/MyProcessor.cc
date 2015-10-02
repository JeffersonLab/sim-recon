// Author: David Lawrence  Sat Jan 29 09:37:37 EST 2011
//
//
// MyProcessor.cc
//

#include <iostream>
#include <cmath>
using namespace std;

#include <strings.h>

#include "MyProcessor.h"

#include <JANA/JEvent.h>

#include <HDDM/DEventSourceHDDM.h>
#include <TRACKING/DMCThrown.h>

extern void Smear(hddm_s::HDDM *record);
extern char *OUTFILENAME;

static pthread_mutex_t output_file_mutex;
static pthread_t output_file_mutex_last_owner;

#include <JANA/JCalibrationFile.h>
static JCalibration *jcalib=NULL;

void mcsmear_thread_HUP_sighandler(int sig)
{
   jerr<<" Caught HUP signal for thread 0x"<<hex<<pthread_self()<<dec<<" thread exiting..."<<endl;

   // We use output_file_mutex_owner to keep track (sort of)
   // of which thread has the mutex locked. This mutex is only
   // locked at the end of the evnt method. Once the lock is
   // obtained, this value is set to hold the id of the thread
   // that locked it. It may help in debugging to know the last
   // known thread to have locked the mutex when the signal
   // handler was called
   jerr<<endl;
   jerr<<" Last thread to lock output file mutex: 0x"<<hex<<pthread_self()<<dec<<endl;
   jerr<<" Attempting to unlock mutex to avoid deadlock." <<endl;
   jerr<<" However, the output file is likely corrupt at" <<endl;
   jerr<<" this point and the process should be restarted ..." <<endl;
   jerr<<endl;
   pthread_mutex_unlock(&output_file_mutex);
   pthread_exit(NULL);
}


//------------------------------------------------------------------
// init   -Open output file 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
   // open HDDM file
   ofs = new ofstream(OUTFILENAME);
   if (!ofs->is_open()){
      cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
      exit(-1);
   }
   fout = new hddm_s::ostream(*ofs);
   Nevents_written = 0;

   HDDM_USE_COMPRESSION = 2;
   gPARMS->SetDefaultParameter("HDDM:USE_COMPRESSION", HDDM_USE_COMPRESSION,
                          "Turn on/off compression of the output HDDM stream."
                          " \"0\"=no compression, \"1\"=bz2 compression, \"2\"=z compression (default)");
   HDDM_USE_INTEGRITY_CHECKS = true;
   gPARMS->SetDefaultParameter("HDDM:USE_INTEGRITY_CHECKS",
                                HDDM_USE_INTEGRITY_CHECKS,
                          "Turn on/off automatic integrity checking on the"
                          " output HDDM stream."
                          " Set to \"0\" to turn off (it's on by default)");

   // enable on-the-fly bzip2 compression on output stream
   if (HDDM_USE_COMPRESSION == 0) {
      jout << " HDDM compression disabled" << std::endl;
   } else if (HDDM_USE_COMPRESSION == 1) {
      jout << " Enabling bz2 compression of output HDDM file stream" 
           << std::endl;
      fout->setCompression(hddm_s::k_bz2_compression);
   } else {
      jout << " Enabling z compression of output HDDM file stream (default)" 
           << std::endl;
      fout->setCompression(hddm_s::k_z_compression);
   }

   // enable a CRC data integrity check at the end of each event record
   if (HDDM_USE_INTEGRITY_CHECKS) {
      jout << " Enabling CRC data integrity check in output HDDM file stream"
           << std::endl;
      fout->setIntegrityChecks(hddm_s::k_crc32_integrity);
   }
   else {
      jout << " HDDM integrity checks disabled" << std::endl;
   }

   // We set the mutex type to "ERRORCHECK" so that if the
   // signal handler is called, we can unlock the mutex
   // safely whether we have it locked or not.
   pthread_mutexattr_t attr;
   pthread_mutexattr_init(&attr);
   pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK);
   pthread_mutex_init(&output_file_mutex, NULL);
   
   // pthreads does not provide an "invalid" value for 
   // a pthread_t that we can initialize with. Furthermore,
   // the pthread_t may be a simple as an integer or as
   // a complicated structure. Hence, to make this portable
   // we clear it with bzero.
   bzero(&output_file_mutex_last_owner, sizeof(pthread_t));
   
   return NOERROR;
}

jerror_t MyProcessor::brun(JEventLoop *loop, int locRunNumber)
{
   DApplication* locDApp = dynamic_cast<DApplication*>(japp);
   jcalib = locDApp->GetJCalibration(locRunNumber);
   DGeometry *dgeom=locDApp->GetDGeometry(locRunNumber);

   // Make sure jcalib is set
   if(!jcalib){
     _DBG_<<"ERROR - jcalib not set!"<<endl;
     _DBG_<<"ERROR - Exiting ..."<<endl;
     return 0;
   }
   
   // get the TOF parameters
   {
     cout<<"get TOF/tof_parms parameters from calibDB"<<endl;
     map<string, double> tofparms;
     jcalib->Get("TOF/tof_parms", tofparms);
     TOF_SIGMA =  tofparms["TOF_SIGMA"];
     TOF_PHOTONS_PERMEV =  tofparms["TOF_PHOTONS_PERMEV"];
   }

   // get the BCAL parameters
   {
     cout<<"get BCAL/bcal_parms parameters from calibDB"<<endl;
     map<string, double> bcalparms;
     jcalib->Get("BCAL/bcal_parms", bcalparms);
     BCAL_DARKRATE_GHZ         =  bcalparms["BCAL_DARKRATE_GHZ"];
     BCAL_SIGMA_SIG_RELATIVE   = bcalparms["BCAL_SIGMA_SIG_RELATIVE"];
     BCAL_SIGMA_PED_RELATIVE   = bcalparms["BCAL_SIGMA_PED_RELATIVE"];
     BCAL_SIPM_GAIN_VARIATION   = bcalparms["BCAL_SIPM_GAIN_VARIATION"];
     BCAL_XTALK_FRACT         = bcalparms["BCAL_XTALK_FRACT"];
     BCAL_INTWINDOW_NS         = bcalparms["BCAL_INTWINDOW_NS"];
     BCAL_DEVICEPDE         = bcalparms["BCAL_DEVICEPDE"];
     BCAL_SAMPLING_FRACT      = bcalparms["BCAL_SAMPLING_FRACT"];
     BCAL_AVG_DARK_DIGI_VALS_PER_EVENT      = bcalparms["BCAL_AVG_DARK_DIGI_VALS_PER_EVENT"];
     BCAL_PHOTONSPERSIDEPERMEV_INFIBER = bcalparms["BCAL_PHOTONSPERSIDEPERMEV_INFIBER"];
     BCAL_SAMPLINGCOEFA = bcalparms["BCAL_SAMPLINGCOEFA"];
     BCAL_SAMPLINGCOEFB = bcalparms["BCAL_SAMPLINGCOEFB"];
     BCAL_TIMEDIFFCOEFA = bcalparms["BCAL_TIMEDIFFCOEFA"];
     BCAL_TIMEDIFFCOEFB = bcalparms["BCAL_TIMEDIFFCOEFB"];
     BCAL_TWO_HIT_RESOL = bcalparms["BCAL_TWO_HIT_RESOL"];
   }

   {
     cout<<"get BCAL/attenuation_parameters from calibDB"<<endl;
     vector< vector<double> > in_atten_parameters;
     jcalib->Get("BCAL/attenuation_parameters", in_atten_parameters);
     attenuation_parameters.clear();
     int channel = 0;
     for (int module=1; module<=BCAL_NUM_MODULES; module++) {
   	  for (int layer=1; layer<=BCAL_NUM_LAYERS; layer++) {
   		for (int sector=1; sector<=BCAL_NUM_SECTORS; sector++) {
			//int cell_id = GetCalibIndex(module,layer,sector);

			vector<double> new_params(3,0.);
			//attenuation_parameters[cell_id][0] = in_atten_parameters[channel][0];
			//attenuation_parameters[cell_id][1] = in_atten_parameters[channel][1];
			//attenuation_parameters[cell_id][2] = in_atten_parameters[channel][2];
			// hack to workaround odd CCDB behavior
			//attenuation_parameters[cell_id][0] = in_atten_parameters[channel][1];
			//attenuation_parameters[cell_id][1] = in_atten_parameters[channel][2];
			//attenuation_parameters[cell_id][2] = in_atten_parameters[channel][0];
			 
			new_params[0] = in_atten_parameters[channel][0];
			new_params[1] = in_atten_parameters[channel][1];
			new_params[2] = in_atten_parameters[channel][2];
			attenuation_parameters.push_back( new_params );

			channel++;
		}
	  }
     }
   }

   {
     cout<<"get BCAL/effective_velocities parameters from calibDB"<<endl;
     vector <double> effective_velocities_temp;
     jcalib->Get("BCAL/effective_velocities", effective_velocities_temp);
     for (unsigned int i = 0; i < effective_velocities_temp.size(); i++){
       effective_velocities.push_back(effective_velocities_temp.at(i));
     }
   }

   {
     cout<<"get BCAL/digi_scales parameters from calibDB"<<endl;
     map<string, double> bcaldigiscales;
     jcalib->Get("BCAL/digi_scales", bcaldigiscales);
     BCAL_NS_PER_ADC_COUNT = bcaldigiscales["BCAL_ADC_TSCALE"];
     BCAL_NS_PER_TDC_COUNT = bcaldigiscales["BCAL_TDC_SCALE"];
   }

   {
     cout<<"get BCAL/base_time_offset parameters from calibDB"<<endl;
     map<string, double> bcaltimeoffsets;
     jcalib->Get("BCAL/base_time_offset", bcaltimeoffsets);
     BCAL_BASE_TIME_OFFSET = bcaltimeoffsets["BCAL_BASE_TIME_OFFSET"];
     BCAL_TDC_BASE_TIME_OFFSET = bcaltimeoffsets["BCAL_TDC_BASE_TIME_OFFSET"];
   }

   {
     cout<<"get FCAL/fcal_parms parameters from calibDB"<<endl;
     map<string, double> fcalparms;
     jcalib->Get("FCAL/fcal_parms", fcalparms);
     if (FCAL_PHOT_STAT_COEF == 0.0)
       FCAL_PHOT_STAT_COEF   = fcalparms["FCAL_PHOT_STAT_COEF"]; 
     if (FCAL_BLOCK_THRESHOLD == 0.0)
       FCAL_BLOCK_THRESHOLD  = fcalparms["FCAL_BLOCK_THRESHOLD"];
   }
   {
     cout<<"get CDC/cdc_parms parameters from calibDB"<<endl;
     map<string, double> cdcparms;
     jcalib->Get("CDC/cdc_parms", cdcparms);
     if (CDC_TDRIFT_SIGMA == 0.0)
       CDC_TDRIFT_SIGMA   = cdcparms["CDC_TDRIFT_SIGMA"]; 
     if (CDC_TIME_WINDOW == 0.0)
       CDC_TIME_WINDOW    = cdcparms["CDC_TIME_WINDOW"];
     if (CDC_PEDESTAL_SIGMA == 0.0)
       CDC_PEDESTAL_SIGMA = cdcparms["CDC_PEDESTAL_SIGMA"]; 
     if (CDC_THRESHOLD_FACTOR == 0.0)
       CDC_THRESHOLD_FACTOR = cdcparms["CDC_THRESHOLD_FACTOR"];
   }

   {
     cout<<"get FDC/fdc_parms parameters from calibDB"<<endl;
     map<string, double> fdcparms;
     jcalib->Get("FDC/fdc_parms", fdcparms);

     if (FDC_TDRIFT_SIGMA == 0.0)
       FDC_TDRIFT_SIGMA      = fdcparms["FDC_TDRIFT_SIGMA"];
     if (FDC_CATHODE_SIGMA ==0.0)
       FDC_CATHODE_SIGMA     = fdcparms["FDC_CATHODE_SIGMA"];
     if (FDC_THRESHOLD_FACTOR == 0.0)
       FDC_THRESHOLD_FACTOR = fdcparms["FDC_THRESHOLD_FACTOR"];
     FDC_PED_NOISE         = fdcparms["FDC_PED_NOISE"];

     if (FDC_TIME_WINDOW == 0.0)
       FDC_TIME_WINDOW       = fdcparms["FDC_TIME_WINDOW"];

     if (FDC_HIT_DROP_FRACTION == 0.0)
       FDC_HIT_DROP_FRACTION = fdcparms["FDC_HIT_DROP_FRACTION"];  
     if (FDC_THRESH_KEV == 0.0)
       FDC_THRESH_KEV = fdcparms["FDC_THRESH_KEV"]; 
   }

   {
     cout<<"get START_COUNTER/start_parms parameters from calibDB"<<endl;
     map<string, double> startparms;
     jcalib->Get("START_COUNTER/start_parms", startparms);

     START_SIGMA = startparms["START_SIGMA"] ;
     START_PHOTONS_PERMEV = startparms["START_PHOTONS_PERMEV"];

   }

   // hist file
   fdc_drift_time_smear_hist=new TH2F("fdc_drift_time_smear_hist","Drift time smearing for FDC",
                  300,0.0,0.6,400,-200,200);
   fdc_drift_dist_smear_hist=new TH2F("fdc_drift_dist_smear_hist","Drift distance smearing for FDC",
                  100,0.0,0.6,400,-0.5,0.5);
   double tmax=TRIGGER_LOOKBACK_TIME+FDC_TIME_WINDOW;
   int num_time_bins=int(FDC_TIME_WINDOW);
   fdc_drift_time=new TH2F("fdc_drift_time","FDC drift distance vs. time",num_time_bins,TRIGGER_LOOKBACK_TIME,tmax,100,0,1.);
   
   fdc_anode_mult = new TH1F("fdc_anode_mult","wire hit multiplicity",20,-0.5,19.5);
   fdc_cathode_charge = new TH1F("fdc_cathode_charge","charge on strips",1000,0,1000);

   tmax=TRIGGER_LOOKBACK_TIME+CDC_TIME_WINDOW;
   num_time_bins=int(CDC_TIME_WINDOW);
   cdc_drift_time = new TH2F("cdc_drift_time","CDC drift distance vs time",num_time_bins,TRIGGER_LOOKBACK_TIME,tmax,80,0.,0.8);

   cdc_drift_smear = new TH2F("cdc_drift_smear","CDC drift smearing",
               100,0.0,800.0,100,-0.1,0.1);
   
   cdc_charge  = new TH1F("cdc_charge","Measured charge in straw",1000,-10e3,40e3);

   // Get number of cdc wires per ring and the radii of each ring
   vector<vector<DCDCWire *> >cdcwires;
   dgeom->GetCDCWires(cdcwires);
   for (unsigned int i=0;i<cdcwires.size();i++) {
      NCDC_STRAWS.push_back(cdcwires[i].size());
      CDC_RING_RADIUS.push_back(cdcwires[i][0]->origin.Perp());
   }  
   // Get the FDC z positions for the wire planes
   dgeom->GetFDCZ(FDC_LAYER_Z);

   // Coefficient used to calculate FDCsingle wire rate. We calculate
   // it once here just to save calculating it for every wire in every event
   FDC_RATE_COEFFICIENT = exp(-log(4.0)/23.0)/2.0/log(24.0)*FDC_TIME_WINDOW/1000.0E-9;
   
   // Something is a little off in my calculation above so I scale it down via
   // an emprical factor:
   FDC_RATE_COEFFICIENT *= 0.353;

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *loop, int eventnumber)
{
   JEvent& event = loop->GetJEvent();
   JEventSource *source = event.GetJEventSource();
   DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
   if (!hddm_source) {
      cerr << " This program MUST be used with an HDDM file as input!" << endl;
      exit(-1);
   }
   hddm_s::HDDM *record = (hddm_s::HDDM*)event.GetRef();
   if (!record)
      return NOERROR;
   
   // Smear values and add noise hits
   Smear(record);
   
   // Write event to output file
   pthread_mutex_lock(&output_file_mutex);
   output_file_mutex_last_owner = pthread_self();
   *fout << *record;
   Nevents_written++;
   pthread_mutex_unlock(&output_file_mutex);

   return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
   if (fout)
      delete fout;
   if (ofs) {
      ofs->close();
      cout << endl << "Closed HDDM file" << endl;
   }
   cout << " " << Nevents_written << " event written to " << OUTFILENAME
        << endl;
   
   return NOERROR;
}

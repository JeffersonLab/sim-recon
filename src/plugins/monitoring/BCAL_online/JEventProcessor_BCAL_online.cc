// $Id$
//
//    File: JEventProcessor_FCAL_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#include <stdint.h>
#include <vector>

#include "JEventProcessor_BCAL_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "BCAL/DBCALDigiHit.h"
#include "BCAL/DBCALTDCDigiHit.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALTDCHit.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALGeometry.h"

#include "DAQ/DF1TDCHit.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250WindowRawData.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile2D.h>
#include <TStyle.h>


// root hist pointers
static TH1I *bcal_fadc_digi_integral = NULL;
static TH1I *bcal_fadc_digi_pedestal = NULL;
static TH1I *bcal_fadc_digi_good_pedestal = NULL;
static TH1I *bcal_fadc_digi_QF = NULL;
static TH1I *bcal_fadc_digi_time = NULL;
static TH2I *bcal_fadc_digi_occ = NULL;
static TProfile2D *bcal_fadc_digi_pedestal_ave = NULL;
static TH1I *bcal_fadc_digi_nsamples_integral = NULL;
static TH1I *bcal_fadc_digi_nsamples_pedestal = NULL;
static TH1I *bcal_fadc_digi_occ_layer1 = NULL;
static TH1I *bcal_fadc_digi_occ_layer2 = NULL;
static TH1I *bcal_fadc_digi_occ_layer3 = NULL;
static TH1I *bcal_fadc_digi_occ_layer4 = NULL;
// static TH1I *bcal_fadc_digi_nhits_chan = NULL;
static TH1I *bcal_fadc_digi_nhits_evnt = NULL;

static TH1I *bcal_tdc_digi_time = NULL;
static TH2I *bcal_tdc_digi_reltime = NULL;
static TH2I *bcal_tdc_digi_occ = NULL;
static TH1I *bcal_tdc_digi_occ_layer1 = NULL;
static TH1I *bcal_tdc_digi_occ_layer2 = NULL;
static TH1I *bcal_tdc_digi_occ_layer3 = NULL;
// static TH1I *bcal_tdc_digi_nhits_chan = NULL;
static TH1I *bcal_tdc_digi_nhits_evnt = NULL;

static TH1I *bcal_fadc_E = NULL;
static TH1I *bcal_fadc_t = NULL;
static TH2I *bcal_fadc_occ = NULL;
// keep track of the integrated energy in each channel so that
// we can calculate the average energy per hit
static TH2D *bcal_fadc_avgE = NULL;
static TH2I *bcal_fadc_saturated = NULL;

static TH1I *bcal_tdc_t = NULL;
static TH2I *bcal_tdc_occ = NULL;

static TH1I *bcal_Uhit_E = NULL;
static TH1I *bcal_Uhit_t = NULL;
static TH1I *bcal_Uhit_t_ADC = NULL;
static TH1I *bcal_Uhit_t_TDC = NULL;
static TH1I *bcal_Uhit_tdiff = NULL;
static TH1I *bcal_Uhit_tTDC_twalk = NULL;
static TH1I *bcal_Uhit_noTDC_E = NULL;
static TH2I *bcal_Uhit_tTDC_tADC = NULL;
static TH2I *bcal_Uhit_tTDC_E = NULL;
static TH2I *bcal_Uhit_tADC_E = NULL;
static TProfile2D *bcal_Uhit_tdiff_ave = NULL;

static TH1I *bcal_point_E = NULL;
static TH1I *bcal_point_t = NULL;
static TH1I *bcal_point_rho = NULL;
static TH1I *bcal_point_sigRho = NULL;
static TH1I *bcal_point_theta = NULL;
static TH1I *bcal_point_sigTheta = NULL;
static TH1I *bcal_point_phi = NULL;
static TH1I *bcal_point_sigPhi = NULL;
static TH1I *bcal_point_z = NULL;
static TH1I *bcal_point_sigZ = NULL;
static TProfile2D *bcal_point_z_dist = NULL;
static TH2I *bcal_point_z_sector = NULL;
static TH2I *bcal_point_E_sector = NULL;
static TH1I *bcal_point_E_layer1 = NULL;
static TH1I *bcal_point_E_layer2 = NULL;
static TH1I *bcal_point_E_layer3 = NULL;
static TH1I *bcal_point_E_layer4 = NULL;
static TProfile *bcal_point_aveE_sector_layer1 = NULL;
static TProfile *bcal_point_aveE_sector_layer2 = NULL;
static TProfile *bcal_point_aveE_sector_layer3 = NULL;
static TProfile *bcal_point_aveE_sector_layer4 = NULL;


static TH1I *bcal_cluster_nCells = NULL;
static TH1I *bcal_cluster_E = NULL;
static TH1I *bcal_cluster_t = NULL;
static TH1I *bcal_cluster_sigT = NULL;
static TH1I *bcal_cluster_t0 = NULL;
static TH1I *bcal_cluster_rho = NULL;
static TH1I *bcal_cluster_sigRho = NULL;
static TH1I *bcal_cluster_theta = NULL;
static TH1I *bcal_cluster_sigTheta = NULL;
static TH1I *bcal_cluster_phi = NULL;
static TH1I *bcal_cluster_sigPhi = NULL;
static TH2I *bcal_cluster_rho_theta = NULL;

static TH1I *bcal_shower_N_cell = NULL;
static TH1I *bcal_shower_E = NULL;
static TH1I *bcal_shower_E_raw = NULL;
static TH1I *bcal_shower_x = NULL;
static TH1I *bcal_shower_y = NULL;
static TH1I *bcal_shower_z = NULL;
static TH1I *bcal_shower_t = NULL;
static TH1I *bcal_shower_xErr = NULL;
static TH1I *bcal_shower_yErr = NULL;
static TH1I *bcal_shower_zErr = NULL;
static TH1I *bcal_shower_tErr = NULL;
static TH2I *bcal_shower_plane = NULL;

static TH1I *bcal_num_events;

//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_BCAL_online());
	}
}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_online::JEventProcessor_BCAL_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_online::~JEventProcessor_BCAL_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_BCAL_online::init(void) {
	
	// lock all root operations
	japp->RootWriteLock();
	
	// First thread to get here makes all histograms. If one pointer is
	// already not NULL, assume all histograms are defined and return now
	if(bcal_fadc_digi_integral != NULL){
		japp->RootUnLock();
		return NOERROR;
	}
	
	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("bcal")->cd();
	//gStyle->SetOptStat(111110);

	// book hists
	int timemin_ns = -200;
	int timemax_ns = 400;
	int timeminTDC_ns = -300;
	int timemaxTDC_ns = 900;
	float Ehit_min = -0.2;
	float Ehit_max = 2.0;
	float Eclust_min = -0.2;
	float Eclust_max = 6.0;

	int numphibins = 192;
	float pi=3.14159265358;
	float phibinsize=2*pi/numphibins;
	float phimin=-phibinsize;
	float phimax=2*pi+phibinsize;
	numphibins+=2;

	gStyle->SetTitleOffset(1, "Y");
  	gStyle->SetTitleSize(0.05,"xyz");
	gStyle->SetTitleSize(0.08,"h");
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetTitleX(0);
	gStyle->SetTitleAlign(13);
	gStyle->SetNdivisions(505,"xy");

	gStyle->SetStatH(0.20);
	gStyle->SetStatW(0.30);
	gStyle->SetStatX(0.99);
	gStyle->SetStatY(0.99);

	bcal_num_events = new TH1I("bcal_num_events","BCAL Number of events",1, 0.5, 1.5);

	bcal_fadc_digi_integral = new TH1I("bcal_fadc_digi_integral","BCAL Integral (DBCALDigiHit);Integral (fADC counts)", 500, 0, 40000);
	bcal_fadc_digi_pedestal = new TH1I("bcal_fadc_digi_pedestal","BCAL Pedestal (DBCALDigiHit);Pedestal (fADC counts)", 150, -10, 140);
	bcal_fadc_digi_good_pedestal = new TH1I("bcal_fadc_digi_good_pedestal","BCAL Good Pedestal (DBCALDigiHit);Pedestal (fADC counts)", 150, -10, 140);
	bcal_fadc_digi_QF = new TH1I("bcal_fadc_digi_QF","Qualtiy Factor (DBCALDigiHit);Qualtiy Factor", 20, 0, 20);
	bcal_fadc_digi_time = new TH1I("bcal_fadc_digi_time","ADC Time (DBCALDigiHit);Time (fADC time/62.5 ps)", 410, -160, 6400);
	bcal_fadc_digi_occ = new TH2I("bcal_fadc_digi_occ","ADC occupancy (DBCALDigiHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	bcal_fadc_digi_pedestal_ave = new TProfile2D("bcal_fadc_digi_pedestal_ave",
						     "Mean pedestal per cell (DBCALDigiHit);Module", 
						     48, 0.5, 48.5, 33, 0.5, 33.5);
	bcal_fadc_digi_nsamples_integral = new TH1I("bcal_fadc_digi_nsamples_integral","Number of samples: Integral (DBCALDigiHit);Number of samples", 
						    100, 0, 100);
	bcal_fadc_digi_nsamples_pedestal = new TH1I("bcal_fadc_digi_nsamples_pedestal","Number of samples: Pedestal (DBCALDigiHit);Number of samples", 
						    10, 0, 10);
	bcal_fadc_digi_occ_layer1 = new TH1I("bcal_fadc_digi_occ_layer1","Occupancy in layer 1 (DBCALDigiHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_fadc_digi_occ_layer2 = new TH1I("bcal_fadc_digi_occ_layer2","Occupancy in layer 2 (DBCALDigiHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_fadc_digi_occ_layer3 = new TH1I("bcal_fadc_digi_occ_layer3","Occupancy in layer 3 (DBCALDigiHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_fadc_digi_occ_layer4 = new TH1I("bcal_fadc_digi_occ_layer4","Occupancy in layer 4 (DBCALDigiHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	// bcal_fadc_digi_nhits_chan = new TH1I("bcal_fadc_digi_nhits_chan","ADC hits per channel;hits per channel",5,-0.5,4.5);
	bcal_fadc_digi_nhits_evnt = new TH1I("bcal_fadc_digi_nhits_evnt","ADC hits per event;hits per event",125,-0.5,124.5);


	bcal_tdc_digi_time = new TH1I("bcal_tdc_digi_time","TDC Time (DBCALDigiTDCHit);Time (F1TDC counts)", 500, 0, 66000);
	bcal_tdc_digi_reltime = new TH2I("bcal_tdc_digi_reltime","Relative TDC Time (DBCALDigiTDCHit);Time (F1TDC counts); TDC trig time", 
					 100, 0, 70000, 100, 0, 600);
	bcal_tdc_digi_occ = new TH2I("bcal_tdc_digi_occ","TDC occupancy (DBCALDigiTDCHit);Module", 48, 0.5, 48.5, 25, 0.5, 25.5);

	bcal_tdc_digi_occ_layer1 = new TH1I("bcal_tdc_digi_occ_layer1","Occupancy in layer 1 (DBCALDigiTDCHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_tdc_digi_occ_layer2 = new TH1I("bcal_tdc_digi_occ_layer2","Occupancy in layer 2 (DBCALDigiTDCHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_tdc_digi_occ_layer3 = new TH1I("bcal_tdc_digi_occ_layer3","Occupancy in layer 3 (DBCALDigiTDCHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	// bcal_tdc_digi_nhits_chan = new TH1I("bcal_tdc_digi_nhits_chan","TDC hits per channel;hits per channel",5,-0.5,4.5);
	bcal_tdc_digi_nhits_evnt = new TH1I("bcal_tdc_digi_nhits_evnt","TDC hits per event;hits per event",125,-0.5,124.5);

	bcal_fadc_E = new TH1I("bcal_fadc_E","Uncorrected Energy (DBCALHit);Energy, not atten. corr. (GeV)", 500, Ehit_min, Ehit_max);
	bcal_fadc_t = new TH1I("bcal_fadc_t","ADC Time (DBCALHit);Time (ns)", 1200, timemin_ns, timemax_ns);
	bcal_fadc_occ = new TH2I("bcal_fadc_occ","ADC occupancy (DBCALHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	bcal_fadc_avgE = new TH2D("bcal_fadc_avgE","Average Energy x Occupancy (DBCALHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	bcal_fadc_saturated = new TH2I("bcal_fadc_saturated","Saturated Occupancy (DBCALHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	
	bcal_tdc_t = new TH1I("bcal_tdc_t","TDC Time (DBCALTDCHit);Time (ns)", 1200, timeminTDC_ns, timemaxTDC_ns);
	bcal_tdc_occ = new TH2I("bcal_tdc_occ","TDC occupancy (DBCALTDCHit);Module", 48, 0.5, 48.5, 25, 0.5, 25.5);

	bcal_Uhit_E = new TH1I("bcal_Uhit_E","Uncorrected Energy (DBCALUnifiedHit);Energy, not atten. corr. (GeV)", 500, Ehit_min, Ehit_max);
	bcal_Uhit_t = new TH1I("bcal_Uhit_t","Unified time (DBCALUnifiedHit);time  (ns)", 1200, timeminTDC_ns, timemaxTDC_ns);
	bcal_Uhit_t_TDC = new TH1I("bcal_Uhit_t_TDC","TDC time (DBCALUnifiedHit);Timewalk corrected TDC time (ns)", 1200, timeminTDC_ns, timemaxTDC_ns);
	bcal_Uhit_t_ADC = new TH1I("bcal_Uhit_t_ADC","ADC time (DBCALUnifiedHit);ADC time (ns)", 1000, timemin_ns, timemax_ns);
	bcal_Uhit_noTDC_E = new TH1I("bcal_Uhit_noTDC_E","Energy for no TDC  (DBCALUnifiedHit);Energy", 500, Ehit_min, Ehit_max);
	bcal_Uhit_tTDC_tADC = new TH2I("bcal_Uhit_tTDC_tADC","ADC vs TDC time (DBCALUnifiedHit);TDC time;ADC time (ns)", 
				     100, timemin_ns, timemax_ns, 100, timemin_ns, timemax_ns);
	bcal_Uhit_tdiff = new TH1I("bcal_Uhit_tdiff","time diff. (ADC-TDC)  (DBCALUnifiedHit);(TDC - ADC) time (ns)", 
				   200, -300, 300);
	bcal_Uhit_tTDC_E = new TH2I("bcal_Uhit_tTDC_E","TDC time vs Energy (DBCALUnifiedHit);Energy, not atten. corr. (GeV);Timewalk corrected TDC time", 
				 100, Ehit_min, Ehit_max, 100, timeminTDC_ns, timemaxTDC_ns);
	bcal_Uhit_tADC_E = new TH2I("bcal_Uhit_tADC_E","ADC time vs Energy (DBCALUnifiedHit);Energy, not atten. corr. (GeV);ADC time (ns)", 
				     100, Ehit_min, Ehit_max, 100, timemin_ns, timemax_ns);
	bcal_Uhit_tTDC_twalk = new TH1I("bcal_Uhit_tTDC_twalk","TDC timewalk correction (DBCALUnifiedHit);Timewalk correction  (ns)", 
				 150, -140, 10);
	bcal_Uhit_tdiff_ave = new TProfile2D("bcal_Uhit_tdiff_ave", "Mean time diff. (TDC-ADC) (DBCALDigiHit);Module", 
					     48, 0.5, 48.5, 33, 0.5, 33.5);

	bcal_point_E = new TH1I("bcal_point_E","Energy (DBCALPoint);Energy  (GeV)", 500, Ehit_min, Ehit_max);
	bcal_point_t = new TH1I("bcal_point_t","time (DBCALPoint);Time (ns)", 500, timemin_ns, timemax_ns);
	bcal_point_rho = new TH1I("bcal_point_rho","rho (DBCALPoint);rho", 450, 0, 450);
	bcal_point_sigRho = new TH1I("bcal_point_sigRho","sigRho (DBCALPoint)", 100, 0, 30);
	bcal_point_theta = new TH1I("bcal_point_theta","theta (DBCALPoint)", 500, 0, 3);
	bcal_point_sigTheta = new TH1I("bcal_point_sigTheta","sigTheta (DBCALPoint)", 100, 0, 1);
	bcal_point_phi = new TH1I("bcal_point_phi","phi (DBCALPoint)",numphibins,phimin,phimax);
	bcal_point_sigPhi = new TH1I("bcal_point_sigPhi","sigPhi (DBCALPoint)", 100, 0, 0.1);
	bcal_point_z = new TH1I("bcal_point_z","Z wrt target center (DBCALPoint);Z wrt target center (cm)", 600, -100, 500);
	bcal_point_sigZ = new TH1I("bcal_point_sigZ","sigZ (DBCALPoint)", 100, 0, 35);
	bcal_point_z_dist = new TProfile2D("bcal_point_z_dist","Mean Z positions per cell (DBCALPoint);Module;Sector, Layer", 
					   48, 0.5, 48.5, 16, 0.5, 16.5);
	bcal_point_z_sector = new TH2I("bcal_point_z_sector","Z vs global sector (DBCALPoint);global sector  (4 x module + sector)",
				       192, 0.5, 192.5, 300, -100, 500);
	bcal_point_E_sector = new TH2I("bcal_point_E_sector","Energy vs global sector (DBCALPoint);global sector  (4 x module + sector);Energy  (GeV)",
				       192, 0.5, 192.5, 200, Ehit_min, Ehit_max);
	bcal_point_E_layer1 = new TH1I("bcal_point_E_layer1","Energy in layer 1 (DBCALPoint);Energy  (GeV)", 225, Ehit_min, Ehit_max);
	bcal_point_E_layer2 = new TH1I("bcal_point_E_layer2","Energy in layer 2 (DBCALPoint);Energy  (GeV)", 225, Ehit_min, Ehit_max);
	bcal_point_E_layer3 = new TH1I("bcal_point_E_layer3","Energy in layer 3 (DBCALPoint);Energy  (GeV)", 225, Ehit_min, Ehit_max);
	bcal_point_E_layer4 = new TH1I("bcal_point_E_layer4","Energy in layer 4 (DBCALPoint);Energy  (GeV)", 225, Ehit_min, Ehit_max);
	bcal_point_aveE_sector_layer1 = new TProfile("bcal_point_aveE_sector_layer1",
						     "Mean energy in layer 1 (DBCALPoint);global sector  (4 x module + sector);Energy  (GeV)", 
						     192, 0.5, 192.5);
	bcal_point_aveE_sector_layer2 = new TProfile("bcal_point_aveE_sector_layer2",
						     "Mean energy in layer 2 (DBCALPoint);global sector  (4 x module + sector);Energy  (GeV)", 
						     192, 0.5, 192.5);
	bcal_point_aveE_sector_layer3 = new TProfile("bcal_point_aveE_sector_layer3",
						     "Mean energy in layer 3 (DBCALPoint);global sector  (4 x module + sector);Energy  (GeV)", 
						     192, 0.5, 192.5);
	bcal_point_aveE_sector_layer4 = new TProfile("bcal_point_aveE_sector_layer4",
						     "Mean energy in layer 4 (DBCALPoint);global sector  (4 x module + sector);Energy  (GeV)", 
						     192, 0.5, 192.5);

	bcal_cluster_nCells = new TH1I("bcal_cluster_nCells","Number of Cells (DBCALCluster)", 50, 0, 50);
	bcal_cluster_E = new TH1I("bcal_cluster_E","Energy (DBCALCluster);Energy  (GeV)", 450, Eclust_min, Eclust_max);
	bcal_cluster_t = new TH1I("bcal_cluster_t","time (DBCALCluster);Time (ns)", 500, timemin_ns, timemax_ns);
	bcal_cluster_sigT = new TH1I("bcal_cluster_sigT","sigT (DBCALCluster)", 100, 0, 50);
	bcal_cluster_t0 = new TH1I("bcal_cluster_t0","t0 (DBCALCluster);Inner raduis time (ns)", 500, timemin_ns, timemax_ns);
	bcal_cluster_rho = new TH1I("bcal_cluster_rho","rho (DBCALCluster)", 450, 0, 450);
	bcal_cluster_sigRho = new TH1I("bcal_cluster_sigRho","sigRho (DBCALCluster)", 100, 0, 50);
	bcal_cluster_theta = new TH1I("bcal_cluster_theta","theta (DBCALCluster)", 500, 0, 3);
	bcal_cluster_sigTheta = new TH1I("bcal_cluster_sigTheta","sigTheta (DBCALCluster)", 100, 0, 0.5);

	numphibins = 384;
	phibinsize=2*pi/numphibins;
	phimin=-phibinsize;
	phimax=2*pi+phibinsize;
	numphibins+=2;
	bcal_cluster_phi = new TH1I("bcal_cluster_phi","phi (DBCALCluster)",numphibins,phimin,phimax);
	bcal_cluster_sigPhi = new TH1I("bcal_cluster_sigPhi","sigPhi (DBCALCluster)", 100, 0, 0.1);
	bcal_cluster_rho_theta = new TH2I("bcal_cluster_rho_theta","theta vs rho (DBCALCluster);rho;theta", 100, 0, 360, 100, 0, 2.2);

	bcal_shower_N_cell = new TH1I("bcal_shower_N_cell","Number of Cells (DBCALShower)", 50, 0, 50);
	bcal_shower_E = new TH1I("bcal_shower_E","Energy (DBCALShower);Energy  (GeV)", 500, Eclust_min, Eclust_max);
	bcal_shower_E_raw = new TH1I("bcal_shower_E_raw","Raw Energy (DBCALShower)", 500, Eclust_min, Eclust_max);
	bcal_shower_x = new TH1I("bcal_shower_x","x (DBCALShower);X position  (cm)", 500, -100, 100);
	bcal_shower_y = new TH1I("bcal_shower_y","y (DBCALShower);Y position  (cm)", 500, -100, 100);
	bcal_shower_z = new TH1I("bcal_shower_z","z (DBCALShower);Z position  (cm)", 600, -100, 500);
	bcal_shower_t = new TH1I("bcal_shower_t","Time (DBCALShower);Time (ns)", 500, timemin_ns, timemax_ns);
	bcal_shower_xErr = new TH1I("bcal_shower_xErr","xErr (DBCALShower)", 100, 0, 30);
	bcal_shower_yErr = new TH1I("bcal_shower_yErr","yErr (DBCALShower)", 100, 0, 30);
	bcal_shower_zErr = new TH1I("bcal_shower_zErr","zErr (DBCALShower)", 100, 0, 40);
	bcal_shower_tErr = new TH1I("bcal_shower_tErr","tErr (DBCALShower)", 100, 0, 50);
	bcal_shower_plane = new TH2I("bcal_shower_plane","Shower position (DBCALShower);X position  (cm);Y position  (cm)", 100, -100, 100, 100, -100, 100);

	// Turn off stats window for occupancy plots
	bcal_fadc_digi_occ->SetStats(0);
	bcal_tdc_digi_occ->SetStats(0);
	bcal_tdc_digi_reltime->SetStats(0);
	bcal_fadc_occ->SetStats(0);
	bcal_fadc_avgE->SetStats(0);
	bcal_fadc_saturated->SetStats(0);
	bcal_tdc_occ->SetStats(0);
	bcal_Uhit_tTDC_tADC->SetStats(0);
	bcal_Uhit_tTDC_E->SetStats(0);
	bcal_Uhit_tADC_E->SetStats(0);
	bcal_Uhit_tTDC_twalk->SetStats(0);
	bcal_point_z_dist->SetStats(0);
	bcal_point_z_sector->SetStats(0);
	bcal_fadc_digi_pedestal_ave->SetStats(0);
	bcal_Uhit_tdiff_ave->SetStats(0);
	bcal_shower_plane->SetStats(0);
	bcal_cluster_rho_theta->SetStats(0);

	// Set y-axis labels for occupancy plots
	for(int ibin=1; ibin<=16; ibin++){
		int idy = ibin-1; // convenient to use index that starts from zero!
		
		int layer  = 1 + (idy%4);
		int sector = 1 + idy/4;
		
		stringstream ss;
		ss<<"D  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_fadc_digi_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_fadc_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_fadc_avgE->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_fadc_saturated->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_fadc_digi_pedestal_ave->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_Uhit_tdiff_ave->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());

		ss.str("");
		ss.clear();
		ss<<"U  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_fadc_digi_occ->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_fadc_occ->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_fadc_avgE->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_fadc_saturated->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_fadc_digi_pedestal_ave->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_Uhit_tdiff_ave->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
	}
	
	// Occupancy plots for TDC (without layer 4)
	// Set y-axis labels for occupancy plots
	for(int ibin=1; ibin<=12; ibin++){
		int idy = ibin-1; // convenient to use index that starts from zero!
		
		int layer  = 1 + (idy%3);
		int sector = 1 + idy/3;
		
		stringstream ss;
		ss<<"D  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_tdc_digi_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_tdc_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());

		ss.str("");
		ss.clear();
		ss<<"U  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_tdc_digi_occ->GetYaxis()->SetBinLabel(ibin+13, ss.str().c_str());
		bcal_tdc_occ->GetYaxis()->SetBinLabel(ibin+13, ss.str().c_str());
	}

	// Occupancy plots WITHOUT 2 ends
	// Set y-axis labels for occupancy plots
	for(int ibin=1; ibin<=16; ibin++){
		int idy = ibin-1; // convenient to use index that starts from zero!
		
		int layer  = 1 + (idy%4);
		int sector = 1 + idy/4;
		
		stringstream ss;
		ss<<"  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_point_z_dist->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_point_z_dist->GetYaxis()->SetTitleOffset(1.4);
	}

	
	// back to main dir
	main->cd();
	
	// unlock
	japp->RootUnLock();
	
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_online::brun(JEventLoop *eventLoop, int runnumber) {
	// This is called whenever the run number changes
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_online::evnt(JEventLoop *loop, int eventnumber) {
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop-Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	vector<const DBCALDigiHit*> dbcaldigihits;
	vector<const DBCALTDCDigiHit*> dbcaltdcdigihits;
	vector<const DBCALHit*> dbcalhits;
	vector<const DBCALTDCHit*> dbcaltdchits;
	vector<const DBCALUnifiedHit*> dbcaluhits;
	vector<const DBCALPoint*> dbcalpoints;
	vector<const DBCALCluster*> dbcalclusters;
	vector<const DBCALShower*> dbcalshowers;
	
	loop->Get(dbcaldigihits);
	loop->Get(dbcaltdcdigihits);
	loop->Get(dbcalhits);
	loop->Get(dbcaltdchits);
	loop->Get(dbcaluhits);
	loop->Get(dbcalpoints);
	loop->Get(dbcalclusters);
	loop->Get(dbcalshowers);
	
	// Lock ROOT
	japp->RootWriteLock();

	if( (dbcaldigihits.size() > 0) || (dbcaltdcdigihits.size() > 0) )
		bcal_num_events->Fill(1);

	bcal_fadc_digi_nhits_evnt->Fill(dbcaldigihits.size());
	bcal_tdc_digi_nhits_evnt->Fill(dbcaltdcdigihits.size());

	// The map is created to get the nunmber of hits per cell
	// map<readout_channel, cellHits> cellHitMap;
	// for( vector<const DBCALDigiHit*>::const_iterator hitPtr = dbcaldigihits.begin();
	// 	hitPtr != dbcaldigihits.end();
	// 	++hitPtr ){
	  
	// 	const DBCALDigiHit& hit = (**hitPtr);
	  
	// 	int id = DBCALGeometry::cellId( hit.module, hit.layer, hit.sector );
	// 	readout_channel chan(id, hit.end);
	  
	// 	//this will create cellHitMap[chan] if it doesn't already exist
	// 	cellHitMap[chan].hits.push_back(*hitPtr);
	// }


	// Digitized fADC hits for bcal
	for(unsigned int i=0; i<dbcaldigihits.size(); i++) {
		const DBCALDigiHit *hit = dbcaldigihits[i];
		
		bcal_fadc_digi_integral->Fill(hit->pulse_integral);
		bcal_fadc_digi_pedestal->Fill(hit->pedestal);
		if (hit->pedestal > 0) bcal_fadc_digi_good_pedestal->Fill(hit->pedestal);
		bcal_fadc_digi_QF->Fill(hit->QF);
		bcal_fadc_digi_time->Fill(hit->pulse_time);
		bcal_fadc_digi_nsamples_integral->Fill(hit->nsamples_integral);
		bcal_fadc_digi_nsamples_pedestal->Fill(hit->nsamples_pedestal);
		int layer = hit->layer;
		int glosect = DBCALGeometry::getglobalsector(hit->module, hit->sector);
		if (layer==1) bcal_fadc_digi_occ_layer1->Fill(glosect);
		if (layer==2) bcal_fadc_digi_occ_layer2->Fill(glosect);
		if (layer==3) bcal_fadc_digi_occ_layer3->Fill(glosect);
		if (layer==4) bcal_fadc_digi_occ_layer4->Fill(glosect);

		// Occupancy histogram defined to give better aspect ratio
		int ix = hit->module;
		int iy = (hit->sector-1)*4 + hit->layer;
		if(hit->end == DBCALGeometry::kUpstream) {
		  bcal_fadc_digi_occ->Fill(ix, iy+17);
		  if ( hit->pedestal > 0 )
		    bcal_fadc_digi_pedestal_ave->Fill(ix, iy+17, hit->pedestal);
		}
		if(hit->end == DBCALGeometry::kDownstream) {
		  bcal_fadc_digi_occ->Fill(ix, iy);
		  if ( hit->pedestal > 0 )
		    bcal_fadc_digi_pedestal_ave->Fill(ix, iy, hit->pedestal);
		}

	}
	
	// Digitized TDC hits for bcal
	for(unsigned int i=0; i<dbcaltdcdigihits.size(); i++) {
		const DBCALTDCDigiHit *hit = dbcaltdcdigihits[i];
		
		bcal_tdc_digi_time->Fill(hit->time);
		vector<const DF1TDCHit*> f1tdchits;
		hit->Get(f1tdchits);
		if (f1tdchits.size()>0) bcal_tdc_digi_reltime->Fill(f1tdchits[0]->time,f1tdchits[0]->trig_time);

		int layer = hit->layer;
		int glosect = DBCALGeometry::getglobalsector(hit->module, hit->sector);
		if (layer==1) bcal_tdc_digi_occ_layer1->Fill(glosect);
		if (layer==2) bcal_tdc_digi_occ_layer2->Fill(glosect);
		if (layer==3) bcal_tdc_digi_occ_layer3->Fill(glosect);

		// Occupancy histogram defined to give better aspect ratio
		int ix = hit->module;
		int iy = (hit->sector-1)*3 + hit->layer; // TDC has 3 layers per sector
		if(hit->end == DBCALGeometry::kUpstream) bcal_tdc_digi_occ->Fill(ix, iy+13);
		if(hit->end == DBCALGeometry::kDownstream) bcal_tdc_digi_occ->Fill(ix, iy);
	}

	// Calibrated fADC hits for bcal
	for(unsigned int i=0; i<dbcalhits.size(); i++) {
		const DBCALHit *hit = dbcalhits[i];

		bool saturated=0;
		vector<const DBCALDigiHit*> assoc_BCALDigiHit;
		hit->Get(assoc_BCALDigiHit);
		if (assoc_BCALDigiHit.size()>0) {
			vector<const Df250PulseIntegral*> f250PulseIntegral;
			assoc_BCALDigiHit[0]->Get(f250PulseIntegral);
			//printf("got %i DBCALDigiHit %i f250PulseIntegral\n", assoc_BCALDigiHit.size(),f250PulseIntegral.size());
			if (f250PulseIntegral.size()>0) {
				vector<const Df250WindowRawData*> f250WindowRawData;
				f250PulseIntegral[0]->Get(f250WindowRawData);
				if (f250WindowRawData.size()>0) {
					//printf("got Df250WindowRawData\n");
					if (f250WindowRawData[0]->overflow==1) {
						saturated=1;
					}
				}
			}
		}

		// Occupancy histogram defined to give better aspect ratio
		int ix = hit->module;
		int iy = (hit->sector-1)*4 + hit->layer;
		if (!saturated) {
			bcal_fadc_E->Fill(hit->E);
			bcal_fadc_t->Fill(hit->t);
			if(hit->end == DBCALGeometry::kUpstream) {
				bcal_fadc_occ->Fill(ix, iy+17);
				bcal_fadc_avgE->Fill(ix, iy+17, hit->E);
			}
			if(hit->end == DBCALGeometry::kDownstream) {
				bcal_fadc_occ->Fill(ix, iy);
				bcal_fadc_avgE->Fill(ix, iy, hit->E);
			}
		} else {
			if(hit->end == DBCALGeometry::kUpstream) {
				bcal_fadc_saturated->Fill(ix, iy+17);
			}
			if(hit->end == DBCALGeometry::kDownstream) {
				bcal_fadc_saturated->Fill(ix, iy);
			}
		}
	}
	
	// Calibrated TDC hits for bcal
	for(unsigned int i=0; i<dbcaltdchits.size(); i++) {
		const DBCALTDCHit *hit = dbcaltdchits[i];
		
		bcal_tdc_t->Fill(hit->t);
		
		// Occupancy histogram defined to give better aspect ratio
		int ix = hit->module;
		int iy = (hit->sector-1)*3 + hit->layer; // TDC has 3 layers per sector
		if(hit->end == DBCALGeometry::kUpstream) bcal_tdc_occ->Fill(ix, iy+13);
		if(hit->end == DBCALGeometry::kDownstream) bcal_tdc_occ->Fill(ix, iy);
	}

	// BCAL Unified Hit (combined ADC and TDC info)
	for(unsigned int i=0; i<dbcaluhits.size(); i++) {
		const DBCALUnifiedHit *Uhit = dbcaluhits[i];
		bcal_Uhit_E->Fill(Uhit->E);
		bcal_Uhit_t->Fill(Uhit->t);
		bcal_Uhit_t_ADC->Fill(Uhit->t_ADC);
		bcal_Uhit_tADC_E->Fill(Uhit->E,Uhit->t_ADC);
		if (Uhit->t_TDC != 0) {
			bcal_Uhit_t_TDC->Fill(Uhit->t_TDC);
			float t_diff = Uhit->t_TDC - Uhit->t_ADC;
			bcal_Uhit_tdiff->Fill(t_diff);
			bcal_Uhit_tTDC_E->Fill(Uhit->E,Uhit->t_TDC);
			bcal_Uhit_tTDC_tADC->Fill(Uhit->t_TDC,Uhit->t_ADC);
			vector<const DBCALTDCHit*> tdchits;
			Uhit->Get(tdchits);
			bcal_Uhit_tTDC_twalk->Fill(Uhit->t_TDC - tdchits[0]->t);

			int ix = Uhit->module;
			int iy = (Uhit->sector-1)*4 + Uhit->layer;
			if(Uhit->end == DBCALGeometry::kUpstream) {
			  bcal_Uhit_tdiff_ave->Fill(ix, iy+17, t_diff);
			}
			if(Uhit->end == DBCALGeometry::kDownstream) {
			  bcal_Uhit_tdiff_ave->Fill(ix, iy, t_diff);
			}
		} else {
			bcal_Uhit_noTDC_E->Fill(Uhit->E);
		}
	}
	
	// BCAL Points (combined hits from each end)
	for(unsigned int i=0; i<dbcalpoints.size(); i++) {
		const DBCALPoint *point = dbcalpoints[i];
		bcal_point_E->Fill(point->E());
		bcal_point_t->Fill(point->t());
		bcal_point_rho->Fill(point->rho());
		bcal_point_sigRho->Fill(point->sigRho());
		bcal_point_theta->Fill(point->theta());
		bcal_point_sigTheta->Fill(point->sigTheta());
		bcal_point_phi->Fill(point->phi());
		bcal_point_sigPhi->Fill(point->sigPhi());
		bcal_point_z->Fill(point->z());
		bcal_point_sigZ->Fill(point->sigZ());

		int layer = point->layer();
		if (layer==1) bcal_point_E_layer1->Fill(point->E());
		if (layer==2) bcal_point_E_layer2->Fill(point->E());
		if (layer==3) bcal_point_E_layer3->Fill(point->E());
		if (layer==4) bcal_point_E_layer4->Fill(point->E());

		vector<const DBCALUnifiedHit*> endhits;
		point->Get(endhits);

		int glosect = DBCALGeometry::getglobalsector(endhits[0]->module, endhits[0]->sector);
		bcal_point_z_sector->Fill(glosect,point->z());
		bcal_point_E_sector->Fill(glosect,point->E());
		if (layer==1) bcal_point_aveE_sector_layer1->Fill(glosect,point->E());
		if (layer==2) bcal_point_aveE_sector_layer2->Fill(glosect,point->E());
		if (layer==3) bcal_point_aveE_sector_layer3->Fill(glosect,point->E());
		if (layer==4) bcal_point_aveE_sector_layer4->Fill(glosect,point->E());

		int ix = endhits[0]->module;
		int iy = (endhits[0]->sector-1)*4 + endhits[0]->layer;
		bcal_point_z_dist->Fill(ix, iy, point->z());
	}

	// BCAL Clusters (combined neighboring points)
	for(unsigned int i=0; i<dbcalclusters.size(); i++) {
		const DBCALCluster *cluster = dbcalclusters[i];
		bcal_cluster_nCells->Fill(cluster->nCells());
		bcal_cluster_E->Fill(cluster->E());
		bcal_cluster_t->Fill(cluster->t());
		bcal_cluster_sigT->Fill(cluster->sigT());
		bcal_cluster_t0->Fill(cluster->t0());
		bcal_cluster_rho->Fill(cluster->rho());
		bcal_cluster_sigRho->Fill(cluster->sigRho());
		bcal_cluster_theta->Fill(cluster->theta());
		bcal_cluster_sigTheta->Fill(cluster->sigTheta());
		bcal_cluster_phi->Fill(cluster->phi());
		bcal_cluster_sigPhi->Fill(cluster->sigPhi());
		bcal_cluster_rho_theta->Fill(cluster->rho(),cluster->theta());
	}

	// BCAL Showers (combined cluster with track)
	for(unsigned int i=0; i<dbcalshowers.size(); i++) {
		const DBCALShower *shower = dbcalshowers[i];
		bcal_shower_N_cell->Fill(shower->N_cell);
		bcal_shower_E->Fill(shower->E);
		bcal_shower_E_raw->Fill(shower->E_raw);
		bcal_shower_x->Fill(shower->x);
		bcal_shower_y->Fill(shower->y);
		bcal_shower_z->Fill(shower->z);
		bcal_shower_t->Fill(shower->t);
		bcal_shower_xErr->Fill(shower->xErr);
		bcal_shower_yErr->Fill(shower->yErr);
		bcal_shower_zErr->Fill(shower->zErr);
		bcal_shower_tErr->Fill(shower->tErr);
		bcal_shower_plane->Fill(shower->x,shower->y);
	}


	// Unlock ROOT
	japp->RootUnLock();
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_online::erun(void) {
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_online::fini(void) {
	// Called before program exit after event processing is finished.
	return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

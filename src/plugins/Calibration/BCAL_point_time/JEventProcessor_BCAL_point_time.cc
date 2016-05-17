// $Id$
//
//    File: JEventProcessor_BCAL_point_time.cc
// Created: Fri Apr  8 12:59:18 EDT 2016
// Creator: dalton (on Linux gluon109.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#include <stdint.h>
#include <vector>

#include "JEventProcessor_BCAL_point_time.h"
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

using namespace std;
using namespace jana;

#include "BCAL/DBCALPoint.h"
#include "TRACKING/DMCThrown.h"
//#include <DAQ/DF1TDCHit.h>
//#include <DAQ/Df250PulseIntegral.h>
//#include <DAQ/JEventSource_EVIO.h>
//#include <TTAB/DTranslationTable.h>

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <TStyle.h>
#include <TROOT.h>

static const double degperrad = 180/3.14159265;

static const int nummodule=48;
static const int numlayer=4;
static const int numsector=4;
// root hist pointers
static TH2I *hitus_TimeVsZ_layer[numlayer];
static TH2I *hitds_TimeVsZ_layer[numlayer];
static TH2I *hit_TimediffVsZ_layer[numlayer];
static TH2I *hit_TimesumVsZ_layer[numlayer];
static TH2I *point_TimeVsZ_chan[nummodule][numlayer][numsector];
static TH2I *point_TimeVsZ_layer[numlayer];
static TH1I *point_NVsZ_layer[numlayer];
static TH1I *point_NVsTheta_layer[numlayer];
static TH1I *thrown_NVsZ;
static TH1I *thrown_NVsTheta;
static TH1F *point_NormVsZ_layer[numlayer];
static TH1F *point_NormVsTheta_layer[numlayer];

static TH2I *test_coords;

static TProfile *hitus_TimeVsZ_layer_prof[numlayer];
static TProfile *hitds_TimeVsZ_layer_prof[numlayer];
static TProfile *hit_TimediffVsZ_layer_prof[numlayer];
static TProfile *hit_TimesumVsZ_layer_prof[numlayer];
static TProfile *point_TimeVsZ_layer_prof[numlayer];

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_BCAL_point_time());
}
} // "C"


//------------------
// JEventProcessor_BCAL_point_time (Constructor)
//------------------
JEventProcessor_BCAL_point_time::JEventProcessor_BCAL_point_time()
{
	VERBOSE = 0;

	if(gPARMS){
		gPARMS->SetDefaultParameter("BCAL_point_time:VERBOSE", VERBOSE, "BCAL_point_time verbosity level");
	}
}

//------------------
// ~JEventProcessor_BCAL_point_time (Destructor)
//------------------
JEventProcessor_BCAL_point_time::~JEventProcessor_BCAL_point_time()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_point_time::init(void)
{
	if (VERBOSE>=1) printf("JEventProcessor_pedestal_online::init()\n");

	// lock all root operations
	japp->RootWriteLock();
		
	gStyle->SetTitleSize(0.06,"xyz");
	gStyle->SetTitleSize(0.07,"h");
	gStyle->SetLabelSize(0.06,"xyz");


	// create root folder for DAQ and cd to it, store main dir
	maindir = gDirectory;
	peddir = maindir->mkdir("bcalpointime");
	peddir->cd();

	int timebins = 200;
	float mintime = -5;
	float maxtime = 45;

	int zbins = 100;
	float minz = -80;
	float maxz = 350;

	int thetabins = 120;
	float mintheta = 10;
	float maxtheta = 130;

	int timediffbins = 150;
	float mintimediff = -30;
	float maxtimediff = 30;

	int timesumbins = 150;
	float mintimesum = 20;
	float maxtimesum = 60;

	test_coords = new TH2I("test_coords","test_coords",thetabins,mintheta,maxtheta,zbins,minz,maxz);

	char histname[255], histtitle[255];
	for (int j=0; j<numlayer; j++) {
		if (hitus_TimeVsZ_layer_prof[j] == NULL) {
			sprintf(histname,"hitus_TimeVsZ_layer%i_prof",j+1);
			sprintf(histtitle,"arrival time;Z  (cm);time  (ns)");
			hitus_TimeVsZ_layer_prof[j] = new TProfile(histname,histtitle,zbins,minz,maxz,mintime,maxtime);
			//int color = j+1; if (color>2) color++;
			hitus_TimeVsZ_layer_prof[j]->SetLineColor(j+1);
			hitus_TimeVsZ_layer_prof[j]->SetMarkerColor(j+1);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (hitds_TimeVsZ_layer_prof[j] == NULL) {
			sprintf(histname,"hitds_TimeVsZ_layer%i_prof",j+1);
			sprintf(histtitle,"arrival time;Z  (cm);time  (ns)");
			hitds_TimeVsZ_layer_prof[j] = new TProfile(histname,histtitle,zbins,minz,maxz,mintime,maxtime);
			hitds_TimeVsZ_layer_prof[j]->SetLineColor(j+1);
			hitds_TimeVsZ_layer_prof[j]->SetMarkerColor(j+1);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (hit_TimediffVsZ_layer_prof[j] == NULL) {
			sprintf(histname,"hit_TimediffVsZ_layer%i_prof",j+1);
			sprintf(histtitle,"end time difference;Z  (cm);#Delta t (ns)");
			hit_TimediffVsZ_layer_prof[j] = new TProfile(histname,histtitle,zbins,minz,maxz,mintimediff,maxtimediff);
			hit_TimediffVsZ_layer_prof[j]->SetLineColor(j+1);
			hit_TimediffVsZ_layer_prof[j]->SetMarkerColor(j+1);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (hit_TimesumVsZ_layer_prof[j] == NULL) {
			sprintf(histname,"hit_TimesumVsZ_layer%i_prof",j+1);
			sprintf(histtitle,"end time sum;Z  (cm);#Sigma t (ns)");
			hit_TimesumVsZ_layer_prof[j] = new TProfile(histname,histtitle,zbins,minz,maxz,mintimesum,maxtimesum);
			hit_TimesumVsZ_layer_prof[j]->SetLineColor(j+1);
			hit_TimesumVsZ_layer_prof[j]->SetMarkerColor(j+1);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (point_TimeVsZ_layer_prof[j] == NULL) {
			sprintf(histname,"point_TimeVsZ_layer%i_prof",j+1);
			sprintf(histtitle,"calculated point time;Z  (cm);time  (ns)");
			point_TimeVsZ_layer_prof[j] = new TProfile(histname,histtitle,zbins,minz,maxz,mintime,maxtime);
			point_TimeVsZ_layer_prof[j]->SetLineColor(j+1);
			point_TimeVsZ_layer_prof[j]->SetMarkerColor(j+1);
		} 
	}

	for (int j=0; j<numlayer; j++) {
		if (hitus_TimeVsZ_layer[j] == NULL) {
			sprintf(histname,"hitus_TimeVsZ_layer%i",j+1);
			sprintf(histtitle,"arrival time upstream [layer %i];Z  (cm);time  (ns)",j+1);
			hitus_TimeVsZ_layer[j] = new TH2I(histname,histtitle,260,-120,400,timebins,mintime,maxtime);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (hitds_TimeVsZ_layer[j] == NULL) {
			sprintf(histname,"hitds_TimeVsZ_layer%i",j+1);
			sprintf(histtitle,"arrival time downstream [layer %i];Z  (cm);time  (ns)",j+1);
			hitds_TimeVsZ_layer[j] = new TH2I(histname,histtitle,260,-120,400,timebins,mintime,maxtime);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (hit_TimediffVsZ_layer[j] == NULL) {
			sprintf(histname,"hit_TimediffVsZ_layer%i",j+1);
			sprintf(histtitle,"end time difference [layer %i];Z  (cm);#Delta t (ns)",j+1);
			hit_TimediffVsZ_layer[j] = new TH2I(histname,histtitle,260,-120,400,timediffbins,mintimediff,maxtimediff);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (hit_TimesumVsZ_layer[j] == NULL) {
			sprintf(histname,"hit_TimesumVsZ_layer%i",j+1);
			sprintf(histtitle,"end time sum [layer %i];Z  (cm);#Delta t (ns)",j+1);
			hit_TimesumVsZ_layer[j] = new TH2I(histname,histtitle,260,-120,400,timesumbins,mintimesum,maxtimesum);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (point_TimeVsZ_layer[j] == NULL) {
			sprintf(histname,"point_TimeVsZ_layer%i",j+1);
			sprintf(histtitle,"calculated point time [layer %i];Z  (cm);time  (ns)",j+1);
			point_TimeVsZ_layer[j] = new TH2I(histname,histtitle,260,-120,400,timebins,mintime,maxtime);
		} 
	}
	sprintf(histname,"thrown_NVsZ");
	sprintf(histtitle,"number of points;Z  (cm)");
	thrown_NVsZ = new TH1I(histname,histtitle,zbins,minz,maxz);
	sprintf(histname,"thrown_NVsTheta");
	sprintf(histtitle,"number of points;#theta_{p}  (deg)");
	thrown_NVsTheta  = new TH1I(histname,histtitle,thetabins,mintheta,maxtheta);

	for (int j=0; j<numlayer; j++) {
		if (point_NormVsZ_layer[j] == NULL) {
			sprintf(histname,"point_NormVsZ_layer%i",j+1);
			sprintf(histtitle,"number of points per thrown particle;Z  (cm);points/thrown");
			point_NormVsZ_layer[j] = new TH1F(histname,histtitle,zbins,minz,maxz);
			point_NormVsZ_layer[j]->SetLineColor(j+1);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (point_NormVsTheta_layer[j] == NULL) {
			sprintf(histname,"point_NormVsTheta_layer%i",j+1);
			sprintf(histtitle,"number of points per thrown particle;#theta_{p}  (cm);points/thrown");
			point_NormVsTheta_layer[j] = new TH1F(histname,histtitle,thetabins,mintheta,maxtheta);
			point_NormVsTheta_layer[j]->SetLineColor(j+1);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (point_NVsZ_layer[j] == NULL) {
			sprintf(histname,"point_NVsZ_layer%i",j+1);
			sprintf(histtitle,"number of points;Z  (cm)");
			point_NVsZ_layer[j] = new TH1I(histname,histtitle,zbins,minz,maxz);
			point_NVsZ_layer[j]->SetLineColor(j+1);
		} 
	}
	for (int j=0; j<numlayer; j++) {
		if (point_NVsTheta_layer[j] == NULL) {
			sprintf(histname,"point_NVsTheta_layer%i",j+1);
			sprintf(histtitle,"number of points;#theta_{p}  (deg)");
			point_NVsTheta_layer[j] = new TH1I(histname,histtitle,thetabins,mintheta,maxtheta);
			point_NVsTheta_layer[j]->SetLineColor(j+1);
		} 
	}

	// Initialise histograms and variables
	for (int i=0; i<nummodule; i++) {
		for (int j=0; j<numlayer; j++) {
			for (int k=0; k<numsector; k++) {
				sprintf(histname,"point_TimeVsZ_chan_%02i%02i%02i",i,j,k);
				sprintf(histtitle,"points [%2i,%2i,%2i];Z  (cm);time  (ns)",i,j,k);
				point_TimeVsZ_chan[i][j][k] = new TH2I(histname,histtitle,130,-120,400,timebins/2,mintime,maxtime);
			}
		}
	}
	
	// back to main dir
	maindir->cd();
	
	// unlock
	japp->RootUnLock();

	return NOERROR;

}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_point_time::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_point_time::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	vector<const DBCALPoint*> points;
	vector<const DBCALHit*> hits;
	vector<const DMCThrown*> thrown;

	loop->Get(points);
	loop->Get(thrown);

	
    // FILL HISTOGRAMS
    // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	if (peddir!=NULL) peddir->cd();

	float z_coord = 0;
	float theta_thrown = 0;
	unsigned int numthrown = thrown.size();
	if (VERBOSE>=2) printf("event %i, thrown %i\n",(int)eventnumber,numthrown);
	if (numthrown==1) {
		float pz = thrown[0]->pz();
		float pt = sqrt(thrown[0]->px()*thrown[0]->px() + thrown[0]->py()*thrown[0]->py());
		float z_p = pz/pt * DBCALGeometry::m_radius[0];
		theta_thrown = degperrad*atan2(pt,pz);
		//float z_targ = thrown[0]->z();
		z_coord = z_p;
		test_coords->Fill(theta_thrown,z_coord);
		thrown_NVsZ->Fill(z_coord);
		thrown_NVsTheta->Fill(theta_thrown);

	} else {
		printf("Event %6lu numthrown %2i   \n",eventnumber,numthrown);
		japp->RootUnLock();
		return NOERROR;
	}

	unsigned int numpoints = points.size();
	if (VERBOSE>=2) printf("event %i, points %i\n",(int)eventnumber,numpoints);
	// Loop over all  objects in this event
	for(unsigned int pointn=0; pointn<numpoints; pointn++){
		const DBCALPoint *point = points[pointn];
		int module = point->module();
		int layer = point->layer();
		int sector = point->sector();
		float time = point->t();
		//		float zpos = point->z();
		float zpos = z_coord;

		if (module > nummodule || layer > numlayer || sector > numsector) {
			printf ("(%i,%i,%i) is greater than (%i,%i,%i)\n",module,layer,sector,nummodule,numlayer,numsector);
			continue;
		}
		if (VERBOSE>=3) printf ("      point %2i  (%i,%i,%i)  %8.2f  %8.3f\n",pointn,module,layer,sector,zpos,time);
		point_TimeVsZ_chan[module-1][layer-1][sector-1]->Fill(zpos,time);
		point_TimeVsZ_layer[layer-1]->Fill(zpos,time);
		point_TimeVsZ_layer_prof[layer-1]->Fill(zpos,time);
		point_NVsZ_layer[layer-1]->Fill(zpos);
		point_NVsTheta_layer[layer-1]->Fill(theta_thrown);

		point->Get(hits);
		unsigned int numhits = hits.size();
		if (VERBOSE>=2) printf("event %i, hits %i\n",(int)eventnumber,numhits);

		float ustime=0;
		float dstime=0;
		// Loop over all hits 
		for(unsigned int hitn=0; hitn<numhits; hitn++){
			const DBCALHit *hit = hits[hitn];
			int module = hit->module;
			int layer = hit->layer;
			int sector = hit->sector;
			int end = hit->end;  // 0=US
			float hittime = hit->t;

			if (module > nummodule || layer > numlayer || sector > numsector) {
				printf ("(%i,%i,%i) is greater than (%i,%i,%i)\n",module,layer,sector,nummodule,numlayer,numsector);
				continue;
			}
			if (VERBOSE>=3) printf ("      hit %2i  (%i,%i,%i,%i)   %8.3f\n",hitn,module,layer,sector,end,hittime);
			if (end==0)	{
				hitus_TimeVsZ_layer[layer-1]->Fill(zpos,hittime);
				hitus_TimeVsZ_layer_prof[layer-1]->Fill(zpos,hittime);
				ustime = hittime;
			}
			if (end==1)	{
				hitds_TimeVsZ_layer[layer-1]->Fill(zpos,hittime);
				hitds_TimeVsZ_layer_prof[layer-1]->Fill(zpos,hittime);
				dstime = hittime;
			}
			if (ustime!=0 && dstime!=0) {
				float timediff = ustime - dstime;
				float timesum = ustime + dstime;
				hit_TimediffVsZ_layer[layer-1]->Fill(zpos,timediff);
				hit_TimediffVsZ_layer_prof[layer-1]->Fill(zpos,timediff);
				hit_TimesumVsZ_layer[layer-1]->Fill(zpos,timesum);
				hit_TimesumVsZ_layer_prof[layer-1]->Fill(zpos,timesum);
			}
		}
	}



	maindir->cd();

    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_point_time::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

    // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	
	int nbins = thrown_NVsZ->GetXaxis()->GetNbins();
	for (int j=0; j<numlayer; j++) {
		for (int i=1; i<=nbins; i++) {
			float points = point_NVsZ_layer[j]->GetBinContent(i);
			float thrown = thrown_NVsZ->GetBinContent(i);
			//printf("%i  %i  %i  %i  %i\n",nbins,i,j,points, thrown);
			if (thrown>0) point_NormVsZ_layer[j]->SetBinContent(i,points/thrown);
			 
		}
	}
	nbins = thrown_NVsTheta->GetXaxis()->GetNbins();
	for (int j=0; j<numlayer; j++) {
		for (int i=1; i<=nbins; i++) {
			float points = point_NVsTheta_layer[j]->GetBinContent(i);
			float thrown = thrown_NVsTheta->GetBinContent(i);
			if (thrown>0) point_NormVsTheta_layer[j]->SetBinContent(i,points/thrown);
		}
	}

    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_point_time::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}


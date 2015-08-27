// $Id$
//
//    File: JEventProcessor_BCAL_attenlength_gainratio.cc
// Created: Mon Aug 10 10:17:48 EDT 2015
// Creator: dalton (on Linux gluon02.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_BCAL_attenlength_gainratio.h"
using namespace jana;

#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALDigiHit.h"
#include "DANA/DApplication.h"

#include <TDirectory.h>
#include <TStyle.h>
#include <TF1.h>


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_BCAL_attenlength_gainratio());
}
} // "C"


//------------------
// JEventProcessor_BCAL_attenlength_gainratio (Constructor)
//------------------
JEventProcessor_BCAL_attenlength_gainratio::JEventProcessor_BCAL_attenlength_gainratio()
{
	VERBOSE = 0;

	if(gPARMS){
		gPARMS->SetDefaultParameter("BCAL_ALGR:VERBOSE", VERBOSE, "Verbosity level");
	}

}

//------------------
// ~JEventProcessor_BCAL_attenlength_gainratio (Destructor)
//------------------
JEventProcessor_BCAL_attenlength_gainratio::~JEventProcessor_BCAL_attenlength_gainratio()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_attenlength_gainratio::init(void)
{
	// lock all root operations
	japp->RootWriteLock();

	// Set style
	gStyle->SetTitleOffset(1, "Y");
	gStyle->SetTitleOffset(1.3, "z");
  	gStyle->SetTitleSize(0.06,"xyz");
	gStyle->SetTitleSize(0.07,"h");
	gStyle->SetLabelSize(0.06,"xy");
	gStyle->SetLabelSize(0.04,"z");
	gStyle->SetTitleX(0);
	gStyle->SetTitleAlign(13);
	gStyle->SetNdivisions(505,"xy");

	gStyle->SetOptStat(11);
	gStyle->SetOptFit(1);

	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;  // save current directory
	gDirectory->mkdir("bcalgainratio")->cd();
	char histname[255], modtitle[255], histtitle[255];

	sprintf(histtitle,"Atten. length from integ.;-2 / slope of [ln(A_{U}/A_{D}) vs position];BCAL cells");
	hist_attenlength = new TH1I("hist_attenlength",histtitle,100,200,900);
	sprintf(histtitle,"Gain ratio from integ.;G_{U}/G_{D} from [ln(A_{U}/A_{D}) vs position];BCAL cells");
	hist_gainratio = new TH1I("hist_gainratio",histtitle,100,0.3,1.7);

	sprintf(histtitle,"Atten. length from integ. (Uncertainty);fit uncertainty;BCAL cells");
	hist_attenlength_err = new TH1I("hist_attenlength_err",histtitle,100,0,0);
	sprintf(histtitle,"Gain ratio from integ. (Uncertainty);fit uncertainty;BCAL cells");
	hist_gainratio_err = new TH1I("hist_gainratio_err",histtitle,100,0,0);

	sprintf(histtitle,"Atten. length from integ. (Uncertainty);rel. fit uncertainty;BCAL cells");
	hist_attenlength_relerr = new TH1I("hist_attenlength_relerr",histtitle,100,0,0);
	sprintf(histtitle,"Gain ratio from integ. (Uncertainty);rel. fit uncertainty;BCAL cells");
	hist_gainratio_relerr = new TH1I("hist_gainratio_relerr",histtitle,100,0,0);

	sprintf(histtitle,"Atten. length from integ.;Module;Layer and Sector");
	hist2D_attenlength = new TH2F("hist2D_attenlength",histtitle,48,0.5,48.5,16,0.5,16.5);
	sprintf(histtitle,"Gain ratio from integ.;Module;Layer and Sector;G_{U}/G_{D}");
	hist2D_gainratio = new TH2F("hist2D_gainratio",histtitle,48,0.5,48.5,16,0.5,16.5);
	EvsZ_all = new TH2I("EvsZ_all","E vs Z;Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);
	EvsZ_layer[0] = new TH2I("EvsZ_layer1","E vs Z (layer 1);Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);
	EvsZ_layer[1] = new TH2I("EvsZ_layer2","E vs Z (layer 2);Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);
	EvsZ_layer[2] = new TH2I("EvsZ_layer3","E vs Z (layer 3);Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);
	EvsZ_layer[3] = new TH2I("EvsZ_layer4","E vs Z (layer 4);Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);



	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
 
	gDirectory->mkdir("channels")->cd();
	// Create histograms
	for (int module=0; module<nummodule; module++) {
		for (int layer=0; layer<numlayer; layer++) {
			for (int sector=0; sector<numsector; sector++) {
				// sprintf(histname,"logpeakratiovsZ_%02i%i%i",module+1,layer+1,sector+1);
				// sprintf(modtitle,"Channel (M%i,L%i,S%i)",module+1,layer+1,sector+1);
				// sprintf(histtitle,"%s;Position (cm) [from ADC time diff [US-DS]];log of pulse height ratio US/DS",modtitle);
				// logpeakratiovsZ[module][layer][sector] = new TH2I(histname,histtitle,500,-225.0,225.0,500,-3,3);
				sprintf(histname,"logintratiovsZ_%02i%i%i",module+1,layer+1,sector+1);
				sprintf(modtitle,"Channel (M%i,L%i,S%i)",module+1,layer+1,sector+1);
				//sprintf(histtitle,"%s;Position (cm) [from ADC time diff [US-DS]];log of pulse height ratio US/DS",modtitle);
				sprintf(histtitle,"%s;Z Position (cm);log of integral ratio US/DS",modtitle);
				logintratiovsZ[module][layer][sector] = new TH2I(histname,histtitle,500,-250.0,250.0,500,-3,3);
				sprintf(histname,"EvsZ_%02i%i%i",module+1,layer+1,sector+1);
				sprintf(modtitle,"Channel (M%i,L%i,S%i)",module+1,layer+1,sector+1);
				sprintf(histtitle,"%s;Z Position (cm);Energy",modtitle);
				EvsZ[module][layer][sector] = new TH2I(histname,histtitle,100,-250.0,250.0,200,0,0.2);
			}
		}
	}



	// back to main dir
	main->cd();

	japp->RootUnLock();
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_attenlength_gainratio::brun(JEventLoop *eventLoop, int runnumber)
{
	// This is called whenever the run number changes

	DApplication* app = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
	DGeometry* geom = app->GetDGeometry(runnumber);
	geom->GetTargetZ(z_target_center);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_attenlength_gainratio::evnt(JEventLoop *loop, int eventnumber)
{

	// Start with matched points
	vector<const DBCALPoint*> dbcalpoints;
	loop->Get(dbcalpoints);

	japp->RootWriteLock();

	for(unsigned int i=0; i<dbcalpoints.size(); i++) {
		const DBCALPoint *point = dbcalpoints[i];
		int module = point->module();
		int layer = point->layer();
		int sector = point->sector();
		float Energy = point->E();

		// get the associated digi hits
		vector<const DBCALDigiHit*> digihits;
		point->Get(digihits);
		if (digihits.size()!=2) {
			printf("Warning: BCAL_attenlength_gainratio: event %i: wrong number of BCALDigiHit objects found %i\n",
				   eventnumber,(int)digihits.size());
			continue;
		}
		if (digihits[0]->end==digihits[1]->end) {
			printf("Warning: BCAL_attenlength_gainratio: event %i: two hits in same end of point\n",eventnumber);
			continue;
		}
		float integralUS, integralDS;
		// end 0=upstream, 1=downstream
		if (digihits[0]->end==0) {
			integralUS = digihits[0]->pulse_integral - ((float)digihits[0]->nsamples_integral*digihits[0]->pedestal)/
				digihits[0]->nsamples_pedestal;
			integralDS = digihits[1]->pulse_integral - ((float)digihits[1]->nsamples_integral*digihits[1]->pedestal)/
				digihits[1]->nsamples_pedestal;
		} else { 
			integralDS = digihits[0]->pulse_integral - ((float)digihits[0]->nsamples_integral*digihits[0]->pedestal)/
				digihits[0]->nsamples_pedestal;
			integralUS = digihits[1]->pulse_integral - ((float)digihits[1]->nsamples_integral*digihits[1]->pedestal)/
				digihits[1]->nsamples_pedestal;
		}

		//float timediff = t_ADCus_vec[0]-t_ADCds_vec[0];
		//float zpos = (timediff)*17./2;
		float zpos = point->z() - DBCALGeometry::GLOBAL_CENTER + z_target_center;
		float intratio = (float)integralUS/(float)integralDS;
		float logintratio = log(intratio);
		if (VERBOSE>4) printf("%5i  %2i %i %i  %8.1f  %8.1f  %8.3f  %8.3f  %8.3f\n", 
							  eventnumber,module,layer,sector,integralUS,integralDS,intratio,logintratio,zpos);
		if (Energy > 0.01) {  // 10 MeV cut to remove bias due to attenuation
			logintratiovsZ[module-1][layer-1][sector-1]->Fill(zpos, logintratio);
		} 
		EvsZ[module-1][layer-1][sector-1]->Fill(zpos, Energy);
		EvsZ_all->Fill(zpos, Energy);
		EvsZ_layer[layer-1]->Fill(zpos, Energy);
	}
	japp->RootUnLock();


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_attenlength_gainratio::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_attenlength_gainratio::fini(void)
{
	// Called before program exit after event processing is finished.


	
	// for (int module=0; module<nummodule; module++) {
	// 	if (VERBOSE>0) printf("M%i ",module);
	// 	for (int layer=0; layer<numlayer; layer++) {
	// 		for (int sector=0; sector<numsector; sector++) {
	// 			int entries = logintratiovsZ[module][layer][sector]->GetEntries();
	// 			if (VERBOSE>0) printf("(%i,%i) %3i  ", layer,sector,entries);
	// 		}
	// 	}
	// 	if (VERBOSE>0) printf("\n");
	// }

	printf("\n");
	for (int module=0; module<nummodule; module++) {
		for (int layer=0; layer<numlayer; layer++) {
			for (int sector=0; sector<numsector; sector++) {
				int layersect = (layer)*4 + sector + 1;
				int entries = logintratiovsZ[module][layer][sector]->GetEntries();
				if (entries>10) {
					logintratiovsZ[module][layer][sector]->Fit("pol1","q");
					TF1 *myfit = (TF1*)logintratiovsZ[module][layer][sector]->GetFunction("pol1");
					float p0 = myfit->GetParameter(0);
					float p1 = myfit->GetParameter(1);
					float p0err = myfit->GetParError(0);
					float p1err = myfit->GetParError(1);
					
					float attenlength = -2./p1;
					float gainratio = exp(p0);
					float attenlengtherr = 2/p1/p1*p1err;
					float gainratioerr = exp(p0)*p0err;
					if (VERBOSE>0) printf("(%2i,%i,%i) %3i %8.3f %8.3f   ", module, layer,sector,entries,attenlength,gainratio);
					hist_attenlength->Fill(attenlength);
					hist_gainratio->Fill(gainratio);
					hist_attenlength_err->Fill(attenlengtherr);
					hist_gainratio_err->Fill(gainratioerr);
					hist_attenlength_relerr->Fill(attenlengtherr/attenlength);
					hist_gainratio_relerr->Fill(gainratioerr/gainratio);
					hist2D_attenlength->SetBinContent(module+1,layersect,attenlength);
					hist2D_gainratio->SetBinContent(module+1,layersect,gainratio);
				}
			}
		}
	}
	return NOERROR;
}


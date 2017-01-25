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

	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	if (logintratiovsZ_all != NULL){
		return NOERROR;
	}

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
	TDirectory *bcalgainratio = main->mkdir("bcalgainratio");
    bcalgainratio->cd();

	char histname[255], modtitle[255], histtitle[255];

	sprintf(histtitle,"All channels;Z Position (cm);log of integral ratio US/DS");
	logintratiovsZ_all = new TH2I("logintratiovsZ_all",histtitle,500,-250.0,250.0,500,-3,3);
	sprintf(histtitle,"All channels;Z Position (cm);log of peak ratio US/DS");
	logpeakratiovsZ_all = new TH2I("logpeakratiovsZ_all",histtitle,500,-250.0,250.0,500,-3,3);

	sprintf(histtitle,"Atten. length from integ.;-2 / slope of [ln(A_{U}/A_{D}) vs position];BCAL cells");
	hist_attenlength = new TH1I("hist_attenlength",histtitle,100,200,900);
	sprintf(histtitle,"Gain ratio from integ.;G_{U}/G_{D} from [ln(A_{U}/A_{D}) vs position];BCAL cells");
	hist_gainratio = new TH1I("hist_gainratio",histtitle,100,0.3,1.7);

	sprintf(histtitle,"Atten. length from integ. (Uncertainty);fit uncertainty;BCAL cells");
	hist_attenlength_err = new TH1I("hist_attenlength_err",histtitle,100,0,0);
	sprintf(histtitle,"Gain ratio from integ. (Uncertainty);fit uncertainty;BCAL cells");
	hist_gainratio_err = new TH1I("hist_gainratio_err",histtitle,100,0,0);

	sprintf(histtitle,"Atten. length from integ. (Uncertainty);rel. fit uncertainty  (%%);BCAL cells");
	hist_attenlength_relerr = new TH1I("hist_attenlength_relerr",histtitle,100,0,0);
	sprintf(histtitle,"Gain ratio from integ. (Uncertainty);rel. fit uncertainty  (%%);BCAL cells");
	hist_gainratio_relerr = new TH1I("hist_gainratio_relerr",histtitle,100,0,0);

	sprintf(histtitle,"Atten. length from integ.;Module;Layer and Sector");
	hist2D_intattenlength = new TH2F("hist2D_intattenlength",histtitle,48,0.5,48.5,16,0.5,16.5);
	sprintf(histtitle,"Gain ratio from integ.;Module;Layer and Sector;G_{U}/G_{D}");
	hist2D_intgainratio = new TH2F("hist2D_intgainratio",histtitle,48,0.5,48.5,16,0.5,16.5);

	sprintf(histtitle,"Atten. length from peak;Module;Layer and Sector");
	hist2D_peakattenlength = new TH2F("hist2D_peakattenlength",histtitle,48,0.5,48.5,16,0.5,16.5);
	sprintf(histtitle,"Gain ratio from peak;Module;Layer and Sector;G_{U}/G_{D}");
	hist2D_peakgainratio = new TH2F("hist2D_peakgainratio",histtitle,48,0.5,48.5,16,0.5,16.5);


	EvsZ_all = new TH2I("EvsZ_all","E vs Z;Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);
	EvsZ_layer[0] = new TH2I("EvsZ_layer1","E vs Z (layer 1);Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);
	EvsZ_layer[1] = new TH2I("EvsZ_layer2","E vs Z (layer 2);Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);
	EvsZ_layer[2] = new TH2I("EvsZ_layer3","E vs Z (layer 3);Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);
	EvsZ_layer[3] = new TH2I("EvsZ_layer4","E vs Z (layer 4);Z Position (cm);Energy",100,-250.0,250.0,200,0,0.2);

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TDirectory *dirlogpeakratiovsZ = bcalgainratio->mkdir("logpeakratiovsZ");
	TDirectory *dirlogintratiovsZ = bcalgainratio->mkdir("logintratiovsZ");
	TDirectory *dirEvsZ = bcalgainratio->mkdir("EvsZ");

	// Create histograms
    dirlogpeakratiovsZ->cd();
	for (int module=0; module<nummodule; module++) {
		for (int layer=0; layer<numlayer; layer++) {
			for (int sector=0; sector<numsector; sector++) {
				sprintf(histname,"logpeakratiovsZ_%02i%i%i",module+1,layer+1,sector+1);
				sprintf(modtitle,"Channel (M%i,L%i,S%i)",module+1,layer+1,sector+1);
				sprintf(histtitle,"%s;Z Position (cm);log of pulse height ratio US/DS",modtitle);
				logpeakratiovsZ[module][layer][sector] = new TH2I(histname,histtitle,500,-225.0,225.0,500,-3,3);
			}
		}
	}
    dirlogintratiovsZ->cd();
	for (int module=0; module<nummodule; module++) {
		for (int layer=0; layer<numlayer; layer++) {
			for (int sector=0; sector<numsector; sector++) {
				sprintf(histname,"logintratiovsZ_%02i%i%i",module+1,layer+1,sector+1);
				sprintf(modtitle,"Channel (M%i,L%i,S%i)",module+1,layer+1,sector+1);
				sprintf(histtitle,"%s;Z Position (cm);log of integral ratio US/DS",modtitle);
				logintratiovsZ[module][layer][sector] = new TH2I(histname,histtitle,500,-250.0,250.0,500,-3,3);
			}
		}
	}
    dirEvsZ->cd();
	for (int module=0; module<nummodule; module++) {
		for (int layer=0; layer<numlayer; layer++) {
			for (int sector=0; sector<numsector; sector++) {
				sprintf(histname,"EvsZ_%02i%i%i",module+1,layer+1,sector+1);
				sprintf(modtitle,"Channel (M%i,L%i,S%i)",module+1,layer+1,sector+1);
				sprintf(histtitle,"%s;Z Position (cm);Energy",modtitle);
				EvsZ[module][layer][sector] = new TH2I(histname,histtitle,100,-250.0,250.0,200,0,0.2);
			}
		}
	}

	// back to main dir
	main->cd();

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_attenlength_gainratio::brun(JEventLoop *eventLoop, int32_t runnumber)
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
jerror_t JEventProcessor_BCAL_attenlength_gainratio::evnt(JEventLoop *loop, uint64_t eventnumber)
{

	// Start with matched points
	vector<const DBCALPoint*> dbcalpoints;
	loop->Get(dbcalpoints);

    // pull out the associated digihits
    vector< vector<const DBCALDigiHit*> > digihits_vec;
	for(unsigned int i=0; i<dbcalpoints.size(); i++) {
        const DBCALPoint *point = dbcalpoints[i];
        vector<const DBCALDigiHit*> digihits;
        point->Get(digihits);
        
        digihits_vec.push_back(digihits);
    }

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	for(unsigned int i=0; i<dbcalpoints.size(); i++) {
		const DBCALPoint *point = dbcalpoints[i];
		int module = point->module();
		int layer = point->layer();
		int sector = point->sector();
		float Energy = point->E();

		// get the associated digi hits
		vector<const DBCALDigiHit*> &digihits = digihits_vec[i];
		//point->Get(digihits);
		if (digihits.size()!=2) {
			printf("Warning: BCAL_attenlength_gainratio: event %llu: wrong number of BCALDigiHit objects found %i\n",
				   (long long unsigned int)eventnumber,(int)digihits.size());
			continue;
		}
		if (digihits[0]->end==digihits[1]->end) {
			printf("Warning: BCAL_attenlength_gainratio: event %llu: two hits in same end of point\n",(long long unsigned int)eventnumber);
			continue;
		}
		int Vmid0 = (digihits[0]->pulse_peak+digihits[0]->pedestal)/2;
		int Vmid1 = (digihits[1]->pulse_peak+digihits[1]->pedestal)/2;
		if (Vmid0 <= 105 || Vmid1 <= 105 || digihits[0]->pulse_time > 2880 || digihits[1]->pulse_time > 2880)  { // 2880 = 45 samples x 64 subsamples
			continue;
		}

		float peakUS, peakDS;
		float integralUS, integralDS;
		// end 0=upstream, 1=downstream
		if (digihits[0]->end==0) {
			integralUS = digihits[0]->pulse_integral - ((float)digihits[0]->nsamples_integral*digihits[0]->pedestal)/
				digihits[0]->nsamples_pedestal;
			integralDS = digihits[1]->pulse_integral - ((float)digihits[1]->nsamples_integral*digihits[1]->pedestal)/
				digihits[1]->nsamples_pedestal;
			peakUS = digihits[0]->pulse_peak - digihits[0]->pedestal;
			peakDS = digihits[1]->pulse_peak - digihits[1]->pedestal;
		} else { 
			integralDS = digihits[0]->pulse_integral - ((float)digihits[0]->nsamples_integral*digihits[0]->pedestal)/
				digihits[0]->nsamples_pedestal;
			integralUS = digihits[1]->pulse_integral - ((float)digihits[1]->nsamples_integral*digihits[1]->pedestal)/
				digihits[1]->nsamples_pedestal;
			peakDS = digihits[0]->pulse_peak - digihits[0]->pedestal;
			peakUS = digihits[1]->pulse_peak - digihits[1]->pedestal;
		}

		//float timediff = t_ADCus_vec[0]-t_ADCds_vec[0];
		//float zpos = (timediff)*17./2;
		float zpos = point->z() - DBCALGeometry::GetBCAL_center() + z_target_center;
		float intratio = (float)integralUS/(float)integralDS;
		float logintratio = log(intratio);
		float peakratio = (float)peakUS/(float)peakDS;
		float logpeakratio = log(peakratio);
		if (VERBOSE>4) printf("%5llu  %2i %i %i  %8.1f  %8.1f  %8.3f  %8.3f  %8.3f\n", 
							  (long long unsigned int)eventnumber,module,layer,sector,integralUS,integralDS,intratio,logintratio,zpos);
		if (Energy > 0.01) {  // 10 MeV cut to remove bias due to attenuation
			logintratiovsZ[module-1][layer-1][sector-1]->Fill(zpos, logintratio);
			logintratiovsZ_all->Fill(zpos, logintratio);
			logpeakratiovsZ[module-1][layer-1][sector-1]->Fill(zpos, logpeakratio);
			logpeakratiovsZ_all->Fill(zpos, logpeakratio);
		} 
		EvsZ[module-1][layer-1][sector-1]->Fill(zpos, Energy);
		EvsZ_all->Fill(zpos, Energy);
		EvsZ_layer[layer-1]->Fill(zpos, Energy);
	}

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK


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


	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

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
	printf("BCAL_attenlength_gainratio::fini >>  Fitting all histograms\n");

	for (int module=0; module<nummodule; module++) {
		for (int layer=0; layer<numlayer; layer++) {
			for (int sector=0; sector<numsector; sector++) {
				int layersect = (layer)*4 + sector + 1;
				int entries = logintratiovsZ[module][layer][sector]->GetEntries();
				if (entries>10) {
					logintratiovsZ[module][layer][sector]->Fit("pol1","q");
					TF1 *intfit = (TF1*)logintratiovsZ[module][layer][sector]->GetFunction("pol1");
					float p0 = intfit->GetParameter(0);
					float p1 = intfit->GetParameter(1);
					float p0err = intfit->GetParError(0);
					float p1err = intfit->GetParError(1);
					
					float attenlength = -2./p1;
					float gainratio = exp(p0);
					float attenlengtherr = 2/p1/p1*p1err;
					float gainratioerr = exp(p0)*p0err;
					if (VERBOSE>0) printf("(%2i,%i,%i) %3i %8.3f %8.3f   ", module, layer,sector,entries,attenlength,gainratio);
					hist_attenlength->Fill(attenlength);
					hist_gainratio->Fill(gainratio);
					hist_attenlength_err->Fill(attenlengtherr);
					hist_gainratio_err->Fill(gainratioerr);
					hist_attenlength_relerr->Fill(attenlengtherr/attenlength*100);
					hist_gainratio_relerr->Fill(gainratioerr/gainratio*100);
					hist2D_intattenlength->SetBinContent(module+1,layersect,attenlength);
					hist2D_intattenlength->SetBinError(module+1,layersect,attenlengtherr);
					hist2D_intgainratio->SetBinContent(module+1,layersect,gainratio);
					hist2D_intgainratio->SetBinError(module+1,layersect,gainratioerr);

					logpeakratiovsZ[module][layer][sector]->Fit("pol1","q");
					TF1 *peakfit = (TF1*)logpeakratiovsZ[module][layer][sector]->GetFunction("pol1");
					p0 = peakfit->GetParameter(0);
					p1 = peakfit->GetParameter(1);
					p0err = peakfit->GetParError(0);
					p1err = peakfit->GetParError(1);
					attenlength = -2./p1;
					gainratio = exp(p0);
					attenlengtherr = 2/p1/p1*p1err;
					gainratioerr = exp(p0)*p0err;
					hist2D_peakattenlength->SetBinContent(module+1,layersect,attenlength);
					hist2D_peakattenlength->SetBinError(module+1,layersect,attenlengtherr);
					hist2D_peakgainratio->SetBinContent(module+1,layersect,gainratio);
					hist2D_peakgainratio->SetBinError(module+1,layersect,gainratioerr);

				}
			}
		}
	}
	hist2D_peakattenlength->SetBinContent(0,0,1);
	hist2D_peakgainratio->SetBinContent(0,0,1);
	hist2D_intattenlength->SetBinContent(0,0,1);
	hist2D_intgainratio->SetBinContent(0,0,1);

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	return NOERROR;
}


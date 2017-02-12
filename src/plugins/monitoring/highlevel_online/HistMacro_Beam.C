// hnamepath: /highlevel/RFBeamBunchPeriod
// hnamepath: /highlevel/RFBeamBunchPeriod_DFT
// hnamepath: /highlevel/BeamEnergy
//
// e-mail: davidl@jlab.org
// e-mail: pmatt@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
//

{
	// This taken from the bin contents of the BeamEnergy histogram from the
	// online monitoring root file for run 30437, an amorphous target run.
	// The macro used is in /gluonwork1/Users/davidl/2017.02.08.amorphous_norm
	Double_t amorphous_data[] = {
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,    134126.0,    139064.0, 
		   126794.0,         0.0,    121863.0,    115792.0,    122810.0,    117894.0,    116636.0,    115590.0,         0.0,    112548.0, 
		   105316.0,    106372.0,    116118.0,    107304.0,    106213.0,         0.0,    101547.0,     95310.0,     91782.0,     91995.0, 
		    90393.0,     89470.0,         0.0,     86933.0,         0.0,     87181.0,     81565.0,     83063.0,     78357.0,     76459.0, 
		    77421.0,         0.0,     74627.0,     73528.0,     75292.0,     69292.0,         0.0,     74155.0,     65897.0,         0.0, 
		    68477.0,     68845.0,     67523.0,     62534.0,     61955.0,     59223.0,         0.0,     63448.0,     59471.0,     58280.0, 
		        7.0,     55112.0,     56541.0,         0.0,     54530.0,     53913.0,     52699.0,     49423.0,     50351.0,     50088.0, 
		    52160.0,         0.0,     48647.0,     47294.0,     45889.0,     60040.0,     56780.0,     71158.0,         0.0,     54132.0, 
		    69983.0,     53413.0,     54215.0,     63994.0,     65048.0,         0.0,     61921.0,     61022.0,     59527.0,     58432.0, 
		    48333.0,         0.0,     48505.0,     48674.0,     46976.0,     35317.0,     65127.0,     44354.0,         0.0,     45218.0, 
		    42768.0,     39849.0,     42103.0,     95483.0,     90990.0,     98062.0,     90283.0,    101098.0,     82158.0,     72890.0, 
		    70762.0,     83613.0,     80477.0,     82991.0,     60528.0,     79189.0,     78583.0,     74303.0,     75870.0,     73930.0, 
		    72032.0,     72951.0,     75030.0,    102318.0,     62335.0,     61197.0,     88373.0,     53576.0,    103091.0,     60006.0, 
		    54252.0,     83595.0,     60955.0,     75724.0,     50129.0,     80061.0,     80477.0,     76996.0,     52265.0,     77370.0, 
		    76922.0,     72300.0,     72957.0,     68846.0,     68733.0,     57080.0,     68629.0,     66338.0,     98329.0,     64495.0, 
		    62137.0,     59695.0,     84638.0,     70312.0,     61333.0,     65529.0,     63363.0,     60876.0,     58171.0,     57680.0, 
		    82559.0,     52618.0,     74584.0,     47807.0,     74043.0,     69587.0,     65811.0,     58158.0,     56437.0,     81881.0, 
		    50398.0,     68305.0,     47781.0,     73504.0,     43212.0,     55776.0,     62530.0,     50128.0,     41531.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
	0.0};



	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	TH1* locHist_RFBeamBunchPeriod = (TH1*)gDirectory->Get("RFBeamBunchPeriod");
	TH1* locHist_RFBeamBunchPeriod_DFT = (TH1*)gDirectory->Get("RFBeamBunchPeriod_DFT");
	TH1* locHist_BeamEnergy = (TH1*)gDirectory->Get("BeamEnergy");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Beam", "Beam", 1200, 600); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 1);

	TLatex latex;
	latex.SetTextSize(0.04);
	char str[256];

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if( (locHist_RFBeamBunchPeriod_DFT != NULL) && (locHist_RFBeamBunchPeriod != NULL))
	{
		locHist_RFBeamBunchPeriod_DFT->GetXaxis()->SetTitleSize(0.05);
		locHist_RFBeamBunchPeriod_DFT->GetXaxis()->SetLabelSize(0.05);
		locHist_RFBeamBunchPeriod_DFT->GetYaxis()->SetLabelSize(0.03);
		locHist_RFBeamBunchPeriod_DFT->SetFillStyle(3001);
		locHist_RFBeamBunchPeriod_DFT->SetFillColor(kGreen);
		locHist_RFBeamBunchPeriod_DFT->SetLineColor(kGreen-2);
		locHist_RFBeamBunchPeriod_DFT->SetLineWidth(2);
		locHist_RFBeamBunchPeriod_DFT->SetStats(0);
		locHist_RFBeamBunchPeriod_DFT->Draw();

		sprintf(str, "%d entries", (uint32_t)locHist_RFBeamBunchPeriod->GetEntries());
		latex.DrawLatex(300.0, 1.06*locHist_RFBeamBunchPeriod_DFT->GetMaximum(), str);
		
		// Plot RFBeamBunchPeriod as inset
		TPad *p = (TPad*)gDirectory->FindObjectAny("RF_DFT");
		if(!p) p = new TPad("RF_DFT", "insert", 0.1, 0.6, 0.5, 0.9);
		p->Draw();
		p->cd();

		locHist_RFBeamBunchPeriod->GetXaxis()->SetRangeUser(30.0, 90.0);
		locHist_RFBeamBunchPeriod->GetXaxis()->SetTitleSize(0.05);
		locHist_RFBeamBunchPeriod->GetXaxis()->SetLabelSize(0.05);
		locHist_RFBeamBunchPeriod->GetYaxis()->SetLabelSize(0.03);
		locHist_RFBeamBunchPeriod->SetFillStyle(3001);
		locHist_RFBeamBunchPeriod->SetFillColor(kMagenta);
		locHist_RFBeamBunchPeriod->SetStats(0);
		locHist_RFBeamBunchPeriod->Draw();
		
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BeamEnergy != NULL)
	{
		// Create normalized histogram
		TH1D* locHist_BeamEnergy_norm = (TH1D*)gDirectory->Get("BeamEnergy_norm");
		if(!locHist_BeamEnergy_norm){
			locHist_BeamEnergy_norm = new TH1D("BeamEnergy_norm", "Reconstructed photon beam energy normalized to amorphous run 30437;Beam #gamma energy (GeV)", 240, 0.0, 12.0);

			locHist_BeamEnergy_norm->GetXaxis()->SetTitleSize(0.05);
			locHist_BeamEnergy_norm->GetXaxis()->SetLabelSize(0.05);
			locHist_BeamEnergy_norm->GetYaxis()->SetLabelSize(0.03);
			locHist_BeamEnergy_norm->SetFillStyle(3001);
			locHist_BeamEnergy_norm->SetFillColor(kOrange);
			locHist_BeamEnergy_norm->SetLineColor(kRed-2);
			locHist_BeamEnergy_norm->SetLineWidth(2);
			locHist_BeamEnergy_norm->SetStats(0);
		}
		if(locHist_BeamEnergy_norm){
			for(int ibin=1; ibin<=locHist_BeamEnergy_norm->GetNbinsX(); ibin++){
				Double_t norm = amorphous_data[ibin-1];
				if( norm < 1000.0) continue;

				Double_t v = (Double_t)locHist_BeamEnergy->GetBinContent(ibin);
				locHist_BeamEnergy_norm->SetBinContent(ibin, v/norm);
			}
			locHist_BeamEnergy_norm->Draw();
		}

		TPad *beamenergypad = (TPad*)gDirectory->FindObjectAny("beamenergypad");
		if(!beamenergypad) beamenergypad = new TPad("beamenergypad", "", 0.105, 0.65, 0.5, 0.895);
		beamenergypad->SetTicks();
		beamenergypad->Draw();
		beamenergypad->cd();
	
		//locHist_BeamEnergy->Rebin(5);
		locHist_BeamEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_BeamEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_BeamEnergy->GetYaxis()->SetLabelSize(0.03);
		locHist_BeamEnergy->SetFillStyle(3001);
		locHist_BeamEnergy->SetFillColor(kCyan);
		locHist_BeamEnergy->SetLineColor(kBlue);
		locHist_BeamEnergy->SetLineWidth(2);
		locHist_BeamEnergy->SetStats(0);
		locHist_BeamEnergy->Draw();

		sprintf(str, "%d entries", (uint32_t)locHist_BeamEnergy->GetEntries());
		latex.DrawLatex(2.0, 1.06*locHist_BeamEnergy->GetMaximum(), str);
	}
}


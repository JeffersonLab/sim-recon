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
	// online monitoring root file for run 30618, an amorphous target run.
	// The macro used is in /gluonwork1/Users/davidl/2017.02.08.amorphous_norm
	string amorphous_label = "Normalized to Amorphous run 30635";

	Double_t amorphous_data[] = {
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,   6480495.0,   6690711.0, 
		  6072670.0,     72706.0,   5900728.0,   5549488.0,   5960327.0,   5603308.0,   5616157.0,   5536458.0,     60999.0,   5424051.0, 
		  5067497.0,   5128762.0,   5547658.0,   5192341.0,   5110649.0,     62981.0,   4864402.0,   4609476.0,   4450636.0,   4375820.0, 
		  4337742.0,   4324010.0,     33510.0,   4172194.0,     37347.0,   4186537.0,   3919215.0,   3977801.0,   3743500.0,   3670956.0, 
		  3723197.0,     30552.0,   3589869.0,   3537469.0,   3610762.0,   3310315.0,         0.0,   3567404.0,   3184722.0,     22259.0, 
		  3305376.0,   3301878.0,   3205892.0,   3020151.0,   2985414.0,   2853633.0,     18794.0,   3039094.0,   2862473.0,   2768211.0, 
		      184.0,   2664602.0,   2724465.0,     18551.0,   2624968.0,   2597269.0,   2523356.0,   2368681.0,   2409156.0,   2407717.0, 
		  2504065.0,     13381.0,   2336774.0,   2262871.0,   2198572.0,   2875947.0,   2740131.0,   3431027.0,     20435.0,   2608577.0, 
		  3371560.0,   2554849.0,   2593550.0,   3070361.0,   3117995.0,     20972.0,   2977981.0,   2921731.0,   2860211.0,   2775698.0, 
		  2319696.0,     16369.0,   2322039.0,   2324831.0,   2237517.0,   1697384.0,   3152409.0,   2131754.0,     10109.0,   2166119.0, 
		  2056859.0,   1913953.0,   2009933.0,   4581584.0,   4350164.0,   4716292.0,   4333120.0,   4822289.0,   3660088.0,   3018872.0, 
		  2949291.0,   3697225.0,   3277885.0,   3546454.0,   2724227.0,   3374876.0,   3367514.0,   3259825.0,   3165758.0,   3151004.0, 
		  3044984.0,   3210430.0,   3199895.0,   3943359.0,   2811515.0,   2354102.0,   3826929.0,   2404213.0,   3945339.0,   2879126.0, 
		  3013577.0,   3569501.0,   2916109.0,   3606024.0,   2414077.0,   3823164.0,   3842668.0,   3698215.0,   2628706.0,   3557183.0, 
		  3671198.0,   3462811.0,   3471662.0,   3292499.0,   3271814.0,   2783141.0,   3264741.0,   3187159.0,   4672615.0,   3089888.0, 
		  2974877.0,   2874173.0,   4056566.0,   3367096.0,   2930448.0,   3135150.0,   3033929.0,   2920002.0,   2844815.0,   2726648.0, 
		  3946867.0,   2555125.0,   3565399.0,   2311862.0,   3496942.0,   3338569.0,   3144022.0,   2796863.0,   2701768.0,   3942299.0, 
		  2429586.0,   3289743.0,   2349074.0,   3589586.0,   2139634.0,   2753187.0,   3093805.0,   2510955.0,   2040040.0,         0.0, 
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
		double max = locHist_RFBeamBunchPeriod_DFT->GetMaximum()*2.0;
		locHist_RFBeamBunchPeriod_DFT->GetXaxis()->SetTitleSize(0.05);
		locHist_RFBeamBunchPeriod_DFT->GetXaxis()->SetLabelSize(0.05);
		locHist_RFBeamBunchPeriod_DFT->GetYaxis()->SetLabelSize(0.03);
		locHist_RFBeamBunchPeriod_DFT->SetFillStyle(3001);
		locHist_RFBeamBunchPeriod_DFT->SetFillColor(kGreen);
		locHist_RFBeamBunchPeriod_DFT->SetLineColor(kGreen-2);
		locHist_RFBeamBunchPeriod_DFT->SetLineWidth(2);
		locHist_RFBeamBunchPeriod_DFT->SetStats(0);
		locHist_RFBeamBunchPeriod_DFT->GetYaxis()->SetRangeUser(0.0, max);
		locHist_RFBeamBunchPeriod_DFT->Draw();

		sprintf(str, "%g entries", (double)locHist_RFBeamBunchPeriod->GetEntries());
		latex.DrawLatex(300.0, 1.02*max, str);
		
		// Plot RFBeamBunchPeriod as inset
		TPad *p = (TPad*)gDirectory->FindObjectAny("RF_DFT");
		if(!p) p = new TPad("RF_DFT", "insert", 0.16, 0.6, 0.89, 0.9);
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
			locHist_BeamEnergy_norm = new TH1D("BeamEnergy_norm", "Reconstructed Photon Beam Energy;Beam #gamma energy (GeV)", 240, 0.0, 12.0);

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
		
			// Find leftmost non-zero bin for scaling whole hist
			double scale = 0.0;
			for(int ibin=1; ibin<=locHist_BeamEnergy->GetNbinsX(); ibin++){
				if( amorphous_data[ibin-1] < 1000.0) continue;
				Double_t v = (Double_t)locHist_BeamEnergy->GetBinContent(ibin);
				if(v>0.0){
					scale = v/amorphous_data[ibin-1];
					break;
				}
			}

			// Normalize to amorphous baseline including scaling factor
			// above so that leftmost non-zero bin is always 1
			double max = 0.0;
			for(int ibin=1; ibin<=locHist_BeamEnergy_norm->GetNbinsX(); ibin++){
				Double_t norm = amorphous_data[ibin-1];
				if( norm < 1000.0) continue;
				norm *= scale;

				Double_t v = (Double_t)locHist_BeamEnergy->GetBinContent(ibin);
				locHist_BeamEnergy_norm->SetBinContent(ibin, v/norm);
				
				if(v/norm  > max) max = v/norm;
			}
			max = (max-0.75)*1.05 + 0.75;
			if(max > 3.0) max = 3.0;
			if(max < 1.4) max = 1.4;
			locHist_BeamEnergy_norm->GetYaxis()->SetRangeUser(0.75, max);
			locHist_BeamEnergy_norm->Draw();
			
			// If maximum is > 1.5 then assume this is not an amorphous run
			// and draw a label of the peak energy
			if(max > 1.5){
				double Epeak = locHist_BeamEnergy_norm->GetBinCenter(locHist_BeamEnergy_norm->GetMaximumBin());
				sprintf(str, "Epeak: %3.2f GeV", Epeak);
				latex.SetTextAlign(12);
				latex.DrawLatex(1.0, 0.5*(max-0.75)+0.75, str);
			}

			sprintf(str, "%g entries", (double)locHist_BeamEnergy->GetEntries());
			latex.SetTextAlign(22);
			latex.DrawLatex(6.0, 1.035*(max-0.75)+0.75, str);

			latex.SetTextAngle(270);
			latex.DrawLatex(12.5, 0.5*(max-0.75)+0.75, amorphous_label.c_str());
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
	}
}


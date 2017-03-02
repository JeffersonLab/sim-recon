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
	string amorphous_label = "Normalized to Amorphous run 30839";

	Double_t amorphous_data[] = {
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,     59002.0,     55554.0,     54664.0,         0.0,     51564.0,     52808.0,     51092.0,     50718.0, 
		    49471.0,     46579.0,     46933.0,     46995.0,     48366.0,     45802.0,     44626.0,   8788738.0,   8304443.0,   8180075.0, 
		    66517.0,   7739883.0,   7943276.0,     49450.0,   7646915.0,   7580332.0,   7405749.0,   7008976.0,   7095583.0,   7051528.0, 
		  7307102.0,         0.0,   6872107.0,   6660319.0,   6489457.0,   8375291.0,   7964083.0,   9935822.0,     42136.0,   7615902.0, 
		  9715414.0,   7460525.0,   7654250.0,   8942988.0,   9108202.0,     99239.0,   8697543.0,   8482876.0,   8328684.0,   8214948.0, 
		  6854912.0,     76483.0,   6842285.0,   6852213.0,   6654311.0,   5067410.0,   9154980.0,   6318726.0,     64752.0,   6385929.0, 
		  6116479.0,   5690496.0,   5969774.0,  13616461.0,  12930485.0,  11712460.0,  12890584.0,  14397880.0,   8436157.0,   7935417.0, 
		  8128062.0,  10089952.0,  10093994.0,   9516037.0,   7281862.0,   9468792.0,   9400915.0,   9379595.0,   8798702.0,   8645830.0, 
		  8546285.0,   8052767.0,   7998508.0,   9574554.0,   6495475.0,   7683286.0,   8535819.0,   7198713.0,  12682343.0,   9512809.0, 
		  7752982.0,  11895880.0,   8727576.0,  10791075.0,   7178164.0,  11520077.0,  11521476.0,  11085626.0,   7463967.0,  11118820.0, 
		 11022009.0,  10392706.0,  10417068.0,   9915180.0,   9846869.0,   8214725.0,   9768965.0,   9471939.0,  13930227.0,   8850982.0, 
		  8798099.0,   8486296.0,  12094598.0,   9984665.0,   8702569.0,   9285678.0,   8991655.0,   8663223.0,   8286451.0,   8188784.0, 
		 11802725.0,   7493393.0,  10670477.0,   6835822.0,  10394698.0,   9840321.0,   9334801.0,   8278408.0,   7967535.0,  11610292.0, 
		  7072116.0,   9654340.0,   6766942.0,  10398622.0,   6094841.0,   7851124.0,   8808280.0,   7058197.0,   5778494.0,         0.0, 
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
		

			// Normalize to amorphous baseline 
			for(int ibin=1; ibin<=locHist_BeamEnergy_norm->GetNbinsX(); ibin++){
				Double_t norm = amorphous_data[ibin-1];
				if( norm < 100000.0) continue;

				Double_t v = (Double_t)locHist_BeamEnergy->GetBinContent(ibin);
				locHist_BeamEnergy_norm->SetBinContent(ibin, v/norm);
			}

			// Find maximum and scale histogram so peak is at 1
			double scale = (Double_t)locHist_BeamEnergy_norm->GetMaximum();
			locHist_BeamEnergy_norm->Scale(1.0/scale);		

			// Find leftmost non-zero bin 
			double left_scale = 0.0;
			for(int ibin=1; ibin<=locHist_BeamEnergy_norm->GetNbinsX(); ibin++){
				if( amorphous_data[ibin-1] < 1000.0) continue;
				Double_t v = (Double_t)locHist_BeamEnergy_norm->GetBinContent(ibin);
				if(v>0.1){
					left_scale = v;
					break;
				}
			}
			
			// Set min and max of plot to show range
			double diff = 1.0 - left_scale;
			double min = left_scale - diff*0.1;
			double max = left_scale + diff*1.2;
			if( fabs(max-min) < 0.4 ){
				min = 0.5;
				max =1.5;
			}
			
			locHist_BeamEnergy_norm->GetYaxis()->SetRangeUser(min, max);
			locHist_BeamEnergy_norm->Draw();
			
			// If diff is > 0.3 then assume this is not an amorphous run
			// and draw a label of the peak energy
			if(diff > 0.3){
				double Epeak = locHist_BeamEnergy_norm->GetBinCenter(locHist_BeamEnergy_norm->GetMaximumBin());
				sprintf(str, "Epeak: %3.2f GeV", Epeak);
				latex.SetTextAlign(12);
				latex.DrawLatex(1.0, (min+max)/2.0, str);
			}

			sprintf(str, "%g entries", (double)locHist_BeamEnergy->GetEntries());
			latex.SetTextAlign(22);
			latex.DrawLatex(6.0, 1.035*(max-min)+min, str);

			latex.SetTextAngle(270);
			latex.DrawLatex(12.5, (min+max)/2.0, amorphous_label.c_str());
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


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
	string amorphous_label = "Normalized to Amorphous run 30808";

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
		        0.0,         0.0,    100768.0,     95650.0,     94150.0,         0.0,     88692.0,     91703.0,     87631.0,     87937.0, 
		    85061.0,     80428.0,     81562.0,     80657.0,     83799.0,     78764.0,     76626.0,   8792911.0,   8317217.0,   8201228.0, 
		   114649.0,   7752784.0,   7966827.0,     85925.0,   7663143.0,   7597263.0,   7428433.0,   7023925.0,   7109020.0,   7076531.0, 
		  7326945.0,         0.0,   6887713.0,   6676535.0,   6503579.0,   8375514.0,   7963187.0,   9935152.0,     71453.0,   7623876.0, 
		  9708973.0,   7463279.0,   7672715.0,   8967789.0,   9131950.0,    170299.0,   8756174.0,   8471753.0,   8354280.0,   8227098.0, 
		  6875793.0,    132330.0,   6863070.0,   6888314.0,   6677740.0,   5104907.0,   9164354.0,   6343633.0,    111025.0,   6405366.0, 
		  6143878.0,   5709977.0,   5982726.0,  13606937.0,  12944968.0,  11738335.0,  12914732.0,  14409128.0,   8140100.0,   7934019.0, 
		  8156563.0,  10123315.0,  10118132.0,   9535642.0,   7301814.0,   9500738.0,   9435191.0,   9403256.0,   8835431.0,   8683418.0, 
		  8557447.0,   8067092.0,   8032815.0,   9529378.0,   6536625.0,   7713771.0,   8542687.0,   7204053.0,  12720393.0,   9492931.0, 
		  7792337.0,  11883632.0,   8746341.0,  10808774.0,   7226866.0,  11542299.0,  11531197.0,  11096566.0,   7476257.0,  11142268.0, 
		 11046503.0,  10421249.0,  10418394.0,   9934460.0,   9869122.0,   8247500.0,   9784511.0,   9562820.0,  13904267.0,   8830215.0, 
		  8777431.0,   8462338.0,  12066954.0,   9954624.0,   8677448.0,   9264007.0,   8969777.0,   8645879.0,   8267051.0,   8164448.0, 
		 11791013.0,   7482360.0,  10661033.0,   6834552.0,  10380241.0,   9840795.0,   9335792.0,   8281151.0,   7994163.0,  11674755.0, 
		  7133660.0,   9775747.0,   6873756.0,  10611020.0,   6247860.0,   8085751.0,   9124534.0,   7335488.0,   6031261.0,         0.0, 
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
				if( norm < 100000.0) continue;
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


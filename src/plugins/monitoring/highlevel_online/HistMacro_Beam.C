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
	string amorphous_label = "Normalized to Amorphous run 30618";

	Double_t amorphous_data[] = {
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,  13329101.0,  13767317.0, 
		 12500182.0,    146010.0,  12127205.0,  11405040.0,  12254831.0,  11519598.0,  11548892.0,  11381739.0,    123400.0,  11137576.0, 
		 10413430.0,  10546938.0,  11404778.0,  10675637.0,  10508055.0,    127267.0,  10000171.0,   9476094.0,   9155157.0,   8999014.0, 
		  8913742.0,   8885251.0,     66640.0,   8577218.0,     75579.0,   8598786.0,   8048925.0,   8187400.0,   7685723.0,   7544269.0, 
		  7638135.0,     61565.0,   7378249.0,   7262264.0,   7414243.0,   6794648.0,         0.0,   7327027.0,   6537204.0,     44564.0, 
		  6791211.0,   6783999.0,   6580332.0,   6199766.0,   6134268.0,   5856378.0,     38034.0,   6239810.0,   5883325.0,   5692989.0, 
		      378.0,   5470929.0,   5594835.0,     37553.0,   5388025.0,   5334477.0,   5182329.0,   4860160.0,   4950248.0,   4940420.0, 
		  5143916.0,     27042.0,   4798379.0,   4651533.0,   4517473.0,   5911847.0,   5630595.0,   7049546.0,     41236.0,   5358871.0, 
		  6927380.0,   5246739.0,   5326394.0,   6315068.0,   6404606.0,     41933.0,   6117280.0,   5998597.0,   5875528.0,   5701562.0, 
		  4766096.0,     32367.0,   4763668.0,   4777561.0,   4597048.0,   3489429.0,   6472842.0,   4384097.0,     20294.0,   4449566.0, 
		  4224787.0,   3932106.0,   4128773.0,   9405183.0,   8936691.0,   9668481.0,   8902766.0,   9893552.0,   7531781.0,   6194120.0, 
		  6063309.0,   7590929.0,   6733585.0,   7295465.0,   5604939.0,   6921966.0,   6915215.0,   6702086.0,   6498747.0,   6471440.0, 
		  6269830.0,   6612158.0,   6595438.0,   8118514.0,   5803073.0,   4861235.0,   7843359.0,   4941311.0,   8105759.0,   5916816.0, 
		  6181315.0,   7338504.0,   5985119.0,   7402252.0,   4967463.0,   7854303.0,   7878154.0,   7603787.0,   5393027.0,   7301684.0, 
		  7540077.0,   7106116.0,   7132262.0,   6768025.0,   6720524.0,   5715467.0,   6702646.0,   6582716.0,   9601140.0,   6339619.0, 
		  6101534.0,   5895715.0,   8328724.0,   6910694.0,   6019151.0,   6431977.0,   6228568.0,   5994324.0,   5842111.0,   5595672.0, 
		  8109059.0,   5247465.0,   7321161.0,   4742971.0,   7175381.0,   6854878.0,   6462803.0,   5751180.0,   5556189.0,   8099901.0, 
		  4996609.0,   6776256.0,   4833678.0,   7394648.0,   4409897.0,   5677410.0,   6393813.0,   5188749.0,   4207139.0,         0.0, 
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


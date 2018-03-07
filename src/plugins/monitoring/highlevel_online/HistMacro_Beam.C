// hnamepath: /highlevel/RFBeamBunchPeriod
// hnamepath: /highlevel/RFBeamBunchPeriod_DFT
// hnamepath: /highlevel/BeamEnergy
//
// e-mail: davidl@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
// e-mail: tbritton@jlab.org
//

{
// This taken from the bin contents of the BeamEnergy histogram from the
// online monitoring root file for run 41180, an amorphous target run.
// The macro used is in /gluonwork1/Users/davidl/2017.02.08.amorphous_norm
string amorphous_label = "Normalized to Amorphous run 41180";

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
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,    209738.0,    198138.0,         0.0,    193285.0, 
		   182928.0,    189803.0,    188623.0,    182261.0,    180752.0,    176329.0,         0.0,    167626.0,    169955.0,    168081.0, 
		   173281.0,    163542.0,    158071.0,         0.0,    154996.0,    200971.0,    190962.0,    238371.0,    182560.0,    233500.0, 
		        0.0,    179005.0,    183107.0,    212746.0,    217870.0,    206988.0,    203560.0,    199817.0,         0.0,    196465.0, 
		   163990.0,    163596.0,    163323.0,    158890.0,         0.0,    120653.0,    219450.0,    150392.0,    152018.0,    146696.0, 
		   135199.0,         0.0,    242118.0,    311081.0,    343239.0,    397432.0,    303272.0,    321248.0,    217471.0,    179834.0, 
		   223478.0,    204684.0,    175483.0,    212535.0,    199347.0,    204616.0,    206993.0,    200380.0,    192812.0,    187961.0, 
		   186785.0,    192361.0,    185509.0,    221479.0,    155732.0,    165711.0,    211286.0,    213044.0,    162172.0,    218470.0, 
		   326773.0,    173718.0,    260586.0,    175565.0,    283390.0,    190320.0,    277235.0,    271381.0,    270542.0,    258743.0, 
		   254601.0,    253674.0,    242940.0,    239718.0,    271785.0,    233869.0,    227794.0,    220098.0,    222537.0,    211732.0, 
		   202200.0,    302707.0,    216286.0,    216766.0,    227206.0,    218474.0,    214509.0,    207579.0,    287877.0,    191310.0, 
		   278819.0,    175601.0,    245601.0,    252577.0,    124981.0,    231491.0,    212180.0,    299163.0,    190505.0,    175994.0, 
		   250525.0,    153673.0,    258291.0,    154869.0,    207703.0,    219492.0,    165417.0,    124383.0,     24692.0,         0.0, 
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
		}
		if(locHist_BeamEnergy_norm){

			locHist_BeamEnergy_norm->Reset();

			// Normalize to amorphous baseline 
			double scale = 0.0;
			for(int ibin=1; ibin<=locHist_BeamEnergy_norm->GetNbinsX(); ibin++){
				Double_t norm = amorphous_data[ibin-1];
				if( norm < 10000.0) continue;

				Double_t v = (Double_t)locHist_BeamEnergy->GetBinContent(ibin);
				v /= norm;
				locHist_BeamEnergy_norm->SetBinContent(ibin, v);
				if(v > scale) scale = v;
			}

			// Find maximum and scale histogram so peak is at 1
			if(scale != 0.0) locHist_BeamEnergy_norm->Scale(1.0/scale);		

			// Find leftmost non-zero bin 
			double left_scale = 0.0;
			for(int ibin=1; ibin<=locHist_BeamEnergy_norm->GetNbinsX(); ibin++){
				if( amorphous_data[ibin-1] < 10000.0) continue;
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
			if( (max-min) < 0.4 ){
				min = 0.25;
				if(max<1.50) max =1.50;
			}
			
			// max = 1.50;
			// min = 0.25;
			
			locHist_BeamEnergy_norm->GetXaxis()->SetTitleSize(0.05);
			locHist_BeamEnergy_norm->GetXaxis()->SetLabelSize(0.05);
			locHist_BeamEnergy_norm->GetYaxis()->SetLabelSize(0.03);
			locHist_BeamEnergy_norm->SetFillStyle(3001);
			locHist_BeamEnergy_norm->SetFillColor(kOrange);
			locHist_BeamEnergy_norm->SetLineColor(kRed-2);
			locHist_BeamEnergy_norm->SetLineWidth(2);
			locHist_BeamEnergy_norm->SetStats(0);		
			locHist_BeamEnergy_norm->GetYaxis()->SetRangeUser(min, max);
			locHist_BeamEnergy_norm->Draw("HIST");
			
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


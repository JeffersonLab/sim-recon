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
	// online monitoring root file for run 30407, an amorphous target run.
	// The macro used is in /gluonwork1/Users/davidl/2017.02.08.amorphous_norm
	Double_t amorphous_data[] = {
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		   268815.0,    277165.0,    253090.0,    244280.0,    231710.0,    245631.0,    234140.0,    232976.0,   8271512.0,   8542142.0, 
		  7742462.0,    299650.0,   7543624.0,   7096752.0,   7601820.0,   7139019.0,   7153091.0,   7035370.0,    256485.0,   6888551.0, 
		  6448155.0,   6349710.0,   7059659.0,   6608391.0,   6495564.0,    241405.0,   6181591.0,   5862059.0,   5670413.0,   5567122.0, 
		  5526595.0,   5505217.0,    176622.0,   5171378.0,    193276.0,   5325682.0,   4985749.0,   5066447.0,   4762809.0,   4668980.0, 
		  4730295.0,    155103.0,   4570354.0,   4492606.0,   4583049.0,   4090891.0,    110200.0,   4534973.0,   4048191.0,    134425.0, 
		  4199169.0,   4186735.0,   4074072.0,   3842270.0,   3793906.0,   3624450.0,    116349.0,   3855018.0,   3664349.0,   3544608.0, 
		   142743.0,   3406249.0,   3508911.0,    128214.0,   3355178.0,   3343800.0,   3249953.0,   3049285.0,   3101735.0,   3095700.0, 
		  3215569.0,     16683.0,   2987397.0,   2896420.0,   2815198.0,   3653910.0,   3459118.0,   4374118.0,    111709.0,   3317103.0, 
		  4260374.0,   3241622.0,   3352044.0,   4007611.0,   4065627.0,    238653.0,   3939288.0,   3925749.0,   4255606.0,   3961364.0, 
		  3634160.0,    648586.0,   3326504.0,   3491000.0,   3325419.0,   2652889.0,   4513419.0,   3210662.0,    378791.0,   3056241.0, 
		  3004919.0,   2801270.0,   2931685.0,   6121998.0,   5769367.0,   6012010.0,   5544758.0,   6171668.0,  11266127.0,  14535415.0, 
		 13910871.0,  16323468.0,  22401178.0,  15267750.0,  10666054.0,  16501947.0,  15467434.0,  13475669.0,  15574225.0,  15366381.0, 
		 12578323.0,   9430888.0,   9435130.0,  13386909.0,  11086492.0,  11059303.0,  14301954.0,  10089998.0,   9532201.0,   2981357.0, 
		  3906894.0,   4547482.0,   3786136.0,   4636059.0,   3158099.0,   4932486.0,   4895401.0,   4723181.0,   3382190.0,   4569131.0, 
		  4702079.0,   4436688.0,   4398409.0,   4233152.0,   4193743.0,   3600813.0,   4160292.0,   4069417.0,   5828748.0,   3820512.0, 
		  3691527.0,   3569672.0,   5026853.0,   4175317.0,   3639197.0,   3883132.0,   3757802.0,   3616470.0,   3530535.0,   3372157.0, 
		  4883899.0,   3171457.0,   4408503.0,   2855637.0,   4325779.0,   4129167.0,   3883291.0,   3450540.0,   3324912.0,   4830021.0, 
		  2961071.0,   3995855.0,   2835169.0,   4312612.0,   2552515.0,   3265914.0,   3663745.0,   2978067.0,   2433008.0,         0.0, 
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
			locHist_BeamEnergy_norm = new TH1D("BeamEnergy_norm", "Reconstructed photon beam energy normalized to amorphous run 30407;Beam #gamma energy (GeV)", 240, 0.0, 12.0);

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
				if( norm == 0.0) continue;

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


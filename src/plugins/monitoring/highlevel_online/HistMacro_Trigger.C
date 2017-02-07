//
// NOTE: The L1GTPRate histo now comes from the 
// occupancy_online plugin since it should see many
// more sync events.
//
//
//
// hnamepath: /occupancy/L1GTPRate
// hnamepath: /occupancy/L1livetime
// hnamepath: /highlevel/BCALVsFCAL_TrigBit1
// hnamepath: /highlevel/L1bits_gtp
// hnamepath: /highlevel/L1bits_fp
//
// e-mail: davidl@jlab.org
// e-mail: pmatt@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
//

{
	TDirectory *locTopDirectory = gDirectory;


	// Grab remaining histos from highlevel directory
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	TH2* locHist_BCALVsFCAL_TrigBit1 = (TH2*)gDirectory->Get("BCALVsFCAL_TrigBit1");
	TH1* locHist_L1bits_gtp = (TH1*)gDirectory->Get("L1bits_gtp");
	TH1* locHist_L1bits_fp = (TH1*)gDirectory->Get("L1bits_fp");

	// Grab a couple of histos from occupancy directory
	locDirectory = (TDirectory*)locTopDirectory->FindObjectAny("occupancy");
	TH2* locHist_L1GTPRate  = NULL;
	TH1* locHist_L1livetime = NULL;
	if(locDirectory){
		locHist_L1GTPRate = (TH2*)locDirectory->Get("L1GTPRate");
		locHist_L1livetime = (TH1*)locDirectory->Get("L1livetime");
	}

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Kinematics", "Kinematics", 1200, 400); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 1);
	
	TLatex latex;
	latex.SetTextSize(0.04);
	char str[256];

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_L1GTPRate != NULL)
	{
		locHist_L1GTPRate->GetXaxis()->SetTitleSize(0.05);
		locHist_L1GTPRate->GetYaxis()->SetTitleSize(0.04);
		locHist_L1GTPRate->GetXaxis()->SetLabelSize(0.05);
		locHist_L1GTPRate->GetYaxis()->SetLabelSize(0.05);
		locHist_L1GTPRate->SetStats(0);
		locHist_L1GTPRate->Draw("colz");
		
		sprintf(str, "from %d sync events", (uint32_t)locHist_L1GTPRate->GetEntries()/8);
		latex.DrawLatex(1.0, 101.0, str);
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALVsFCAL_TrigBit1 != NULL)
	{
		locHist_BCALVsFCAL_TrigBit1->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALVsFCAL_TrigBit1->GetYaxis()->SetTitleSize(0.04);
		locHist_BCALVsFCAL_TrigBit1->SetStats(0);
		locHist_BCALVsFCAL_TrigBit1->Draw("colz");

		sprintf(str, "%d entries", (uint32_t)locHist_BCALVsFCAL_TrigBit1->GetEntries());
		latex.DrawLatex(500.0, 50000.0*1.01, str);

		gPad->SetLogz();
		gPad->Update();
	}

	// We want to draw two plots in the right 1/3 of canvas. The
	// Divide function will clear the canvas wiping out he above
	// histos we just drew. Instead, explicitly make a TPad being
	// sure to recycle any existing one to avoid memory leaks or
	// error messages.
	locCanvas->cd(0);
	TPad *trigpad1 = (TPad*)gDirectory->FindObjectAny("trigpad1");
	if(!trigpad1) trigpad1 = new TPad("trigpad1", "", 0.66, 0.5, 1.0, 1.0);
	trigpad1->SetTicks();
	//trigpad1->SetLeftMargin(0.10);
	//trigpad1->SetRightMargin(0.15);
	trigpad1->Draw();
	trigpad1->cd();

	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_L1bits_gtp!=NULL && locHist_L1bits_fp!=NULL)
	{
		double max_gtp = locHist_L1bits_gtp->GetMaximum();
		double max_fp  = locHist_L1bits_fp->GetMaximum();
		double max = (max_gtp>max_fp) ? max_gtp:max_fp;
		
		locHist_L1bits_gtp->GetYaxis()->SetRangeUser(1.0, max*2.0);
		locHist_L1bits_gtp->GetXaxis()->SetRangeUser(0.0, 17.0);
	
		locHist_L1bits_gtp->GetXaxis()->SetTitleSize(0.05);
		locHist_L1bits_gtp->GetYaxis()->SetTitleSize(0.04);
		locHist_L1bits_gtp->SetStats(0);

		locHist_L1bits_gtp->SetLineColor(kBlack);
		locHist_L1bits_gtp->SetBarOffset(0.15);
		locHist_L1bits_gtp->SetBarWidth(0.35);
		locHist_L1bits_fp->SetBarOffset(0.5);
		locHist_L1bits_fp->SetBarWidth(0.35);

		locHist_L1bits_gtp->SetFillColor(kOrange);
		locHist_L1bits_fp->SetFillColor(kRed-4);

		locHist_L1bits_gtp->Draw("bar");
		locHist_L1bits_fp->Draw("bar same");

		TLegend *legend_gtp = new TLegend(0.5,0.85,0.7,0.9);
		TLegend *legend_fp  = new TLegend(0.7,0.85,0.9,0.9);
		legend_gtp->AddEntry(locHist_L1bits_gtp,"GTP","f");
		legend_fp->AddEntry(locHist_L1bits_fp,"FP","f");
		legend_gtp->Draw();
		legend_fp->Draw();

		gPad->SetLogy();
		gPad->Update();
	}

	locCanvas->cd(0);
	TPad *trigpad2 = (TPad*)gDirectory->FindObjectAny("trigpad2");
	if(!trigpad2) trigpad2 = new TPad("trigpad2", "", 0.66, 0.0, 1.0, 0.5);
	trigpad2->SetTicks();
	trigpad2->Draw();
	trigpad2->cd();

	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_L1livetime!=NULL)
	{
		locHist_L1livetime->SetLineColor(kGreen-3);
		locHist_L1livetime->SetFillColor(kGreen);
		locHist_L1livetime->SetFillStyle(3001);
		locHist_L1livetime->Draw();
		gPad->SetLogy();
		gPad->Update();
	}
}

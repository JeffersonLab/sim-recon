// hnamepath: /p3pi_preco_2FCAL/Hist_InvariantMass_Pi0/InvariantMass
// hnamepath: /p3pi_preco_2FCAL/Hist_InvariantMass_Omega_PostKinFitCut/InvariantMass
// hnamepath: /p3pi_preco_2BCAL/Hist_InvariantMass_Pi0/InvariantMass
// hnamepath: /p3pi_preco_2BCAL/Hist_InvariantMass_Omega_PostKinFitCut/InvariantMass
// hnamepath: /p3pi_preco_FCAL-BCAL/Hist_InvariantMass_Pi0/InvariantMass
// hnamepath: /p3pi_preco_FCAL-BCAL/Hist_InvariantMass_Omega_PostKinFitCut/InvariantMass
// hnamepath: /p3pi_preco_any_kinfit/Hist_MissingMassSquared/MissingMassSquared
// hnamepath: /p3pi_preco_any_kinfit/Hist_MissingMassSquared_PostKinFitCut/MissingMassSquared
// hnamepath: /p3pi_preco_any_kinfit/Hist_InvariantMass_Pi0/InvariantMass
// hnamepath: /p3pi_preco_any_kinfit/Hist_InvariantMass_Pi0_PostKinFitCut/InvariantMass
// hnamepath: /p3pi_preco_any_kinfit/Hist_KinFitResults/ConfidenceLevel
// hnamepath: /p3pi_preco_any_kinfit/Hist_InvariantMass_Omega_PostKinFitCut/InvariantMass
// hnamepath: /p3pi_preco_any_kinfit/Hist_InvariantMass_Omega_KinFit_PostKinFitCut/InvariantMass


{	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locReactionDirectory_2FCAL = (TDirectory*)locInitDirectory->FindObjectAny("p3pi_preco_2FCAL");
	TDirectory *locReactionDirectory_2BCAL = (TDirectory*)locInitDirectory->FindObjectAny("p3pi_preco_2BCAL");
	TDirectory *locReactionDirectory_Both = (TDirectory*)locInitDirectory->FindObjectAny("p3pi_preco_FCAL-BCAL");
	TDirectory *locReactionDirectory_Any = (TDirectory*)locInitDirectory->FindObjectAny("p3pi_preco_any");
	TDirectory *locReactionDirectory_KinFit = (TDirectory*)locInitDirectory->FindObjectAny("p3pi_preco_any_kinfit");
	if((locReactionDirectory_2FCAL == NULL) || (locReactionDirectory_2BCAL == NULL) || (locReactionDirectory_Both == NULL))
		return;
	if((locReactionDirectory_Any == NULL) || (locReactionDirectory_KinFit == NULL))
		return;

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("p3pi", "p3pi", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	// Get overall normalization to number of triggers
	TH1D* locHist_NumEvents = (TH1D*) locReactionDirectory_KinFit->Get("NumEventsSurvivedAction");
	double n_triggers = locHist_NumEvents->GetBinContent(1);

	double n_omega_kinfit;
	double omega_mass = 0;
	double omega_width= 0;

	TH1I* locHist_Pi0_2FCAL = (TH1I*)locReactionDirectory_2FCAL->Get("Hist_InvariantMass_Pi0/InvariantMass");
	TH1I* locHist_Omega_2FCAL = (TH1I*)locReactionDirectory_2FCAL->Get("Hist_InvariantMass_Omega_PostKinFitCut/InvariantMass");

	TH1I* locHist_Pi0_2BCAL = (TH1I*)locReactionDirectory_2BCAL->Get("Hist_InvariantMass_Pi0/InvariantMass");
	TH1I* locHist_Omega_2BCAL = (TH1I*)locReactionDirectory_2BCAL->Get("Hist_InvariantMass_Omega_PostKinFitCut/InvariantMass");

	TH1I* locHist_Pi0_Both = (TH1I*)locReactionDirectory_Both->Get("Hist_InvariantMass_Pi0/InvariantMass");
	TH1I* locHist_Omega_Both = (TH1I*)locReactionDirectory_Both->Get("Hist_InvariantMass_Omega_PostKinFitCut/InvariantMass");

	TH1I* locHist_MM2 = (TH1I*)locReactionDirectory_KinFit->Get("Hist_MissingMassSquared/MissingMassSquared");
	TH1I* locHist_MM2_KinFitCut = (TH1I*)locReactionDirectory_KinFit->Get("Hist_MissingMassSquared_PostKinFitCut/MissingMassSquared");

	TH1I* locHist_Pi0 = (TH1I*)locReactionDirectory_KinFit->Get("Hist_InvariantMass_Pi0/InvariantMass");
	TH1I* locHist_Pi0_KinFitCut = (TH1I*)locReactionDirectory_KinFit->Get("Hist_InvariantMass_Pi0_PostKinFitCut/InvariantMass");

	TH1I* locHist_KinFitConLev = (TH1I*)locReactionDirectory_KinFit->Get("Hist_KinFitResults/ConfidenceLevel");

	TH1I* locHist_Omega_KinFitCut = (TH1I*)locReactionDirectory_KinFit->Get("Hist_InvariantMass_Omega_PostKinFitCut/InvariantMass");
	TH1I* locHist_KinFitOmega_KinFitCut = (TH1I*)locReactionDirectory_KinFit->Get("Hist_InvariantMass_Omega_KinFit_PostKinFitCut/InvariantMass");

	//Draw
	int locNumRebin = 5;

	//pi0 for 3 cases (all on top of each other)
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_Pi0_2FCAL != NULL) && (locHist_Pi0_2BCAL != NULL) && (locHist_Pi0_Both != NULL))
	{
		locHist_Pi0_2FCAL->Rebin(locNumRebin);
		locHist_Pi0_2BCAL->Rebin(locNumRebin);
		locHist_Pi0_Both->Rebin(locNumRebin);

		double locMaxHeight = locHist_Pi0_2FCAL->GetBinContent(locHist_Pi0_2FCAL->GetMaximumBin());
		double locOtherHeight = locHist_Pi0_2BCAL->GetBinContent(locHist_Pi0_2BCAL->GetMaximumBin());
		if(locOtherHeight > locMaxHeight)
			locMaxHeight = locOtherHeight;
		locOtherHeight = locHist_Pi0_Both->GetBinContent(locHist_Pi0_Both->GetMaximumBin());
		if(locOtherHeight > locMaxHeight)
			locMaxHeight = locOtherHeight;

		locHist_Pi0_2FCAL->SetTitle("Measured");
		locHist_Pi0_2FCAL->GetXaxis()->SetTitleSize(0.05);
		locHist_Pi0_2FCAL->GetYaxis()->SetTitle("");
		locHist_Pi0_2FCAL->GetXaxis()->SetLabelSize(0.05);
		locHist_Pi0_2FCAL->GetYaxis()->SetLabelSize(0.05);
		locHist_Pi0_2FCAL->SetLineColor(kRed);
		locHist_Pi0_2FCAL->SetLineWidth(2);
		locHist_Pi0_2FCAL->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_Pi0_2FCAL->Draw();

		locHist_Pi0_2BCAL->GetXaxis()->SetTitleSize(0.05);
		locHist_Pi0_2BCAL->GetYaxis()->SetTitle("");
		locHist_Pi0_2BCAL->GetXaxis()->SetLabelSize(0.05);
		locHist_Pi0_2BCAL->GetYaxis()->SetLabelSize(0.05);
		locHist_Pi0_2BCAL->SetLineColor(kBlue);
		locHist_Pi0_2BCAL->SetLineWidth(2);
		locHist_Pi0_2BCAL->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_Pi0_2BCAL->Draw("SAME");

		locHist_Pi0_Both->GetXaxis()->SetTitleSize(0.05);
		locHist_Pi0_Both->GetYaxis()->SetTitle("");
		locHist_Pi0_Both->GetXaxis()->SetLabelSize(0.05);
		locHist_Pi0_Both->GetYaxis()->SetLabelSize(0.05);
		locHist_Pi0_Both->SetLineColor(kBlack);
		locHist_Pi0_Both->SetLineWidth(2);
		locHist_Pi0_Both->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_Pi0_Both->Draw("SAME");

		TLegend *locLegend = new TLegend(0.14, 0.70, 0.39, 0.86); //botleft x/y, topright x/y
		locLegend->SetHeader("Legend");
		locLegend->AddEntry(locHist_Pi0_2FCAL, "2#gamma in FCAL", "F");
		locLegend->AddEntry(locHist_Pi0_2BCAL, "2#gamma in BCAL", "F");
		locLegend->AddEntry(locHist_Pi0_Both, "1 #gamma in Each", "F");
		locLegend->Draw();
	}

	//con lev
	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitConLev != NULL)
	{
		locHist_KinFitConLev->GetXaxis()->SetTitle("Confidence Level");
		locHist_KinFitConLev->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitConLev->GetYaxis()->SetTitle("");
		locHist_KinFitConLev->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitConLev->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitConLev->SetFillColor(kTeal + 1);
		locHist_KinFitConLev->Draw();
		gPad->SetLogy();
	}

	//mm2
	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_MM2 != NULL) && (locHist_MM2_KinFitCut != NULL))
	{
		locHist_MM2->Rebin(locNumRebin);
		locHist_MM2->GetXaxis()->SetRangeUser(-0.05, 0.05);
		locHist_MM2->GetXaxis()->SetTitleSize(0.045);
		locHist_MM2->GetYaxis()->SetTitle("");
		locHist_MM2->GetXaxis()->SetLabelSize(0.045);
		locHist_MM2->GetYaxis()->SetLabelSize(0.05);
		locHist_MM2->SetFillColor(kRed - 7);
		locHist_MM2->GetYaxis()->SetRangeUser(0.0, 1.05*locHist_MM2->GetBinContent(locHist_MM2->GetMaximumBin()));
		locHist_MM2->Draw();

		locHist_MM2_KinFitCut->Rebin(locNumRebin);
		locHist_MM2_KinFitCut->GetXaxis()->SetTitleSize(0.05);
		locHist_MM2_KinFitCut->GetYaxis()->SetTitle("");
		locHist_MM2_KinFitCut->GetXaxis()->SetLabelSize(0.05);
		locHist_MM2_KinFitCut->GetYaxis()->SetLabelSize(0.05);
		locHist_MM2_KinFitCut->SetFillColor(kAzure + 1);
		locHist_MM2_KinFitCut->Draw("SAME");

		TLegend *locLegend = new TLegend(0.14, 0.74, 0.39, 0.86); //botleft x/y, topright x/y
		locLegend->SetHeader("Legend");
		locLegend->AddEntry(locHist_MM2, "No KinFit Cut", "F");
		locLegend->AddEntry(locHist_MM2_KinFitCut, "KinFit Cut", "F");
		locLegend->Draw();
	}

	//pi0 before/after kinfit
	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_Pi0 != NULL) && (locHist_Pi0_KinFitCut != NULL))
	{
		locHist_Pi0->Rebin(locNumRebin);
		locHist_Pi0->GetXaxis()->SetRangeUser(-0.05, 0.05);
		locHist_Pi0->GetXaxis()->SetTitleSize(0.045);
		locHist_Pi0->GetYaxis()->SetTitle("");
		locHist_Pi0->GetXaxis()->SetLabelSize(0.045);
		locHist_Pi0->GetYaxis()->SetLabelSize(0.05);
		locHist_Pi0->SetFillColor(kRed - 7);
		locHist_Pi0->GetYaxis()->SetRangeUser(0.0, 1.05*locHist_Pi0->GetBinContent(locHist_Pi0->GetMaximumBin()));
		locHist_Pi0->Draw();

		locHist_Pi0_KinFitCut->Rebin(locNumRebin);
		locHist_Pi0_KinFitCut->GetXaxis()->SetTitleSize(0.05);
		locHist_Pi0_KinFitCut->GetYaxis()->SetTitle("");
		locHist_Pi0_KinFitCut->GetXaxis()->SetLabelSize(0.05);
		locHist_Pi0_KinFitCut->GetYaxis()->SetLabelSize(0.05);
		locHist_Pi0_KinFitCut->SetFillColor(kAzure + 1);
		locHist_Pi0_KinFitCut->Draw("SAME");

		TLegend *locLegend = new TLegend(0.14, 0.74, 0.39, 0.86); //botleft x/y, topright x/y
		locLegend->SetHeader("Legend");
		locLegend->AddEntry(locHist_Pi0, "No KinFit Cut", "F");
		locLegend->AddEntry(locHist_Pi0_KinFitCut, "KinFit Cut", "F");
		locLegend->Draw();
	}

	//compare omega kinfit / measured after kinfit cut
	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_KinFitOmega_KinFitCut != NULL) && (locHist_Omega_KinFitCut != NULL))
	{
		locHist_Omega_KinFitCut->Rebin(2*locNumRebin);
		locHist_KinFitOmega_KinFitCut->Rebin(2*locNumRebin);

		double locMaxHeight = locHist_Omega_KinFitCut->GetBinContent(locHist_Omega_KinFitCut->GetMaximumBin());
		double locOtherHeight = locHist_KinFitOmega_KinFitCut->GetBinContent(locHist_KinFitOmega_KinFitCut->GetMaximumBin());
		if(locOtherHeight > locMaxHeight)
			locMaxHeight = locOtherHeight;

		locHist_Omega_KinFitCut->SetTitle("Post KinFit Cut");
		locHist_Omega_KinFitCut->GetXaxis()->SetTitleSize(0.05);
		locHist_Omega_KinFitCut->GetYaxis()->SetTitle("");
		locHist_Omega_KinFitCut->GetXaxis()->SetLabelSize(0.05);
		locHist_Omega_KinFitCut->GetYaxis()->SetLabelSize(0.05);
		locHist_Omega_KinFitCut->SetLineColor(kRed);
		locHist_Omega_KinFitCut->SetLineWidth(2);
		locHist_Omega_KinFitCut->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_Omega_KinFitCut->Draw("SAME");

		locHist_KinFitOmega_KinFitCut->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitOmega_KinFitCut->GetYaxis()->SetTitle("");
		locHist_KinFitOmega_KinFitCut->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitOmega_KinFitCut->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitOmega_KinFitCut->SetLineColor(kBlue);
		locHist_KinFitOmega_KinFitCut->SetLineWidth(2);
		locHist_KinFitOmega_KinFitCut->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_KinFitOmega_KinFitCut->Draw("SAME");

		TLegend *locLegend = new TLegend(0.14, 0.74, 0.39, 0.86); //botleft x/y, topright x/y
		locLegend->SetHeader("Legend");
		locLegend->AddEntry(locHist_Omega_KinFitCut, "Measured", "F");
		locLegend->AddEntry(locHist_KinFitOmega_KinFitCut, "KinFit", "F");
		locLegend->Draw();

		n_omega_kinfit = locHist_Omega_KinFitCut->Integral(100./locNumRebin, 400./locNumRebin);
		
		// Determine Mass and Width
		TF1 *fomega = new TF1("fomega", "gaus", 0.6, 0.9);
		fomega->SetParameter(1,0.782);
		fomega->SetParameter(2,0.03);
		locHist_KinFitOmega_KinFitCut->Fit("fomega", "RQ0");
		omega_mass = fomega->GetParameter(1);
		omega_width = fomega->GetParameter(2);

	}

	//measured omega for 3 cases (on top of each other)
	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_Omega_2FCAL != NULL) && (locHist_Omega_2BCAL != NULL) && (locHist_Omega_Both != NULL))
	{
		locHist_Omega_2FCAL->Rebin(2*locNumRebin);
		locHist_Omega_2BCAL->Rebin(2*locNumRebin);
		locHist_Omega_Both->Rebin(2*locNumRebin);

		double locMaxHeight = locHist_Omega_2FCAL->GetBinContent(locHist_Omega_2FCAL->GetMaximumBin());
		double locOtherHeight = locHist_Omega_2BCAL->GetBinContent(locHist_Omega_2BCAL->GetMaximumBin());
		if(locOtherHeight > locMaxHeight)
			locMaxHeight = locOtherHeight;
		locOtherHeight = locHist_Omega_Both->GetBinContent(locHist_Omega_Both->GetMaximumBin());
		if(locOtherHeight > locMaxHeight)
			locMaxHeight = locOtherHeight;

		locHist_Omega_2FCAL->SetTitle("Post KinFit Cut");
		locHist_Omega_2FCAL->GetXaxis()->SetTitleSize(0.05);
		locHist_Omega_2FCAL->GetYaxis()->SetTitle("");
		locHist_Omega_2FCAL->GetXaxis()->SetLabelSize(0.05);
		locHist_Omega_2FCAL->GetYaxis()->SetLabelSize(0.05);
		locHist_Omega_2FCAL->SetLineColor(kRed);
		locHist_Omega_2FCAL->SetLineWidth(2);
		locHist_Omega_2FCAL->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_Omega_2FCAL->Draw();

		locHist_Omega_2BCAL->GetXaxis()->SetTitleSize(0.05);
		locHist_Omega_2BCAL->GetYaxis()->SetTitle("");
		locHist_Omega_2BCAL->GetXaxis()->SetLabelSize(0.05);
		locHist_Omega_2BCAL->GetYaxis()->SetLabelSize(0.05);
		locHist_Omega_2BCAL->SetLineColor(kBlue);
		locHist_Omega_2BCAL->SetLineWidth(2);
		locHist_Omega_2BCAL->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_Omega_2BCAL->Draw("SAME");

		locHist_Omega_Both->GetXaxis()->SetTitleSize(0.05);
		locHist_Omega_Both->GetYaxis()->SetTitle("");
		locHist_Omega_Both->GetXaxis()->SetLabelSize(0.05);
		locHist_Omega_Both->GetYaxis()->SetLabelSize(0.05);
		locHist_Omega_Both->SetLineColor(kBlack);
		locHist_Omega_Both->SetLineWidth(2);
		locHist_Omega_Both->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_Omega_Both->Draw("SAME");

		TLegend *locLegend = new TLegend(0.14, 0.70, 0.39, 0.86); //botleft x/y, topright x/y
		locLegend->SetHeader("Legend");
		locLegend->AddEntry(locHist_Omega_2FCAL, "2#gamma in FCAL", "F");
		locLegend->AddEntry(locHist_Omega_2BCAL, "2#gamma in BCAL", "F");
		locLegend->AddEntry(locHist_Omega_Both, "1 #gamma in Each", "F");
		locLegend->Draw();
	}

	// Print the mass, width and number of reconstructed omegas per trigger (Fitted values)
	locCanvas->cd(2);
	if(locHist_KinFitConLev != NULL){
	  TLatex tx;
	  tx.SetTextAlign(11);
	  tx.SetTextSize(0.06);
	  char text[100];
	  sprintf(text, "E_{#gamma} > 7 GeV");
	  tx.DrawLatex(0.1, locHist_KinFitConLev->GetMaximum()/4., text);
	  sprintf(text, "Post KinFit");
	  tx.DrawLatex(0.1, locHist_KinFitConLev->GetMaximum()/16., text);
	  sprintf(text, "M(#omega) = %0.3f GeV/c^{2}", omega_mass);
	  tx.DrawLatex(0.1, locHist_KinFitConLev->GetMaximum()/64., text);
	  sprintf(text, "#Gamma(#omega) = %0.3f GeV/c^{2}", omega_width);
	  tx.DrawLatex(0.1, locHist_KinFitConLev->GetMaximum()/256., text);
	  sprintf(text, "N(#omega) = %0.2f / 1k Trigger", n_omega_kinfit/n_triggers*1000);
	  tx.DrawLatex(0.1,  locHist_KinFitConLev->GetMaximum()/1024., text);
	}
}


// hnamepath: /p2pi_preco_kinfit/Hist_MissingMassSquared/MissingMassSquared
// hnamepath: /p2pi_preco_kinfit/Hist_MissingMassSquared_PostKinFitCut/MissingMassSquared
// hnamepath: /p2pi_preco_kinfit/HHist_KinFitResults/ConfidenceLevel
// hnamepath: /p2pi_preco/Hist_InvariantMass_Rho/InvariantMass
// hnamepath: /p2pi_preco_kinfit/Hist_InvariantMass_Rho_PostKinFitCut/InvariantMass
// hnamepath: /p2pi_preco_kinfit/Hist_InvariantMass_Rho_KinFit_PostKinFitCut/InvariantMass
// hnamepath: /p2pi_preco/Custom_p2pi_hists/PiPlusPsi_t

{
	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locReactionDirectory_NoKinFit = (TDirectory*)locInitDirectory->FindObjectAny("p2pi_preco");
	TDirectory *locReactionDirectory_KinFit = (TDirectory*)locInitDirectory->FindObjectAny("p2pi_preco_kinfit");
	if((locReactionDirectory_NoKinFit == NULL) || (locReactionDirectory_KinFit == NULL))
		return;

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("p2pi", "p2pi", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);


	TH1I* locHist_MM2 = (TH1I*)locReactionDirectory_KinFit->Get("Hist_MissingMassSquared/MissingMassSquared");
	TH1I* locHist_MM2_KinFitCut = (TH1I*)locReactionDirectory_KinFit->Get("Hist_MissingMassSquared_PostKinFitCut/MissingMassSquared");

	TH1I* locHist_KinFitConLev = (TH1I*)locReactionDirectory_KinFit->Get("Hist_KinFitResults/ConfidenceLevel");

	TH1I* locHist_RhoMass_MMCut = (TH1I*)locReactionDirectory_NoKinFit->Get("Hist_InvariantMass_Rho/InvariantMass");
	TH1I* locHist_RhoMass_KinFitCut = (TH1I*)locReactionDirectory_KinFit->Get("Hist_InvariantMass_Rho_PostKinFitCut/InvariantMass");
	TH1I* locHist_KinFitRhoMass_KinFitCut = (TH1I*)locReactionDirectory_KinFit->Get("Hist_InvariantMass_Rho_KinFit_PostKinFitCut/InvariantMass");

	TH2I* locHist_PiPlusPsi_t = (TH2I*)locReactionDirectory_NoKinFit->Get("Custom_p2pi_hists/PiPlusPsi_t");

	//Draw
	int locNumRebin = 4;

	//rho: mm cut
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RhoMass_MMCut != NULL)
	{
		locHist_RhoMass_MMCut->SetTitle("With MM Cut");
		locHist_RhoMass_MMCut->Rebin(locNumRebin);
		locHist_RhoMass_MMCut->GetYaxis()->SetRangeUser(0.0, 1.05*locHist_RhoMass_MMCut->GetBinContent(locHist_RhoMass_MMCut->GetMaximumBin()));
		locHist_RhoMass_MMCut->GetXaxis()->SetTitleSize(0.05);
		locHist_RhoMass_MMCut->GetYaxis()->SetTitle("");
		locHist_RhoMass_MMCut->GetXaxis()->SetLabelSize(0.05);
		locHist_RhoMass_MMCut->GetYaxis()->SetLabelSize(0.05);
		locHist_RhoMass_MMCut->SetFillColor(kAzure + 1);
		locHist_RhoMass_MMCut->Draw("");
	}

	//psi
	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PiPlusPsi_t != NULL)
	{
		locHist_PiPlusPsi_t->SetTitle("#psi_{#pi^{+}} vs E_{#gamma}: Proton dE/dx > 2.2; E_{#gamma}; #psi_{#pi^{+}}");
		TH1D *locHist_TimingCut_PiPlusPsi = (TH1D*)locHist_PiPlusPsi_t->ProjectionY();
		locHist_TimingCut_PiPlusPsi->Rebin(4);
		locHist_TimingCut_PiPlusPsi->SetMinimum(0);
		locHist_TimingCut_PiPlusPsi->GetXaxis()->SetTitleSize(0.05);
		locHist_TimingCut_PiPlusPsi->GetXaxis()->SetLabelSize(0.05);
		locHist_TimingCut_PiPlusPsi->Draw();

		// fit 1+cos(2*phi) distribution
		TF1* fit = new TF1("psiFit","[0]*(1.0 + [1]*cos(2*(x + [2])/180.*3.14159))");
		locHist_TimingCut_PiPlusPsi->Fit(fit, "Q", "");
		locHist_TimingCut_PiPlusPsi->Draw("e");

		// print fit parameters to canvas
		Double_t PSigma = fit->GetParameter(1);
		Double_t PSigmaErr = fit->GetParError(1);
		Double_t Phi0 = fit->GetParameter(2);
		Double_t Phi0Err = fit->GetParError(2);
		TLatex tx;
		tx.SetTextAlign(21);
		tx.SetTextSize(0.05);
		char text[100];
		sprintf(text, "P#Sigma=%0.2f#pm%0.2f", PSigma, PSigmaErr);
		tx.DrawLatex(0., locHist_TimingCut_PiPlusPsi->GetMaximum()*0.15, text);
		sprintf(text, "#phi_{0}=%0.2f#pm%0.2f", Phi0, Phi0Err);
		tx.DrawLatex(0., locHist_TimingCut_PiPlusPsi->GetMaximum()*0.06, text);

		fit->Delete();
	}

	//con lev
	locCanvas->cd(3);
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

	//mm
	locCanvas->cd(4);
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

	//compare rho kinfit / measured after kinfit cut
	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_KinFitRhoMass_KinFitCut != NULL) && (locHist_RhoMass_KinFitCut != NULL))
	{
		locHist_RhoMass_KinFitCut->Rebin(locNumRebin);
		locHist_KinFitRhoMass_KinFitCut->Rebin(locNumRebin);

		double locMaxHeight = locHist_RhoMass_KinFitCut->GetBinContent(locHist_RhoMass_KinFitCut->GetMaximumBin());
		double locOtherHeight = locHist_KinFitRhoMass_KinFitCut->GetBinContent(locHist_KinFitRhoMass_KinFitCut->GetMaximumBin());
		if(locOtherHeight > locMaxHeight)
			locMaxHeight = locOtherHeight;

		locHist_RhoMass_KinFitCut->SetTitle("Post KinFit Cut");
		locHist_RhoMass_KinFitCut->GetXaxis()->SetTitleSize(0.05);
		locHist_RhoMass_KinFitCut->GetYaxis()->SetTitle("");
		locHist_RhoMass_KinFitCut->GetXaxis()->SetLabelSize(0.05);
		locHist_RhoMass_KinFitCut->GetYaxis()->SetLabelSize(0.05);
		locHist_RhoMass_KinFitCut->SetLineColor(kRed);
		locHist_RhoMass_KinFitCut->SetLineWidth(2);
		locHist_RhoMass_KinFitCut->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_RhoMass_KinFitCut->Draw("SAME");

		locHist_KinFitRhoMass_KinFitCut->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitRhoMass_KinFitCut->GetYaxis()->SetTitle("");
		locHist_KinFitRhoMass_KinFitCut->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitRhoMass_KinFitCut->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitRhoMass_KinFitCut->SetLineColor(kBlue);
		locHist_KinFitRhoMass_KinFitCut->SetLineWidth(2);
		locHist_KinFitRhoMass_KinFitCut->GetYaxis()->SetRangeUser(0.0, 1.05*locMaxHeight);
		locHist_KinFitRhoMass_KinFitCut->Draw("SAME");

		TLegend *locLegend = new TLegend(0.14, 0.74, 0.39, 0.86); //botleft x/y, topright x/y
		locLegend->SetHeader("Legend");
		locLegend->AddEntry(locHist_RhoMass_KinFitCut, "Measured", "F");
		locLegend->AddEntry(locHist_KinFitRhoMass_KinFitCut, "KinFit", "F");
		locLegend->Draw();
	}
}


// hnamepath: /p2pi_pmiss/Custom_p2pi_hists_NoCut_Measured/MM_M2pi
// hnamepath: /p2pi_pmiss/Custom_p2pi_hists_KinConverge_Measured/MM_M2pi
// hnamepath: /p2pi_pmiss/Custom_p2pi_hists_KinConverge_Measured/Egamma

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory;
	if((TDirectory*)locTopDirectory->FindObjectAny("p2pi_pmiss") != 0)
	  locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("p2pi_pmiss");
	else
	  return;

	//Go to NoCut directory
	locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_NoCut_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_NoCut_MM_M2pi = (TH2I*)gDirectory->Get("MM_M2pi");

	//Go to KinFitConverge directory
	locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_KinFitConverge_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_KinConverge_MM_M2pi = (TH2I*)gDirectory->Get("MM_M2pi");
	TH1I* locHist_KinConverge_Egamma = (TH1I*)gDirectory->Get("Egamma");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("p2pi_pmiss", "p2pi_pmiss", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoCut_MM_M2pi != NULL) {	
		locHist_NoCut_MM_M2pi->Rebin2D();
		locHist_NoCut_MM_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoCut_MM_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_NoCut_MM_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoCut_MM_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoCut_MM_M2pi->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinConverge_MM_M2pi != NULL) {
		locHist_KinConverge_MM_M2pi->Rebin2D();
		locHist_KinConverge_MM_M2pi->SetTitle("MM off #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}: Timing cut and KinFit converges");
		locHist_KinConverge_MM_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinConverge_MM_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_KinConverge_MM_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinConverge_MM_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinConverge_MM_M2pi->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinConverge_MM_M2pi != NULL) {
		
		TH1I* locHist_KinConverge_M2pi = (TH1I*)locHist_KinConverge_MM_M2pi->ProjectionX()->Clone();
		locHist_KinConverge_M2pi->SetTitle("M_{#pi^{+}#pi^{-}}: Timing cut and KinFit converges");
		locHist_KinConverge_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinConverge_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinConverge_M2pi->GetYaxis()->SetTitleSize(0.05);
                locHist_KinConverge_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinConverge_M2pi->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinConverge_Egamma != NULL) {
		
		locHist_KinConverge_Egamma->GetXaxis()->SetTitleSize(0.05);
		locHist_KinConverge_Egamma->GetXaxis()->SetLabelSize(0.05);
		locHist_KinConverge_Egamma->GetYaxis()->SetTitleSize(0.05);
                locHist_KinConverge_Egamma->GetYaxis()->SetLabelSize(0.05);
		locHist_KinConverge_Egamma->Draw();
	}

}


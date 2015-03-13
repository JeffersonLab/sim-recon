// hnamepath: /p2pi_pmiss/Custom_p2pi_hists_NoKinFit_Measured/MM_M2pi
// hnamepath: /p2pi_pmiss/Custom_p2pi_hists_KinFitCut10_Measured/MM_M2pi
// hnamepath: /p2pi_pmiss/Custom_p2pi_hists_KinFitCut10_Measured/Egamma

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("p2pi_pmiss");

	//Go to NoKinFit directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_NoKinFit_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_NoKinFit_MM_M2pi = (TH2I*)gDirectory->Get("MM_M2pi");

	//Go to KinFitCut10 directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_KinFitCut10_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_KinFitCut10_MM_M2pi = (TH2I*)gDirectory->Get("MM_M2pi");
	TH1I* locHist_KinFitCut10_Egamma = (TH1I*)gDirectory->Get("Egamma");

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
	if(locHist_NoKinFit_MM_M2pi != NULL) {	
		locHist_NoKinFit_MM_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MM_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MM_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MM_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MM_M2pi->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_MM_M2pi != NULL) {
		locHist_KinFitCut10_MM_M2pi->SetTitle("MM off #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}: KinFit CL > 0.1");
		locHist_KinFitCut10_MM_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_MM_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_MM_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_MM_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_MM_M2pi->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_MM_M2pi != NULL) {
		
		TH1I* locHist_KinFitCut10_M2pi = (TH1I*)locHist_KinFitCut10_MM_M2pi->ProjectionX()->Clone();
		locHist_KinFitCut10_M2pi->SetTitle("M_{#pi^{+}#pi^{-}}: KinFit CL > 0.1");
		locHist_KinFitCut10_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_M2pi->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_Egamma != NULL) {
		
		locHist_KinFitCut10_Egamma->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_Egamma->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_Egamma->Draw();
	}

}


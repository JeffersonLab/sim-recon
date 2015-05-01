// hnamepath: /p2pi_preco/Custom_p2pi_hists_NoKinFit_Measured/MM2_M2pi
// hnamepath: /p2pi_preco/Custom_p2pi_hists_KinFitConverge_Measured/MM2_M2pi
// hnamepath: /p2pi_preco/Custom_p2pi_hists_KinFitConverge_Measured/Egamma

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory;
	if((TDirectory*)locTopDirectory->FindObjectAny("p2pi_preco") != 0)
	  locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("p2pi_preco");
	else
	  return;

	//Go to NoKinFit directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_NoKinFit_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_NoKinFit_MM2_M2pi = (TH2I*)gDirectory->Get("MM2_M2pi");

	//Go to KinFitConverge directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_KinFitConverge_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_KinFitConverge_MM2_M2pi = (TH2I*)gDirectory->Get("MM2_M2pi");
	TH1I* locHist_KinFitConverge_M2pi = (TH1I*)locHist_KinFitConverge_MM2_M2pi->ProjectionX("temp",96,105)->Clone();
	TH1I* locHist_KinFitConverge_Egamma = (TH1I*)gDirectory->Get("Egamma");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("p2pi_preco1", "p2pi_preco1", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoKinFit_MM2_M2pi != NULL) {	
		locHist_NoKinFit_MM2_M2pi->Rebin2D();
		locHist_NoKinFit_MM2_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MM2_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MM2_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MM2_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MM2_M2pi->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitConverge_MM2_M2pi != NULL) {
		locHist_KinFitConverge_MM2_M2pi->Rebin2D();
		locHist_KinFitConverge_MM2_M2pi->SetTitle("MM^{2} off p#pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}: KinFit Converge");
		locHist_KinFitConverge_MM2_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitConverge_MM2_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_KinFitConverge_MM2_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitConverge_MM2_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitConverge_MM2_M2pi->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitConverge_M2pi != NULL) {
		
		locHist_KinFitConverge_MM2_M2pi->SetTitle("M_{#pi^{+}#pi^{-}}: KinFit Converge |MM^2| < 0.05");
		locHist_KinFitConverge_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitConverge_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitConverge_M2pi->GetYaxis()->SetTitleSize(0.05);
                locHist_KinFitConverge_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitConverge_M2pi->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitConverge_Egamma != NULL) {
		
		locHist_KinFitConverge_Egamma->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitConverge_Egamma->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitConverge_Egamma->GetYaxis()->SetTitleSize(0.05);
                locHist_KinFitConverge_Egamma->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitConverge_Egamma->Draw();
	}

}


// hnamepath: /p3pi_preco_FCAL-BCAL/Custom_p3pi_hists_NoKinFit_Measured/MM2_M3pi
// hnamepath: /p3pi_preco_FCAL-BCAL/Custom_p3pi_hists_KinFitCut10_Measured/MM2_M3pi
// hnamepath: /p3pi_preco_FCAL-BCAL/Hist_InvariantMass_NoKinFit_Measured/InvariantMass

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("p3pi_preco_FCAL-BCAL");

	//Go to Pi0 mass directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Hist_InvariantMass_NoKinFit_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH1I* locHist_NoKinFit_MPi0 = (TH1I*)gDirectory->Get("InvariantMass");

	//Go to NoKinFit directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p3pi_hists_NoKinFit_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_NoKinFit_MM2_M3pi = (TH2I*)gDirectory->Get("MM2_M3pi");

	//Go to KinFitCut10 directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p3pi_hists_KinFitCut10_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_KinFitCut10_MM2_M3pi = (TH2I*)gDirectory->Get("MM2_M3pi");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("p3pi_preco", "p3pi_preco", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoKinFit_MPi0 != NULL) {
		locHist_NoKinFit_MPi0->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MPi0->GetYaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MPi0->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MPi0->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MPi0->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoKinFit_MM2_M3pi != NULL) {	
		locHist_NoKinFit_MM2_M3pi->Rebin2D();
		locHist_NoKinFit_MM2_M3pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MM2_M3pi->GetYaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MM2_M3pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MM2_M3pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MM2_M3pi->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_MM2_M3pi != NULL) {
		locHist_KinFitCut10_MM2_M3pi->Rebin2D();
		locHist_KinFitCut10_MM2_M3pi->SetTitle("MM^{2} off #pi^{+}#pi^{-}#pi^{0} vs M_{#pi^{+}#pi^{-}#pi^{0}}: KinFit CL > 0.1");
		locHist_KinFitCut10_MM2_M3pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_MM2_M3pi->GetYaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_MM2_M3pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_MM2_M3pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_MM2_M3pi->Draw("colz");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_MM2_M3pi != NULL) {
		
		TH1I* locHist_KinFitCut10_M3pi = (TH1I*)locHist_KinFitCut10_MM2_M3pi->ProjectionX()->Clone();
		locHist_KinFitCut10_M3pi->SetTitle("M_{#pi^{+}#pi^{-}#pi^{0}}: KinFit CL > 0.1");
		locHist_KinFitCut10_M3pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_M3pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_M3pi->GetYaxis()->SetTitleSize(0.05);
                locHist_KinFitCut10_M3pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_M3pi->Draw();
	}
}


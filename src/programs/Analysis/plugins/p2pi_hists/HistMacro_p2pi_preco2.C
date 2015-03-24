// hnamepath: /p2pi_preco/Custom_p2pi_hists_NoKinFit_Measured/MM2_M2pi
// hnamepath: /p2pi_preco/Custom_p2pi_hists_KinFitCut10_Measured/MM2_M2pi
// hnamepath: /p2pi_preco/Custom_p2pi_hists_KinFitCut10_Measured/Egamma

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory;
	if((TDirectory*)locTopDirectory->FindObjectAny("p2pi_preco") != 0)
	  locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("p2pi_preco");
	else
	  return;

	//Go to KinFitCut10 directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_KinFitCut10_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_KinFitCut10_Proton_dEdx_P = (TH2I*)gDirectory->Get("Proton_dEdx_P");
	TH2I* locHist_KinFitCut10_Proton_P_Theta = (TH2I*)gDirectory->Get("Proton_P_Theta");
	TH2I* locHist_KinFitCut10_DeltaE_M2pi = (TH2I*)gDirectory->Get("DeltaE_M2pi");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("p2pi_preco2", "p2pi_preco2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_Proton_dEdx_P != NULL) {
		locHist_KinFitCut10_Proton_dEdx_P->Rebin2D();	
		locHist_KinFitCut10_Proton_dEdx_P->SetTitle("dE/dx vs p: KinFit CL > 0.1");
		locHist_KinFitCut10_Proton_dEdx_P->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_Proton_dEdx_P->GetYaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_Proton_dEdx_P->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_Proton_dEdx_P->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_Proton_dEdx_P->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_Proton_P_Theta != NULL) {
		locHist_KinFitCut10_Proton_P_Theta->SetTitle("p vs. #theta: KinFit CL > 0.1");
		locHist_KinFitCut10_Proton_P_Theta->GetYaxis()->SetRangeUser(0.,5.);
		locHist_KinFitCut10_Proton_P_Theta->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_Proton_P_Theta->GetYaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_Proton_P_Theta->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_Proton_P_Theta->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_Proton_P_Theta->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_DeltaE_M2pi != NULL) {
		locHist_KinFitCut10_DeltaE_M2pi->Rebin2D();
		locHist_KinFitCut10_DeltaE_M2pi->SetTitle("#DeltaE vs M_{#pi^{+}#pi^{-}}: KinFit CL > 0.1");
		locHist_KinFitCut10_DeltaE_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_DeltaE_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_DeltaE_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_DeltaE_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_DeltaE_M2pi->Draw("colz");	
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_DeltaE_M2pi != NULL) {
		
		TH1I* locHist_KinFitCut10_DeltaE = (TH1I*)locHist_KinFitCut10_DeltaE_M2pi->ProjectionY();
		locHist_KinFitCut10_DeltaE->SetTitle("#DeltaE: KinFit CL > 0.1");
		locHist_KinFitCut10_DeltaE->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_DeltaE->GetYaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_DeltaE->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_DeltaE->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_DeltaE->Draw();	
	}

}


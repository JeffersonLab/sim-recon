// hnamepath: /p3pi_preco_2FCAL/Custom_p3pi_hists_NoKinFit_Measured/MM2_M3pi
// hnamepath: /p3pi_preco_2FCAL/Custom_p3pi_hists_KinFitCut10_Measured/MM2_M3pi
// hnamepath: /p3pi_preco_2FCAL/Hist_InvariantMass_NoKinFit_Measured/InvariantMass

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory;
	if((TDirectory*)locTopDirectory->FindObjectAny("p3pi_preco_2FCAL") != 0)
	  locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("p3pi_preco_2FCAL");
	else
	  return;

	//Go to Pi0 mass directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Hist_InvariantMass_NoKinFit_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH1I* locHist_NoKinFit_MPi0 = (TH1I*)gDirectory->Get("InvariantMass");

	//Go to NoKinFit directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p3pi_hists_CutPi0_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_NoKinFit_MM2_M3pi = (TH2I*)gDirectory->Get("MM2_M3pi");
	TH2I* locHist_NoKinFit_Proton_dEdx_P = (TH2I*)gDirectory->Get("Proton_dEdx_P");
	TH2I* locHist_NoKinFit_Egamma_M3pi = (TH2I*)gDirectory->Get("dEgamma_M3pi_ProtonTag");
	TH1I* locHist_NoKinFit_M3pi = (TH1I*)locHist_NoKinFit_Egamma_M3pi->ProjectionX("M3pi");

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
		locHist_NoKinFit_MPi0->GetXaxis()->SetRangeUser(0.0, 0.6);
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
		locHist_NoKinFit_MM2_M3pi->SetTitle("MM^{2} off #pi^{+}#pi^{-}#pi^{0} vs M_{#pi^{+}#pi^{-}#pi^{0}}");
		locHist_NoKinFit_MM2_M3pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MM2_M3pi->GetYaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_MM2_M3pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MM2_M3pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_MM2_M3pi->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoKinFit_Proton_dEdx_P != NULL) {
		locHist_NoKinFit_Proton_dEdx_P->Rebin2D();
		locHist_NoKinFit_Proton_dEdx_P->SetTitle("Proton dE/dx vs momentum: |MM^{2}|<0.05 and |Missing Energy| < 0.5");
		locHist_NoKinFit_Proton_dEdx_P->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_Proton_dEdx_P->GetYaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_Proton_dEdx_P->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_Proton_dEdx_P->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_Proton_dEdx_P->Draw("colz");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoKinFit_M3pi != NULL) {
		
		locHist_NoKinFit_M3pi->SetTitle("M_{#pi^{+}#pi^{-}#pi^{0}}: Proton dE/dx > 2.2");
		locHist_NoKinFit_M3pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_M3pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_M3pi->GetYaxis()->SetTitleSize(0.05);
                locHist_NoKinFit_M3pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_M3pi->Draw();
	}
}


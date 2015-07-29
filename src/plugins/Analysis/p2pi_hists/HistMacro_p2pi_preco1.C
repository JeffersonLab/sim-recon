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
	TH2I* locHist_NoKinFit_Proton_dEdx_P = (TH2I*)gDirectory->Get("Proton_dEdx_P");
	TH2I* locHist_NoKinFit_Egamma_M2pi = (TH2I*)gDirectory->Get("Egamma_M2pi");
	TH2I* locHist_NoKinFit_PiPlusPsi_Egamma = (TH2I*)gDirectory->Get("PiPlusPsi_Egamma");

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
	if(locHist_NoKinFit_Proton_dEdx_P != NULL) {
		locHist_NoKinFit_Proton_dEdx_P->Rebin2D();
		locHist_NoKinFit_Proton_dEdx_P->SetTitle("Proton dE/dx vs momentum: |MM^{2}|<0.02 and |Missing Energy| < 0.2");
		locHist_NoKinFit_Proton_dEdx_P->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_Proton_dEdx_P->GetYaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_Proton_dEdx_P->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_Proton_dEdx_P->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_Proton_dEdx_P->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoKinFit_Egamma_M2pi != NULL) {
		
		locHist_NoKinFit_Egamma_M2pi->SetTitle("E_{#gamma} vs M_{#pi^{+}#pi^{-}}: Proton dE/dx > 2.2");
		locHist_NoKinFit_Egamma_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_Egamma_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_Egamma_M2pi->GetYaxis()->SetTitleSize(0.05);
                locHist_NoKinFit_Egamma_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_Egamma_M2pi->Draw("colz");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoKinFit_PiPlusPsi_Egamma != NULL) {

		locHist_NoKinFit_PiPlusPsi_Egamma->SetTitle("#psi_{#pi^{+}} vs E_{#gamma}: Proton dE/dx > 2.2; E_{#gamma}; #psi_{#pi^{+}}");
		locHist_NoKinFit_PiPlusPsi_Egamma->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_PiPlusPsi_Egamma->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_PiPlusPsi_Egamma->GetYaxis()->SetTitleSize(0.05);
                locHist_NoKinFit_PiPlusPsi_Egamma->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_PiPlusPsi_Egamma->Draw("colz");
	}
}


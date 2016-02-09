// hnamepath: /p2pi_preco/Custom_p2pi_hists_TimingCut_Measured/MM2_M2pi

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory;
	if((TDirectory*)locTopDirectory->FindObjectAny("p2pi_preco") != 0)
	  locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("p2pi_preco");
	else
	  return;

	//Go to TimingCut directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_TimingCut_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_TimingCut_MM2_M2pi = (TH2I*)gDirectory->Get("MM2_M2pi");
	TH2I* locHist_TimingCut_Proton_dEdx_P = (TH2I*)gDirectory->Get("Proton_dEdx_P");
	TH2I* locHist_TimingCut_Egamma_M2pi = (TH2I*)gDirectory->Get("Egamma_M2pi");
	TH2I* locHist_TimingCut_PiPlusPsi_t = (TH2I*)gDirectory->Get("PiPlusPsi_t");

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
	if(locHist_TimingCut_MM2_M2pi != NULL) {	
		locHist_TimingCut_MM2_M2pi->Rebin2D();
		locHist_TimingCut_MM2_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_TimingCut_MM2_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_TimingCut_MM2_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_TimingCut_MM2_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_TimingCut_MM2_M2pi->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TimingCut_Proton_dEdx_P != NULL) {
		locHist_TimingCut_Proton_dEdx_P->Rebin2D();
		locHist_TimingCut_Proton_dEdx_P->SetTitle("Proton dE/dx vs momentum: |MM^{2}|<0.02 and |Missing Energy| < 0.2");
		locHist_TimingCut_Proton_dEdx_P->GetXaxis()->SetTitleSize(0.05);
		locHist_TimingCut_Proton_dEdx_P->GetYaxis()->SetTitleSize(0.05);
		locHist_TimingCut_Proton_dEdx_P->GetXaxis()->SetLabelSize(0.05);
		locHist_TimingCut_Proton_dEdx_P->GetYaxis()->SetLabelSize(0.05);
		locHist_TimingCut_Proton_dEdx_P->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TimingCut_Egamma_M2pi != NULL) {
		
		locHist_TimingCut_Egamma_M2pi->SetTitle("E_{#gamma} vs M_{#pi^{+}#pi^{-}}: Proton dE/dx > 2.2");
		locHist_TimingCut_Egamma_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_TimingCut_Egamma_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_TimingCut_Egamma_M2pi->GetYaxis()->SetTitleSize(0.05);
                locHist_TimingCut_Egamma_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_TimingCut_Egamma_M2pi->Draw("colz");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TimingCut_PiPlusPsi_t != NULL) {

		locHist_TimingCut_PiPlusPsi_t->SetTitle("#psi_{#pi^{+}} vs E_{#gamma}: Proton dE/dx > 2.2; E_{#gamma}; #psi_{#pi^{+}}");
		TH1D *locHist_TimingCut_PiPlusPsi = (TH1D*)locHist_TimingCut_PiPlusPsi_t->ProjectionY();
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
}


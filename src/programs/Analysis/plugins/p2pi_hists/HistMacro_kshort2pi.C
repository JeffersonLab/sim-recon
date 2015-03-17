// hnamepath: /kshort2pi/Custom_p2pi_hists_NoKinFit_Measured/MM_M2pi
// hnamepath: /kshort2pi/Custom_p2pi_hists_KinFitCut10_Measured/MM_M2pi

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("kshort2pi");

	//Go to NoKinFit directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_NoKinFit_Measured");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_NoKinFit_DeltaZ_M2pi = (TH2I*)gDirectory->Get("DeltaZ_M2pi");

	//Go to KinFitCut10 directory
	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Custom_p2pi_hists_KinFitCut10_KinFit");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH2I* locHist_KinFitCut10_DeltaZ_M2pi = (TH2I*)gDirectory->Get("DeltaZ_M2pi");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("kshort2pi", "kshort2pi", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();
	if(locHist_NoKinFit_DeltaZ_M2pi != NULL) {
		locHist_NoKinFit_DeltaZ_M2pi->GetXaxis()->SetRangeUser(0.,1.);
		locHist_NoKinFit_DeltaZ_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_DeltaZ_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_DeltaZ_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_DeltaZ_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_DeltaZ_M2pi->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();
	if(locHist_KinFitCut10_DeltaZ_M2pi != NULL) {
		locHist_KinFitCut10_DeltaZ_M2pi->SetTitle("Vertex #DeltaZ of #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}: KinFit CL > 0.1");
		locHist_KinFitCut10_DeltaZ_M2pi->GetXaxis()->SetRangeUser(0.,1.);
		locHist_KinFitCut10_DeltaZ_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_DeltaZ_M2pi->GetYaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_DeltaZ_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_DeltaZ_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_DeltaZ_M2pi->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NoKinFit_DeltaZ_M2pi != NULL) {
		
		int binMinZ = locHist_NoKinFit_DeltaZ_M2pi->GetYaxis()->FindBin(5.);

		TH1I* locHist_NoKinFit_M2pi = (TH1I*)locHist_NoKinFit_DeltaZ_M2pi->ProjectionX("NoKinFit_DeltaZ_GT5",binMinZ,1000)->Clone();
		locHist_NoKinFit_M2pi->SetTitle("M_{#pi^{+}#pi^{-}}: #DeltaZ > 5 cm");
		locHist_NoKinFit_M2pi->GetXaxis()->SetRangeUser(0.,1.);
		locHist_NoKinFit_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_NoKinFit_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_M2pi->GetYaxis()->SetTitleSize(0.05);
                locHist_NoKinFit_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_NoKinFit_M2pi->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_KinFitCut10_DeltaZ_M2pi != NULL) {
		
		int binMinZ = locHist_KinFitCut10_DeltaZ_M2pi->GetYaxis()->FindBin(5.);

		TH1I* locHist_KinFitCut10_M2pi = (TH1I*)locHist_KinFitCut10_DeltaZ_M2pi->ProjectionX("KinFitCut10_DeltaZ_GT5",binMinZ,1000)->Clone();
		locHist_KinFitCut10_M2pi->SetTitle("M_{#pi^{+}#pi^{-}}: KinFit CL > 0.1 and #DeltaZ > 5 cm");
		locHist_KinFitCut10_M2pi->GetXaxis()->SetRangeUser(0.,1.);
		locHist_KinFitCut10_M2pi->GetXaxis()->SetTitleSize(0.05);
		locHist_KinFitCut10_M2pi->GetXaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_M2pi->GetYaxis()->SetTitleSize(0.05);
                locHist_KinFitCut10_M2pi->GetYaxis()->SetLabelSize(0.05);
		locHist_KinFitCut10_M2pi->Draw();
	}

}


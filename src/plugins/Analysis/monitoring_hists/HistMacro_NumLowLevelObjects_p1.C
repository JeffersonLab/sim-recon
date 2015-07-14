// hnamepath: /Independent/Hist_NumReconstructedObjects/NumRFSignals
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTAGMHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTAGHHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumSCHits

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_NumReconstructedObjects");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH1I* locHist_NumRFSignals = (TH1I*)gDirectory->Get("NumRFSignals");
	TH1I* locHist_NumTAGMHits = (TH1I*)gDirectory->Get("NumTAGMHits");
	TH1I* locHist_NumTAGHHits = (TH1I*)gDirectory->Get("NumTAGHHits");
	TH1I* locHist_NumSCHits = (TH1I*)gDirectory->Get("NumSCHits");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("NumLowLevelObjects_p1", "NumLowLevelObjects_p1", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumRFSignals != NULL)
	{
		locHist_NumRFSignals->GetXaxis()->SetTitleSize(0.05);
		locHist_NumRFSignals->GetXaxis()->SetLabelSize(0.05);
		locHist_NumRFSignals->GetYaxis()->SetLabelSize(0.05);
		locHist_NumRFSignals->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumTAGMHits != NULL)
	{
		locHist_NumTAGMHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumTAGMHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumTAGMHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumTAGMHits->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumTAGHHits != NULL)
	{
		locHist_NumTAGHHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumTAGHHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumTAGHHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumTAGHHits->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumSCHits != NULL)
	{
		locHist_NumSCHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumSCHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumSCHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumSCHits->Draw();
	}
	gPad->SetLogy();
	gPad->Update();

}


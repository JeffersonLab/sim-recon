// hnamepath: /Independent/Hist_NumReconstructedObjects/NumHighLevelObjects

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_NumReconstructedObjects");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH2D* locHist_NumHighLevel = (TH2D*)gDirectory->Get("NumHighLevelObjects");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("NumHighLevelObjects", "NumHighLevelObjects", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumHighLevel != NULL)
		locHist_NumHighLevel->Draw("COLZ");
	gPad->SetLogz();
	gPad->Update();
}


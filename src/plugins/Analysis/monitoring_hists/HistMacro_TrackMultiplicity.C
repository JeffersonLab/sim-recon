// hnamepath: /Independent/Hist_TrackMultiplicity/NumGoodReconstructedParticles

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_TrackMultiplicity");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH2D* locHist = (TH2D*)gDirectory->Get("NumGoodReconstructedParticles");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TrackMultiplicity", "TrackMultiplicity", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist != NULL)
		locHist->Draw("COLZ");
	gPad->SetLogz();
	gPad->Update();
}


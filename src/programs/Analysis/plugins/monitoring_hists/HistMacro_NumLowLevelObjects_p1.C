// hnamepath: /Independent/Hist_NumReconstructedObjects/NumSCHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumCDCHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumFDCWireHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumFDCCathodeHits

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_NumReconstructedObjects");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH1I* locHist_NumSCHits = (TH1I*)gDirectory->Get("NumSCHits");
	TH1I* locHist_NumCDCHits = (TH1I*)gDirectory->Get("NumCDCHits");
	TH1I* locHist_NumFDCWireHits = (TH1I*)gDirectory->Get("NumFDCWireHits");
	TH1I* locHist_NumFDCCathodeHits = (TH1I*)gDirectory->Get("NumFDCCathodeHits");

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
	if(locHist_NumSCHits != NULL)
	{
		locHist_NumSCHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumSCHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumSCHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumSCHits->Draw();
	}
	gPad->SetLogy();
	gPad->Update();

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumCDCHits != NULL)
	{
		locHist_NumCDCHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumCDCHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumCDCHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumCDCHits->Draw();
	}
	gPad->SetLogy();
	gPad->Update();

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumFDCWireHits != NULL)
	{
		locHist_NumFDCWireHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumFDCWireHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumFDCWireHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumFDCWireHits->Draw();
	}
	gPad->SetLogy();
	gPad->Update();

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumFDCCathodeHits != NULL)
	{
		locHist_NumFDCCathodeHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumFDCCathodeHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumFDCCathodeHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumFDCCathodeHits->Draw();
	}
	gPad->SetLogy();
	gPad->Update();
}


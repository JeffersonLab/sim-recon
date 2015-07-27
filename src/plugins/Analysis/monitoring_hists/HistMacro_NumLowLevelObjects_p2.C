// hnamepath: /Independent/Hist_NumReconstructedObjects/NumCDCHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumFDCWireHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumFDCCathodeHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTOFHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumBCALHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumFCALHits

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_NumReconstructedObjects");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH1I* locHist_NumCDCHits = (TH1I*)gDirectory->Get("NumCDCHits");
	TH1I* locHist_NumFDCWireHits = (TH1I*)gDirectory->Get("NumFDCWireHits");
	TH1I* locHist_NumFDCCathodeHits = (TH1I*)gDirectory->Get("NumFDCCathodeHits");
	TH1I* locHist_NumTOFHits = (TH1I*)gDirectory->Get("NumTOFHits");
	TH1I* locHist_NumBCALHits = (TH1I*)gDirectory->Get("NumBCALHits");
	TH1I* locHist_NumFCALHits = (TH1I*)gDirectory->Get("NumFCALHits");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("NumLowLevelObjects_p2", "NumLowLevelObjects_p2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
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

	locCanvas->cd(2);
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

	locCanvas->cd(3);
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

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumTOFHits != NULL)
	{
	        locHist_NumTOFHits->GetXaxis()->SetRangeUser(0.0, 200.0);
		locHist_NumTOFHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumTOFHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumTOFHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumTOFHits->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumBCALHits != NULL)
	{
		locHist_NumBCALHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumBCALHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumBCALHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumBCALHits->Draw();
	}
	gPad->SetLogy();
	gPad->Update();

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumFCALHits != NULL)
	{
	        locHist_NumFCALHits->GetXaxis()->SetRangeUser(0.0, 200.0);
		locHist_NumFCALHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumFCALHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumFCALHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumFCALHits->Draw();
	}
}


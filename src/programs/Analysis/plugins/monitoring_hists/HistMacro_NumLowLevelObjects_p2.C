// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTAGMHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTAGHHits
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
	TH1I* locHist_NumTAGMHits = (TH1I*)gDirectory->Get("NumTAGMHits");
	TH1I* locHist_NumTAGHHits = (TH1I*)gDirectory->Get("NumTAGHHits");
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
	if(locHist_NumTAGMHits != NULL)
	{
		locHist_NumTAGMHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumTAGMHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumTAGMHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumTAGMHits->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumTAGHHits != NULL)
	{
		locHist_NumTAGHHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumTAGHHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumTAGHHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumTAGHHits->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumTOFHits != NULL)
	{
		locHist_NumTOFHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumTOFHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumTOFHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumTOFHits->Draw();
	}

	locCanvas->cd(4);
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

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumFCALHits != NULL)
	{
		locHist_NumFCALHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumFCALHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumFCALHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumFCALHits->Draw();
	}
}


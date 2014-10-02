// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTAGMHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTAGHHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumSCHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTOFHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumBCALHits
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumFCALHits
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
	TH1I* locHist_NumTAGMHits = (TH1I*)gDirectory->Get("NumTAGMHits");
	TH1I* locHist_NumTAGHHits = (TH1I*)gDirectory->Get("NumTAGHHits");
	TH1I* locHist_NumSCHits = (TH1I*)gDirectory->Get("NumSCHits");
	TH1I* locHist_NumTOFHits = (TH1I*)gDirectory->Get("NumTOFHits");
	TH1I* locHist_NumBCALHits = (TH1I*)gDirectory->Get("NumBCALHits");
	TH1I* locHist_NumFCALHits = (TH1I*)gDirectory->Get("NumFCALHits");
	TH1I* locHist_NumCDCHits = (TH1I*)gDirectory->Get("NumCDCHits");
	TH1I* locHist_NumFDCWireHits = (TH1I*)gDirectory->Get("NumFDCWireHits");
	TH1I* locHist_NumFDCCathodeHits = (TH1I*)gDirectory->Get("NumFDCCathodeHits");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("NumLowLevelObjects"); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 3);

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
	if(locHist_NumSCHits != NULL)
	{
		locHist_NumSCHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumSCHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumSCHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumSCHits->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumTOFHits != NULL)
	{
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

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumFCALHits != NULL)
	{
		locHist_NumFCALHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumFCALHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumFCALHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumFCALHits->Draw();
	}

	locCanvas->cd(7);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumCDCHits != NULL)
	{
		locHist_NumCDCHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumCDCHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumCDCHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumCDCHits->Draw();
	}

	locCanvas->cd(8);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumFDCWireHits != NULL)
	{
		locHist_NumFDCWireHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumFDCWireHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumFDCWireHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumFDCWireHits->Draw();
	}

	locCanvas->cd(9);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumFDCCathodeHits != NULL)
	{
		locHist_NumFDCCathodeHits->GetXaxis()->SetTitleSize(0.05);
		locHist_NumFDCCathodeHits->GetXaxis()->SetLabelSize(0.05);
		locHist_NumFDCCathodeHits->GetYaxis()->SetLabelSize(0.05);
		locHist_NumFDCCathodeHits->Draw();
	}
}


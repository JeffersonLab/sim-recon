// hnamepath: /RF/DeltaT_RF_OtherRFs/RFDeltaT_FDC_TOF
// hnamepath: /RF/DeltaT_RF_OtherRFs/RFDeltaT_FDC_TAGH
// hnamepath: /RF/DeltaT_RF_OtherRFs/RFDeltaT_FDC_PSC
// hnamepath: /RF/DeltaT_RF_OtherRFs/RFDeltaT_PSC_TAGH
// hnamepath: /RF/DeltaT_RF_OtherRFs/RFDeltaT_PSC_TOF
// hnamepath: /RF/DeltaT_RF_OtherRFs/RFDeltaT_TAGH_TOF

{
	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("DeltaT_RF_OtherRFs");
	TH1I* locHist_RFDeltaT_FDC_TOF = (TH1I*)gDirectory->Get("RFDeltaT_FDC_TOF");
	TH1I* locHist_RFDeltaT_FDC_TAGH = (TH1I*)gDirectory->Get("RFDeltaT_FDC_TAGH");
	TH1I* locHist_RFDeltaT_FDC_PSC = (TH1I*)gDirectory->Get("RFDeltaT_FDC_PSC");
	TH1I* locHist_RFDeltaT_PSC_TAGH = (TH1I*)gDirectory->Get("RFDeltaT_PSC_TAGH");
	TH1I* locHist_RFDeltaT_PSC_TOF = (TH1I*)gDirectory->Get("RFDeltaT_PSC_TOF");
	TH1I* locHist_RFDeltaT_TAGH_TOF = (TH1I*)gDirectory->Get("RFDeltaT_TAGH_TOF");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("RF_p2", "RF_p2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_TAGH_TOF != NULL)
	{
		locHist_RFDeltaT_TAGH_TOF->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_TAGH_TOF->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_TAGH_TOF->GetXaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_TAGH_TOF->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_TAGH_TOF->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_PSC_TOF != NULL)
	{
		locHist_RFDeltaT_PSC_TOF->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_PSC_TOF->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_PSC_TOF->GetXaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_PSC_TOF->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_PSC_TOF->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_FDC_TOF != NULL)
	{
		locHist_RFDeltaT_FDC_TOF->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_TOF->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_TOF->GetXaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_TOF->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_TOF->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_FDC_TAGH != NULL)
	{
		locHist_RFDeltaT_FDC_TAGH->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_TAGH->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_TAGH->GetXaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_TAGH->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_TAGH->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_FDC_PSC != NULL)
	{
		locHist_RFDeltaT_FDC_PSC->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_PSC->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_PSC->GetXaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_PSC->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_PSC->Draw();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_PSC_TAGH != NULL)
	{
		locHist_RFDeltaT_PSC_TAGH->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_PSC_TAGH->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_PSC_TAGH->GetXaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_PSC_TAGH->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_PSC_TAGH->Draw();
	}
}


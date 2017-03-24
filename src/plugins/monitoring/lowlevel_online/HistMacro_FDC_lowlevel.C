// hnamepath: /lowlevel_online/FDC/fdc_adc_multi
// hnamepath: /lowlevel_online/FDC/fdc_tdc_multi
// hnamepath: /lowlevel_online/FDC/fdc_tdc_time
// hnamepath: /lowlevel_online/FDC/fdc_adc_integral

{
	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locDirectory = (TDirectory*)locInitDirectory->FindObjectAny("lowlevel_online");
	if(locDirectory == NULL)
		return;

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("FDC", "FDC", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);


	TH1I* locHist_ADCmulti = (TH1I*)locDirectory->Get("FDC/fdc_adc_multi");
	TH1I* locHist_TDCmulti = (TH1I*)locDirectory->Get("FDC/fdc_tdc_multi");
	TH1I* locHist_ADCtime  = (TH1I*)locDirectory->Get("FDC/fdc_adc_time");
	TH1I* locHist_TDCtime  = (TH1I*)locDirectory->Get("FDC/fdc_tdc_time");

	//# ADC hits
	locCanvas->cd(1);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_ADCmulti != NULL)
	{
		locHist_ADCmulti->SetTitle("FDC ADC Multiplicity");
		//locHist_ADCmulti->Rebin(2);
		//locHist_ADCmulti->GetYaxis()->SetRangeUser(0.0, 1.05*locHist_ADCmulti->GetBinContent(locHist_ADCmulti->GetMaximumBin()));
		locHist_ADCmulti->GetXaxis()->SetTitleSize(0.05);
		locHist_ADCmulti->GetXaxis()->SetTitle("# ADC Hits");
		//locHist_ADCmulti->GetYaxis()->SetTitle("");
		locHist_ADCmulti->GetXaxis()->SetLabelSize(0.05);
		locHist_ADCmulti->GetYaxis()->SetLabelSize(0.05);
		locHist_ADCmulti->SetFillColor(kYellow);
		locHist_ADCmulti->Draw("");
	}

	//# TDC hits
	locCanvas->cd(2);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_TDCmulti != NULL)
	{
		locHist_TDCmulti->SetTitle("FDC TDC Multiplicity");
		//locHist_TDCmulti->Rebin(2);
		//locHist_TDCmulti->GetYaxis()->SetRangeUser(0.0, 1.05*locHist_TDCmulti->GetBinContent(locHist_TDCmulti->GetMaximumBin()));
		locHist_TDCmulti->GetXaxis()->SetTitleSize(0.05);
		locHist_TDCmulti->GetXaxis()->SetTitle("# TDC Hits");
		//locHist_TDCmulti->GetYaxis()->SetTitle("");
		locHist_TDCmulti->GetXaxis()->SetLabelSize(0.05);
		locHist_TDCmulti->GetYaxis()->SetLabelSize(0.05);
		locHist_TDCmulti->SetFillColor(kYellow);
		locHist_TDCmulti->Draw("");
	}

	//ADC time
	locCanvas->cd(3);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_ADCtime != NULL)
	{
		locHist_ADCtime->SetTitle("FDC ADC Time");
		//locHist_ADCtime->Rebin(2);
		locHist_ADCtime->GetXaxis()->SetRangeUser(100, 600);
		locHist_ADCtime->GetXaxis()->SetTitleSize(0.05);
		//locHist_ADCtime->GetYaxis()->SetTitle("");
		locHist_ADCtime->GetXaxis()->SetLabelSize(0.05);
		locHist_ADCtime->GetYaxis()->SetLabelSize(0.05);
		locHist_ADCtime->SetFillColor(kYellow);
		locHist_ADCtime->Draw("");
	}
	
	//TDC time
	locCanvas->cd(4);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_TDCtime != NULL)
	{
		locHist_TDCtime->SetTitle("FDC TDC Time");
		//locHist_TDCtime->Rebin(2);
		locHist_TDCtime->GetXaxis()->SetRangeUser(-500, 8000);
		locHist_TDCtime->GetXaxis()->SetTitleSize(0.05);
		//locHist_TDCtime->GetYaxis()->SetTitle("");
		locHist_TDCtime->GetXaxis()->SetLabelSize(0.05);
		locHist_TDCtime->GetYaxis()->SetLabelSize(0.05);
		locHist_TDCtime->SetFillColor(kYellow);
		locHist_TDCtime->Draw("");
	}
}

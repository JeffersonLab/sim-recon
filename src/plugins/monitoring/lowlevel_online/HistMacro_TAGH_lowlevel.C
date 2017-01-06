// hnamepath: /lowlevel_online/TAGH/tagh_adc_multi
// hnamepath: /lowlevel_online/TAGH/tagh_tdc_multi
// hnamepath: /lowlevel_online/TAGH/tagh_tdc_time
// hnamepath: /lowlevel_online/TAGH/tagh_adc_integral

{
	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locDirectory = (TDirectory*)locInitDirectory->FindObjectAny("lowlevel_online");
	if(locDirectory == NULL)
		return;

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TAGH", "TAGH", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);


	TH1I* locHist_ADCmulti = (TH1I*)locDirectory->Get("TAGH/tagh_adc_multi");
	TH1I* locHist_TDCmulti = (TH1I*)locDirectory->Get("TAGH/tagh_tdc_multi");
	TH1I* locHist_TDCtime  = (TH1I*)locDirectory->Get("TAGH/tagh_tdc_time");
	TH1I* locHist_ADCintegral = (TH1I*)locDirectory->Get("TAGH/tagh_adc_integral");

	//# ADC hits
	locCanvas->cd(1);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_ADCmulti != NULL)
	{
		locHist_ADCmulti->SetTitle("TAGH ADC Multiplicity");
		locHist_ADCmulti->Rebin(2);
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
		locHist_TDCmulti->SetTitle("TAGH TDC Multiplicity");
		locHist_TDCmulti->Rebin(2);
		//locHist_TDCmulti->GetYaxis()->SetRangeUser(0.0, 1.05*locHist_TDCmulti->GetBinContent(locHist_TDCmulti->GetMaximumBin()));
		locHist_TDCmulti->GetXaxis()->SetTitleSize(0.05);
		locHist_TDCmulti->GetXaxis()->SetTitle("# TDC Hits");
		//locHist_TDCmulti->GetYaxis()->SetTitle("");
		locHist_TDCmulti->GetXaxis()->SetLabelSize(0.05);
		locHist_TDCmulti->GetYaxis()->SetLabelSize(0.05);
		locHist_TDCmulti->SetFillColor(kYellow);
		locHist_TDCmulti->Draw("");
	}

	//ADC integral
	locCanvas->cd(3);
	gPad->SetTicks();
	//gPad->SetLogy();
	//gPad->SetGrid();
	if(locHist_ADCintegral != NULL)
	{
		locHist_ADCintegral->SetTitle("TAGH ADC Integral");
		//locHist_ADCintegral->Rebin(2);
		locHist_ADCintegral->GetXaxis()->SetRangeUser(1000.0, 4000.);
		locHist_ADCintegral->GetXaxis()->SetTitleSize(0.05);
		//locHist_ADCintegral->GetYaxis()->SetTitle("");
		locHist_ADCintegral->GetXaxis()->SetLabelSize(0.05);
		locHist_ADCintegral->GetYaxis()->SetLabelSize(0.05);
		locHist_ADCintegral->SetFillColor(kYellow);
		locHist_ADCintegral->Draw("");
	}

	//TDC time
	locCanvas->cd(4);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_TDCtime != NULL)
	{
		locHist_TDCtime->SetTitle("TAGH TDC Time");
		locHist_TDCtime->Rebin(2);
		locHist_TDCtime->GetXaxis()->SetRangeUser(0, 1500);
		locHist_TDCtime->GetXaxis()->SetTitleSize(0.05);
		//locHist_TDCtime->GetYaxis()->SetTitle("");
		locHist_TDCtime->GetXaxis()->SetLabelSize(0.05);
		locHist_TDCtime->GetYaxis()->SetLabelSize(0.05);
		locHist_TDCtime->SetFillColor(kYellow);
		locHist_TDCtime->Draw("");
	}

}
// hnamepath: /lowlevel_online/CDC/cdc_adc_multi
// hnamepath: /lowlevel_online/CDC/cdc_adc_pedestal
// hnamepath: /lowlevel_online/CDC/cdc_adc_time
// hnamepath: /lowlevel_online/CDC/cdc_adc_integral

{
	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locDirectory = (TDirectory*)locInitDirectory->FindObjectAny("lowlevel_online");
	if(locDirectory == NULL)
		return;

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("CDC", "CDC", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);


	TH1I* locHist_ADCmulti = (TH1I*)locDirectory->Get("CDC/cdc_adc_multi");
	TH1I* locHist_ADCpedestal = (TH1I*)locDirectory->Get("CDC/cdc_adc_pedestal");
	TH1I* locHist_ADCtime  = (TH1I*)locDirectory->Get("CDC/cdc_adc_time");
	TH1I* locHist_ADCintegral = (TH1I*)locDirectory->Get("CDC/cdc_adc_integral");

	//# ADC hits
	locCanvas->cd(1);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_ADCmulti != NULL)
	{
		locHist_ADCmulti->SetTitle("CDC ADC Multiplicity");
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

	//# ADC hits
	locCanvas->cd(2);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_ADCpedestal != NULL)
	{
		locHist_ADCpedestal->SetTitle("CDC ADC Pedestal");
		//locHist_ADCpedestal->Rebin(2);
		//locHist_ADCpedestal->GetYaxis()->SetRangeUser(0.0, 1.05*locHist_ADCpedestal->GetBinContent(locHist_ADCpedestal->GetMaximumBin()));
		locHist_ADCpedestal->GetXaxis()->SetRangeUser(0.0, 500.);
		locHist_ADCpedestal->GetXaxis()->SetTitleSize(0.05);
		locHist_ADCpedestal->GetXaxis()->SetTitle("# ADC Hits");
		//locHist_ADCpedestal->GetYaxis()->SetTitle("");
		locHist_ADCpedestal->GetXaxis()->SetLabelSize(0.05);
		locHist_ADCpedestal->GetYaxis()->SetLabelSize(0.05);
		locHist_ADCpedestal->SetFillColor(kYellow);
		locHist_ADCpedestal->Draw("");
	}

	//ADC integral
	locCanvas->cd(3);
	gPad->SetTicks();
	//gPad->SetLogy();
	//gPad->SetGrid();
	if(locHist_ADCintegral != NULL)
	{
		locHist_ADCintegral->SetTitle("CDC ADC Integral");
		//locHist_ADCintegral->Rebin(2);
		locHist_ADCintegral->GetXaxis()->SetRangeUser(0.0, 5000.);
		locHist_ADCintegral->GetXaxis()->SetTitleSize(0.05);
		//locHist_ADCintegral->GetYaxis()->SetTitle("");
		locHist_ADCintegral->GetXaxis()->SetLabelSize(0.05);
		locHist_ADCintegral->GetYaxis()->SetLabelSize(0.05);
		locHist_ADCintegral->SetFillColor(kYellow);
		locHist_ADCintegral->Draw("");
	}

	//ADC time
	locCanvas->cd(4);
	gPad->SetTicks();
	//gPad->SetGrid();
	if(locHist_ADCtime != NULL)
	{
		locHist_ADCtime->SetTitle("CDC ADC Time");
		locHist_ADCtime->Rebin(5);
		locHist_ADCtime->GetXaxis()->SetRangeUser(0, 1800);
		locHist_ADCtime->GetXaxis()->SetTitleSize(0.05);
		//locHist_ADCtime->GetYaxis()->SetTitle("");
		locHist_ADCtime->GetXaxis()->SetLabelSize(0.05);
		locHist_ADCtime->GetYaxis()->SetLabelSize(0.05);
		locHist_ADCtime->SetFillColor(kYellow);
		locHist_ADCtime->Draw("");
	}

}
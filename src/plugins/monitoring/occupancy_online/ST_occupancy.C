
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/st_num_events
// hnamepath: /occupancy/st_adc_occ
// hnamepath: /occupancy/st_tdc_occ

{
	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	double Nevents = 1.0;
	TH1I *st_num_events = (TH1I*)gDirectory->FindObjectAny("st_num_events");
	if(st_num_events) Nevents = (double)st_num_events->GetBinContent(1);

	TH1I *st_adc_occ = (TH1I*)gDirectory->FindObjectAny("st_adc_occ");
	TH1I *st_tdc_occ = (TH1I*)gDirectory->FindObjectAny("st_tdc_occ");

	TLegend *legend_sa = new TLegend(0.1,0.85,0.3,0.9);
	TLegend *legend_na = new TLegend(0.3,0.85,0.5,0.9);

	if(st_tdc_occ){
		st_tdc_occ->SetBarWidth(0.5);
		st_tdc_occ->SetBarOffset(0);
		st_tdc_occ->SetFillColor(kGreen);
		st_tdc_occ->SetStats(0);
		st_tdc_occ->SetXTitle("Channel number");
		st_tdc_occ->SetYTitle("fADC/TDC hit count");
		st_tdc_occ->SetTitleSize(0.05,"X");
		st_tdc_occ->GetXaxis()->CenterTitle();
		st_tdc_occ->SetTitleSize(0.05,"Y");
		st_tdc_occ->GetYaxis()->CenterTitle();
		st_tdc_occ->GetYaxis()->SetRangeUser(0.0, st_tdc_occ->GetMaximum());
		if(st_adc_occ)st_adc_occ->GetYaxis()->SetRangeUser(0.0, 1.05*st_tdc_occ->GetMaximum());
	}
	
	if(st_adc_occ){
		st_adc_occ->SetBarWidth(0.5);
		st_adc_occ->SetBarOffset(0.5);
		st_adc_occ->SetFillColor(kRed);
		st_adc_occ->SetStats(0);
	}

	legend_sa->AddEntry(st_adc_occ,"fADC","f");
	legend_na->AddEntry(st_tdc_occ,"TDC","f");

	// Just for testing
	if(gPad == NULL){
		TCanvas *c1 = new TCanvas("c1");
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}
	if(!gPad) {savedir->cd(); return;}

	TCanvas *c1 = gPad->GetCanvas();
	c1->cd(0);
	c1->Clear();

	gPad->SetTicks();
	gPad->SetGridy();
	if(st_adc_occ) st_adc_occ->DrawCopy("BAR");
	if(st_tdc_occ) st_tdc_occ->DrawCopy("BAR same");

	legend_sa->Draw();
	legend_na->Draw();
	
	if(st_tdc_occ){
		char str[256];
		sprintf(str,"%0.0f events", Nevents);
		TLatex lat(15.0, 1.07*st_tdc_occ->GetMaximum(), str);
		lat.SetTextAlign(22);
		lat.SetTextSize(0.035);
		lat.Draw();
	}
}

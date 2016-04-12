
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/st_adc_occ
// hnamepath: /occupancy/st_tdc_occ

{
	TDirectory *savedir = gDirectory;

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	TH1I *st_adc_occ = (TH1I*)gDirectory->FindObjectAny("st_adc_occ");
	TH1I *st_tdc_occ = (TH1I*)gDirectory->FindObjectAny("st_tdc_occ");

	TLegend *legend_sa = new TLegend(0.1,0.85,0.3,0.9);
	TLegend *legend_na = new TLegend(0.3,0.85,0.5,0.9);

//	TLegend *legend_s = new TLegend(0.1,0.85,0.3,0.9);
//	TLegend *legend_n = new TLegend(0.3,0.85,0.5,0.9);

	st_tdc_occ->SetBarWidth(0.5);
	st_tdc_occ->SetBarOffset(0);
	st_tdc_occ->SetFillColor(2);
	st_tdc_occ->SetStats(0);
	st_tdc_occ->SetXTitle("Channel number");
	st_tdc_occ->SetYTitle("fADC/TDC occupancy");
	st_tdc_occ->SetTitleSize(0.05,"X");
	st_tdc_occ->GetXaxis()->CenterTitle();
	st_tdc_occ->SetTitleSize(0.05,"Y");
	st_tdc_occ->GetYaxis()->CenterTitle();
	st_tdc_occ->GetYaxis()->SetRangeUser(0.0, st_tdc_occ->GetMaximum());

	st_adc_occ->SetBarWidth(0.5);
	st_adc_occ->SetBarOffset(0.5);
	st_adc_occ->SetFillColor(3);
	st_adc_occ->SetStats(0);

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
	if(st_adc_occ) st_tdc_occ->Draw("BAR");
	if(st_tdc_occ) st_adc_occ->Draw("BAR sames");

	legend_sa->Draw();
	legend_na->Draw();

	savedir->cd();
}

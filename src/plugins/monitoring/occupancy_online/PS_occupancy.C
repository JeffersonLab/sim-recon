
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/ps_num_events
// hnamepath: /occupancy/psc_adc_left_occ
// hnamepath: /occupancy/psc_adc_right_occ
// hnamepath: /occupancy/psc_tdc_left_occ
// hnamepath: /occupancy/psc_tdc_right_occ
// hnamepath: /occupancy/ps_left_occ
// hnamepath: /occupancy/ps_right_occ
//
// e-mail: davidl@jlab.org
// e-mail: somov@jlab.org
//

{
	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	double Nevents = 1.0;
	TH1I *ps_num_events = (TH1I*)gDirectory->FindObjectAny("ps_num_events");
	if(ps_num_events) Nevents = (double)ps_num_events->GetBinContent(1);

	TH1I *psc_adc_left_occ  = (TH1I*)gDirectory->FindObjectAny("psc_adc_left_occ");
	TH1I *psc_adc_right_occ = (TH1I*)gDirectory->FindObjectAny("psc_adc_right_occ");
	TH1I *psc_tdc_left_occ  = (TH1I*)gDirectory->FindObjectAny("psc_tdc_left_occ");
	TH1I *psc_tdc_right_occ = (TH1I*)gDirectory->FindObjectAny("psc_tdc_right_occ");
	TH1I *ps_left_occ       = (TH1I*)gDirectory->FindObjectAny("ps_left_occ");
	TH1I *ps_right_occ      = (TH1I*)gDirectory->FindObjectAny("ps_right_occ");

	TLegend *legend_sa = new TLegend(0.1,0.85,0.3,0.9);
	TLegend *legend_na = new TLegend(0.3,0.85,0.5,0.9);

	if(psc_tdc_left_occ){
		psc_tdc_left_occ->SetBarWidth(0.5);
		psc_tdc_left_occ->SetBarOffset(0);
		psc_tdc_left_occ->SetFillColor(kGreen);
		psc_tdc_left_occ->SetStats(0);
		psc_tdc_left_occ->SetXTitle("Module number");
		psc_tdc_left_occ->SetYTitle("fADC/TDC occupancy");
		psc_tdc_left_occ->SetTitleSize(0.05,"X");
		psc_tdc_left_occ->GetXaxis()->CenterTitle();
		psc_tdc_left_occ->SetTitleSize(0.05,"Y");
		psc_tdc_left_occ->GetYaxis()->CenterTitle();
	}
	
	if(psc_adc_left_occ){
		psc_adc_left_occ->SetBarWidth(0.5);
		psc_adc_left_occ->SetBarOffset(0.5);
		psc_adc_left_occ->SetFillColor(kRed);
		psc_adc_left_occ->SetStats(0);
	}
	
	if(psc_tdc_right_occ){
		psc_tdc_right_occ->SetBarWidth(0.5);
		psc_tdc_right_occ->SetBarOffset(0);
		psc_tdc_right_occ->SetFillColor(kGreen);
		psc_tdc_right_occ->SetStats(0);
		psc_tdc_right_occ->SetXTitle("Module number");
		psc_tdc_right_occ->SetYTitle("fADC/TDC occupancy");
		psc_tdc_right_occ->SetTitleSize(0.05,"X");
		psc_tdc_right_occ->GetXaxis()->CenterTitle();
		psc_tdc_right_occ->SetTitleSize(0.05,"Y");
		psc_tdc_right_occ->GetYaxis()->CenterTitle();
	}

	if(psc_adc_right_occ){
		psc_adc_right_occ->SetBarWidth(0.5);
		psc_adc_right_occ->SetBarOffset(0.5);
		psc_adc_right_occ->SetFillColor(kRed);
		psc_adc_right_occ->SetStats(0);
	}

	if(ps_left_occ){
		ps_left_occ->SetBarWidth(0.5);
		ps_left_occ->SetBarOffset(0);
		ps_left_occ->SetFillColor(kRed);
		ps_left_occ->SetStats(0);
		ps_left_occ->SetXTitle("Column number");
		ps_left_occ->SetYTitle("fADC occupancy");
		ps_left_occ->SetTitleSize(0.05,"X");
		ps_left_occ->GetXaxis()->CenterTitle();
		ps_left_occ->SetTitleSize(0.05,"Y");
		ps_left_occ->GetYaxis()->CenterTitle();
	}
	
	if(ps_right_occ){
		ps_right_occ->SetBarWidth(0.5);
		ps_right_occ->SetBarOffset(0);
		ps_right_occ->SetFillColor(kRed);
		ps_right_occ->SetStats(0);
		ps_right_occ->SetXTitle("Column number");
		ps_right_occ->SetYTitle("fADC occupancy");
		ps_right_occ->SetTitleSize(0.05,"X");
		ps_right_occ->GetXaxis()->CenterTitle();
		ps_right_occ->SetTitleSize(0.05,"Y");
		ps_right_occ->GetYaxis()->CenterTitle();
	}

	legend_sa->AddEntry(psc_tdc_left_occ,"TDC","f");
	legend_na->AddEntry(psc_adc_left_occ,"fADC","f");


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

	c1->Divide(2,2);
	TVirtualPad *pad1 = c1->cd(1);
	pad1->SetTicks();
	pad1->SetGridy();
	if(psc_tdc_left_occ) psc_tdc_left_occ->DrawCopy("BAR");
	if(psc_adc_left_occ) psc_adc_left_occ->DrawCopy("BAR same");

	legend_sa->Draw();
	legend_na->Draw();

	TVirtualPad *pad2 = c1->cd(2);
	pad2->SetTicks();
	pad2->SetGridy();
	if(psc_tdc_right_occ) psc_tdc_right_occ->DrawCopy("BAR");
	if(psc_adc_right_occ) psc_adc_right_occ->DrawCopy("BAR same");

	legend_sa->Draw();
	legend_na->Draw();

	TVirtualPad *pad3 = c1->cd(3);
	pad3->SetTicks();
	pad3->SetGridy();
	if(ps_left_occ) ps_left_occ->DrawCopy();

	legend_na->Draw("BAR");


	TVirtualPad *pad4 = c1->cd(4);
	pad4->SetTicks();
	pad4->SetGridy();
	if(ps_right_occ) ps_right_occ->DrawCopy();

	legend_na->Draw("BAR");

}

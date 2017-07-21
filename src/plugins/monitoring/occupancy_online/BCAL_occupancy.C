
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/bcal_num_events
// hnamepath: /occupancy/bcal_adc_occ
// hnamepath: /occupancy/bcal_tdc_occ
//
// e-mail: davidl@jlab.org
// e-mail: elton@jlab.org
// e-mail: dalton@jlab.org
// e-mail: zisis@uregina.ca
//

{
	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.
        TDirectory *savedir = gDirectory;
	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	TH2I *bcal_adc_occ = (TH2I*)gDirectory->FindObjectAny("bcal_adc_occ");
	TH2I *bcal_tdc_occ = (TH2I*)gDirectory->FindObjectAny("bcal_tdc_occ");

	double Nevents = 1.0;
	TH1I *bcal_num_events = (TH1I*)gDirectory->FindObjectAny("bcal_num_events");
	if(bcal_num_events) Nevents = (double)bcal_num_events->GetBinContent(1);

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
	
	c1->Divide(2,1);
	
	TVirtualPad *pad1 = c1->cd(1);
	pad1->SetTicks();
	pad1->SetLeftMargin(0.15);
	pad1->SetRightMargin(0.15);
	if(bcal_adc_occ){
		bcal_adc_occ->SetStats(0);
		bcal_adc_occ->Draw("colz");
	}

	TVirtualPad *pad2 = c1->cd(2);
	pad2->SetTicks();
	pad2->SetLeftMargin(0.15);
	pad2->SetRightMargin(0.15);
	if(bcal_tdc_occ){
		bcal_tdc_occ->SetStats(0);
		bcal_tdc_occ->Draw("colz");
	}

	char str[256];
	sprintf(str,"%0.0f events", Nevents);
	TLatex lat(50.0, 26.5, str);
	lat.SetTextAlign(22);
	lat.SetTextSize(0.035);
	lat.Draw();

}



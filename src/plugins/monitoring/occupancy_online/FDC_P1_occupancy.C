
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/fdc_num_events
// hnamepath: /occupancy/fdc_cathode_occ_1
// hnamepath: /occupancy/fdc_wire__occ_1

{
	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	TH2F *fdc_cathode_occ = (TH2F*)gDirectory->FindObjectAny("fdc_cathode_occ_1");
	TH2F *fdc_wire_occ = (TH2F*)gDirectory->FindObjectAny("fdc_wire_occ_1");

	double Nevents = 1.0;
	TH1I *fdc_num_events = (TH1I*)gDirectory->FindObjectAny("fdc_num_events");
	if(fdc_num_events) Nevents = (double)fdc_num_events->GetBinContent(1);

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
	
	c1->Divide(1,2);
	
	TVirtualPad *pad1 = c1->cd(1);
	pad1->SetTicks();
    pad1->SetLogz();
	pad1->SetLeftMargin(0.15);
	pad1->SetRightMargin(0.15);
	if(fdc_cathode_occ){
		fdc_cathode_occ->SetStats(0);
        fdc_cathode_occ->Scale(1./Nevents);
		fdc_cathode_occ->Draw("colz");
	}

	TVirtualPad *pad2 = c1->cd(2);
	pad2->SetTicks();
    pad2->SetLogz();
	pad2->SetLeftMargin(0.15);
	pad2->SetRightMargin(0.15);
	if(fdc_wire_occ){
		fdc_wire_occ->SetStats(0);
        fdc_wire_occ->Scale(1./Nevents);
		fdc_wire_occ->Draw("colz");
	}

}



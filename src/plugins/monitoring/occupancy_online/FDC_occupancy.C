
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/fdc_num_events
// hnamepath: /occupancy/fdc_cathode_occ
// hnamepath: /occupancy/fdc_wire_occ

{
	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	TH2I *fdc_cathode_occ = (TH2I*)gDirectory->FindObjectAny("fdc_cathode_occ");
	TH2I *fdc_wire_occ = (TH2I*)gDirectory->FindObjectAny("fdc_wire_occ");

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
	c1->Clear();
	
	TPad *pad1 = new TPad("pad1", "", 0.0, 0.0, 0.66, 1.0);
	pad1->SetTicks();
	pad1->SetLogz();
	pad1->SetLeftMargin(0.10);
	pad1->SetRightMargin(0.15);
	pad1->Draw();
	pad1->cd();
	if(fdc_cathode_occ){
		fdc_cathode_occ->SetStats(0);
		fdc_cathode_occ->SetYTitle("strip");
		fdc_cathode_occ->SetXTitle("cathode plane");
		TH1* h = fdc_cathode_occ->DrawCopy("colz");
		h->Scale(1./Nevents);
		h->GetZaxis()->SetRangeUser(0.001, 1.0);
	}

	c1->cd(0);
	TPad *pad2 = new TPad("pad2", "", 0.66, 0.0, 1.0, 1.0);
	pad2->SetTicks();
	pad2->SetLogz();
	pad2->SetLeftMargin(0.10);
	pad2->SetRightMargin(0.12);
	pad2->Draw();
	pad2->cd();
	if(fdc_wire_occ){
		fdc_wire_occ->SetStats(0);
		fdc_wire_occ->SetYTitle("wire");
		fdc_wire_occ->SetXTitle("wire plane");
		TH1* h = fdc_wire_occ->DrawCopy("colz");
		h->Scale(1./Nevents);
		h->GetZaxis()->SetRangeUser(0.001, 1.0);
	}

}



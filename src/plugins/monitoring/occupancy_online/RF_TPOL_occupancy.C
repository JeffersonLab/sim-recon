
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/rf_occ
// hnamepath: /occupancy/tpol_occ
// hnamepath: /occupancy/rf_num_events
//
// e-mail: davidl@jlab.org
// e-mail: nsparks@jlab.org
//

{
	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	TH2I *rf_occ = (TH2I*)gDirectory->FindObjectAny("rf_occ");
	TH2I *tpol_occ = (TH2I*)gDirectory->FindObjectAny("tpol_occ");
	TH1I *rf_num_events = (TH1I*)gDirectory->FindObjectAny("rf_num_events");

	double Nevents = 1.0;
	if(rf_num_events) Nevents = (double)rf_num_events->GetBinContent(1);

	TLegend *legend_sa = new TLegend(0.1,0.85,0.3,0.9);
	legend_sa->AddEntry(rf_occ, "TDC","f");

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

	c1->Divide(1,2);
	
	TVirtualPad *pad1 = c1->cd(1);
	pad1->SetTicks();
	pad1->SetGrid();
	if(rf_occ){
		
		// Draw axes
		TH1D *rf_axes = (TH1D *)dir->Get("rf_axes");
		if(!rf_axes) rf_axes = (TH1D*)rf_occ->Clone("rf_axes");

		rf_axes->SetStats(0);
		rf_axes->GetYaxis()->SetRangeUser(0.0, 1.05*rf_occ->GetMaximum());
		rf_axes->Draw();
		rf_occ->SetFillColor(kBlue);
		rf_occ->SetLineColor(kBlack);
		rf_occ->SetLineWidth(5);
		rf_occ->Draw("same");
		
		char str[256];
		sprintf(str,"%0.0f events", Nevents);
		TLatex lat(2.5, 1.075*rf_occ->GetMaximum(), str);
		lat.SetTextAlign(22);
		lat.SetTextSize(0.035);
		lat.Draw();

		legend_sa->Draw();
	}

	TVirtualPad *pad2 = c1->cd(2);
	pad2->SetTicks();
	pad2->SetGrid();
	if(tpol_occ){
		tpol_occ->SetFillColor(kOrange);
		tpol_occ->SetLineWidth(5);
		tpol_occ->SetLineColor(kBlack);
		tpol_occ->Draw();
	}

}

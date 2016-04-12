
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/fcal_occ
// hnamepath: /occupancy/fcal_num_events

{
	TDirectory *savedir = gDirectory;

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	TH2I *fcal_occ = (TH2I*)gDirectory->FindObjectAny("fcal_occ");
	TH1I *fcal_num_events = (TH1I*)gDirectory->FindObjectAny("fcal_num_events");

	double Nevents = 1.0;
	if(fcal_num_events) Nevents = (double)fcal_num_events->GetBinContent(1);

	TLegend *legend_sa = new TLegend(0.1,0.85,0.3,0.9);
	legend_sa->AddEntry(fcal_occ, "fADC","f");

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
	gPad->SetGrid();
	gPad->SetRightMargin(2.0);
	gPad->SetLeftMargin(2.0);
	if(fcal_occ){
		fcal_occ->SetStats(0);
		fcal_occ->Draw("colz");
		
		char str[256];
		sprintf(str,"%0.0f events", Nevents);
		TLatex lat(30.0, 61.75, str);
		lat.SetTextAlign(22);
		lat.SetTextSize(0.035);
		lat.Draw();
	}

	legend_sa->Draw();

	savedir->cd();
}

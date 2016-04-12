
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/bcal_adc_occ
// hnamepath: /occupancy/bcal_tdc_occ

{
	TDirectory *savedir = gDirectory;

	TH2I *bcal_adc_occ = (TH2I*)gDirectory->FindObjectAny("bcal_adc_occ");
	TH2I *bcal_tdc_occ = (TH2I*)gDirectory->FindObjectAny("bcal_tdc_occ");

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

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

	savedir->cd();
}



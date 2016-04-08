


// hnamepath: /bcal/bcal_tdc_t
// hnamepath: /bcal/bcal_fadc_t
// hnamepath: /bcal/bcal_Uhit_t_TDC
// hnamepath: /bcal/bcal_Uhit_tTDC_tADC
// hnamepath: /bcal/bcal_Uhit_tdiff
// hnamepath: /bcal/bcal_Uhit_tdiff_ave

{
	gStyle->SetOptStat(111110);	
	gStyle->SetStatX(0.99);
	gStyle->SetStatY(0.99);
	gStyle->SetStatW(0.25);

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal");
	if(dir) dir->cd();

	TH1I *h1 = (TH1I*)gDirectory->FindObjectAny("bcal_fadc_t");
	TH1I *h2 = (TH1I*)gDirectory->FindObjectAny("bcal_tdc_t");
	TH1I *h3 = (TH1I*)gDirectory->FindObjectAny("bcal_Uhit_t_TDC");
	TH1I *h5 = (TH1I*)gDirectory->FindObjectAny("bcal_Uhit_tTDC_tADC");
	TH1I *h4 = (TH1I*)gDirectory->FindObjectAny("bcal_Uhit_tdiff");
	TH1I *h6 = (TH1I*)gDirectory->FindObjectAny("bcal_Uhit_tdiff_ave");

	// Just for testing
	if(gPad == NULL){
		TCanvas *c1 = new TCanvas("c1");
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}
	
	if(!gPad) return;
	TCanvas *c1 = gPad->GetCanvas();
	c1->Divide(3, 2);

	c1->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogy();
	if(h1) h1->Draw();

	c1->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogy();
	if(h2) h2->Draw();

	c1->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogy();
	if(h3) h3->Draw();

	c1->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogy();
	if(h4) h4->Draw();

	c1->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();
	if(h5) h5->Draw("colz");

	c1->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	//gPad->SetLogy();
	if(h6) h6->Draw("colz");
}



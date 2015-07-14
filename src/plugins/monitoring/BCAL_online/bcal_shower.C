

// hnamepath: /bcal/bcal_shower_plane
// hnamepath: /bcal/bcal_shower_t
// hnamepath: /bcal/bcal_shower_z   


{
	gStyle->SetOptStat(111110);	
	gStyle->SetStatX(0.99);
	gStyle->SetStatY(0.99);
	gStyle->SetStatW(0.25);

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal");
	if(dir) dir->cd();

	TH1I *h1 = (TH1I*)gDirectory->FindObjectAny("bcal_shower_plane");
	TH1I *h2 = (TH1I*)gDirectory->FindObjectAny("bcal_shower_z");
	TH1I *h3 = (TH1I*)gDirectory->FindObjectAny("bcal_shower_t");

	// Just for testing
	if(gPad == NULL){
		TCanvas *c1 = new TCanvas("c1");
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}
	
	if(!gPad) return;
	TCanvas *c1 = gPad->GetCanvas();
	c1->Divide(2, 2);
	
	c1->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(h1) h1->Draw("colz");

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

}



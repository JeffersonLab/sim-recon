

// hnamepath: /bcal/bcal_point_E
// hnamepath: /bcal/bcal_cluster_E
// hnamepath: /bcal/bcal_point_z
// hnamepath: /bcal/bcal_cluster_nCells

{
	gStyle->SetOptStat(111110);	
	gStyle->SetStatX(0.99);
	gStyle->SetStatY(0.99);
	gStyle->SetStatW(0.25);

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal");
	if(dir) dir->cd();

	TH1I *h1 = (TH1I*)gDirectory->FindObjectAny("bcal_point_E");
	TH1I *h2 = (TH1I*)gDirectory->FindObjectAny("bcal_cluster_E");
	TH1I *h3 = (TH1I*)gDirectory->FindObjectAny("bcal_point_z");
	TH1I *h4 = (TH1I*)gDirectory->FindObjectAny("bcal_cluster_nCells");
	
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
}



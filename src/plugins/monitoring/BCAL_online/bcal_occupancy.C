

// hnamepath: /bcal/bcal_fadc_digi_occ
// hnamepath: /bcal/bcal_fadc_digi_occ_layer1
// hnamepath: /bcal/bcal_fadc_digi_occ_layer2
// hnamepath: /bcal/bcal_fadc_digi_occ_layer3
// hnamepath: /bcal/bcal_fadc_digi_occ_layer4
// hnamepath: /bcal/bcal_tdc_digi_occ
// hnamepath: /bcal/bcal_tdc_digi_occ_layer1
// hnamepath: /bcal/bcal_tdc_digi_occ_layer2
// hnamepath: /bcal/bcal_tdc_digi_occ_layer3



{
	gStyle->SetOptStat(111110);	
	gStyle->SetStatX(0.99);
	gStyle->SetStatY(0.99);
	gStyle->SetStatW(0.25);

	gStyle->SetTitleOffset(1, "Y");
  	gStyle->SetTitleSize(0.05,"xyz");
	gStyle->SetTitleSize(0.08,"h");
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetTitleX(0);
	gStyle->SetTitleAlign(13);
	gStyle->SetNdivisions(505,"xy");

	int col_layer1 = 1;
	int col_layer2 = 2;
	int col_layer3 = 3;
	int col_layer4 = 4;

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal");
	if(dir) dir->cd();

	TH1I *bcal_fadc_digi_occ = (TH1I*)gDirectory->FindObjectAny("bcal_fadc_digi_occ");
	TH1I *bcal_tdc_digi_occ = (TH1I*)gDirectory->FindObjectAny("bcal_tdc_digi_occ");
	TH1I *bcal_fadc_digi_occ_layer1 = (TH1I*)gDirectory->FindObjectAny("bcal_fadc_digi_occ_layer1");
	TH1I *bcal_fadc_digi_occ_layer2 = (TH1I*)gDirectory->FindObjectAny("bcal_fadc_digi_occ_layer2");
	TH1I *bcal_fadc_digi_occ_layer3 = (TH1I*)gDirectory->FindObjectAny("bcal_fadc_digi_occ_layer3");
	TH1I *bcal_fadc_digi_occ_layer4 = (TH1I*)gDirectory->FindObjectAny("bcal_fadc_digi_occ_layer4");
	
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
	if(bcal_fadc_digi_occ) bcal_fadc_digi_occ->Draw("colz");

	c1->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(bcal_tdc_digi_occ) bcal_tdc_digi_occ->Draw("colz");

	c1->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	THStack *occ_stack = new THStack("occ_stack","Occupancy by layer (DBCALDigiHit);global sector  (4 x module + sector);hits");
	TLegend *occ_stack_legend = new TLegend(0.85,0.75,0.99,0.99);
 	if (bcal_fadc_digi_occ_layer1 != NULL) {
		bcal_fadc_digi_occ_layer1->SetLineColor(col_layer1);
		bcal_fadc_digi_occ_layer1->SetMarkerColor(col_layer1);
		bcal_fadc_digi_occ_layer1->SetMarkerStyle(21);
		occ_stack->Add(bcal_fadc_digi_occ_layer1);
		occ_stack_legend->AddEntry(bcal_fadc_digi_occ_layer1,"Layer 1","p");  
	}
	if (bcal_fadc_digi_occ_layer2 != NULL) {
		bcal_fadc_digi_occ_layer2->SetLineColor(col_layer2);
		bcal_fadc_digi_occ_layer2->SetMarkerColor(col_layer2);
		bcal_fadc_digi_occ_layer2->SetMarkerStyle(21);
		occ_stack->Add(bcal_fadc_digi_occ_layer2);
		occ_stack_legend->AddEntry(bcal_fadc_digi_occ_layer2,"Layer 2","p");  
	}
	if (bcal_fadc_digi_occ_layer3 != NULL) {
		bcal_fadc_digi_occ_layer3->SetLineColor(col_layer3);
		bcal_fadc_digi_occ_layer3->SetMarkerColor(col_layer3);
		bcal_fadc_digi_occ_layer3->SetMarkerStyle(21);
		occ_stack->Add(bcal_fadc_digi_occ_layer3);
		occ_stack_legend->AddEntry(bcal_fadc_digi_occ_layer3,"Layer 3","p");  
	}
	if (bcal_fadc_digi_occ_layer4 != NULL) {
		bcal_fadc_digi_occ_layer4->SetLineColor(col_layer4);
		bcal_fadc_digi_occ_layer4->SetMarkerColor(col_layer4);
		bcal_fadc_digi_occ_layer4->SetMarkerStyle(21);
		occ_stack->Add(bcal_fadc_digi_occ_layer4);
		occ_stack_legend->AddEntry(bcal_fadc_digi_occ_layer4,"Layer 4","p");  
	}
	occ_stack->Draw("nostack");
	occ_stack_legend->Draw();




	c1->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	TH1I *bcal_tdc_digi_occ_layer1 = (TH1I*)gROOT->FindObject("bcal_tdc_digi_occ_layer1");
	TH1I *bcal_tdc_digi_occ_layer2 = (TH1I*)gROOT->FindObject("bcal_tdc_digi_occ_layer2");
	TH1I *bcal_tdc_digi_occ_layer3 = (TH1I*)gROOT->FindObject("bcal_tdc_digi_occ_layer3");
	THStack *tdc_occ_stack = new THStack("tdc_occ_stack","TDC Occupancy by layer (DBCALDigiHit);global sector  (4 x module + sector);hits");
	TLegend *tdc_occ_stack_legend = new TLegend(0.85,0.75,0.99,0.99);
 	if (bcal_tdc_digi_occ_layer1 != NULL) {
		bcal_tdc_digi_occ_layer1->SetLineColor(col_layer1);
		bcal_tdc_digi_occ_layer1->SetMarkerColor(col_layer1);
		bcal_tdc_digi_occ_layer1->SetMarkerStyle(21);
		tdc_occ_stack->Add(bcal_tdc_digi_occ_layer1);    
		tdc_occ_stack_legend->AddEntry(bcal_tdc_digi_occ_layer1,"Layer 1","p");  
	}
	if (bcal_tdc_digi_occ_layer2 != NULL) {
		bcal_tdc_digi_occ_layer2->SetLineColor(col_layer2);
		bcal_tdc_digi_occ_layer2->SetMarkerColor(col_layer2);
		bcal_tdc_digi_occ_layer2->SetMarkerStyle(21);
		tdc_occ_stack->Add(bcal_tdc_digi_occ_layer2);
		tdc_occ_stack_legend->AddEntry(bcal_tdc_digi_occ_layer2,"Layer 2","p");  
	}
	if (bcal_tdc_digi_occ_layer3 != NULL) {
		bcal_tdc_digi_occ_layer3->SetLineColor(col_layer3);
		bcal_tdc_digi_occ_layer3->SetMarkerColor(col_layer3);
		bcal_tdc_digi_occ_layer3->SetMarkerStyle(21);
		tdc_occ_stack->Add(bcal_tdc_digi_occ_layer3);
		tdc_occ_stack_legend->AddEntry(bcal_tdc_digi_occ_layer3,"Layer 3","p");  
	}
	tdc_occ_stack->Draw("nostack");
	tdc_occ_stack_legend->Draw();

}



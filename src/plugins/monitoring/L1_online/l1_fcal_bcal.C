

{

  //    gROOT->Reset();
  //    gROOT->SetStyle("Plain");
  //    gROOT->ForceStyle();


        gStyle->SetLineWidth(1.5);
        gStyle->SetTextSize(1.5);

        gStyle->SetTitleFont(132,"xy");
        gStyle->SetLabelFont(62,"xy");


        gStyle->SetOptStat( 0);	
	//      gStyle->SetOptStat(111110);	
	//   	gStyle->SetStatX(0.99);
	//	gStyle->SetStatY(0.99);
	//	gStyle->SetStatW(0.25);

	gStyle->SetTitleSize(0.065,"h");

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("L1");
	if(dir) dir->cd();

	TH2I *fcal_bcal_en1 = (TH2I*)gDirectory->FindObjectAny("fcal_bcal_en1");
	TH1I *fcal_en1      = (TH1I*)gDirectory->FindObjectAny("fcal_en1");
	TH2I *fcal_bcal_en6 = (TH2I*)gDirectory->FindObjectAny("fcal_bcal_en6");
	TH1I *fcal_en6      = (TH1I*)gDirectory->FindObjectAny("fcal_en6");

	
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

	gPad->SetLeftMargin(0.12);
	
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();


	if(fcal_bcal_en1){	 
	  fcal_bcal_en1->SetTitle("TRIG BIT 1");
	  fcal_bcal_en1->GetXaxis()->SetRangeUser(-1,8000.);
	  fcal_bcal_en1->GetYaxis()->SetRangeUser(0.,50000.);
	  fcal_bcal_en1->GetXaxis()->SetTitle("E (FCAL)  (count)");
	  fcal_bcal_en1->GetYaxis()->SetTitle("E (BCAL)  (count)");

	  fcal_bcal_en1->SetTitleSize(0.05,"xy");
	  fcal_bcal_en1->SetTitleOffset(1.2,"y");
	  fcal_bcal_en1->SetTitleOffset(1.,"x");

	  fcal_bcal_en1->Draw("zcol");

	}


	c1->cd(2);

	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogy();
	if(fcal_en1){
	  fcal_en1->SetTitle("TRIG BIT 1");
	  fcal_en1->GetXaxis()->SetRangeUser(0.,5000.);
	  fcal_en1->GetXaxis()->SetTitle("E (FCAL)  (count)");

	  fcal_en1->SetTitleSize(0.05,"xy");
	  fcal_en1->SetTitleOffset(1.2,"y");
	  fcal_en1->SetTitleOffset(1.,"x");

	  fcal_en1->Draw();
	}


	c1->cd(3);

	gPad->SetLeftMargin(0.12);

	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();
	if(fcal_bcal_en6){

	  fcal_bcal_en6->SetTitle("TRIG BIT 6");

	  fcal_bcal_en6->GetXaxis()->SetRangeUser(-1,8000.);
	  fcal_bcal_en6->GetYaxis()->SetRangeUser(0.,50000.);
	  fcal_bcal_en6->GetXaxis()->SetTitle("E (FCAL)  (count)");
	  fcal_bcal_en6->GetYaxis()->SetTitle("E (BCAL)  (count)");

	  fcal_bcal_en6->SetTitleSize(0.05,"xy");
	  fcal_bcal_en6->SetTitleOffset(1.2,"y");
	  fcal_bcal_en6->SetTitleOffset(1.,"x");

	  fcal_bcal_en6->Draw("zcol");

	}


	c1->cd(4);

	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogy();
	if(fcal_en6){

	  fcal_en6->SetTitle("TRIG BIT 6");
	  fcal_en6->GetXaxis()->SetRangeUser(0.,5000.);
	  fcal_en6->GetXaxis()->SetTitle("E (FCAL)  (count)");

	  fcal_en6->SetTitleSize(0.05,"xy");
	  fcal_en6->SetTitleOffset(1.2,"y");
	  fcal_en6->SetTitleOffset(1.,"x");

	  fcal_en6->Draw();

	}

}



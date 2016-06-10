// hnamepath:  /L1/rate_gtp_0
// hnamepath:  /L1/rate_gtp_1
// hnamepath:  /L1/rate_gtp_2
// hnamepath:  /L1/rate_gtp_3
// hnamepath:  /L1/rate_gtp_5
// hnamepath:  /L1/rate_gtp_6

{
 
  //    gROOT->Reset();
  //    gROOT->SetStyle("Plain");
  //    gROOT->ForceStyle();


        gStyle->SetLineWidth(1.5);
        gStyle->SetTextSize(1.5);

        gStyle->SetTitleFont(132,"xy");
        gStyle->SetLabelFont(62,"xy");


        gStyle->SetOptStat( 0);	
	gStyle->SetOptStat(111110);	
	gStyle->SetStatX(0.99);
	gStyle->SetStatY(0.99);
	gStyle->SetStatW(0.25);

	gStyle->SetTitleSize(0.075,"h");


	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("L1");
	if(dir) dir->cd();

	TH1I *rate_gtp_0 = (TH1I*)gDirectory->FindObjectAny("rate_gtp_0");
	TH1I *rate_gtp_1 = (TH1I*)gDirectory->FindObjectAny("rate_gtp_1");
	TH1I *rate_gtp_2 = (TH1I*)gDirectory->FindObjectAny("rate_gtp_2");
	TH1I *rate_gtp_3 = (TH1I*)gDirectory->FindObjectAny("rate_gtp_3");
	TH1I *rate_gtp_5 = (TH1I*)gDirectory->FindObjectAny("rate_gtp_5");
	TH1I *rate_gtp_6 = (TH1I*)gDirectory->FindObjectAny("rate_gtp_6");

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
	gPad->SetLogz();

	if(rate_gtp_0){	 
	  rate_gtp_0->SetTitle("TRIG BIT 1");
	  rate_gtp_0->GetXaxis()->SetRangeUser(0,50);

	  rate_gtp_0->SetLabelSize(0.05,"xy");
	  rate_gtp_0->SetTitleSize(0.06,"xy");
	  rate_gtp_0->SetTitleOffset(0.8,"x");

	  rate_gtp_0->Draw();
	}


	c1->cd(2);

	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();

	if(rate_gtp_1){	 
	  rate_gtp_1->SetTitle("TRIG BIT 2");
	  rate_gtp_1->GetXaxis()->SetRangeUser(0,50);

	  rate_gtp_1->SetLabelSize(0.05,"xy");
	  rate_gtp_1->SetTitleSize(0.06,"xy");
	  rate_gtp_1->SetTitleOffset(0.8,"x");

	  rate_gtp_1->Draw();
	}


	c1->cd(3);

	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();

	if(rate_gtp_2){	 
	  rate_gtp_2->SetTitle("TRIG BIT 3");
	  rate_gtp_2->GetXaxis()->SetRangeUser(0,50);

	  rate_gtp_2->SetLabelSize(0.05,"xy");
	  rate_gtp_2->SetTitleSize(0.06,"xy");
	  rate_gtp_2->SetTitleOffset(0.8,"x");

	  rate_gtp_2->Draw();
	}

	c1->cd(4);

	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();

	if(rate_gtp_3){	 
	  rate_gtp_3->SetTitle("TRIG BIT 4");
	  rate_gtp_3->GetXaxis()->SetRangeUser(0,6);
	  rate_gtp_3->GetXaxis()->SetTitle("Rate  (kHz)");

	  rate_gtp_3->SetLabelSize(0.05,"xy");
	  rate_gtp_3->SetTitleSize(0.06,"xy");
	  rate_gtp_3->SetTitleOffset(0.8,"x");

	  rate_gtp_3->Draw();
	}


	c1->cd(5);

	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();

	if(rate_gtp_5){	 
	  rate_gtp_5->SetTitle("TRIG BIT 6");
	  rate_gtp_5->GetXaxis()->SetRangeUser(0,50);
	  rate_gtp_5->GetXaxis()->SetTitle("Rate  (kHz)");

	  rate_gtp_5->SetLabelSize(0.05,"xy");
	  rate_gtp_5->SetTitleSize(0.06,"xy");
	  rate_gtp_5->SetTitleOffset(0.8,"x");

	  rate_gtp_5->Draw();
	}


	c1->cd(6);

	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogz();

	if(rate_gtp_6){	 
	  rate_gtp_6->SetTitle("TRIG BIT 7");
	  rate_gtp_6->GetXaxis()->SetRangeUser(0,50);
	  rate_gtp_6->GetXaxis()->SetTitle("Rate  (kHz)");

	  rate_gtp_6->SetLabelSize(0.05,"xy");
	  rate_gtp_6->SetTitleSize(0.06,"xy");
	  rate_gtp_6->SetTitleOffset(0.8,"x");

	  rate_gtp_6->Draw();
	}



}



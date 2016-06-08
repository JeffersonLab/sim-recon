

{

  //    gROOT->Reset();
  //    gROOT->SetStyle("Plain");
  //    gROOT->ForceStyle();


        gStyle->SetLineWidth(1.5);
        gStyle->SetTextSize(1.5);

        gStyle->SetTitleFont(132,"xy");
        gStyle->SetLabelFont(62,"xy");


        gStyle->SetOptStat( 0);	
	gStyle->SetOptStat(110);	
	gStyle->SetStatX(0.99);
	gStyle->SetStatY(0.99);
	gStyle->SetStatW(0.25);

	gStyle->SetTitleSize(0.08,"h");

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("L1");
	if(dir) dir->cd();

	TH1I *fcal_en7     = (TH1I*)gDirectory->FindObjectAny("fcal_en7");
	TH1I *st_hit7      = (TH1I*)gDirectory->FindObjectAny("st_hit7");

	TH1I *bcal_en3     = (TH1I*)gDirectory->FindObjectAny("bcal_en3");
	TH1I *fcal_en2     = (TH1I*)gDirectory->FindObjectAny("fcal_en2");

	TH1I *tagh_occup2  = (TH1I*)gDirectory->FindObjectAny("tagh_occup2");
	TH1I *st_hit2      = (TH1I*)gDirectory->FindObjectAny("st_hit2");

	
	// Just for testing
	if(gPad == NULL){
		TCanvas *c1 = new TCanvas("c1");
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}
	

	if(!gPad) return;


	TCanvas *c1 = gPad->GetCanvas();

	c1->Divide(2, 3);

	c1->cd(1);

	gPad->SetBottomMargin(0.14);
	
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogy();


	if(fcal_en7){	 
	  fcal_en7->SetTitle("TRIG BIT 7");

	  fcal_en7->GetXaxis()->SetRangeUser(0,5000.);
	  fcal_en7->GetXaxis()->SetTitle("E (FCAL)  (count)");

	  fcal_en7->SetLabelSize(0.06,"xy");
	  fcal_en7->SetTitleSize(0.07,"xy");
	  fcal_en7->SetTitleOffset(1.0,"x");

	  fcal_en7->Draw();

	}


	c1->cd(2);

	gPad->SetBottomMargin(0.14);
	
	gPad->SetTicks();
	gPad->SetGrid();
	//	gPad->SetLogy();


	if(st_hit7){	 
	  st_hit7->SetTitle("TRIG BIT 7");
	  st_hit7->GetXaxis()->SetRangeUser(0,15);
	  st_hit7->GetXaxis()->SetTitle("ST hits ");

	  st_hit7->SetLabelSize(0.06,"xy");
	  st_hit7->SetTitleSize(0.07,"xy");
	  st_hit7->SetTitleOffset(1.0,"x");

	  st_hit7->Draw();

	}


	c1->cd(3);

	gPad->SetBottomMargin(0.14);
	
	gPad->SetTicks();
	gPad->SetGrid();
	gPad->SetLogy();


	if(bcal_en3){	 

	  bcal_en3->SetTitle("TRIG BIT 3");
	  bcal_en3->GetXaxis()->SetRangeUser(0,45000.);
	  bcal_en3->GetXaxis()->SetTitle("E (BCAL)  (count)");

	  bcal_en3->SetLabelSize(0.06,"xy");
	  bcal_en3->SetTitleSize(0.07,"xy");
	  bcal_en3->SetTitleOffset(1.0,"x");

	  bcal_en3->Draw();

	}


	c1->cd(4);

	gPad->SetBottomMargin(0.14);
	
	gPad->SetTicks();
	gPad->SetGrid();
	//	gPad->SetLogy();


	if(fcal_en2){	 

	  fcal_en2->SetTitle("TRIG BIT 2 (FCAL)");
	  fcal_en2->GetXaxis()->SetRangeUser(0,5000.);
	  fcal_en2->GetXaxis()->SetTitle("E (FCAL)  (count)");

	  fcal_en2->SetLabelSize(0.06,"xy");
	  fcal_en2->SetTitleSize(0.07,"xy");
	  fcal_en2->SetTitleOffset(1.0,"x");

	  fcal_en2->Draw();

	}


	c1->cd(5);

	gPad->SetBottomMargin(0.14);
	
	gPad->SetTicks();
	gPad->SetGrid();
	//	gPad->SetLogy();


	if(tagh_occup2){	 

	  tagh_occup2->SetTitle("TRIG BIT 2 (TAGH & ST)");
	  tagh_occup2->GetXaxis()->SetRangeUser(0,200.);
	  tagh_occup2->GetXaxis()->SetTitle("TAGH counter");

	  tagh_occup2->SetLabelSize(0.06,"xy");
	  tagh_occup2->SetTitleSize(0.07,"xy");
	  tagh_occup2->SetTitleOffset(1.0,"x");

	  tagh_occup2->Draw();

	}


	c1->cd(6);

	gPad->SetBottomMargin(0.14);
	
	gPad->SetTicks();
	gPad->SetGrid();
	//	gPad->SetLogy();


	if(st_hit2){	 

	  st_hit2->SetTitle("TRIG BIT 2 (TAGH & ST)");
	  st_hit2->GetXaxis()->SetRangeUser(0,15.);
	  st_hit2->GetXaxis()->SetTitle("ST hits");

	  st_hit2->SetLabelSize(0.06,"xy");
	  st_hit2->SetTitleSize(0.07,"xy");
	  st_hit2->SetTitleOffset(1.0,"x");

	  st_hit2->Draw();

	}

}

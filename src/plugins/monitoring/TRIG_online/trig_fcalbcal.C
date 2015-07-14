
// File: trig_fcalbcal.C
// Created: 05/01/15
// Creator: Elton Smith
// Purpose: Display basic energy plots in fcal and bcal for trigger
// Uses output of TRIG_online/JEventProcessor_TRIG_online.cc

// hnamepath: /trig/h1trig_fcal
// hnamepath: /trig/h1trig_fcalN
// hnamepath: /trig/h1trig_bcal
// hnamepath: /trig/h1trig_bcalN
// hnamepath: /trig/h2trig_fcalVSbcal

{
	gStyle->SetPalette(1,0);
	gStyle->SetOptStat(kTRUE);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);

	// get histograms from trig subdirectory

   TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("trig");
   if(!dir){
		cout << "Can't find trig TDirectory!" << endl;
		return;
	}

   TH1F *h1trig_fcal         = (TH1F*)dir->Get("h1trig_fcal"         );
   TH1F *h1trig_fcalN        = (TH1F*)dir->Get("h1trig_fcalN"         );
   TH1F *h1trig_bcal         = (TH1F*)dir->Get("h1trig_bcal"         );
   TH1F *h1trig_bcalN        = (TH1F*)dir->Get("h1trig_bcalN"         );
   TH1F *h2trig_fcalVSbcal   = (TH1F*)dir->Get("h2trig_fcalVSbcal"         );
	
   if( !h1trig_fcal         ) { cout << "Can't find h1trig_fcal!"         << endl; return; }
   if( !h1trig_fcalN        ) { cout << "Can't find h1trig_fcalN!"         << endl; return; }
   if( !h1trig_bcal         ) { cout << "Can't find h1trig_bcal!"         << endl; return; }
   if( !h1trig_bcalN        ) { cout << "Can't find h1trig_bcalN!"         << endl; return; }
   if( !h2trig_fcalVSbcal   ) { cout << "Can't find h2trig_fcalVSbcal!"         << endl; return; }

   //
   if(gPad == NULL){
	  TCanvas *c1 = new TCanvas("c1","c1 trig fcal/bcal ",200,10,1200,800);
	  c1->cd(0);
	  c1->Draw();
	  c1->Update();
	}
	
   if(!gPad) return;
   TCanvas *c1 = gPad->GetCanvas();

   c1->SetBorderMode(0);
   c1->SetFillColor(0);

   /*Double_t ymin=0.5;
   Double_t ymin=0.;
   Double_t ymax=1.1;*/

   c1->SetGridx();
   c1->SetGridy();
   c1->SetBorderMode(0);
   c1->SetFillColor(0);

   c1->Divide(2,2);
   c1->cd(1);

   h1trig_fcal->SetTitle("");
   // h1trig_fcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig_fcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig_fcal->GetXaxis()->SetTitleSize(0.05);
   h1trig_fcal->GetYaxis()->SetTitleSize(0.05);
   h1trig_fcal->GetXaxis()->SetTitle("Fcal Energy (GeV)");
   h1trig_fcal->SetLineColor(2);
   h1trig_fcal->Draw("");

   c1->cd(2);

   h2trig_fcalVSbcal->SetTitle("");
   // h2trig_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig_fcalVSbcal->SetLineColor(2);
   h2trig_fcalVSbcal->Draw("colz");

   c1->cd(3);

   h1trig_bcal->SetTitle("");
   // h1trig_bcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig_bcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig_bcal->GetXaxis()->SetTitleSize(0.05);
   h1trig_bcal->GetYaxis()->SetTitleSize(0.05);
   h1trig_bcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h1trig_bcal->SetLineColor(2);
   h1trig_bcal->Draw("");

   c1->cd(4); 
   c1_4->SetLogz();

   h2trig_fcalVSbcal->SetTitle("");
   // h2trig_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig_fcalVSbcal->SetLineColor(2);
   h2trig_fcalVSbcal->Draw("colz");

}



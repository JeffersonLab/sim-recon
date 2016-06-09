
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

  char string[256];
    
  TString filename = "hd_rawdata_010270_000";
  TFile* f = new TFile(filename+".root");

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

   TH1F *h1trig_trgbits         = (TH1F*)dir->Get("h1trig_trgbits"         );

   TH1F *h1trig1_fcal         = (TH1F*)dir->Get("h1trig1_fcal"         );
   TH1F *h1trig1_fcalN        = (TH1F*)dir->Get("h1trig1_fcalN"         );
   TH1F *h1trig1_bcal         = (TH1F*)dir->Get("h1trig1_bcal"         );
   TH1F *h1trig1_bcalN        = (TH1F*)dir->Get("h1trig1_bcalN"         );
   TH1F *h2trig1_fcalVSbcal   = (TH1F*)dir->Get("h2trig1_fcalVSbcal"         );

   TH1F *h1trig3_fcal         = (TH1F*)dir->Get("h1trig3_fcal"         );
   TH1F *h1trig3_fcalN        = (TH1F*)dir->Get("h1trig3_fcalN"         );
   TH1F *h1trig3_bcal         = (TH1F*)dir->Get("h1trig3_bcal"         );
   TH1F *h1trig3_bcalN        = (TH1F*)dir->Get("h1trig3_bcalN"         );
   TH1F *h2trig3_fcalVSbcal   = (TH1F*)dir->Get("h2trig3_fcalVSbcal"         );

   TH1F *h1trig5_fcal         = (TH1F*)dir->Get("h1trig5_fcal"         );
   TH1F *h1trig5_fcalN        = (TH1F*)dir->Get("h1trig5_fcalN"         );
   TH1F *h1trig5_bcal         = (TH1F*)dir->Get("h1trig5_bcal"         );
   TH1F *h1trig5_bcalN        = (TH1F*)dir->Get("h1trig5_bcalN"         );
   TH1F *h2trig5_fcalVSbcal   = (TH1F*)dir->Get("h2trig5_fcalVSbcal"         );

   TH1F *h1trig7_fcal         = (TH1F*)dir->Get("h1trig7_fcal"         );
   TH1F *h1trig7_fcalN        = (TH1F*)dir->Get("h1trig7_fcalN"         );
   TH1F *h1trig7_bcal         = (TH1F*)dir->Get("h1trig7_bcal"         );
   TH1F *h1trig7_bcalN        = (TH1F*)dir->Get("h1trig7_bcalN"         );
   TH1F *h2trig7_fcalVSbcal   = (TH1F*)dir->Get("h2trig7_fcalVSbcal"         );

   double fcal_max = 2;
   double bcal_max = 1;
   double thf = 0.7;
   double thb = 0.2;

   TCanvas *c0 = new TCanvas("c0", "c0",200,10,700,700);

   c0->Divide(2,2);
   c0->cd(1);

   h1trig_trgbits->SetTitle(filename);
   // h1trig_trgbits->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig_trgbits->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig_trgbits->GetXaxis()->SetTitleSize(0.05);
   h1trig_trgbits->GetYaxis()->SetTitleSize(0.05);
   h1trig_trgbits->GetXaxis()->SetTitle("trig_mask || (10+fp_trig_mask)");
   h1trig_trgbits->SetLineColor(2);
   h1trig_trgbits->Draw("");

   c0->cd(2);
   gPad->SetLogy();

   // h1trig_trgbits->SetTitle("");
   // h1trig_trgbits->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig_trgbits->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig_trgbits->GetXaxis()->SetTitleSize(0.05);
   h1trig_trgbits->GetYaxis()->SetTitleSize(0.05);
   h1trig_trgbits->GetXaxis()->SetTitle("trig_mask || (10+fp_trig_mask)");
   h1trig_trgbits->SetLineColor(2);
   h1trig_trgbits->Draw("");
    
    c0->cd(3);

    
    TH2F *temp = new TH2F("temp","Template",10,0,1,10,0,2);
    temp->GetXaxis()->SetTitleSize(0.05);
    temp->GetYaxis()->SetTitleSize(0.05);
    temp->GetXaxis()->SetTitle("Bcal Energy (GeV)");
    temp->GetYaxis()->SetTitle("Fcal Energy (GeV)");
    temp->Draw();
    
    TLine *thresh_fb = new TLine(0,0.5,0.14,0);
    thresh_fb->SetLineWidth(2);
    thresh_fb->Draw();
    
    TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
    thresh_f->SetLineWidth(2);
    thresh_f->SetLineStyle(2);
    thresh_f->Draw();
    
    TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
    thresh_b->SetLineWidth(2);
    thresh_b->SetLineStyle(2);
    thresh_b->Draw();
    
    
    sprintf (string,"trg 5 = 101\n");
    printf("string=%s",string);
    TLatex *t1 = new TLatex(0.5,0.4,string);
    // t1->SetNDC();
    t1->SetTextSize(0.04);
    t1->Draw();
    
    sprintf (string,"trg 7 = 111\n");
    t1->DrawLatex(0.5,1.4,string);
    
    sprintf (string,"trg 3=\n");
    t1->DrawLatex(0.02,1.4,string);
    
    sprintf (string,"011\n");
    t1->DrawLatex(0.01,1.2,string);
    
    sprintf (string,"trg 1=\n");
    t1->DrawLatex(0.01,0.6,string);
    
    sprintf (string,"001\n");
    t1->DrawLatex(0.1,0.4,string);
    
    c0->cd(4);
    
    
    sprintf (string,"Trigger configuration\n");
    t1->DrawLatex(0.1,0.9,string);
    sprintf (string,"001 : 10Fcal +3Bcal > 18000 counts\n");
    t1->DrawLatex(0.1,0.8,string);
    sprintf (string,"\t\t: Fcal + 3.6 Bcal > 0.49 GeV \n");
    t1->DrawLatex(0.1,0.75,string);
    sprintf (string,"010 : Fcal > 1800 counts : Fcal > 0.49 GeV \n");
    t1->DrawLatex(0.1,0.65,string);
    sprintf (string,"100 : Bcal > 6000 counts: Bcal > 0.14 GeV\n");
    t1->DrawLatex(0.1,0.6,string);
    
    
   TCanvas *c1 = new TCanvas("c1", "c1",200,10,700,700);

   c1->Divide(2,2);
   c1->cd(1);

   // h1trig_fcal->SetTitle("");
   // h1trig_fcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig_fcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig_fcal->GetXaxis()->SetTitleSize(0.05);
   h1trig_fcal->GetYaxis()->SetTitleSize(0.05);
   h1trig_fcal->GetXaxis()->SetTitle("Fcal Energy (GeV)");
   h1trig_fcal->SetLineColor(2);
    h1trig_fcal->Draw("");
    
    double ymax;
    ymax = h1trig_fcal->GetMaximum();
    TLine *thrf = new TLine(thf,0,thf,ymax);
    thrf->SetLineWidth(2);
    thrf->SetLineStyle(2);
    thrf->DrawLine(thf,0,thf,ymax);
    
    
   c1->cd(2);

   temp->SetTitle("");
   // h2trig_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig_fcalVSbcal->SetLineColor(2);
   h2trig_fcalVSbcal->Draw("colz");

   thresh_fb->Draw();
   thresh_f->Draw();
   thresh_b->Draw();

   c1->cd(3);

   h1trig_bcal->SetTitle("");
   // h1trig_bcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig_bcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig_bcal->GetXaxis()->SetTitleSize(0.05);
   h1trig_bcal->GetYaxis()->SetTitleSize(0.05);
   h1trig_bcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h1trig_bcal->SetLineColor(2);
    h1trig_bcal->Draw("");
    
    ymax = h1trig_bcal->GetMaximum();
    TLine *thrb = new TLine(thb,0,thb,ymax);
    thrb->SetLineWidth(2);
    thrb->SetLineStyle(2);
    thrb->DrawLine(thb,0,thb,ymax);

   c1->cd(4); 
   gPad->SetLogz();

   h2trig_fcalVSbcal->SetTitle("");
   // h2trig_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig_fcalVSbcal->SetLineColor(2);
   h2trig_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();

   TCanvas *c2 = new TCanvas("c2", "c2",200,10,700,700);

   c2->Divide(2,2);
   c2->cd(1);

   // h1trig1_fcal->SetTitle("");
   // h1trig1_fcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig1_fcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig1_fcal->GetXaxis()->SetTitleSize(0.05);
   h1trig1_fcal->GetYaxis()->SetTitleSize(0.05);
   h1trig1_fcal->GetXaxis()->SetTitle("Fcal Energy (GeV)");
   h1trig1_fcal->SetLineColor(2);
    h1trig1_fcal->Draw("");
    
    ymax = h1trig1_fcal->GetMaximum();
    thrf->DrawLine(thf,0,thf,ymax);

   c2->cd(2);

   h2trig1_fcalVSbcal->SetTitle("");
   // h2trig1_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig1_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig1_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig1_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig1_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig1_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig1_fcalVSbcal->SetLineColor(2);
   h2trig1_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();

   c2->cd(3);

   h1trig1_bcal->SetTitle("");
   // h1trig1_bcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig1_bcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig1_bcal->GetXaxis()->SetTitleSize(0.05);
   h1trig1_bcal->GetYaxis()->SetTitleSize(0.05);
   h1trig1_bcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h1trig1_bcal->SetLineColor(2);
   h1trig1_bcal->Draw("");
    
    ymax = h1trig1_bcal->GetMaximum();
    thrb->DrawLine(thb,0,thb,ymax);
    
   c2->cd(4); 
   gPad->SetLogz();

   h2trig1_fcalVSbcal->SetTitle("");
   // h2trig1_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig1_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig1_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig1_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig1_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig1_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig1_fcalVSbcal->SetLineColor(2);
   h2trig1_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();

   TCanvas *c3 = new TCanvas("c3", "c3",200,10,700,700);

   c3->Divide(2,2);
   c3->cd(1);

   // h1trig3_fcal->SetTitle("");
   // h1trig3_fcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig3_fcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig3_fcal->GetXaxis()->SetTitleSize(0.05);
   h1trig3_fcal->GetYaxis()->SetTitleSize(0.05);
   h1trig3_fcal->GetXaxis()->SetTitle("Fcal Energy (GeV)");
   h1trig3_fcal->SetLineColor(2);
    h1trig3_fcal->Draw("");
    
    ymax = h1trig3_fcal->GetMaximum();
    thrf->DrawLine(thf,0,thf,ymax);
    
   c3->cd(2);

   h2trig3_fcalVSbcal->SetTitle("");
   // h2trig3_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig3_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig3_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig3_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig3_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig3_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig3_fcalVSbcal->SetLineColor(2);
   h2trig3_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();

   c3->cd(3);

   h1trig3_bcal->SetTitle("");
   // h1trig3_bcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig3_bcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig3_bcal->GetXaxis()->SetTitleSize(0.05);
   h1trig3_bcal->GetYaxis()->SetTitleSize(0.05);
   h1trig3_bcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h1trig3_bcal->SetLineColor(2);
    h1trig3_bcal->Draw("");
    
    ymax = h1trig3_bcal->GetMaximum();
    thrb->DrawLine(thb,0,thb,ymax);
    

   c3->cd(4); 
   gPad->SetLogz();

   h2trig3_fcalVSbcal->SetTitle("");
   // h2trig3_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig3_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig3_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig3_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig3_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig3_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig3_fcalVSbcal->SetLineColor(2);
   h2trig3_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();

   TCanvas *c4 = new TCanvas("c4", "c4",200,10,700,700);

   c4->Divide(2,2);
   c4->cd(1);

   // h1trig5_fcal->SetTitle("");
   // h1trig5_fcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig5_fcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig5_fcal->GetXaxis()->SetTitleSize(0.05);
   h1trig5_fcal->GetYaxis()->SetTitleSize(0.05);
   h1trig5_fcal->GetXaxis()->SetTitle("Fcal Energy (GeV)");
   h1trig5_fcal->SetLineColor(2);
   h1trig5_fcal->Draw("");
    
    ymax = h1trig5_fcal->GetMaximum();
    thrf->DrawLine(thf,0,thf,ymax);
    
   c4->cd(2);

   h2trig5_fcalVSbcal->SetTitle("");
   // h2trig5_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig5_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig5_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig5_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig5_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig5_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig5_fcalVSbcal->SetLineColor(2);
   h2trig5_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();

   c4->cd(3);

   h1trig5_bcal->SetTitle("");
   // h1trig5_bcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig5_bcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig5_bcal->GetXaxis()->SetTitleSize(0.05);
   h1trig5_bcal->GetYaxis()->SetTitleSize(0.05);
   h1trig5_bcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h1trig5_bcal->SetLineColor(2);
    h1trig5_bcal->Draw("");
    
    ymax = h1trig5_bcal->GetMaximum();
    thrb->DrawLine(thb,0,thb,ymax);
    

   c4->cd(4); 
   gPad->SetLogz();

   h2trig5_fcalVSbcal->SetTitle("");
   // h2trig5_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig5_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig5_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig5_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig5_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig5_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig5_fcalVSbcal->SetLineColor(2);
   h2trig5_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();



   TCanvas *c5 = new TCanvas("c5", "c5",200,10,700,700);

   c5->Divide(2,2);
   c5->cd(1);

   // h1trig7_fcal->SetTitle("");
   // h1trig7_fcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig7_fcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig7_fcal->GetXaxis()->SetTitleSize(0.05);
   h1trig7_fcal->GetYaxis()->SetTitleSize(0.05);
   h1trig7_fcal->GetXaxis()->SetTitle("Fcal Energy (GeV)");
   h1trig7_fcal->SetLineColor(2);
   h1trig7_fcal->Draw("");
    
    
    ymax = h1trig7_fcal->GetMaximum();
    thrf->DrawLine(thf,0,thf,ymax);
    
   c5->cd(2);

   h2trig7_fcalVSbcal->SetTitle("");
   // h2trig7_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig7_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig7_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig7_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig7_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig7_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig7_fcalVSbcal->SetLineColor(2);
   h2trig7_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();

   c5->cd(3);

   h1trig7_bcal->SetTitle("");
   // h1trig7_bcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1trig7_bcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h1trig7_bcal->GetXaxis()->SetTitleSize(0.05);
   h1trig7_bcal->GetYaxis()->SetTitleSize(0.05);
   h1trig7_bcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h1trig7_bcal->SetLineColor(2);
    h1trig7_bcal->Draw("");
    
    ymax = h1trig7_bcal->GetMaximum();
    thrb->DrawLine(thb,0,thb,ymax);
    

   c5->cd(4); 
   gPad->SetLogz();

   h2trig7_fcalVSbcal->SetTitle("");
   // h2trig7_fcalVSbcal->GetXaxis()->SetRangeUser(xmin,xmax);
   // h2trig7_fcalVSbcal->GetYaxis()->SetRangeUser(ymin,ymax);
   h2trig7_fcalVSbcal->GetXaxis()->SetTitleSize(0.05);
   h2trig7_fcalVSbcal->GetYaxis()->SetTitleSize(0.05);
   h2trig7_fcalVSbcal->GetXaxis()->SetTitle("Bcal Energy (GeV)");
   h2trig7_fcalVSbcal->GetYaxis()->SetTitle("Fcal Energy (GeV)");
   h2trig7_fcalVSbcal->SetLineColor(2);
   h2trig7_fcalVSbcal->Draw("colz");

   // TLine *thresh_fb = new TLine(0,0.5,0.14,0);
   thresh_fb->SetLineWidth(2);
   thresh_fb->Draw();

   // TLine *thresh_f = new TLine(0,thf,bcal_max,thf);
   thresh_f->SetLineWidth(2);
   thresh_f->SetLineStyle(2);
   thresh_f->Draw();

   // TLine *thresh_b = new TLine(thb,0,thb,fcal_max);
   thresh_b->SetLineWidth(2);
   thresh_b->SetLineStyle(2);
   thresh_b->Draw();

   sprintf (string,"trig_fcalbcal2.pdf(");
   c0->SaveAs(string);
   sprintf (string,"trig_fcalbcal2.pdf");
   c1->SaveAs(string);
   sprintf (string,"trig_fcalbcal2.pdf");
   c2->SaveAs(string);
   sprintf (string,"trig_fcalbcal2.pdf");
   c3->SaveAs(string);
   sprintf (string,"trig_fcalbcal2.pdf");
   c4->SaveAs(string);
   sprintf (string,"trig_fcalbcal2.pdf)");
   c5->SaveAs(string);



}



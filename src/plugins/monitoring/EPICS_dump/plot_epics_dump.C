// File: plot_epics_dump.C
// Adapt trig_fcalbcal3.C to this application
// Include livetime histograms from synch events

{
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(kTRUE);
    gStyle->SetOptStat(111111);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

    char string[256];
    Int_t const nscalers=32;
    
    TString filename = "outputfile_010707d";
    TFile* f = new TFile(filename+".root");

	// get histograms from trig subdirectory

   TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("EPICS_dump");
   if(!dir){
		cout << "Can't find EPICS_dump TDirectory!" << endl;
		return;
	}

   TH1F *h1epics_trgbits         = (TH1F*)dir->Get("h1epics_trgbits"         );
   TH1F *h1epics_AD00            = (TH1F*)dir->Get("h1epics_AD00"         );
   TH1F *h1epics_pos_inner        = (TH1F*)dir->Get("h1epics_pos_inner"         );
   TH1F *h1epics_pos_outer        = (TH1F*)dir->Get("h1epics_pos_outer"         );
    
   if( !h1epics_trgbits         ) { cout << "Can't find h1epics_trgbits!"         << endl; return; }
   if( !h1epics_AD00        ) { cout << "Can't find h1epics_AD00"         << endl; return; }
   if( !h1epics_pos_inner         ) { cout << "Can't find h1epics_pos_inner!"         << endl; return; }
   if( !h1epics_pos_outer       ) { cout << "Can't find h1epics_pos_outer!"         << endl; return; }
    
    
    TH1F *h1epics_liveinst_VSevent      = (TH1F*)dir->Get("h1epics_liveinst_VSevent");
    TH1F *h1epics_AD00_VSevent       = (TH1F*)dir->Get("h1epics_AD00_VSevent");
    TH1F *h1epics_entries1_VSevent      = (TH1F*)dir->Get("h1epics_entries1_VSevent");
    TH1F *h1epics_entries2_VSevent       = (TH1F*)dir->Get("h1epics_entries2_VSevent");
    TH1F *h1epics_liveinst       = (TH1F*)dir->Get("h1epics_liveinst");
    TH1I *h1_trig_rates[nscalers];
    TH1I *h1_trig_livetimes[nscalers];
    
    for (Int_t j=0; j<nscalers;j++) {
        sprintf (string,"Rates%d(kHz)",j);
        h1_trig_rates[j] = (TH1I*)dir->Get(string);
        sprintf (string,"Livetimes%d",j);
        h1_trig_livetimes[j] = (TH1I*)dir->Get(string);
    }

   TCanvas *c0 = new TCanvas("c0", "c0",200,10,700,700);

   c0->Divide(2,2);
    c0->cd(1);
    gPad->SetLogy();
        
   h1epics_trgbits->SetTitle(filename);
   // h1epics_trgbits->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1epics_trgbits->GetYaxis()->SetRangeUser(ymin,ymax);
   h1epics_trgbits->GetXaxis()->SetTitleSize(0.05);
   h1epics_trgbits->GetYaxis()->SetTitleSize(0.05);
   h1epics_trgbits->GetXaxis()->SetTitle("trig_mask || (20+fp_trig_mask/256)");
   h1epics_trgbits->SetLineColor(2);
    h1epics_trgbits->Draw("");
    
    c0->cd(2);
    
    h1epics_AD00->SetTitle(filename);
    // h1epics_AD00->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1epics_AD00->GetYaxis()->SetRangeUser(ymin,ymax);
    h1epics_AD00->GetXaxis()->SetTitleSize(0.05);
    h1epics_AD00->GetYaxis()->SetTitleSize(0.05);
    h1epics_AD00->GetXaxis()->SetTitle("Electron Current AD00 (nA)");
    h1epics_AD00->SetLineColor(2);
    h1epics_AD00->Draw("");
    
   c0->cd(3);

   // h1epics_pos_inner->SetTitle("");
   // h1epics_pos_inner->GetXaxis()->SetRangeUser(xmin,xmax);
   // h1epics_pos_inner->GetYaxis()->SetRangeUser(ymin,ymax);
   h1epics_pos_inner->GetXaxis()->SetTitleSize(0.05);
   h1epics_pos_inner->GetYaxis()->SetTitleSize(0.05);
   // h1epics_pos_inner->GetXaxis()->SetTitle("Y vs X (Active outer)");
   h1epics_pos_inner->SetLineColor(2);
   h1epics_pos_inner->Draw("colz");
    
    c0->cd(4);
    
    // h1epics_pos_outer->SetTitle("");
    // h1epics_pos_outer->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1epics_pos_outer->GetYaxis()->SetRangeUser(ymin,ymax);
    h1epics_pos_outer->GetXaxis()->SetTitleSize(0.05);
    h1epics_pos_outer->GetYaxis()->SetTitleSize(0.05);
    // h1epics_pos_outer->GetXaxis()->SetTitle("Y vs X (Active outer)");
    h1epics_pos_outer->SetLineColor(2);
    h1epics_pos_outer->Draw("colz");
    
    
    TCanvas *c1 = new TCanvas("c1", "c1",200,10,700,700);
    
    // h1epics_liveinst->SetTitle("");
    // h1epics_liveinst->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1epics_liveinst->GetYaxis()->SetRangeUser(ymin,ymax);
    h1epics_liveinst->GetXaxis()->SetTitleSize(0.05);
    h1epics_liveinst->GetYaxis()->SetTitleSize(0.05);
    h1epics_liveinst->GetXaxis()->SetTitle("Inst live time");
    h1epics_liveinst->SetLineColor(2);
    h1epics_liveinst->Draw("colz");
    
    
    TCanvas *c2 = new TCanvas("c2", "c2",200,10,700,700);
    c2->Divide(2,2);
    
    for (Int_t j=0; j<4;j++) {
        c2->cd(j+1);
        h1_trig_rates[j]->Draw();
        sprintf (string,"Rate (kHz) for bit %d",j+1);
        if (j == 0) {
            h1_trig_rates[j]->SetTitle(filename);
        }
        else {
            h1_trig_rates[j]->SetTitle("");
        }
        h1_trig_rates[j]->GetXaxis()->SetTitleSize(0.05);
        h1_trig_rates[j]->GetYaxis()->SetTitleSize(0.05);
        h1_trig_rates[j]->GetXaxis()->SetTitle(string);
    }
    
    
    TCanvas *c3 = new TCanvas("c3", "c3",200,10,700,700);
    c3->Divide(2,2);
    
    for (Int_t j=0; j<4;j++) {
        c3->cd(j+1);
        h1_trig_livetimes[j]->Draw();
        sprintf (string,"Live time for bit %d",j+1);
        if (j == 0) {
            h1_trig_livetimes[j]->SetTitle(filename);
        }
        else {
            h1_trig_livetimes[j]->SetTitle("");
        }
        h1_trig_livetimes[j]->GetXaxis()->SetTitleSize(0.05);
        h1_trig_livetimes[j]->GetYaxis()->SetTitleSize(0.05);
        h1_trig_livetimes[j]->GetXaxis()->SetTitle(string);
    }
    
    
    TCanvas *c4 = new TCanvas("c4", "c4",200,10,700,700);
    
    c4->Divide(1,2);
    c4->cd(1);
    
    h1epics_liveinst_VSevent->Divide(h1epics_entries1_VSevent);
    h1epics_liveinst_VSevent->SetTitle(filename);
    // h1epics_liveinst_VSevent->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1epics_liveinst_VSevent->GetYaxis()->SetRangeUser(ymin,ymax);
    h1epics_liveinst_VSevent->GetXaxis()->SetTitleSize(0.05);
    h1epics_liveinst_VSevent->GetYaxis()->SetTitleSize(0.05);
    h1epics_liveinst_VSevent->GetXaxis()->SetTitle("Trigger number");
    h1epics_liveinst_VSevent->SetLineColor(2);
    h1epics_liveinst_VSevent->Draw("");
    
    c4->cd(2);
    
    // h1epics_AD00_VSevent->Divide(h1epics_entries2_VSevent);
    h1epics_AD00_VSevent->SetTitle(filename);
    // h1epics_AD00_VSevent->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1epics_AD00_VSevent->GetYaxis()->SetRangeUser(ymin,ymax);
    h1epics_AD00_VSevent->GetXaxis()->SetTitleSize(0.05);
    h1epics_AD00_VSevent->GetYaxis()->SetTitleSize(0.05);
    h1epics_AD00_VSevent->GetXaxis()->SetTitle("Trigger number");
    h1epics_AD00_VSevent->SetLineColor(2);
    h1epics_AD00_VSevent->Draw("");
    
        
   sprintf (string,"plot_epics_dump.pdf(");
    c0->SaveAs(string);
    sprintf (string,"plot_epics_dump.pdf");
    c1->SaveAs(string);
    sprintf (string,"plot_epics_dump.pdf");
    c2->SaveAs(string);
    sprintf (string,"plot_epics_dump.pdf");
    c3->SaveAs(string);
    sprintf (string,"plot_epics_dump.pdf)");
    c4->SaveAs(string);
    
        
        
        }

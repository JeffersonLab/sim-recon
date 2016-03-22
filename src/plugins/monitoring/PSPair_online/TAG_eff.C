// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PSPair/PSC_PS_TAGH/PSTAGH_E
// hnamepath: /PSPair/PSC_PS_TAGM/PSTAGM_E
// hnamepath: /PSPair/PSC_PS/PS_E

{
    TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PSPair");
    if(dir) dir->cd();

    TH1F *E_tagh =(TH1F*)gDirectory->Get("PSC_PS_TAGH/PSTAGH_E");
    TH1F *E_tagm =(TH1F*)gDirectory->Get("PSC_PS_TAGM/PSTAGM_E");
    TH1F *E_ps =(TH1F*)gDirectory->Get("PSC_PS/PS_E");

    TH1F *eff_tagh = (TH1F*)E_tagh->Clone();
    TH1F *eff_tagm = (TH1F*)E_tagm->Clone();
    eff_tagh->Sumw2();
    eff_tagm->Sumw2();
    eff_tagh->Divide(E_tagh,E_ps);
    eff_tagm->Divide(E_tagm,E_ps);

    if(gPad == NULL){
        TCanvas *c1 = new TCanvas("c1","Tagging Efficiency",150,10,990,660);
        c1->cd(0);
        c1->Draw();
        c1->Update();
    }

    if(!gPad) return;
    TCanvas* c1 = gPad->GetCanvas();

    gStyle->SetOptStat("");
    TLegend *tleg = new TLegend(0.75,0.75,0.90,0.90);
    tleg->SetTextSize(0.045);
    tleg->AddEntry(eff_tagm,"TAGM","l");
    tleg->AddEntry(eff_tagh,"TAGH","l");
    eff_tagh->SetTitle("Tagging Efficiency");
    eff_tagh->SetTitleSize(0.045);
    eff_tagh->GetXaxis()->SetTitleSize(0.045);
    eff_tagh->GetYaxis()->SetTitleSize(0.045);
    eff_tagh->GetXaxis()->SetTitle("PS energy [GeV]");
    eff_tagh->GetYaxis()->SetTitle("N(TAGX,PSC,PS) / N(PSC,PS)");
    eff_tagh->GetXaxis()->SetRange(eff_tagh->FindFirstBinAbove(0.001),eff_tagh->FindLastBinAbove(0.001));
    eff_tagh->SetAxisRange(0.0,1.0,"Y");
    eff_tagh->SetLineColor(kGreen);
    eff_tagh->Draw();
    eff_tagm->SetLineColor(kRed);
    eff_tagm->Draw("same");
    tleg->Draw();

}

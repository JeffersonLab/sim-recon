// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PSPair/PSC_PS/PS_PSIDLeftVsIDRight
// hnamepath: /PSPair/PSC_PS_TAGH/PSTAGH_PSIDLeftVsIDRight
// hnamepath: /PSPair/PSC_PS_TAGM/PSTAGM_PSIDLeftVsIDRight
// hnamepath: /PSPair/PSC_PS/PS_ElVsEr
// hnamepath: /PSPair/PSC_PS_TAGH/PSTAGH_ElVsEr
// hnamepath: /PSPair/PSC_PS_TAGM/PSTAGM_ElVsEr

{
    TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PSPair");
    if(dir) dir->cd();

    TH2F *id_LR_D =(TH2F*)gDirectory->Get("PSC_PS/PS_PSIDLeftVsIDRight");
    TH2F *id_LR_N_tagh =(TH2F*)gDirectory->Get("PSC_PS_TAGH/PSTAGH_PSIDLeftVsIDRight");
    TH2F *id_LR_N_tagm =(TH2F*)gDirectory->Get("PSC_PS_TAGM/PSTAGM_PSIDLeftVsIDRight");
    //
    TH2F *E_LR_D =(TH2F*)gDirectory->Get("PSC_PS/PS_ElVsEr");
    TH2F *E_LR_N_tagh =(TH2F*)gDirectory->Get("PSC_PS_TAGH/PSTAGH_ElVsEr");
    TH2F *E_LR_N_tagm =(TH2F*)gDirectory->Get("PSC_PS_TAGM/PSTAGM_ElVsEr");

    TH2F *eff_tagh2d = (TH2F*)id_LR_N_tagh->Clone();
    TH2F *eff_tagm2d = (TH2F*)id_LR_N_tagm->Clone();
    eff_tagh2d->Sumw2();
    eff_tagm2d->Sumw2();
    eff_tagh2d->Divide(id_LR_N_tagh,id_LR_D);
    eff_tagm2d->Divide(id_LR_N_tagm,id_LR_D);
    //
    TH2F *eff_tagh_E = (TH2F*)E_LR_N_tagh->Clone();
    TH2F *eff_tagm_E = (TH2F*)E_LR_N_tagm->Clone();
    eff_tagh_E->Sumw2();
    eff_tagm_E->Sumw2();
    eff_tagh_E->Divide(E_LR_N_tagh,E_LR_D);
    eff_tagm_E->Divide(E_LR_N_tagm,E_LR_D);

    if(gPad == NULL){
        TCanvas *c1 = new TCanvas("c1","Tagger Efficiency",150,10,990,660);
        c1->cd(0);
        c1->Draw();
        c1->Update();
    }

    if(!gPad) return;
    TCanvas* c1 = gPad->GetCanvas();
    c1->Divide(2,2);

    double tsize = 0.0475;
    gStyle->SetOptStat("");
    c1->cd(4);
    eff_tagh2d->SetTitle("TAGH Tagger Efficiency: N(TAGH,PSC,PS) / N(PSC,PS)");
    eff_tagh2d->SetTitleSize(tsize,"xy");
    eff_tagh2d->Draw("colz");
    c1->cd(3);
    eff_tagm2d->SetTitle("TAGM Tagger Efficiency: N(TAGM,PSC,PS) / N(PSC,PS)");
    eff_tagm2d->SetTitleSize(tsize,"xy");
    eff_tagm2d->Draw("colz");
    c1->cd(2);
    eff_tagh_E->SetTitle("TAGH Tagger Efficiency: N(TAGH,PSC,PS) / N(PSC,PS)");
    eff_tagh_E->SetTitleSize(tsize,"xy");
    eff_tagh_E->Draw("colz");
    c1->cd(1);
    eff_tagm_E->SetTitle("TAGM Tagger Efficiency: N(TAGM,PSC,PS) / N(PSC,PS)");
    eff_tagm_E->SetTitleSize(tsize,"xy");
    eff_tagm_E->Draw("colz");

}

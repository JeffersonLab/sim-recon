// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /TAGH/Hit/Hit_tdcTimeVsSlotID
// hnamepath: /TAGH/Hit/Hit_HasTDCvsHasADC
// hnamepath: /TAGH/Hit/Hit_fadcTimeVsSlotID
// hnamepath: /TAGH/Hit/Hit_tdcadcTimeDiffVsSlotID
// hnamepath: /TAGH/Hit/Hit_IntegralVsSlotID
// hnamepath: /TAGH/Hit/Hit_TimeVsIntegral

{
    TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("TAGH/Hit");
    if(dir) dir->cd();

    TH2I* hT_tdc = (TH2I*)gDirectory->FindObjectAny("Hit_tdcTimeVsSlotID");
    TH2F* hStat = (TH2F*)gDirectory->FindObjectAny("Hit_HasTDCvsHasADC");
    TH2I* hT_adc = (TH2I*)gDirectory->FindObjectAny("Hit_fadcTimeVsSlotID");
    TH2I* hTdiff = (TH2I*)gDirectory->FindObjectAny("Hit_tdcadcTimeDiffVsSlotID");
    TH2I* hPI = (TH2I*)gDirectory->FindObjectAny("Hit_IntegralVsSlotID");
    TH2I* hTvsPI = (TH2I*)gDirectory->FindObjectAny("Hit_TimeVsIntegral");

    if(gPad == NULL){
        TCanvas *c1 = new TCanvas("c1","TAGH Hit Monitor II",150,10,990,660);
        c1->cd(0);
        c1->Draw();
        c1->Update();
    }

    if(!gPad) return;
    TCanvas* c1 = gPad->GetCanvas();
    c1->Divide(3,2);

    double tsize = 0.0475;
    gStyle->SetOptStat("emr");
    if(hT_tdc){
        c1->cd(1);
        hT_tdc->SetTitleSize(tsize,"xy");
        hT_tdc->GetYaxis()->SetRange(hT_tdc->FindFirstBinAbove(10.0,2),hT_tdc->FindLastBinAbove(10.0,2));
        hT_tdc->Draw("colz");
    }

    if(hStat){
        c1->cd(2);
        hStat->SetTitleSize(tsize,"xy");
        hStat->GetYaxis()->SetRange(hStat->FindFirstBinAbove(10.0,2),hStat->FindLastBinAbove(10.0,2));
        hStat->SetMarkerColor(kRed);
        hStat->SetMarkerSize(2.0);
        hStat->DrawNormalized("text");
    }

    if(hPI){
        c1->cd(3);
        hPI->SetTitleSize(tsize,"xy");
        hPI->GetYaxis()->SetRange(hPI->FindFirstBinAbove(10.0,2),hPI->FindLastBinAbove(10.0,2));
        hPI->Draw("colz");
    }

    if(hT_adc){
        c1->cd(4);
        hT_adc->SetTitleSize(tsize,"xy");
        hT_adc->GetYaxis()->SetRange(hT_adc->FindFirstBinAbove(10.0,2),hT_adc->FindLastBinAbove(10.0,2));
        hT_adc->Draw("colz");
    }

    if(hTdiff){
        c1->cd(5);
        hTdiff->SetTitleSize(tsize,"xy");
        hTdiff->GetYaxis()->SetRange(hTdiff->FindFirstBinAbove(10.0,2),hTdiff->FindLastBinAbove(10.0,2));
        hTdiff->Draw("colz");
    }

    if(hTvsPI){
        c1->cd(6);
        hTvsPI->SetTitleSize(tsize,"xy");
        hTvsPI->GetXaxis()->SetRange(hTvsPI->FindFirstBinAbove(10.0),hTvsPI->FindLastBinAbove(10.0));
        hTvsPI->GetYaxis()->SetRange(hTvsPI->FindFirstBinAbove(10.0,2),hTvsPI->FindLastBinAbove(10.0,2));
        hTvsPI->Draw("colz");
    }

}

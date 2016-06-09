// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /TPOL/Hit/Hit_NHits

{
    TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("TPOL");
    if(dir) dir->cd();

    TH1I* tpol_num_events = (TH1I*)gDirectory->FindObjectAny("tpol_num_events");
    TH1I* hHit_NHits = (TH1I*)gDirectory->FindObjectAny("hHit_NHits");

    if(gPad == NULL){
        TCanvas *c1 = new TCanvas("c1","TPOL Hit Monitor",150,10,990,660);
        c1->cd(0);
        c1->Draw();
        c1->Update();
    }
    if(!gPad) return;
    TCanvas* c1 = gPad->GetCanvas();
    c1->Divide(1,2);
    double tsize = 0.0475;
    gStyle->SetOptStat("emr");
    if(tpol_num_events){
        tpol_num_events->SetFillColor(kBlue);
        c1->cd(1);
        tpol_num_events->SetTitleSize(tsize,"xy");
        tpol_num_events->Draw();
    }
    if(hHit_NHits){
        hHit_NHits->SetFillColor(kBlue);
        c1->cd(2);
        hHit_NHits->SetTitleSize(tsize,"xy");
        hHit_NHits->GetXaxis()->SetRange(1,hHit_NHits->FindLastBinAbove(1.0));
        hHit_NHits->Draw();
    }

}

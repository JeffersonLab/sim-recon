// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /TAGH/Hit/Hit_NHits
// hnamepath: /TAGH/Hit/Hit_Occupancy
// hnamepath: /TAGH/Hit/Hit_Time
// hnamepath: /TAGH/Hit/Hit_Energy

{
    TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("TAGH/Hit");
    if(dir) dir->cd();

    TH1I* hNHits = (TH1I*)gDirectory->FindObjectAny("Hit_NHits");
    TH1I* hOccupancy = (TH1I*)gDirectory->FindObjectAny("Hit_Occupancy");
    TH1I* hTime = (TH1I*)gDirectory->FindObjectAny("Hit_Time");
    TH1I* hEnergy = (TH1I*)gDirectory->FindObjectAny("Hit_Energy");

    if(gPad == NULL){
        TCanvas *c1 = new TCanvas("c1","TAGH Hit Monitor I",150,10,990,660);
        c1->cd(0);
        c1->Draw();
        c1->Update();
    }
    
    if(!gPad) return;
    TCanvas* c1 = gPad->GetCanvas();
    c1->Divide(2,2);

    double tsize = 0.0475;
    gStyle->SetOptStat("emr");
    if(hOccupancy){
        hOccupancy->SetFillColor(kBlue);
        c1->cd(1);
        hOccupancy->SetTitleSize(tsize,"xy");
        hOccupancy->Draw();
    }

    if(hEnergy){
        hEnergy->SetFillColor(kBlue);
        c1->cd(2);
        hEnergy->SetTitleSize(tsize,"xy");
        hEnergy->Draw();
    }

    if(hTime){
        hTime->SetFillColor(kBlue);
        c1->cd(3);
        hTime->SetTitleSize(tsize,"xy");
        hTime->GetXaxis()->SetRange(hTime->FindFirstBinAbove(1.0),hTime->FindLastBinAbove(1.0));
        hTime->Draw();
    }

    if(hNHits){
        hNHits->SetFillColor(kBlue);
        c1->cd(4);
        hNHits->SetTitleSize(tsize,"xy");
        hNHits->GetXaxis()->SetRange(hNHits->FindFirstBinAbove(1.0),hNHits->FindLastBinAbove(1.0));
        hNHits->Draw();
    }

}

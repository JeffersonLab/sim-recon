// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PSPair/PSC_PS/PS_PSIDLeftVsPSCIDLeft
// hnamepath: /PSPair/PSC_PS/PS_PSIDRightVsPSCIDRight
// hnamepath: /PSPair/PSC_PS/PS_ElVsPSCIDLeft
// hnamepath: /PSPair/PSC_PS/PS_ErVsPSCIDRight

{
    TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PSPair");
    if(dir) dir->cd();

    TH2I* id_LL =(TH2I*)gDirectory->Get("PSC_PS/PS_PSIDLeftVsPSCIDLeft");
    TH2I* id_RR =(TH2I*)gDirectory->Get("PSC_PS/PS_PSIDRightVsPSCIDRight");
    TH2I* Eid_LL =(TH2I*)gDirectory->Get("PSC_PS/PS_ElVsPSCIDLeft");
    TH2I* Eid_RR =(TH2I*)gDirectory->Get("PSC_PS/PS_ErVsPSCIDRight");

    if(gPad == NULL){
        TCanvas *c1 = new TCanvas("c1","PS/PSC Geometrical Coincidence Monitor",150,10,990,660);
        c1->cd(0);
        c1->Draw();
        c1->Update();
    }

    if(!gPad) return;
    TCanvas* c1 = gPad->GetCanvas();
    c1->Divide(2,2);

    double tsize = 0.0475;
    gStyle->SetOptStat("emr");
    if (id_LL){
        c1->cd(3);
        id_LL->SetTitleSize(tsize,"xy");
        id_LL->Draw();
    }
    if (id_RR) {
        c1->cd(4);
        id_RR->SetTitleSize(tsize,"xy");
        id_RR->Draw();
    }
    if (Eid_LL){
        c1->cd(1);
        Eid_LL->SetTitleSize(tsize,"xy");
        Eid_LL->Draw();
    }
    if (Eid_RR) {
        c1->cd(2);
        Eid_RR->SetTitleSize(tsize,"xy");
        Eid_RR->Draw();
    }

}

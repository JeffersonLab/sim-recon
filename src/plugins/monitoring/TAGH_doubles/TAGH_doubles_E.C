// Plot fraction of TAGH hits that are double hits vs. Energy

void TAGH_doubles_E()
{
    TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("TAGH_doubles");
    if(dir) dir->cd();

    TH1F *doubles =(TH1F*)gDirectory->Get("BeforeMergingDoubles/BM2_Energy");
    TH1F *total =(TH1F*)gDirectory->Get("BeforeMergingDoubles/BM1_Energy");

    TH1F *f_doubles = (TH1F*)doubles->Clone();
    f_doubles->Sumw2();
    f_doubles->Divide(doubles, total);

    if(gPad == NULL) {
        TCanvas *c1 = new TCanvas("c1","TAGH double-hit fraction",150,10,990,660);
        c1->cd(0);
        c1->Draw();
        c1->Update();
    }

    if(!gPad) return;
    TCanvas* c1 = gPad->GetCanvas();

    gStyle->SetOptStat("");
    f_doubles->SetTitle("TAGH double-hit fraction vs. energy");
    f_doubles->SetTitleSize(0.045, "XY");
    f_doubles->GetXaxis()->SetTitle("TAGH energy [GeV]");
    f_doubles->GetYaxis()->SetTitle("double-hit fraction");
    f_doubles->SetAxisRange(0.0,1.0,"Y");
    f_doubles->Draw();
}

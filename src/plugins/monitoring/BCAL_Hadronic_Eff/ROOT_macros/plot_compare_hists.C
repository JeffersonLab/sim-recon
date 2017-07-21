void plot_compare_hists(void)
{
// read output of Read_bcal_hadronic_eff2 root files (in folder root) with efficiency histograms from data and MC and plot on same plot.
//

#include <TRandom.h>

gROOT->Reset();
//TTree *Bfield = (TTree *) gROOT->FindObject("Bfield");
gStyle->SetPalette(1,0);
gStyle->SetOptStat(kFALSE);
// gStyle->SetOptStat(11111111);
gStyle->SetPadRightMargin(0.15);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
//

   char string[256];
    
    Int_t const nlayers=4;
    
    Int_t coinc_cut = 3;
    // TString filerun = "011529";
    // TString filerun = "011405";
    // TString prefix = "MC08392Beam";
    TString filerun = "030450";
    TString prefix = "MC01Beam";
    TString titname = "Run "+filerun+" "+prefix;
    
    TFile *datahist[nlayers];
    TFile *MChist[nlayers];
    TString filename;
    
    
    TH1F *h1_eff_up[nlayers];
    TH1F *h1_eff_down[nlayers];
    TH1F *h1_eff_mom_up[nlayers];
    TH1F *h1_eff_mom_down[nlayers];
    TH1F *h1_eff_z_up[nlayers];
    TH1F *h1_eff_z_down[nlayers];
    
    TH1F *h1_MC_eff_up[nlayers];
    TH1F *h1_MC_eff_down[nlayers];
    TH1F *h1_MC_eff_mom_up[nlayers];
    TH1F *h1_MC_eff_mom_down[nlayers];
    TH1F *h1_MC_eff_z_up[nlayers];
    TH1F *h1_MC_eff_z_down[nlayers];
    
    
    for (Int_t j=0; j<nlayers;j++) {
    	filename = "R"+filerun+"_layer"+TString::Itoa(j+1,10)+"_cut"+TString::Itoa(coinc_cut,10);
    
        // open input histogram files
    
    	TString datahist_name = "root/"+filename+"_out.root";
    	datahist[j] = new TFile(datahist_name,"read");
    
    	filename = "R"+prefix+"_"+filerun+"_layer"+TString::Itoa(j+1,10)+"_cut"+TString::Itoa(coinc_cut,10);
    	TString MChist_name = "root/"+filename+"_out.root";
        MChist[j] = new TFile(MChist_name,"read");
        
        h1_eff_up[j] = (TH1F*)datahist[j]->Get("h1_eff_up");
        h1_eff_down[j] = (TH1F*)datahist[j]->Get("h1_eff_down");
        h1_eff_mom_up[j] = (TH1F*)datahist[j]->Get("h1_eff_mom_up");
        h1_eff_mom_down[j] = (TH1F*)datahist[j]->Get("h1_eff_mom_down");
        h1_eff_z_up[j] = (TH1F*)datahist[j]->Get("h1_eff_z_up");
        h1_eff_z_down[j] = (TH1F*)datahist[j]->Get("h1_eff_z_down");
        
        h1_MC_eff_up[j] = (TH1F*)MChist[j]->Get("h1_eff_up");
        h1_MC_eff_down[j] = (TH1F*)MChist[j]->Get("h1_eff_down");
        h1_MC_eff_mom_up[j] = (TH1F*)MChist[j]->Get("h1_eff_mom_up");
        h1_MC_eff_mom_down[j] = (TH1F*)MChist[j]->Get("h1_eff_mom_down");
        h1_MC_eff_z_up[j] = (TH1F*)MChist[j]->Get("h1_eff_z_up");
        h1_MC_eff_z_down[j] = (TH1F*)MChist[j]->Get("h1_eff_z_down");
    }
    
    
    TCanvas *c0 = new TCanvas("c0", "c0",200,10,1000,700);
    gPad->Divide(2,2);
    
    TLegend *leg[nlayers];
    for (Int_t j=0; j<nlayers;j++) {
        c0->cd(j+1);
    	h1_eff_up[j]->SetTitle(titname);
    	// h1_eff_up[j]->GetXaxis()->SetRangeUser(xmin,xmax);
    	h1_eff_up[j]->GetYaxis()->SetRangeUser(0.5,1);
    	h1_eff_up[j]->GetXaxis()->SetTitleSize(0.05);
    	h1_eff_up[j]->GetYaxis()->SetTitleSize(0.05);
    	h1_eff_up[j]->GetYaxis()->SetTitle("Up Efficiency");
    	h1_eff_up[j]->GetXaxis()->SetTitle("Up Channel Number");
    	h1_eff_up[j]->SetMarkerColor(2);
    	h1_eff_up[j]->SetMarkerStyle(20);
    	h1_eff_up[j]->SetMarkerSize(0.3);
    	h1_eff_up[j]->Draw("");
    
    	h1_MC_eff_up[j]->Draw("same");
    	h1_MC_eff_up[j]->SetMarkerColor(4);
    	h1_MC_eff_up[j]->SetMarkerStyle(21);
    	h1_MC_eff_up[j]->SetMarkerSize(0.3);
    
    	leg[j] = new TLegend(0.5,0.5,0.8,0.7);
        TString legend = "Layer "+TString::Itoa(j+1,10)+": Data";
        leg[j]->AddEntry(h1_eff_up[j],legend,"p");
        legend = "Layer "+TString::Itoa(j+1,10)+": MC";
        leg[j]->AddEntry(h1_MC_eff_up[j],legend,"p");
    	leg[j]->Draw();
    }
    
    
    TCanvas *c1 = new TCanvas("c1", "c1",200,10,1000,700);
    gPad->Divide(2,2);
    
    for (Int_t j=0; j<nlayers;j++) {
        c1->cd(j+1);
        h1_eff_down[j]->SetTitle(titname);
        // h1_eff_down[j]->GetXaxis()->SetRangeUser(xmin,xmax);
        h1_eff_down[j]->GetYaxis()->SetRangeUser(0.5,1);
        h1_eff_down[j]->GetXaxis()->SetTitleSize(0.05);
        h1_eff_down[j]->GetYaxis()->SetTitleSize(0.05);
        h1_eff_down[j]->GetYaxis()->SetTitle("Down Efficiency");
        h1_eff_down[j]->GetXaxis()->SetTitle("Down Channel Number");
        h1_eff_down[j]->SetMarkerColor(2);
        h1_eff_down[j]->SetMarkerStyle(20);
        h1_eff_down[j]->SetMarkerSize(0.3);
        h1_eff_down[j]->Draw("");
        
        h1_MC_eff_down[j]->Draw("same");
        h1_MC_eff_down[j]->SetMarkerColor(4);
        h1_MC_eff_down[j]->SetMarkerStyle(21);
        h1_MC_eff_down[j]->SetMarkerSize(0.3);
        
        leg[j] = new TLegend(0.5,0.5,0.8,0.7);
        TString legend = "Layer "+TString::Itoa(j+1,10)+": Data";
        leg[j]->AddEntry(h1_eff_down[j],legend,"p");
        legend = "Layer "+TString::Itoa(j+1,10)+": MC";
        leg[j]->AddEntry(h1_MC_eff_down[j],legend,"p");
        leg[j]->Draw();
    }
    
    TCanvas *c2 = new TCanvas("c2", "c2",200,10,1000,700);
    gPad->Divide(2,2);
    
    for (Int_t j=0; j<nlayers;j++) {
        c2->cd(j+1);
        h1_eff_z_up[j]->SetTitle(titname);
        // h1_eff_z_up[j]->GetXaxis()->SetRangeUser(xmin,xmax);
        h1_eff_z_up[j]->GetYaxis()->SetRangeUser(0.5,1);
        if (j+1 == 4) h1_eff_z_up[j]->GetYaxis()->SetRangeUser(0,1);
        h1_eff_z_up[j]->GetXaxis()->SetTitleSize(0.05);
        h1_eff_z_up[j]->GetYaxis()->SetTitleSize(0.05);
        h1_eff_z_up[j]->GetYaxis()->SetTitle("Up Efficiency");
        h1_eff_z_up[j]->GetXaxis()->SetTitle("z (cm)");
        h1_eff_z_up[j]->SetMarkerColor(2);
        h1_eff_z_up[j]->SetMarkerStyle(20);
        h1_eff_z_up[j]->SetMarkerSize(0.3);
        h1_eff_z_up[j]->Draw("");
        
        h1_MC_eff_z_up[j]->Draw("same");
        h1_MC_eff_z_up[j]->SetMarkerColor(4);
        h1_MC_eff_z_up[j]->SetMarkerStyle(21);
        h1_MC_eff_z_up[j]->SetMarkerSize(0.3);
        
        if (j+1 == 4) {
            leg[j] = new TLegend(0.3,0.2,0.6,0.4);
        }
        else {
            leg[j] = new TLegend(0.3,0.5,0.6,0.7);
        }
        TString legend = "Layer "+TString::Itoa(j+1,10)+": Data";
        leg[j]->AddEntry(h1_eff_z_up[j],legend,"p");
        legend = "Layer "+TString::Itoa(j+1,10)+": MC";
        leg[j]->AddEntry(h1_MC_eff_z_up[j],legend,"p");
        leg[j]->Draw();
    }
    
    
    TCanvas *c3 = new TCanvas("c3", "c3",200,10,1000,700);
    gPad->Divide(2,2);
    
    for (Int_t j=0; j<nlayers;j++) {
        c3->cd(j+1);
        h1_eff_z_down[j]->SetTitle(titname);
        // h1_eff_z_down[j]->GetXaxis()->SetRangeUser(xmin,xmax);
        h1_eff_z_down[j]->GetYaxis()->SetRangeUser(0.5,1);
        if (j+1 == 4) h1_eff_z_down[j]->GetYaxis()->SetRangeUser(0,1);
        h1_eff_z_down[j]->GetXaxis()->SetTitleSize(0.05);
        h1_eff_z_down[j]->GetYaxis()->SetTitleSize(0.05);
        h1_eff_z_down[j]->GetYaxis()->SetTitle("Down Efficiency");
        h1_eff_z_down[j]->GetXaxis()->SetTitle("z (cm)");
        h1_eff_z_down[j]->SetMarkerColor(2);
        h1_eff_z_down[j]->SetMarkerStyle(20);
        h1_eff_z_down[j]->SetMarkerSize(0.3);
        h1_eff_z_down[j]->Draw("");
        
        h1_MC_eff_z_down[j]->Draw("same");
        h1_MC_eff_z_down[j]->SetMarkerColor(4);
        h1_MC_eff_z_down[j]->SetMarkerStyle(21);
        h1_MC_eff_z_down[j]->SetMarkerSize(0.3);
        
        
        if (j+1 == 4) {
            leg[j] = new TLegend(0.3,0.2,0.6,0.4);
        }
        else {
            leg[j] = new TLegend(0.3,0.5,0.6,0.7);
        }
        TString legend = "Layer "+TString::Itoa(j+1,10)+": Data";
        leg[j]->AddEntry(h1_eff_z_down[j],legend,"p");
        legend = "Layer "+TString::Itoa(j+1,10)+": MC";
        leg[j]->AddEntry(h1_MC_eff_z_down[j],legend,"p");
        leg[j]->Draw();
    }
    
    TCanvas *c4 = new TCanvas("c4", "c4",200,10,1000,700);
    gPad->Divide(2,2);
    
    for (Int_t j=0; j<nlayers;j++) {
        c4->cd(j+1);
        h1_eff_mom_up[j]->SetTitle(titname);
        // h1_eff_mom_up[j]->GetXaxis()->SetRangeUser(xmin,xmax);
        h1_eff_mom_up[j]->GetYaxis()->SetRangeUser(0.5,1);
        if (j+1 == 4) h1_eff_mom_up[j]->GetYaxis()->SetRangeUser(0,1);
        h1_eff_mom_up[j]->GetXaxis()->SetTitleSize(0.05);
        h1_eff_mom_up[j]->GetYaxis()->SetTitleSize(0.05);
        h1_eff_mom_up[j]->GetYaxis()->SetTitle("Up Efficiency");
        h1_eff_mom_up[j]->GetXaxis()->SetTitle("Momentum (GeV)");
        h1_eff_mom_up[j]->SetMarkerColor(2);
        h1_eff_mom_up[j]->SetMarkerStyle(20);
        h1_eff_mom_up[j]->SetMarkerSize(0.3);
        h1_eff_mom_up[j]->Draw("");
        
        h1_MC_eff_mom_up[j]->Draw("same");
        h1_MC_eff_mom_up[j]->SetMarkerColor(4);
        h1_MC_eff_mom_up[j]->SetMarkerStyle(21);
        h1_MC_eff_mom_up[j]->SetMarkerSize(0.3);
        
        
        if (j+1 == 4) {
            leg[j] = new TLegend(0.3,0.2,0.6,0.4);
        }
        else {
            leg[j] = new TLegend(0.3,0.5,0.6,0.7);
        }
        TString legend = "Layer "+TString::Itoa(j+1,10)+": Data";
        leg[j]->AddEntry(h1_eff_mom_up[j],legend,"p");
        legend = "Layer "+TString::Itoa(j+1,10)+": MC";
        leg[j]->AddEntry(h1_MC_eff_mom_up[j],legend,"p");
        leg[j]->Draw();
    }
    
    
    TCanvas *c5 = new TCanvas("c5", "c5",200,10,1000,700);
    gPad->Divide(2,2);
    
    for (Int_t j=0; j<nlayers;j++) {
        c5->cd(j+1);
        h1_eff_mom_down[j]->SetTitle(titname);
        // h1_eff_mom_down[j]->GetXaxis()->SetRangeUser(xmin,xmax);
        h1_eff_mom_down[j]->GetYaxis()->SetRangeUser(0.5,1);
        if (j+1 == 4) h1_eff_mom_down[j]->GetYaxis()->SetRangeUser(0,1);
        h1_eff_mom_down[j]->GetXaxis()->SetTitleSize(0.05);
        h1_eff_mom_down[j]->GetYaxis()->SetTitleSize(0.05);
        h1_eff_mom_down[j]->GetYaxis()->SetTitle("Down Efficiency");
        h1_eff_mom_down[j]->GetXaxis()->SetTitle("Momentum (GeV)");
        h1_eff_mom_down[j]->SetMarkerColor(2);
        h1_eff_mom_down[j]->SetMarkerStyle(20);
        h1_eff_mom_down[j]->SetMarkerSize(0.3);
        h1_eff_mom_down[j]->Draw("");
        
        h1_MC_eff_mom_down[j]->Draw("same");
        h1_MC_eff_mom_down[j]->SetMarkerColor(4);
        h1_MC_eff_mom_down[j]->SetMarkerStyle(21);
        h1_MC_eff_mom_down[j]->SetMarkerSize(0.3);
        
        
        if (j+1 == 4) {
            leg[j] = new TLegend(0.3,0.2,0.6,0.4);
        }
        else {
            leg[j] = new TLegend(0.3,0.5,0.6,0.7);
        }
        TString legend = "Layer "+TString::Itoa(j+1,10)+": Data";
        leg[j]->AddEntry(h1_eff_mom_down[j],legend,"p");
        legend = "Layer "+TString::Itoa(j+1,10)+": MC";
        leg[j]->AddEntry(h1_MC_eff_z_down[j],legend,"p");
        leg[j]->Draw();
    }
    
    TCanvas *c6 = new TCanvas("c6", "c6",200,10,1000,500);
    gPad->Divide(2,1);
    
    c6->cd(1);
    
    h1_eff_up[1]->SetTitle("");
    // h1_eff_up[1]->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_eff_up[1]->GetYaxis()->SetRangeUser(0.8,1);
    h1_eff_up[1]->GetXaxis()->SetTitleSize(0.05);
    h1_eff_up[1]->GetYaxis()->SetTitleSize(0.05);
    h1_eff_up[1]->GetYaxis()->SetTitleOffset(1.5);
    h1_eff_up[1]->GetYaxis()->SetTitle("Layer Efficiency");
    h1_eff_up[1]->GetXaxis()->SetTitle("Cell ID");
    h1_eff_up[1]->SetMarkerColor(2);
    h1_eff_up[1]->SetMarkerStyle(20);
    h1_eff_up[1]->SetMarkerSize(0.3);
    h1_eff_up[1]->Draw("");
    
    h1_MC_eff_up[1]->Draw("same");
    h1_MC_eff_up[1]->SetMarkerColor(4);
    h1_MC_eff_up[1]->SetMarkerStyle(20);
    h1_MC_eff_up[1]->SetMarkerSize(0.3);
    
    
    // dummies to enlarge legend
    TH1F *dummy1 = new TH1F("dummy1","dummy1",100,0,100);
    dummy1->SetMarkerColor(2);
    dummy1->SetMarkerStyle(20);
    dummy1->SetMarkerSize(0.7);
    
    TH1F *dummy2 = new TH1F("dummy2","dummy2",100,0,100);
    dummy2->SetMarkerColor(4);
    dummy2->SetMarkerStyle(20);
    dummy2->SetMarkerSize(0.7);
    
    TString legend;
    legend = "Layer 2: Data";
    TLegend *leg1 = new TLegend(0.3,0.2,0.6,0.4);
    leg1->AddEntry(dummy1,legend,"p");
    legend = "Layer 2: MC";
    leg1->AddEntry(dummy2,legend,"p");
    leg1->Draw();

    c6->cd(2);
    
    h1_eff_z_up[1]->SetTitle("");
    // h1_eff_z_up[1]->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_eff_z_up[1]->GetYaxis()->SetRangeUser(0.8,1);
    h1_eff_z_up[1]->GetXaxis()->SetTitleSize(0.05);
    h1_eff_z_up[1]->GetYaxis()->SetTitleSize(0.05);
    h1_eff_z_up[1]->GetYaxis()->SetTitleOffset(1.5);
    h1_eff_z_up[1]->GetYaxis()->SetTitle("Layer Efficiency");
    h1_eff_z_up[1]->GetXaxis()->SetTitle("Position (cm)");
    h1_eff_z_up[1]->SetMarkerColor(2);
    h1_eff_z_up[1]->SetMarkerStyle(20);
    h1_eff_z_up[1]->SetMarkerSize(0.3);
    h1_eff_z_up[1]->Draw("");
    
    h1_MC_eff_z_up[1]->Draw("same");
    h1_MC_eff_z_up[1]->SetMarkerColor(4);
    h1_MC_eff_z_up[1]->SetMarkerStyle(20);
    h1_MC_eff_z_up[1]->SetMarkerSize(0.3);
    
    legend = "Layer 2: Data";
    TLegend *leg2 = new TLegend(0.3,0.2,0.6,0.4);
    leg2->AddEntry(dummy1,legend,"p");
    legend = "Layer 2: MC";
    leg2->AddEntry(dummy2,legend,"p");
    leg2->Draw();
    
    
    c0->SaveAs("plot_compare_hists_"+prefix+".pdf(");
    c1->SaveAs("plot_compare_hists_"+prefix+".pdf");
    c2->SaveAs("plot_compare_hists_"+prefix+".pdf");
    c3->SaveAs("plot_compare_hists_"+prefix+".pdf");
    c4->SaveAs("plot_compare_hists_"+prefix+".pdf");
    c5->SaveAs("plot_compare_hists_"+prefix+".pdf");
    c6->SaveAs("plot_compare_hists_"+prefix+".pdf)");
    
    c6->SaveAs("NIM_hadronic_eff"+prefix+".C");
    c6->SaveAs("NIM_hadronic_eff"+prefix+".pdf");
    
    // datahist[j]->Close();
    // MChist[j]->Close();

}


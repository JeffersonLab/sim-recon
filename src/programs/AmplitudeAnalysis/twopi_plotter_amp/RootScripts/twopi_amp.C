void twopi_amp(void)
{
// File: twopi_amp.C
    // Output histograms and fits generated from amp fitting of parameters.
//

  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  char string[256];
    map<Int_t, TString> sdme;
    sdme[0]="RE g1VM1\t";
    sdme[1]="RE g1VM0\t";
    sdme[2]="IM g1VM0\t\t";
    sdme[3]="RE g1VM-1\t";
    sdme[4]="IM g1VM-1\t";
    sdme[5]="RE g-1VM1\t";
    sdme[6]="IM g-1VM1\t";
    sdme[7]="RE g-1VM0\t";
    sdme[8]="IM g-1VM0\t";
    sdme[9]="RE g-1VM-1\t";
    sdme[10]="IM g-1VM-1\t";
    
    const Int_t nparms=11;
    Double_t parms[nparms];
    Double_t parms_err[nparms];
    
    bool genfile(false);
    
    // TString filename = "twopi_amp";
    TString filename = "twopi_amp_fitPars";
    
    TString infile = filename+".fit2";   // file with parameters
    TFile *f = new TFile(filename+".root","read");
    
    cout << "Opening parameters file: " << infile.Data() << endl;
    cout << "Opening root file: " << (filename+".root").Data() << endl;
    
    
    TH1F *M2pigen = (TH1F*)f->Get("M2pigen");
    TH1F *M2piacc = (TH1F*)f->Get("M2piacc");
    TH1F *M2pidat = (TH1F*)f->Get("M2pidat");
    
    TH1F *cosThetagen = (TH1F*)f->Get("cosThetagen");
    TH1F *cosThetaacc = (TH1F*)f->Get("cosThetaacc");
    TH1F *cosThetadat = (TH1F*)f->Get("cosThetadat");
    
    TH1F *psigen = (TH1F*)f->Get("psigen");
    TH1F *psiacc = (TH1F*)f->Get("psiacc");
    TH1F *psidat = (TH1F*)f->Get("psidat");
    
    TH1F *Phigen = (TH1F*)f->Get("Phigen");
    TH1F *Phiacc = (TH1F*)f->Get("Phiacc");
    TH1F *Phidat = (TH1F*)f->Get("Phidat");
    
    TH1F *phigen = (TH1F*)f->Get("phigen");
    TH1F *phiacc = (TH1F*)f->Get("phiacc");
    TH1F *phidat = (TH1F*)f->Get("phidat");
    
    TH1F *tgen = (TH1F*)f->Get("tgen");
    TH1F *tacc = (TH1F*)f->Get("tacc");
    TH1F *tdat = (TH1F*)f->Get("tdat");
    
    
   TCanvas *c0 = new TCanvas("c0", "c0",200,10,1000,700);

   c0->Divide(3,2);
    c0->cd(1);
    // gPad->SetLogy();
    Double_t xmin = 0;
    Double_t xmax = 2;
    Double_t ymin = 100;
    Double_t ymax = 10000;
        
    M2pigen->SetTitle(filename);
    // M2pigen->GetXaxis()->SetRangeUser(xmin,xmax);
    // M2pigen->GetYaxis()->SetRangeUser(ymin,ymax);
    M2pigen->GetXaxis()->SetTitleSize(0.05);
    M2pigen->GetYaxis()->SetTitleSize(0.05);
    M2pigen->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-})");
    M2pigen->SetMarkerColor(4);
    M2pigen->Draw("p");
    // M2piacc->Draw("samep");
    M2pidat->SetMarkerColor(2);
    M2pidat->SetLineColor(2);
    M2pidat->SetMarkerStyle(20);
    M2pidat->SetMarkerSize(0.1);
    M2pidat->Draw("samep");
    
    TLegend *leg = new TLegend(0.6,0.3,0.8,0.5);
    leg->AddEntry(M2pigen,"Gen","lp");
    leg->AddEntry(M2piacc,"Acc","lp");
    leg->AddEntry(M2pidat,"Data","lp");
    leg->Draw();
    
    c0->cd(2);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 5000;
    
    cosThetagen->SetTitle(filename);
    // cosThetagen->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) cosThetagen->GetYaxis()->SetRangeUser(ymin,ymax);
    cosThetagen->GetXaxis()->SetTitleSize(0.05);
    cosThetagen->GetYaxis()->SetTitleSize(0.05);
    cosThetagen->GetXaxis()->SetTitle("cos(#theta)");
    cosThetagen->SetLineColor(4);
    cosThetagen->Draw("p");
    // cosThetaacc->Draw("samep");
    cosThetadat->SetMarkerColor(2);
    cosThetadat->SetLineColor(2);
    cosThetadat->SetMarkerStyle(20);
    cosThetadat->SetMarkerSize(0.1);
    cosThetadat->Draw("samep");
    
    c0->cd(3);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    TF1 *cos2phi = new TF1("cos2phi","[0]*(1+[1]*cos(2*x))",-3.14159,3.14159);
    
    psigen->SetTitle(filename);
    // psigen->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) psigen->GetYaxis()->SetRangeUser(ymin,ymax);
    psigen->GetXaxis()->SetTitleSize(0.05);
    psigen->GetYaxis()->SetTitleSize(0.05);
    psigen->GetXaxis()->SetTitle("#psi");
    psigen->SetMarkerColor(4);
    psigen->Fit(cos2phi);
    psigen->Draw("p");
    // psiacc->Draw("samep");
    psidat->SetMarkerColor(2);
    psidat->SetLineColor(2);
    psidat->SetMarkerStyle(20);
    psidat->SetMarkerSize(0.1);
    psidat->Draw("samep");
    
    c0->cd(4);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    Phigen->SetTitle(filename);
    // Phigen->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) Phigen->GetYaxis()->SetRangeUser(ymin,ymax);
    Phigen->GetXaxis()->SetTitleSize(0.05);
    Phigen->GetYaxis()->SetTitleSize(0.05);
    Phigen->GetXaxis()->SetTitle("#Phi");
    Phigen->SetMarkerColor(4);
    Phigen->Fit(cos2phi);
    Phigen->Draw("p");
    // Phiacc->Draw("samep");
    Phidat->SetMarkerColor(2);
    Phidat->SetLineColor(2);
    Phidat->SetMarkerStyle(20);
    Phidat->SetMarkerSize(0.1);
    Phidat->Draw("samep");
    
    c0->cd(5);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    phigen->SetTitle(filename);
    // phigen->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) phigen->GetYaxis()->SetRangeUser(ymin,ymax);
    phigen->GetXaxis()->SetTitleSize(0.05);
    phigen->GetYaxis()->SetTitleSize(0.05);
    phigen->GetXaxis()->SetTitle("#phi");
    phigen->SetMarkerColor(4);
    phigen->Fit(cos2phi);
    phigen->Draw("p");
    // phiacc->Draw("samep");
    phidat->SetMarkerColor(2);
    phidat->SetLineColor(2);
    phidat->SetMarkerStyle(20);
    phidat->SetMarkerSize(0.1);
    phidat->Draw("samep");
    
    c0->cd(6);
    gPad->SetLogy();
    xmin = 0;
    xmax = 3;
    
    tgen->SetTitle(filename);
    tgen->GetXaxis()->SetRangeUser(xmin,xmax);
    // tgen->GetYaxis()->SetRangeUser(ymin,ymax);
    tgen->GetXaxis()->SetTitleSize(0.05);
    tgen->GetYaxis()->SetTitleSize(0.05);
    tgen->GetXaxis()->SetTitle("-t");
    tgen->SetMarkerColor(4);
    tgen->Fit("expo","","",0.2,1.3);
    tgen->Draw("p");
    // tacc->Draw("samep");
    tdat->SetMarkerColor(2);
    tdat->SetLineColor(2);
    tdat->SetMarkerStyle(20);
    tdat->SetMarkerSize(0.1);
    tdat->Draw("samep");
    
    
    TCanvas *c2 = new TCanvas("c2", "c2",200,10,1000,700);
    
    c2->Divide(3,2);
    c2->cd(1);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 1.2;
    
    TH1F *M2piAcceptance = (TH1F*)M2piacc->Clone("M2piAcceptance");
    M2piAcceptance->SetTitle("Acceptance");
    // M2piAcceptance->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) M2piAcceptance->GetYaxis()->SetRangeUser(ymin,ymax);
    M2piAcceptance->Divide(M2pigen);
    M2piAcceptance->GetXaxis()->SetTitleSize(0.05);
    M2piAcceptance->GetYaxis()->SetTitleSize(0.05);
    M2piAcceptance->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-})");
    M2piAcceptance->SetMarkerColor(4);
    M2piAcceptance->Draw("p");
    
    c2->cd(2);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 1.2;
    
    TH1F *cosThetaAcceptance = (TH1F*)cosThetaacc->Clone("cosThetaAcceptance");
    cosThetaAcceptance->SetTitle("Acceptance");
    // cosThetaAcceptance->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) cosThetaAcceptance->GetYaxis()->SetRangeUser(ymin,ymax);
    cosThetaAcceptance->Divide(cosThetagen);
    cosThetaAcceptance->GetXaxis()->SetTitleSize(0.05);
    cosThetaAcceptance->GetYaxis()->SetTitleSize(0.05);
    cosThetaAcceptance->GetXaxis()->SetTitle("cos(#theta)");
    cosThetaAcceptance->SetMarkerColor(4);
    cosThetaAcceptance->Draw("p");
    
    c2->cd(3);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 1.2;
    
    TH1F *psiAcceptance = (TH1F*)psiacc->Clone("psiAcceptance");
    psiAcceptance->SetTitle("Acceptance");
    // psiAcceptance->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) psiAcceptance->GetYaxis()->SetRangeUser(ymin,ymax);
    psiAcceptance->Divide(psigen);
    psiAcceptance->GetXaxis()->SetTitleSize(0.05);
    psiAcceptance->GetYaxis()->SetTitleSize(0.05);
    psiAcceptance->GetXaxis()->SetTitle("#psi");
    psiAcceptance->SetMarkerColor(4);
    psiAcceptance->Draw("p");
    
    c2->cd(4);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 1.2;
    
    TH1F *PhiAcceptance = (TH1F*)Phiacc->Clone("PhiAcceptance");
    PhiAcceptance->SetTitle("Acceptance");
    // PhiAcceptance->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) PhiAcceptance->GetYaxis()->SetRangeUser(ymin,ymax);
    PhiAcceptance->Divide(Phigen);
    PhiAcceptance->GetXaxis()->SetTitleSize(0.05);
    PhiAcceptance->GetYaxis()->SetTitleSize(0.05);
    PhiAcceptance->GetXaxis()->SetTitle("#Phi");
    PhiAcceptance->SetMarkerColor(4);
    PhiAcceptance->Draw("p");
    
    c2->cd(5);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 1.2;
    
    TH1F *phiAcceptance = (TH1F*)phiacc->Clone("phiAcceptance");
    phiAcceptance->SetTitle("Acceptance");
    // phiAcceptance->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) phiAcceptance->GetYaxis()->SetRangeUser(ymin,ymax);
    phiAcceptance->Divide(phigen);
    phiAcceptance->GetXaxis()->SetTitleSize(0.05);
    phiAcceptance->GetYaxis()->SetTitleSize(0.05);
    phiAcceptance->GetXaxis()->SetTitle("#phi");
    phiAcceptance->SetMarkerColor(4);
    phiAcceptance->Draw("p");
    
    c2->cd(6);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 1.2;
    xmin = 0;
    xmax = 3;
    
    TH1F *tAcceptance = (TH1F*)tacc->Clone("tAcceptance");
    tAcceptance->SetTitle("Acceptance");
    tAcceptance->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) tAcceptance->GetYaxis()->SetRangeUser(ymin,ymax);
    tAcceptance->Divide(tgen);
    tAcceptance->GetXaxis()->SetTitleSize(0.05);
    tAcceptance->GetYaxis()->SetTitleSize(0.05);
    tAcceptance->GetXaxis()->SetTitle("-t");
    tAcceptance->SetMarkerColor(4);
    tAcceptance->Draw("p");
    
    TCanvas *c1 = new TCanvas("c1", "c1",200,10,1000,700);
    
    c1->Divide(3,2);
    c1->cd(1);
    // gPad->SetLogy();
    ymin = 100;
    ymax = 10000;
    
    M2piacc->SetTitle(filename);
    // M2piacc->GetXaxis()->SetRangeUser(xmin,xmax);
    // M2piacc->GetYaxis()->SetRangeUser(ymin,ymax);
    M2piacc->GetXaxis()->SetTitleSize(0.05);
    M2piacc->GetYaxis()->SetTitleSize(0.05);
    M2piacc->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-})");
    M2piacc->SetMarkerColor(1);
    M2piacc->SetLineColor(1);
    M2piacc->Draw("p");
    // M2piacc->Draw("samep");
    M2pidat->SetMarkerColor(2);
    M2pidat->SetLineColor(2);
    M2pidat->SetMarkerStyle(20);
    M2pidat->SetMarkerSize(0.1);
    M2pidat->Draw("samep");
    
    TLegend *leg1 = new TLegend(0.6,0.3,0.8,0.5);
    leg1->AddEntry(M2pigen,"Gen","lp");
    leg1->AddEntry(M2piacc,"Acc","lp");
    leg1->AddEntry(M2pidat,"Data","lp");
    leg1->Draw();
    
    c1->cd(2);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    cosThetaacc->SetTitle(filename);
    // cosThetaacc->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) cosThetaacc->GetYaxis()->SetRangeUser(ymin,ymax);
    cosThetaacc->GetXaxis()->SetTitleSize(0.05);
    cosThetaacc->GetYaxis()->SetTitleSize(0.05);
    cosThetaacc->GetXaxis()->SetTitle("cos(#theta)");
    cosThetaacc->SetLineColor(1);
    cosThetaacc->SetMarkerColor(1);
    cosThetaacc->Draw("p");
    // cosThetaacc->Draw("samep");
    cosThetadat->SetMarkerColor(2);
    cosThetadat->SetLineColor(2);
    cosThetadat->SetMarkerStyle(20);
    cosThetadat->SetMarkerSize(0.1);
    cosThetadat->Draw("samep");
    
    c1->cd(3);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    psiacc->SetTitle(filename);
    // psiacc->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) psiacc->GetYaxis()->SetRangeUser(ymin,ymax);
    psiacc->GetXaxis()->SetTitleSize(0.05);
    psiacc->GetYaxis()->SetTitleSize(0.05);
    psiacc->GetXaxis()->SetTitle("#psi");
    psiacc->SetMarkerColor(1);
    psiacc->SetLineColor(1);
    psiacc->Fit(cos2phi);
    psiacc->Draw("p");
    // psiacc->Draw("samep");
    psidat->SetMarkerColor(2);
    psidat->SetLineColor(2);
    psidat->SetMarkerStyle(20);
    psidat->SetMarkerSize(0.1);
    psidat->Draw("samep");
    
    c1->cd(4);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    Phiacc->SetTitle(filename);
    // Phiacc->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) Phiacc->GetYaxis()->SetRangeUser(ymin,ymax);
    Phiacc->GetXaxis()->SetTitleSize(0.05);
    Phiacc->GetYaxis()->SetTitleSize(0.05);
    Phiacc->GetXaxis()->SetTitle("#Phi");
    Phiacc->SetMarkerColor(1);
    Phiacc->SetLineColor(1);
    Phiacc->Fit(cos2phi);
    Phiacc->Draw("p");
    // Phiacc->Draw("samep");
    Phidat->SetMarkerColor(2);
    Phidat->SetLineColor(2);
    Phidat->SetMarkerStyle(20);
    Phidat->SetMarkerSize(0.1);
    Phidat->Draw("samep");
    
    
    c1->cd(5);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    phiacc->SetTitle(filename);
    // phiacc->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) phiacc->GetYaxis()->SetRangeUser(ymin,ymax);
    phiacc->GetXaxis()->SetTitleSize(0.05);
    phiacc->GetYaxis()->SetTitleSize(0.05);
    phiacc->GetXaxis()->SetTitle("#phi");
    phiacc->SetMarkerColor(1);
    phiacc->SetLineColor(1);
    phiacc->Fit(cos2phi);
    phiacc->Draw("p");
    // phiacc->Draw("samep");
    phidat->SetMarkerColor(2);
    phidat->SetLineColor(2);
    phidat->SetMarkerStyle(20);
    phidat->SetMarkerSize(0.1);
    phidat->Draw("samep");
    
    c1->cd(6);
    gPad->SetLogy();
    xmin = 0;
    xmax = 3;
    
    tacc->SetTitle(filename);
    tacc->GetXaxis()->SetRangeUser(xmin,xmax);
    // tacc->GetYaxis()->SetRangeUser(ymin,ymax);
    tacc->GetXaxis()->SetTitleSize(0.05);
    tacc->GetYaxis()->SetTitleSize(0.05);
    tacc->GetXaxis()->SetTitle("-t");
    tacc->SetMarkerColor(1);
    tacc->SetLineColor(1);
    tacc->Draw("p");
    // tacc->Draw("samep");
    tdat->SetMarkerColor(2);
    tdat->SetLineColor(2);
    tdat->SetMarkerStyle(20);
    tdat->SetMarkerSize(0.1);
    tdat->Draw("samep");
    
    TCanvas *c3 = new TCanvas("c3", "c3",200,10,700,700);
    
    c3->Divide(2,2);
    c3->cd(1);
    // gPad->SetLogy();
    ymin = 100;
    ymax = 10000;
    
    M2piacc->SetTitle(filename);
    // M2piacc->GetXaxis()->SetRangeUser(xmin,xmax);
    // M2piacc->GetYaxis()->SetRangeUser(ymin,ymax);
    M2piacc->GetXaxis()->SetTitleSize(0.05);
    M2piacc->GetYaxis()->SetTitleSize(0.05);
    M2piacc->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-})");
    M2piacc->SetMarkerColor(1);
    M2piacc->Draw("p");
    // M2piacc->Draw("samep");
    M2pidat->SetMarkerColor(2);
    M2pidat->SetLineColor(2);
    M2pidat->SetMarkerStyle(20);
    M2pidat->SetMarkerSize(0.1);
    M2pidat->Draw("samep");
    
    c3->cd(2);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    cosThetaacc->SetTitle(filename);
    // cosThetaacc->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) cosThetaacc->GetYaxis()->SetRangeUser(ymin,ymax);
    cosThetaacc->GetXaxis()->SetTitleSize(0.05);
    cosThetaacc->GetYaxis()->SetTitleSize(0.05);
    cosThetaacc->GetXaxis()->SetTitle("cos(#theta)");
    cosThetaacc->SetLineColor(1);
    cosThetaacc->Draw("p");
    // cosThetaacc->Draw("samep");
    cosThetadat->SetMarkerColor(2);
    cosThetadat->SetLineColor(2);
    cosThetadat->SetMarkerStyle(20);
    cosThetadat->SetMarkerSize(0.1);
    cosThetadat->Draw("samep");
    
    c3->cd(3);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    // TF1 *cos2phi = new TF1("cos2phi","[0]*(1+[1]*cos(2*x))",-3.14159,3.14159);
    
    psiacc->SetTitle(filename);
    // psiacc->GetXaxis()->SetRangeUser(xmin,xmax);
    if (!genfile) psiacc->GetYaxis()->SetRangeUser(ymin,ymax);
    psiacc->GetXaxis()->SetTitleSize(0.05);
    psiacc->GetYaxis()->SetTitleSize(0.05);
    psiacc->GetXaxis()->SetTitle("#psi");
    psiacc->SetMarkerColor(1);
    psiacc->Fit(cos2phi);
    psiacc->Draw("p");
    // psiacc->Draw("samep");
    psidat->SetMarkerColor(2);
    psidat->SetLineColor(2);
    psidat->SetMarkerStyle(20);
    psidat->SetMarkerSize(0.1);
    psidat->Draw("samep");
    
    
    c3->cd(4);
    
    // now read and print fitted values
    
    ifstream parameters;
    parameters.open (infile.Data());
    if (!parameters) {
        cout << "ERROR: Failed to open data file= " << infile.Data() << endl;
        return;
    }
    
    TString line;
    while (line.ReadLine(parameters)){
        
        TObjArray *tokens = line.Tokenize("\t");
        Int_t ntokens = tokens->GetEntries();
        
        cout << " ntokens=" << ntokens << " line=" << line.Data() << endl;
        Int_t jmax = ntokens/2 > nparms? nparms: ntokens/2;
        for (Int_t j=0; j<jmax; j++){
        	parms[j] = (((TObjString*)tokens->At(2*j))->GetString()).Atof();
        	parms_err[j] = (((TObjString*)tokens->At(2*j+1))->GetString()).Atof();
        }
        
    }   // end loop over lines
    
        sprintf (string,"AmpTool Fit\n");
        printf("string=%s",string);
        TLatex *t1 = new TLatex(0.2,0.95,string);
        // t1->SetNDC();
        t1->SetTextSize(0.04);
        t1->Draw();
        
        for (Int_t j=0; j<nparms; j++) {
            cout << sdme[j] << "=" << parms[j] << "err=" << parms_err[j] << endl;
        
            TString sdmename;
            sdmename = sdme[j];
            sprintf (string,"%s=\t%.3f #pm %.3f\n",sdmename.Data(),parms[j],parms_err[j]);
        	printf("string=%s",string);
        	TLatex *t1 = new TLatex(0.2,0.85 - 0.05*j,string);
        	// t1->SetNDC();
        	t1->SetTextSize(0.04);
        	t1->Draw();
        }
        
    
    parameters.close();
    

    c0->SaveAs(filename+".pdf(");
    c1->SaveAs(filename+".pdf");
    c2->SaveAs(filename+".pdf");
    c3->SaveAs(filename+".pdf)");
}

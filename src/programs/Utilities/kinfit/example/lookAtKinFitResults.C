void lookAtKinFitResults(char *filename)
{

  gStyle->SetOptStat(0);

  gStyle->SetPadBottomMargin (0.18);
  gStyle->SetPadLeftMargin (0.18);

  char name[256];
  char title[256];
  float max = 1.0;

  TCanvas *can[3];
  //for(int i=0;i<3;i++)
  for(int i=0;i<1;i++)
  {
    sprintf(name,"can%d",i);
    if(i==0) sprintf(title,"Fit to nothing missing");
    else if(i==1) sprintf(title,"Fit to missing #pi^{0}");
    else if(i==2) sprintf(title,"Fit to missing K^{0} and constrained K^{0}_{s}");
    can[i] = new TCanvas(name,title,10+10*i,10+10*i,1000,800);
    can[i]->SetFillColor(0);
    can[i]->Divide(4,4);
  }

  TFile *f = new TFile(filename);

  TH1F *hmm[4];
  TH1F *hprob[4];
  TH1F *h2pi[4];
  TH1F *h3mesons[4];
  TH1F *hpulls[4][10];
  TPaveText *text[4][10];

  for(int i=0;i<4;i++) 
  {
    sprintf(name,"hmm%d",i);
    hmm[i] = (TH1F*)f->Get(name);
    hmm[i]->SetNdivisions(8);
    hmm[i]->SetFillColor(8);
    hmm[i]->GetYaxis()->SetTitle("# events");
    hmm[i]->GetYaxis()->SetTitleSize(0.09);
    hmm[i]->GetYaxis()->SetTitleFont(42);
    hmm[i]->GetYaxis()->SetTitleOffset(0.7);
    hmm[i]->GetYaxis()->CenterTitle();
    hmm[i]->GetXaxis()->SetTitle("MM^{2} (GeV^{2}/c^{4})");
    hmm[i]->GetXaxis()->SetLabelSize(0.05);
    hmm[i]->GetXaxis()->SetTitleSize(0.08);
    hmm[i]->GetXaxis()->SetTitleFont(42);
    hmm[i]->GetXaxis()->SetTitleOffset(0.7);
    hmm[i]->GetXaxis()->CenterTitle();

    sprintf(name,"hprob%d",i);
    hprob[i] = (TH1F*)f->Get(name);
    hprob[i]->SetNdivisions(8);
    hprob[i]->SetFillColor(28);
    hprob[i]->GetYaxis()->SetTitle("# events");
    hprob[i]->GetYaxis()->SetTitleSize(0.09);
    hprob[i]->GetYaxis()->SetTitleFont(42);
    hprob[i]->GetYaxis()->SetTitleOffset(0.7);
    hprob[i]->GetYaxis()->CenterTitle();
    hprob[i]->GetXaxis()->SetTitle("Confidence level");
    hprob[i]->GetXaxis()->SetLabelSize(0.05);
    hprob[i]->GetXaxis()->SetTitleSize(0.08);
    hprob[i]->GetXaxis()->SetTitleFont(42);
    hprob[i]->GetXaxis()->SetTitleOffset(0.7);
    hprob[i]->GetXaxis()->CenterTitle();
    hprob[i]->SetMinimum(0);
    max = 5.0*hprob[i]->GetBinContent(50);
    if(max<5.0) max = 5.0;
    hprob[i]->SetMaximum(max);

    sprintf(name,"h2pi%d",i);
    h2pi[i] = (TH1F*)f->Get(name);
    h2pi[i]->SetNdivisions(8);
    h2pi[i]->SetFillColor(39);
    h2pi[i]->GetYaxis()->SetTitle("# events");
    h2pi[i]->GetYaxis()->SetTitleSize(0.09);
    h2pi[i]->GetYaxis()->SetTitleFont(42);
    h2pi[i]->GetYaxis()->SetTitleOffset(0.7);
    h2pi[i]->GetYaxis()->CenterTitle();
    h2pi[i]->GetXaxis()->SetTitle("M(#pi^{+} #pi^{-}) GeV/c^{2}");
    h2pi[i]->GetXaxis()->SetLabelSize(0.05);
    h2pi[i]->GetXaxis()->SetTitleSize(0.08);
    h2pi[i]->GetXaxis()->SetTitleFont(42);
    h2pi[i]->GetXaxis()->SetTitleOffset(0.7);
    h2pi[i]->GetXaxis()->CenterTitle();
    h2pi[i]->SetMinimum(0);

    sprintf(name,"h3mesons%d",i);
    h3mesons[i] = (TH1F*)f->Get(name);
    h3mesons[i]->SetNdivisions(8);
    h3mesons[i]->SetFillColor(2);
    h3mesons[i]->GetYaxis()->SetTitle("# events");
    h3mesons[i]->GetYaxis()->SetTitleSize(0.09);
    h3mesons[i]->GetYaxis()->SetTitleFont(42);
    h3mesons[i]->GetYaxis()->SetTitleOffset(0.7);
    h3mesons[i]->GetYaxis()->CenterTitle();
    h3mesons[i]->GetXaxis()->SetTitle("M(#pi^{+} #pi^{-} MM) GeV/c^{2}");
    h3mesons[i]->GetXaxis()->SetLabelSize(0.05);
    h3mesons[i]->GetXaxis()->SetTitleSize(0.08);
    h3mesons[i]->GetXaxis()->SetTitleFont(42);
    h3mesons[i]->GetXaxis()->SetTitleOffset(0.7);
    h3mesons[i]->GetXaxis()->CenterTitle();
    h3mesons[i]->SetMinimum(0);

    for(int j=0;j<10;j++)
    {
      if(j==0) sprintf(title,"Pull - E_{#gamma}");
      else if(j==1) sprintf(title,"Pull - proton |p|");
      else if(j==2) sprintf(title,"Pull - proton #lambda");
      else if(j==3) sprintf(title,"Pull - proton #phi");
      else if(j==4) sprintf(title,"Pull - #pi^{+} |p|");
      else if(j==5) sprintf(title,"Pull - #pi^{+} #lambda");
      else if(j==6) sprintf(title,"Pull - #pi^{+} #phi");
      else if(j==7) sprintf(title,"Pull - #pi^{-} |p|");
      else if(j==8) sprintf(title,"Pull - #pi^{-} #lambda");
      else if(j==9) sprintf(title,"Pull - #pi^{-} #phi");
      sprintf(name,"hpulls%d_%d",i,j);
      hpulls[i][j] = (TH1F*)f->Get(name);
      hpulls[i][j]->SetNdivisions(8);
      hpulls[i][j]->SetFillColor(33);
      hpulls[i][j]->GetYaxis()->SetTitle("# events");
      hpulls[i][j]->GetYaxis()->SetTitleSize(0.09);
      hpulls[i][j]->GetYaxis()->SetTitleFont(42);
      hpulls[i][j]->GetYaxis()->SetTitleOffset(0.7);
      hpulls[i][j]->GetYaxis()->CenterTitle();
      hpulls[i][j]->GetXaxis()->SetTitle(title);
      hpulls[i][j]->GetXaxis()->SetLabelSize(0.05);
      hpulls[i][j]->GetXaxis()->SetTitleSize(0.08);
      hpulls[i][j]->GetXaxis()->SetTitleFont(42);
      hpulls[i][j]->GetXaxis()->SetTitleOffset(0.7);
      hpulls[i][j]->GetXaxis()->CenterTitle();
      hpulls[i][j]->SetMinimum(0);
    }
  }

  //for(int i=1;i<=3;i++) 
  for(int i=1;i<=1;i++) 
  {
    can[i-1]->cd(1);
    hprob[i]->Draw();

    for(int j=0;j<10;j++)
    {
      can[i-1]->cd(j+2);
      hpulls[i][j]->Draw();
      if(i==1) 
      {
        hpulls[i][j]->Fit("gaus","EQR","",-4,4);
        float mean = hpulls[i][j]->GetFunction("gaus")->GetParameter(1);
        float sigma = hpulls[i][j]->GetFunction("gaus")->GetParameter(2);
        text[i][j] = new TPaveText(0.65,0.7,0.99,0.99,"NDC");
        text[i][j]->SetBorderSize(0);
        text[i][j]->SetFillColor(0);
        sprintf(title,"mean: %3.3f",mean);
        text[i][j]->AddText(title);
        sprintf(title,"sigma: %3.3f",sigma);
        text[i][j]->AddText(title);
        text[i][j]->Draw();
      }
    }

    can[i-1]->cd(13);
    hmm[0]->Draw();

    can[i-1]->cd(14);
    hmm[i]->Draw();

    can[i-1]->cd(15);
    h2pi[i]->Draw();

    can[i-1]->cd(16);
    h3mesons[i]->Draw();
  }

}

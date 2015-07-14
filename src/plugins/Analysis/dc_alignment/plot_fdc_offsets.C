void plot_fdc_offsets(){
  gStyle->SetOptStat(0);

  TCanvas *c1=new TCanvas("c1","Alignment results",1200,800);
  c1->Divide(2,2);
  gStyle->SetTitleYOffset(1.5);

  // Histograms to store offset results
  TH1F *h1=new TH1F("h1","h1",24,0.5,24.5);
  h1->SetMaximum(0.5);
  h1->SetMinimum(-0.5);
  h1->SetXTitle("Wire plane");
  h1->SetYTitle("#Deltau [cm]");
  h1->SetTitle("Wire plane offsets");
  h1->SetMarkerStyle(20);

  TH1F *h1_true=new TH1F("h1_true","h1_true",24,0.5,24.5);
  h1_true->SetMaximum(0.5);
  h1_true->SetMinimum(-0.5);
  h1_true->SetXTitle("Wire plane");
  h1_true->SetYTitle("#Deltau [cm]");
  h1_true->SetTitle("Wire plane offsets");
  h1_true->SetMarkerStyle(20);


  TH1F *h2=new TH1F("h2","h2",24,0.5,24.5);
  h2->SetMaximum(20.0);
  h2->SetMinimum(-20.0);
  h2->SetXTitle("Wire plane ");
  h2->SetYTitle("#Delta#phi [mrad]");
  h2->SetTitle("#phi correction");
  h2->SetMarkerStyle(20);

  TH1F *h2_true=new TH1F("h2_true","h2_true",24,0.5,24.5);
  h2_true->SetMaximum(20.0);
  h2_true->SetMinimum(-20.0);
  h2_true->SetXTitle("Wire plane ");
  h2_true->SetYTitle("#Delta#phi [mrad]");
  h2_true->SetTitle("#phi correction");
  h2_true->SetMarkerStyle(20);


  // Histograms to store offset results
  TH1F *h3=new TH1F("h3","h3",24,0.5,24.5);
  h3->SetMaximum(0.5);
  h3->SetMinimum(-0.5);
  h3->SetXTitle("Wire plane");
  h3->SetYTitle("#Deltau - #Deltau(true) [cm]");
  h3->SetTitle("Wire plane offset residuals");
  h3->SetMarkerStyle(20);

  TH1F *h4=new TH1F("h4","h4",24,0.5,24.5);
  h4->SetMaximum(20.0);
  h4->SetMinimum(-20.0);
  h4->SetXTitle("Wire plane ");
  h4->SetYTitle("#Delta#phi - #Delta#phi(true) [mrad]");
  h4->SetTitle("#phi correction residuals");
  h4->SetMarkerStyle(20);
  
  ifstream fdcfile("fdc_alignment.dat");
  ifstream fdctrue("fdc_alignment_true.dat");

   // Skip first line, used to identify columns in file
  char sdummy[40];
  //  fdcfile.getline(sdummy,40);
  //fdctrue.getline(sdummy,40);
  // loop over remaining entries
  for (unsigned int i=0;i<24;i++){
    unsigned int mylayer=i+1;
    double du,dphi,dz;
    double du_true,dphi_true,dz_true;
    
    fdcfile >> dphi;
    fdcfile >> du;
    fdcfile >> dz;

    fdctrue >> dphi_true;
    fdctrue >> du_true;
    fdctrue >> dz_true;

    h1->Fill(mylayer,du);
    h1_true->Fill(mylayer,-du_true);
    h3->Fill(mylayer,du+du_true);
    //h1->SetBinError(mylayer,sigu);
  
    h2->Fill(mylayer,1000.*dphi);
    h2_true->Fill(mylayer,1000.*dphi_true);
    h4->Fill(mylayer,1000.*(dphi-dphi_true));
    //h2->SetBinError(mylayer,sigphi);

  }
  fdcfile.close();

  c1->cd(1);
  h1->Draw("p");
  h1_true->SetMarkerColor(2);
  h1_true->Draw("p,same");

  c1->cd(2);
  h2->Draw("p");
  h2_true->SetMarkerColor(2);
  h2_true->Draw("p,same");

  c1->cd(3);
  h3->Draw("p");

  c1->cd(4);
  h4->Draw("p");
}

void plot_fdc_offsets(){
  gStyle->SetOptStat(0);

  TCanvas *c1=new TCanvas("c1","Alignment results",1200,800);
  c1->Divide(2,2);
  gStyle->SetTitleYOffset(1.5);

  // Histograms to store offset results
  TH1F *h1=new TH1F("h1","h1",24,0.5,24.5);
  h1->SetMaximum(0.2);
  h1->SetMinimum(-0.2);
  h1->SetXTitle("Wire plane");
  h1->SetYTitle("#Deltau [cm]");
  h1->SetTitle("Wire plane offsets");
  h1->SetMarkerStyle(20);

  TH1F *h2=new TH1F("h2","h2",24,0.5,24.5);
  h2->SetMaximum(2.0);
  h2->SetMinimum(-2.0);
  h2->SetXTitle("Wire plane ");
  h2->SetYTitle("#Delta#phi [degrees]");
  h2->SetTitle("#phi correction");
  h2->SetMarkerStyle(20);
  
  //TF1 *f1=new TF1("f1","-1.*(0.1*sin(0.3*x)+0.05*sin(0.6*x))",1,24);
  TF1 *f1=new TF1("f1","1.*(0.1*sin(1.3*x)+0.05*sin(1.6*x))",1,24);
  TF1 *f2=new TF1("f2","1.*(-0.1*cos(0.3*x)+0.15*cos(0.6*x))",1,24);

  //  TF1 *f3=new TF1("f3","180./3.14159*(-0.01*cos(0.6*x)+0.005*sin(0.3*x))",1,24);
  //TF1 *f3=new TF1("f3","180./3.14159*(-0.01*cos(0.6*x))",1,24);
  TF1 *f3=new TF1("f3","180./3.14159*(-0.01*cos(1.6*x)+0.005*sin(1.3*x))",1,24);


  ifstream fdcfile("fdc_alignment.dat");

   // Skip first line, used to indentify columns in file
  char sdummy[40];
  fdcfile.getline(sdummy,40);
  // loop over remaining entries
  for (unsigned int i=0;i<24;i++){ 
    int mylayer;
    double du,dphi;

    fdcfile >> mylayer;
    fdcfile >> du;
    fdcfile >> dphi;

    h1->Fill(mylayer,du);
    h2->Fill(mylayer,180./3.1415926*dphi);


  }
  fdcfile.close();

  c1->cd(1);
  h1->Draw("p");
  f1->Draw("same");

  c1->cd(2);
  h2->Draw("p");
  f3->Draw("same");
}

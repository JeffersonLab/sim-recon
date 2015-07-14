// Macro to grab offset data for the cdc from the file cdc_alignment.dat and display the results
// for a particular ring.
void plot_cdc_offsets(unsigned int ring){
  unsigned int index=ring-1;
  gStyle->SetOptStat(0);

  TCanvas *c1=new TCanvas("c1","Alignment results",1200,800);
  c1->Divide(2,2);
  gStyle->SetTitleYOffset(1.5);
 
  ifstream cdcfile("cdc_alignment.dat");

  unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
			      135,135,146,146,158,158,170,170,182,182,197,197,
			      209,209};

  // Histograms to store offset results
  TH1F *h1=new TH1F("h1","h1",numstraws[index],0.5,numstraws[index]+0.5);
  h1->SetMaximum(0.1);
  h1->SetMinimum(-0.1);
  h1->SetXTitle("Straw");
  h1->SetYTitle("#Deltax [cm]");
  char mytitle[80];
  sprintf(mytitle,"Offsets at wire center for Ring %d\n",ring);
  h1->SetTitle(mytitle);
  h1->SetMarkerStyle(20);

  TH1F *h2=new TH1F("h2","h2",numstraws[index],0.5,numstraws[index]+0.5);
  h2->SetMaximum(0.1);
  h2->SetMinimum(-0.1);
  h2->SetXTitle("Straw");
  h2->SetYTitle("#Deltay [cm]");
  h2->SetTitle(mytitle);
  h2->SetMarkerStyle(20);

  TH1F *h3=new TH1F("h3","h3",numstraws[index],0.5,numstraws[index]+0.5);
  h3->SetMaximum(0.01);
  h3->SetMinimum(-0.01);
  h3->SetXTitle("Straw");
  h3->SetYTitle("dx/dz");
  char mytitle[80];
  sprintf(mytitle,"Slopes for Ring %d\n",ring);
  h3->SetTitle(mytitle);
  h3->SetMarkerStyle(20);

  TH1F *h4=new TH1F("h4","h4",numstraws[index],0.5,numstraws[index]+0.5);	
  h4->SetMaximum(0.01);
  h4->SetMinimum(-0.01);
  h4->SetXTitle("Straw");
  h4->SetYTitle("dx/dz");
  h4->SetTitle(mytitle);
  h4->SetMarkerStyle(20);
  
  // Skip first line, used to indentify columns in file
  char sdummy[40];
  cdcfile.getline(sdummy,40);
  // loop over straws; a bit clumsy, but...
  for (unsigned int i=0;i<28;i++){    
    for (unsigned int j=0;j<numstraws[i];j++){
      int myring,mystraw;
      double dxd,dyd,dxu,dyu;
      cdcfile >> myring;
      cdcfile >> mystraw;
      cdcfile >> dxu;
      cdcfile >> dyu;
      cdcfile >> dxd;
      cdcfile >> dyd;

      if (myring==ring){
	h1->Fill(mystraw,0.5*(dxu+dxd));
	h2->Fill(mystraw,0.5*(dyu+dyd));
	h3->Fill(mystraw,(dxd-dxu)/150.);
	h4->Fill(mystraw,(dyd-dyu)/150.);
      }
    }
  }
  
  // The following functions were used to mock up the wire deflections
  TF1 *f1=new TF1("f1","[1]*(-0.05*(sin(0.1*(x-1+[0]))+sin(0.2*(x-1+[0]))))",0,numstraws[index]+1);
  TF1 *f2=new TF1("f2","[1]*(-0.05*(cos(0.1*(x-1+[0]))+cos(0.2*(x-1+[0]))))",0,numstraws[index]+1);

  TF1 *f3=new TF1("f3","[1]*(-0.05*(sin(0.1*(x-1+[0]))-sin(0.2*(x-1+[0]))))/75.",0,numstraws[index]+1);
  TF1 *f4=new TF1("f4","[1]*(-0.05*(cos(0.1*(x-1+[0]))-cos(0.2*(x-1+[0]))))/75.",0,numstraws[index]+1);

  c1->cd(1);

  h1->Draw("p");
  f1->SetParameters(ring,0.); 
  f1->Draw("same");
 
  c1->cd(2);

  h2->Draw("p");  
  f2->SetParameters(ring,0);
  f2->Draw("same");

  c1->cd(3);

  h3->Draw("p");
  f3->SetParameters(ring,0);
  f3->Draw("same");

  c1->cd(4);

  h4->Draw("p");
  f4->SetParameters(ring,0);
  f4->Draw("same");

  cdcfile.close();
}

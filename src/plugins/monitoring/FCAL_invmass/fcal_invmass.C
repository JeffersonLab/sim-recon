// hnamepath:  /FCAL_invmass/InvMass1
// hnamepath:  /FCAL_invmass/InvMass2
// hnamepath:  /FCAL_invmass/InvMass3
// hnamepath:  /FCAL_invmass/InvMass4
// hnamepath:  /FCAL_invmass/InvMass5
// hnamepath:  /FCAL_invmass/InvMass6
// hnamepath:  /FCAL_invmass/InvMass7
// hnamepath:  /FCAL_invmass/InvMass8
// hnamepath:  /FCAL_invmass/InvMass9

{
    
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("FCAL_invmass");
  if(dir) dir->cd();

  TH1I* InvMass1 = (TH1I*)gDirectory->FindObjectAny("InvMass1");
  TH1I* InvMass2 = (TH1I*)gDirectory->FindObjectAny("InvMass2");
  TH1I* InvMass3 = (TH1I*)gDirectory->FindObjectAny("InvMass3");
  TH1I* InvMass4 = (TH1I*)gDirectory->FindObjectAny("InvMass4");
  TH1I* InvMass5 = (TH1I*)gDirectory->FindObjectAny("InvMass5");
  TH1I* InvMass6 = (TH1I*)gDirectory->FindObjectAny("InvMass6");
  TH1I* InvMass7 = (TH1I*)gDirectory->FindObjectAny("InvMass7");
  TH1I* InvMass8 = (TH1I*)gDirectory->FindObjectAny("InvMass8");
  TH1I* InvMass9 = (TH1I*)gDirectory->FindObjectAny("InvMass9");
  
  
  float par_1[12]; 
  float par_2[12]; 
  float par_3[12]; 
  float par_4[12]; 
  float par_5[12]; 
  float par_6[12];
  float par_7[12]; 
  float par_8[12]; 
  float par_9[12];  
  
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "FCAL_invmass_plot", 1200, 1200 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide(3,3);

  if( InvMass1 ){

    InvMass1->SetStats(0);
    double max1 = InvMass1->GetMaximum();
    c1->cd(1);
    InvMass1->Draw();
    InvMass1->SetLineWidth(2);
    
   Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass1->Fit(g1,"nR");
   InvMass1->Fit(g2,"nR+");
   InvMass1->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass1->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_1[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_1[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_1[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_1[5]/par_1[4])*100));
    pt_300->Draw();
  }
  
  if( InvMass2 ){

    InvMass2->SetStats(0);
    double max1 = InvMass2->GetMaximum();
    c1->cd(2);
    InvMass2->Draw();
    InvMass2->SetLineWidth(2);
    
  Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass2->Fit(g1,"nR");
   InvMass2->Fit(g2,"nR+");
   InvMass2->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass2->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_2[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_2[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_2[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_2[5]/par_2[4])*100));
    pt_300->Draw();
  }
  
  if( InvMass3 ){

    InvMass3->SetStats(0);
    double max1 = InvMass3->GetMaximum();
    c1->cd(3);
    InvMass3->Draw();
    InvMass3->SetLineWidth(2);
    
   Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass3->Fit(g1,"nR");
   InvMass3->Fit(g2,"nR+");
   InvMass3->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass3->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_3[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_3[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_3[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_3[5]/par_3[4])*100));
    pt_300->Draw();
  }
  
  if( InvMass4 ){

    InvMass4->SetStats(0);
    double max1 = InvMass4->GetMaximum();
    c1->cd(4);
    InvMass4->Draw();
    InvMass4->SetLineWidth(2);
    
   Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass4->Fit(g1,"nR");
   InvMass4->Fit(g2,"nR+");
   InvMass4->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass4->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_4[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_4[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_4[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_4[5]/par_4[4])*100));
    pt_300->Draw();
  }
  
  if( InvMass5 ){

    InvMass5->SetStats(0);
    double max1 = InvMass5->GetMaximum();
    c1->cd(5);
    InvMass5->Draw();
    InvMass5->SetLineWidth(2);
    
   Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass5->Fit(g1,"nR");
   InvMass5->Fit(g2,"nR+");
   InvMass5->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass5->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_5[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_5[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_5[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_5[5]/par_5[4])*100));
    pt_300->Draw();
  }
  
  if( InvMass6 ){

    InvMass6->SetStats(0);
    double max1 = InvMass6->GetMaximum();
    c1->cd(6);
    InvMass6->Draw();
    InvMass6->SetLineWidth(2);
    
   Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass6->Fit(g1,"nR");
   InvMass6->Fit(g2,"nR+");
   InvMass6->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass6->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_6[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_6[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_6[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_6[5]/par_6[4])*100));
    pt_300->Draw();
  }
  
if( InvMass7 ){

    InvMass7->SetStats(0);
    double max1 = InvMass7->GetMaximum();
    c1->cd(7);
    InvMass7->Draw();
    InvMass7->SetLineWidth(2);
    
   Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass7->Fit(g1,"nR");
   InvMass7->Fit(g2,"nR+");
   InvMass7->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass7->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_7[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_7[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_7[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_7[5]/par_7[4])*100));
    pt_300->Draw();
  }
  
 if( InvMass8 ){

    InvMass8->SetStats(0);
    double max1 = InvMass8->GetMaximum();
    c1->cd(8);
    InvMass8->Draw();
    InvMass8->SetLineWidth(2);
    
   Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass8->Fit(g1,"nR");
   InvMass8->Fit(g2,"nR+");
   InvMass8->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass8->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_8[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_8[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_8[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_8[5]/par_8[4])*100));
    pt_300->Draw();
  } 
  
  if( InvMass9 ){

    InvMass9->SetStats(0);
    double max1 = InvMass9->GetMaximum();
    c1->cd(9);
    InvMass9->Draw();
    InvMass9->SetLineWidth(2);
    
   Double_t par[12] = {};  
   TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
   TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
   TF1 *g3    = new TF1("g3","pol2",0.15,0.28);


   TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
   total->SetLineColor(kRed);
   InvMass9->Fit(g1,"nR");
   InvMass9->Fit(g2,"nR+");
   InvMass9->Fit(g3,"nR+");
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   total->SetParameters(par);
   
   total->SetParName(3, "Amplitude");
   total->SetParName(4, "Pi0 Mass");
   total->SetParName(5, "Pi0 Width");

   InvMass9->Fit(total,"R+"); 
   
    for (int i=0; i<12; i++) {
         par_9[i] = total->GetParameter(i);
    }
    
    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_9[4]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_9[5]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_9[5]/par_9[4])*100));
    pt_300->Draw();
  }
}

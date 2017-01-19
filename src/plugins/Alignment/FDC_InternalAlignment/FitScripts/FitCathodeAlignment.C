void FitCathodeAlignment(TString inputFile = "hd_root.root"){

   TFile *file = TFile::Open(inputFile);
   TFile *outFile = new TFile("result.root","RECREATE");

   TH2D *resultHists[24];
   // The old constants are stored in a histogram
   TH1D * constantsHist = (TH1D *) file->Get("FDC_InternalAlignment/CathodeAlignmentConstants");

   // Fit Function
   TF2 *f = new TF2("plane","[0]+[1]*x+[2]*y", -50.,50.,-50.,50.);
   ofstream constantsFile;
   constantsFile.open("CathodeAlignment.txt");

   for(int i=1; i<=24; i++){
      double i_deltaPhiU = constantsHist->GetBinContent((i-1)*4+1);
      double i_deltaU = constantsHist->GetBinContent((i-1)*4+2);
      double i_deltaPhiV = constantsHist->GetBinContent((i-1)*4+3);
      double i_deltaV = constantsHist->GetBinContent((i-1)*4+4);

      double phi_u = 75.0*TMath::DegToRad() + i_deltaPhiU;
      double phi_v = (TMath::Pi() - 75.0*TMath::DegToRad()) + i_deltaPhiV;

      double sinPhiU = sin(phi_u);
      double cosPhiU = cos(phi_u);
      double sinPhiV = sin(phi_v);
      double cosPhiV = cos(phi_v);
      double sinPhiUmPhiV = sin(phi_u-phi_v);
      double sinPhiUmPhiV2 = sinPhiUmPhiV*sinPhiUmPhiV;
      double cosPhiUmPhiV = cos(phi_u-phi_v);

      char name[100];
      sprintf(name,"FDC_InternalAlignment/Plane %.2i Wire Position Vs XY",i);

      TH3I *h = (TH3I *)file->Get(name);

      TF1 *g = new TF1("g", "gaus", -0.5, 0.5);
      TH2D *result = new TH2D(Form("Plane %.2i residual", i), "residual", 100, -50., 50., 100, -50., 50.);
      resultHists[i-1]=result;
      for (int iX=1; iX <= h->GetNbinsX(); iX++){
         for (int iY=1; iY <= h->GetNbinsY(); iY++){
            TString name;
           if(iY==25 && iX%5 ==0) name = Form("Proj bin (%i,%i)",iY,iX);
           else name = "_z"; 
            TH1D * zProj = h->ProjectionZ(name,iX, iX, iY, iY);
            if (zProj->GetEntries() < 5) {
               result->SetBinContent(iX,iY, -1./0.0);
               continue;
            }
            double max = zProj->GetBinCenter(zProj->GetMaximumBin());
            g->SetParameter(1, max);
            TFitResultPtr r = zProj->Fit("g", "Q", "", max - 0.1, max + 0.1 );
            if ((Int_t) r == 0){
               double mean = g->GetParameter(1);
               double err = g->GetParError(1);
               result->SetBinContent(iX,iY,mean);
               result->SetBinError(iX,iY,err);
            }
            else result->SetBinContent(iX,iY, -1./0.0);
         }
      }

      // Use robust fitting to reject up to 25% outliers.
      result->Fit("plane", "+rob=0.75");

      double c0 = f->GetParameter(0);
      double c1 = f->GetParameter(1);
      double c2 = f->GetParameter(2);

      // Calculate the new offsets from the fit parameters
      double deltaU = sinPhiUmPhiV*c0/sinPhiV;
      double deltaV = -sinPhiUmPhiV*c0/sinPhiU;
      double a = -sinPhiU*sinPhiV/sinPhiUmPhiV;
      double b = -a;
      double c = cosPhiU*sinPhiV/sinPhiUmPhiV;
      double d = -cosPhiV*sinPhiU/sinPhiUmPhiV;

      double deltaPhiU = (c1/b-c2/d)/(a/b-c/d);
      double deltaPhiV = (c1/a-c2/c)/(b/a-d/c);

      // Corrections are opposite direction for u and y as delta phi? Implies some sign swapping somewhere...
      cout << " Current Iteration shift plane " << i <<  " deltaU " << -deltaU << " deltaPhiU " << -deltaPhiU << " deltaPhiV " << -deltaPhiV << endl;
      constantsFile << i_deltaPhiU - deltaPhiU << " " << i_deltaU - deltaU << " " << i_deltaPhiV - deltaPhiV << " " << i_deltaV << endl;
   }

   // Draw the results
   TCanvas *cResult = new TCanvas("cResult","cResult",800,600);
   cResult->Divide(6,4);

   for (int i=0; i<24; i++){
      cResult->cd(i+1);
      resultHists[i]->Draw("colz");
      resultHists[i]->SetMinimum(-0.025);
      resultHists[i]->SetMaximum(0.025);
   }

   cResult->SaveAs("Result.png");
   cResult->SaveAs("Result.pdf");
   cResult->Write();

   outFile->Write();
   outFile->Close();
   constantsFile.close();

   return;
}

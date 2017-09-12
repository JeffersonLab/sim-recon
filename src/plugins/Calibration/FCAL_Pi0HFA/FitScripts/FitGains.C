{
   // Script used for fitting the output of the FCAL_Pi0HFA plugin
   // Invoke using root -l -b <rootfile> FitGains.C
   // Force Batch Mode
   gROOT->SetBatch();

   TProfile *hCurrentGainConstants = (TProfile *) gDirectory->Get("FCAL_Pi0HFA/CurrentGainConstants");
   TH1I *hPi0Mass = (TH1I *) gDirectory->Get("FCAL_Pi0HFA/Pi0Mass");
   TH2I *hPi0MassVsChNum = (TH2I *) gDirectory->Get("FCAL_Pi0HFA/Pi0MassVsChNum");
   TH2F *hPi0MassVsChNumWeighted = (TH2F *) gDirectory->Get("FCAL_Pi0HFA/Pi0MassVsChNumWeighted"); 
   TH2F *hPi0MassVsChNumWeightedSquared = (TH2F *) gDirectory->Get("FCAL_Pi0HFA/Pi0MassVsChNumWeightedSquared");
   TH2I *hPi0MassVsE = (TH2I *) gDirectory->Get("FCAL_Pi0HFA/Pi0MassVsE");

   // Make an output file
   TFile *outFile = new TFile("FCALPi0FitResults.root","RECREATE");
   double fitPi0Mean[2800];

   TH1F *hFitPi0Mass = new TH1F("FitPi0Mass", "#pi^{0} Mass fit result", 2800, -0.5, 2799.5);
   TH1F *hFitGain = new TH1F("FitGain", "new gain result", 2800, -0.5, 2799.5);
   // Our fit function
   TF1 *f = new TF1("f","landau(0)+gaus(3)", 0.025, 0.250);
   f->SetParLimits(0,0.0, 1e8);
   f->SetParLimits(1,0.025, 0.080);
   f->SetParLimits(2,0.0, 10.0);
   f->SetParLimits(3,0.0, 1e8);
   f->SetParLimits(4,0.04, 0.25);
   f->SetParLimits(5,0.004, 0.013);
   TF1 *f2 = new TF1("f2","gaus", 0.025, 0.250);
   f2->SetParLimits(0,0.0, 1e6);
   f2->SetParLimits(1,0.05, 0.19);
   f2->SetParLimits(2,0.006, 0.018);
   TCanvas *cFits = new TCanvas("cFits", "cFits", 800, 800);
   cFits->Divide(8,8);
   TCanvas *c = new TCanvas("c", "c", 1200, 1200);
   for (unsigned int i = 1 ; i <= hPi0MassVsChNum->GetNbinsX(); i++){
      c->cd();
      //TH1D *projY = hPi0MassVsChNum->ProjectionY(Form("Channel%.4i", i), i, i);
      TH1D *projY = hPi0MassVsChNumWeightedSquared->ProjectionY(Form("Channel%.4i", i), i, i);

      if (projY->GetEntries() < 400) {
         fitPi0Mean[i-1] = -1.0;
         if (i%64 != 0){
            cFits->cd(i%64);
            projY->Draw();
         }
         else{
            cFits->cd(64);
            projY->Draw();
            cFits->SaveAs(Form("c_%.2i.png",i/64));
            cFits->Clear("D");
         }
         continue;
      }
      // If the maximium bin is near the pi0 mass, switch to gaussian
      double max = projY->GetBinCenter( projY->GetMaximumBin() );
      //cout << max << endl;
      if (i < 1660 && i > 1140 && fabs(max-0.135) > 0.035){
         // Only need to try this fit with the background on the innermost channels
         // Clone function
         TF1 *fClone = (TF1 *) f->Clone(Form("f%.4i", i));
         fClone->SetParLimits(0,0.0, 1e8);
         fClone->SetParLimits(1,0.025, 0.080);
         fClone->SetParLimits(2,0.0, 10.0);
         fClone->SetParLimits(3,0.0, 1e8);
         fClone->SetParLimits(4,0.04, 0.25);
         fClone->SetParLimits(5,0.004, 0.013);
         fClone->SetParameters(projY->GetEntries() / 10.0, 0.05, 0.02, projY->GetEntries() / 100.0, 0.135, 0.01);
         TFitResultPtr r = projY->Fit(Form("f%.4i", i),"SQ+", "", 0.05, 0.3);

         if ((Int_t) r == 0){
            // Fit Converged check chi^2/NDF
            double reducedChi2 = r->Chi2() / r->Ndf();
            // Make some cut here
            if (reducedChi2 < 10.0){
               //Trust the result
               fitPi0Mean[i-1] = r->Parameter(4);
               hFitPi0Mass->SetBinContent(i,r->Parameter(4));
               projY->GetFunction(Form("f%.4i", i))->SetLineColor(kGreen);
               gPad->Update();
            }
            else{
               fitPi0Mean[i-1] = -1.0;
               projY->GetFunction(Form("f%.4i", i))->SetLineColor(kRed);
               // Try plain Gaussian fit
               f2->SetParameters(projY->GetEntries() / 100.0, 0.135, 0.01);
               r = projY->Fit("f2","SQ+", "", 0.05, 0.18);
               if ((Int_t) r == 0){
                  // Fit Converged check chi^2/NDF
                  double reducedChi2 = r->Chi2() / r->Ndf();
                  // Make some cut here
                  if (reducedChi2 < 100.0){
                     //Trust the result
                     fitPi0Mean[i-1] = r->Parameter(1);
                     hFitPi0Mass->SetBinContent(i,r->Parameter(1));
                     projY->GetFunction("f2")->SetLineColor(kGreen);
                     gPad->Update();
                  }
                  else{
                     fitPi0Mean[i-1] = -1.0;
                     projY->GetFunction("f2")->SetLineColor(kRed);
                  }
               }
               else{
                  fitPi0Mean[i-1] = -1.0;
               }
            }
         }
      }
      else{
         // Try plain Gaussian fit
         f2->SetParLimits(0,0.0, 1e6);
         f2->SetParLimits(1,0.05, 0.19);
         f2->SetParLimits(2,0.006, 0.018);
         f2->SetParameters(projY->GetEntries() / 100.0, max, 0.01);
         TFitResultPtr r = projY->Fit("f2","SQ+", "", max - 0.01, max + 0.01);
         if ((Int_t) r == 0){
            // Fit Converged check chi^2/NDF
            double reducedChi2 = r->Chi2() / r->Ndf();
            // Make some cut here
            if (reducedChi2 < 100.0){
               //Trust the result
               fitPi0Mean[i-1] = r->Parameter(1);
               hFitPi0Mass->SetBinContent(i,r->Parameter(1));
               projY->GetFunction("f2")->SetLineColor(kGreen);
               gPad->Update();
            }
            else{
               fitPi0Mean[i-1] = -1.0;
               projY->GetFunction("f2")->SetLineColor(kRed);
            }
         }
         else{
            fitPi0Mean[i-1] = -1.0;
         }
      }
      if (i%64 != 0){
         cFits->cd(i%64);
         projY->Draw();
      }
      else{
         cFits->cd(64);
         projY->Draw();
         cFits->SaveAs(Form("c_%.2i.png",i/64));
         cFits->Clear("D");
      }
   }
   cFits->SaveAs("c_final.png");

   double pdgPi0Mass = 0.1349766;
   double newGains[2800];
   double total = 0., counter = 0.;
   TH1I *hGains = new TH1I("NewGains", "New Gains" , 100, 0.0, 5.0);
   /*
   // Get old gains from file
   double oldGains[2800];
   ifstream textIn;
   textIn.open("oldGains.txt");

   for (unsigned int i=0; i<2800; i++){
   textIn >> oldGains[i];
   }

   textIn.close();
   */
   for (unsigned int i=0; i < 2800; i++){
      if (fitPi0Mean[i] > 0.){
         //newGains[i] = hCurrentGainConstants->GetBinContent(i+1) * (0.5 + 0.5*pdgPi0Mass / fitPi0Mean[i]);  
         newGains[i] = hCurrentGainConstants->GetBinContent(i+1) * (pdgPi0Mass / fitPi0Mean[i]);
         //newGains[i] = oldGains[i] * pdgPi0Mass / fitPi0Mean[i];
         counter += 1.0;
         total += newGains[i];
         hGains->Fill(newGains[i]);
      }
      else{
         newGains[i] = -1.0;
      }
   }

   double averageGain = total / counter;
   for (unsigned int i=0; i < 2800; i++){
      if (newGains[i] < 0.){
         newGains[i] = averageGain;
      }
   }

   // Dump the gains
   ofstream textOut;
   textOut.open("gains.txt");

   counter=1;
   for (auto i: newGains){
      textOut << i << endl;
      hFitGain->SetBinContent(counter,i);
      counter++;
   }

   textOut.close();

   outFile->Write();
   outFile->Close();
}

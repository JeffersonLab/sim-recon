#include <TMath.h>
#include <TH2I.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <fstream>
#include <TMath.h>

void FitCathodeProjections(TString inputFile = "All.root"){

   bool shiftPitch = true;
   bool shiftGap = true;
   bool shiftOffset = false;

   TFile *file = TFile::Open(inputFile);
   TFile *outFile = new TFile("projection_result.root","RECREATE");

   TH2I *hists2DU[24];
   TH2I *hists2DV[24];
   TH1D *resultHistsU[24];
   TH1D *resultHistsV[24];
   TF1  *fLeftU[24];
   TF1  *fLeftV[24];
   TF1  *fMiddleU[24];
   TF1  *fMiddleV[24];
   TF1  *fRightU[24];
   TF1  *fRightV[24];

   // The old constants are stored in a histogram
   TH1D * constantsHist = (TH1D *) file->Get("FDC_InternalAlignment/CathodeAlignmentConstants");

   // Fit Function
   ofstream constantsFile, pitchFile;
   constantsFile.open("CathodeAlignment.txt");
   pitchFile.open("StripPitchesV2.txt");
   for(int i=1; i<=24; i++){
      double i_deltaPhiU = constantsHist->GetBinContent((i-1)*4+1);
      double i_deltaU = constantsHist->GetBinContent((i-1)*4+2);
      double i_deltaPhiV = constantsHist->GetBinContent((i-1)*4+3);
      double i_deltaV = constantsHist->GetBinContent((i-1)*4+4);

      double deltaX=constantsHist->GetBinContent((i-1)*2+401);

      double phi_u = 75.0*TMath::DegToRad() + i_deltaPhiU;
      double phi_v = (TMath::Pi() - 75.0*TMath::DegToRad()) + i_deltaPhiV;

      double sinPhiU = sin(phi_u);
      double cosPhiU = cos(phi_u);
      double sinPhiV = sin(phi_v);
      double cosPhiV = cos(phi_v);
      double sinPhiUmPhiV = sin(phi_u-phi_v);
      double sinPhiUmPhiV2 = sinPhiUmPhiV*sinPhiUmPhiV;
      double cosPhiUmPhiV = cos(phi_u-phi_v);

      double deltaU=0., deltaV=0.;

      for (unsigned int iPlane=1; iPlane<=2; iPlane++){ 
         int offset = (iPlane == 1) ? 0 : 5;
         double i_SP1 = constantsHist->GetBinContent((i-1)*10+101+offset);
         double i_G1  = constantsHist->GetBinContent((i-1)*10+102+offset);
         double i_SP2 = constantsHist->GetBinContent((i-1)*10+103+offset);
         double i_G2  = constantsHist->GetBinContent((i-1)*10+104+offset);
         double i_SP3 = constantsHist->GetBinContent((i-1)*10+105+offset);


         char name[100];
         if(iPlane == 1 ) sprintf(name,"FDC_InternalAlignment/Cathode Projections/Plane %.2i u_res vs u",i);
         else sprintf(name,"FDC_InternalAlignment/Cathode Projections/Plane %.2i v_res vs v",i);
         
         TH2I *h = (TH2I *)file->Get(name);
         if(iPlane == 1) hists2DU[i-1]=h;
         else hists2DV[i-1]=h;

         TF1 *g = new TF1("g", "gaus", -0.5, 0.5);
         if(iPlane == 1) sprintf(name,"Plane %.2i u residual", i);
         else sprintf(name,"Plane %.2i v residual", i);
         TH1D *result = new TH1D(name, "meas-pred", h->GetNbinsX(), -50., 50.);
         if(iPlane == 1) resultHistsU[i-1]=result;
         else resultHistsV[i-1]=result;
         for (int iX=3; iX <= h->GetNbinsX(); iX++){
            TString name;
            if(iX%5 ==0) name = Form("Proj bin (%i)",iX);
            else name = "_z"; 
            TH1D * yProj = h->ProjectionY(name,iX-2, iX);
            if (yProj->GetEntries() < 5) {
               //result->SetBinContent(iX, 0.0);
               continue;
            }
            double max = yProj->GetBinCenter(yProj->GetMaximumBin());
            g->SetParameter(1, max);
            TFitResultPtr r = yProj->Fit("g", "Q", "", max - 0.01, max + 0.01 );
            if ((Int_t) r == 0){
               double mean = g->GetParameter(1);
               double err = g->GetParError(1);
               result->SetBinContent(iX,mean);
               result->SetBinError(iX,err);
            }
            //else result->SetBinContent(iX, 0.0);
         }

         // Use robust fitting to reject up to 10% outliers.

         TString f_name = Form("f%.2i",i);
         TString fLeft_name = Form("fLeft%.2i",i);
         TString fMiddle_name = Form("fMiddle%.2i",i);
         TString fRight_name = Form("fRight%.2i",i);
         if (iPlane == 1) {
            f_name+="u";
            fLeft_name+="u";
            fMiddle_name+="u";
            fRight_name+="u";
         }
         TF1 *f = new TF1(f_name,"[0]+[1]*x", -50.,50.);
         TF1 *fLeft = new TF1(fLeft_name,"[0]+[1]*(x+23.75)", -46.,-25.);
         TF1 *fMiddle = new TF1(fMiddle_name,"[0]+[1]*x", -21.,21.);
         TF1 *fRight = new TF1(fRight_name,"[0]+[1]*(x-23.75)", 25.,46.);

         result->Fit(f_name, "+rob=0.80");
         if (iPlane ==1) deltaU = f->GetParameter(0);
         else deltaV = f->GetParameter(0);
         double slope = f->GetParameter(1);

         result->Fit(fLeft_name, "MR+rob=0.8");
         double deltaSP1=0.5*(fLeft->GetParameter(1)-slope);

         result->Fit(fMiddle_name, "MR+rob=0.85");
         double deltaSP2=0.5*(fMiddle->GetParameter(1)-slope);

         result->Fit(fRight_name, "MR+rob=0.8");
         double deltaSP3=0.5*(fRight->GetParameter(1)-slope);

         if (iPlane ==1) {
            fLeftU[i-1]=fLeft;
            fMiddleU[i-1]=fMiddle;
            fRightU[i-1]=fRight;
         }
         else  {
            fLeftV[i-1]=fLeft;
            fMiddleV[i-1]=fMiddle;
            fRightV[i-1]=fRight;
         }

         double deltaG1 = 0, deltaG2=0.;
         
         if (shiftGap){
            int nBins = h->GetNbinsX();
            TH1D * yProj1 = h->ProjectionY("projFoil1",1, nBins/4);
            TH1D * yProj2 = h->ProjectionY("projFoil2",nBins/4 , 3*nBins/4);
            TH1D * yProj3 = h->ProjectionY("projFoil3",3*nBins/4, nBins);

            double mean1 =0., mean2=0., mean3=0.;
            double max = yProj1->GetBinCenter(yProj1->GetMaximumBin());
            g->SetParameter(1, max);
            TFitResultPtr r = yProj1->Fit("g", "Q", "", max - 0.01, max + 0.01 );
            if ((Int_t) r == 0) mean1 = g->GetParameter(1);

            max = yProj2->GetBinCenter(yProj2->GetMaximumBin());
            g->SetParameter(1, max);
            r = yProj2->Fit("g", "Q", "", max - 0.01, max + 0.01 );
            if ((Int_t) r == 0) mean2 = g->GetParameter(1);

            max = yProj3->GetBinCenter(yProj3->GetMaximumBin());
            g->SetParameter(1, max);
            r = yProj3->Fit("g", "Q", "", max - 0.01, max + 0.01 );
            if ((Int_t) r == 0) mean3 = g->GetParameter(1);

            deltaG1 = mean1 - mean2;
            deltaG2 = mean2 - mean3;
         }

         // Calculate the new offsets from the fit parameters
         if (iPlane ==1 )cout << " Current Iteration shift plane " << i <<  " deltaU " << deltaU << endl;
         else cout << " shift plane " << i <<  " deltaV " << deltaV << endl;
         cout << "dSP1 " << deltaSP1 << " dG1 " << deltaG1 << " dSP2 " << deltaSP2 << " dG2 " << deltaG2 << " dSP3 " << deltaSP3 << endl;

         pitchFile << (shiftPitch ? (i_SP1 - deltaSP1/2.):i_SP1) << " "; 
         pitchFile << (shiftGap   ? (i_G1 + deltaG1)  :i_G1)  << " ";
         pitchFile << (shiftPitch ? (i_SP2 - deltaSP2/2.):i_SP2) << " ";
         pitchFile << (shiftGap   ? (i_G2 + deltaG2)  :i_G2)  << " ";
         pitchFile << (shiftPitch ? (i_SP3 - deltaSP3/2.):i_SP3) << " ";

         // delete f; delete fLeft; delete fMiddle; delete fRight; delete g;

      }

      pitchFile << endl;
      double avgMag = (fabs(deltaU) + fabs(deltaV))/2.;
      double halfAvg = avgMag/2.;

      if(shiftOffset)constantsFile << i_deltaPhiU << " " << i_deltaU - deltaU/fabs(deltaU)*halfAvg << " " << i_deltaPhiV << " " << i_deltaV - deltaV/fabs(deltaV)*halfAvg << endl;
      else constantsFile << i_deltaPhiU << " " << i_deltaU << " " << i_deltaPhiV << " " << i_deltaV << endl;
   }

   // Draw the results
   TCanvas *cResultU = new TCanvas("cResultU","cResultU",1200,900);
   cResultU->Divide(6,4);
   for (unsigned int i=0; i<24; i++){
      cResultU->cd(i+1);
      hists2DU[i]->Draw("colz");
      resultHistsU[i]->Draw("same");
      fLeftU[i]->Draw("same");
      fMiddleU[i]->Draw("same");
      fRightU[i]->Draw("same");
   }

   cResultU->SaveAs("ResultU.png");
   cResultU->SaveAs("ResultU.pdf");
   cResultU->Write();

   TCanvas *cResultV = new TCanvas("cResultV","cResultV",1200,900);
   cResultV->Divide(6,4);
   for (unsigned int i=0; i<24; i++){
      cResultV->cd(i+1);
      hists2DV[i]->Draw("colz");
      resultHistsV[i]->Draw("same");
      fLeftV[i]->Draw("same");
      fMiddleV[i]->Draw("same");
      fRightV[i]->Draw("same");
   }

   cResultV->SaveAs("ResultV.png");
   cResultV->SaveAs("ResultV.pdf");
   cResultV->Write();

   outFile->Write();
   outFile->Close();
   constantsFile.close();

   return;
}

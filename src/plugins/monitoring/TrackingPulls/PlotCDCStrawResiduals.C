{
   TCanvas *cCDCRes = new TCanvas("CDCWireResiduals","CDCWireResiduals",1200,1000);
   cCDCRes->Divide(7,4);

   TH2I *thisCDCRes;
   char histoCDCName[256];
   char histoCDCTitle[256];

   for(auto i=1 ; i <=28; i++){
      cCDCRes->cd(i);
      sprintf(histoCDCName, "/TrackingPulls/CDCPulls_Ring%.2i/Per Straw Residuals", i);
      sprintf(histoCDCTitle, "CDC Residuals Ring %.2i", i);
      gDirectory->GetObject(histoCDCName, thisCDCRes);
      //thisCDCRes->SetMinimum(-0.025);
      //thisCDCRes->SetMaximum(0.025);
      thisCDCRes->SetTitle(histoCDCTitle);
      thisCDCRes->Draw("colz");
   }
}

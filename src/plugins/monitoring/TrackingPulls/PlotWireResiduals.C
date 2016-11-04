{
   TCanvas *cWireRes = new TCanvas("FDCWireResiduals","FDCWireResiduals",1200,1000);
   cWireRes->Divide(6,4);

   TProfile2D *thisWireProfile;
   char histoWName[256];
   char histoWTitle[256];

   for(auto i=1 ; i <=24; i++){
      cWireRes->cd(i);
      sprintf(histoWName, "/TrackingPulls/FDCPulls_Plane%.2i/2D Wire Hit Residuals", i);
      sprintf(histoWTitle, "2D Wire Hit Residuals Plane %.2i", i);
      gDirectory->GetObject(histoWName, thisWireProfile);
      thisWireProfile->SetMinimum(-0.025);
      thisWireProfile->SetMaximum(0.025);
      thisWireProfile->SetTitle(histoWTitle);
      thisWireProfile->Draw("colz");
   }
}

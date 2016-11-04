{
   TCanvas *cCathodeRes = new TCanvas("FDCCathodeResiduals","FDCCathodeResiduals",1200,1000);
   cCathodeRes->Divide(6,4);

   TProfile2D *thisCathodeProfile;
   char histoCName[256];
   char histoCTitle[256];

   for(auto i=1 ; i <=24; i++){
      cCathodeRes->cd(i);
      sprintf(histoCName, "/TrackingPulls/FDCPulls_Plane%.2i/2D Cathode Hit Residuals", i);
      sprintf(histoCTitle, "2D Cathode Hit Residuals Plane %.2i", i);
      gDirectory->GetObject(histoCName, thisCathodeProfile);
      thisCathodeProfile->SetMinimum(-0.025);
      thisCathodeProfile->SetMaximum(0.025);
      thisCathodeProfile->SetTitle(histoCTitle);
      thisCathodeProfile->Draw("colz");
   }
}

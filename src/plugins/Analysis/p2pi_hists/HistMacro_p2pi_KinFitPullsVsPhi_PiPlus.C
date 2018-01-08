{
   TCanvas *cPiPlus = new TCanvas("cPiPlus","cPiPlus", 1200,900);
   cPiPlus->Divide(3,2);

   cPiPlus->cd(1);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton__Proton_Pi+_Pi-/Pi+/Pull_Px_VsPhi");
   h->Draw("colz");

   cPiPlus->cd(2);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton__Proton_Pi+_Pi-/Pi+/Pull_Py_VsPhi");
   h->Draw("colz");

   cPiPlus->cd(3);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton__Proton_Pi+_Pi-/Pi+/Pull_Pz_VsPhi");
   h->Draw("colz");

   cPiPlus->cd(4);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton__Proton_Pi+_Pi-/Pi+/Pull_Xx_VsPhi");
   h->Draw("colz");

   cPiPlus->cd(5);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton__Proton_Pi+_Pi-/Pi+/Pull_Xy_VsPhi");
   h->Draw("colz");

   cPiPlus->cd(6);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton__Proton_Pi+_Pi-/Pi+/Pull_Xz_VsPhi");
   h->Draw("colz");

   cPiPlus->SaveAs("PiPlusKinFitPulls.png");
}



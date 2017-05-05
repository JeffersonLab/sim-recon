{
   TCanvas *cProton = new TCanvas("cProton","cProton", 1200,900);
   cProton->Divide(3,2);

   cProton->cd(1);
   TH2I *h;
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton___Proton_Pi+_Pi-/Proton/Pull_Px_VsPhi");
   h->Draw("colz");

   cProton->cd(2);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton___Proton_Pi+_Pi-/Proton/Pull_Py_VsPhi");
   h->Draw("colz");

   cProton->cd(3);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton___Proton_Pi+_Pi-/Proton/Pull_Pz_VsPhi");
   h->Draw("colz");

   cProton->cd(4);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton___Proton_Pi+_Pi-/Proton/Pull_Xx_VsPhi");
   h->Draw("colz");

   cProton->cd(5);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton___Proton_Pi+_Pi-/Proton/Pull_Xy_VsPhi");
   h->Draw("colz");

   cProton->cd(6);
   h=(TH2I*) gDirectory->Get("/p2pi_preco_kinfit/Hist_KinFitResults/Step0__Photon_Proton___Proton_Pi+_Pi-/Proton/Pull_Xz_VsPhi");
   h->Draw("colz");

   cProton->SaveAs("ProtonKinFitPulls.png");

}



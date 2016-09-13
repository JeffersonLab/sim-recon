{
	TDirectory *mainDir = gDirectory;
	TDirectory *reactionDir;
	TDirectory *customDir;

	if ((TDirectory*)mainDir->FindObjectAny("p2k_preco") != 0)
		reactionDir = (TDirectory*)mainDir->FindObjectAny("p2k_preco");
	else
		return;

	if ((TDirectory*)reactionDir->FindObjectAny("Custom_p2k_hists_NoKinFit_Measured") != 0)
		customDir = (TDirectory*)reactionDir->FindObjectAny("Custom_p2k_hists_NoKinFit_Measured");
	else
		return;

	customDir->cd();

	TH1I* h_InvariantMass = (TH1I*)gROOT->FindObject("InvariantMass");
	TH1I* h_MissingMassSq = (TH1I*)gROOT->FindObject("MissingMassSq");
	TH2I* h_Kp_p_theta_phi = (TH2I*)gROOT->FindObject("Kplus_P_Theta_PhiTag");
	TH2I* h_Km_p_theta_phi = (TH2I*)gROOT->FindObject("Kminus_P_Theta_PhiTag");

	TCanvas *c = NULL;
	if (TVirtualPad::Pad() == NULL)
		c = new TCanvas("p2k_hists","p2k_hists",1200,800);
	else
		c = gPad->GetCanvas();
	c->Divide(2,2);

	c->cd(1);
	if(h_InvariantMass != NULL)
		h_InvariantMass->Draw();

	c->cd(2);
	if(h_MissingMassSq != NULL)
		h_MissingMassSq->Draw();

	c->cd(3);
	if(h_Kp_p_theta_phi != NULL)
		h_Kp_p_theta_phi->Draw("colz");

	c->cd(4);
	if(h_Km_p_theta_phi != NULL)
		h_Km_p_theta_phi->Draw("colz");

}

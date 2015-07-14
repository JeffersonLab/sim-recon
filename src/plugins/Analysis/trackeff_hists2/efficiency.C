

void efficiency(void)
{
	gROOT->Reset();
	gStyle->SetPalette(1,NULL);

	// Open ROOT file and get pointer to tree
	TFile *f = new TFile("hd_root.root");
	f->cd("TRACKING");
	TTree *trkeff = (TTree*)gROOT->FindObject("trkeff");
	
	// It turns out that the resolutions used to determine the pulls
	// could be off which in turn, leads to an inaccurate chisq
	// and an inaccurate efficiency. To avoid this, we first fit the
	// pulls to gaussians in order to determine thier sigmas (which
	// really should be 1). These are then used to scale the event
	// by event pulls such that they are equivalent to sigma=1 and
	// the chisq re-calculated for use in determining the efficiency.
	double sigma_pt_pull = 1.0;
	double sigma_theta_pull = 1.0;
	double sigma_phi_pull = 1.0;
	
	TH1D *h = new TH1D("h","",200, -10.0, 10.0);
	TCut cut1 = "isreconstructable==1";
	
	// sigma pt_pull
	trkeff->Project("h","pt_pull", cut1);
	h->Fit("gaus","0Q");
	sigma_pt_pull = h->GetFunction("gaus")->GetParameter(2);
	
	// sigma sigma_theta_pull
	trkeff->Project("h","theta_pull", cut1);
	h->Fit("gaus","0Q");
	sigma_theta_pull = h->GetFunction("gaus")->GetParameter(2);
	
	// sigma pt_pull
	trkeff->Project("h","phi_pull", cut1);
	h->Fit("gaus","0Q");
	sigma_phi_pull = h->GetFunction("gaus")->GetParameter(2);
	
	cout<<"sigma_pt_pull="<<sigma_pt_pull<<" sigma_theta_pull="<<sigma_theta_pull<<" sigma_phi_pull="<<sigma_phi_pull<<endl;
	
	char chisq_str[256];
	char cut_str[256];
	sprintf(chisq_str, "((pow(pt_pull/%f,2.0) + pow(theta_pull/%f,2.0) + pow(phi_pull/%f,2.0))/3.0)", sigma_pt_pull, sigma_theta_pull, sigma_phi_pull);	
	sprintf(cut_str, "%s<50.0", chisq_str);	
	TCut cut2 = cut_str;

	TH2D *found_vs_pt_vs_theta = new TH2D("found_vs_pt_vs_theta","", 80, 0.0, 160.0, 30, 0.0, 5.0);
	TH2D *thrown_vs_pt_vs_theta = found_vs_pt_vs_theta->Clone("thrown_vs_pt_vs_theta");
	TH2D *eff_vs_pt_vs_theta = found_vs_pt_vs_theta->Clone("eff_vs_pt_vs_theta");
	
	trkeff->Project("found_vs_pt_vs_theta", "E.pthrown.Perp():E.pthrown.Theta()*57.3", cut1 && cut2);
	trkeff->Project("thrown_vs_pt_vs_theta", "E.pthrown.Perp():E.pthrown.Theta()*57.3", cut1);
	eff_vs_pt_vs_theta->Divide(found_vs_pt_vs_theta, thrown_vs_pt_vs_theta);
	
	eff_vs_pt_vs_theta->SetStats(0);
	eff_vs_pt_vs_theta->GetZaxis()->SetRangeUser(0.8, 1.0);
	eff_vs_pt_vs_theta->Draw("colz");
	
	int maxbin = eff_vs_pt_vs_theta->GetMaximumBin();
	cout<<"Maximum efficiency:"<<eff_vs_pt_vs_theta->GetBinContent(maxbin)<<endl;
	
	//trkeff->Draw(chisq_str, cut1&&"chisq<1000000.0");
}


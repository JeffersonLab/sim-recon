// Monitoring macro to check t_PSC - t_RF and the sigma of the distributions

{
	// Define constants
	Int_t NMODULES = 8;

	TDirectory *mainDir = gDirectory;
	TDirectory *pluginDir;

	mainDir->cd();

	TH1F* h_means = new TH1F("means_dt","Mean time difference for each PSC counter;Counter (1-8 left, 9-16 right);#Deltat (PSC - RF) [ns]",16,1.0,17.0);
	TH1F* h_sigmas = new TH1F("sigmas_dt","Sigmas for time difference for each PSC counter;Counter (1-8 left, 9-16);Sigma [ns]",16,1.0,17.0);
	TH2I* h_dt_vs_pp;
	TH1D* h_proj;

	h_means->SetMinimum(-0.3);
	h_means->SetMaximum(0.3);
	h_means->GetYaxis()->SetTitleOffset(1.8);
	h_sigmas->SetMinimum(0.08);
	h_sigmas->SetMaximum(0.22);
	h_sigmas->GetYaxis()->SetTitleOffset(1.8);

	for (Int_t i = 0; i < NMODULES; ++i)
	{
		h_dt_vs_pp = (TH2I*)gROOT->FindObject(Form("h_dt_vs_pp_l_%i",i+1));
		h_proj = h_dt_vs_pp->ProjectionY();
		TFitResultPtr fitptr = h_proj->Fit("gaus","sqn");
		double fitmean = fitptr->Parameters()[1];
		double fitsig = fitptr->Parameters()[2];

		fitptr = h_proj->Fit("gaus","sqn","",fitmean-2*fitsig,fitmean+2*fitsig);
		fitmean = fitptr->Parameters()[1];
		double fitmeanerr = fitptr->ParError(1);
		fitsig = fitptr->Parameters()[2];
		double fitsigerr = fitptr->ParError(2);

		h_means->Fill(i+1,fitmean);
		h_means->SetBinError(i+1,fitmeanerr);
		h_sigmas->Fill(i+1,fitsig);
		h_sigmas->SetBinError(i+1,fitsigerr);
	
		h_dt_vs_pp->Reset();
		h_proj->Reset();
		h_dt_vs_pp = (TH2I*)gROOT->FindObject(Form("h_dt_vs_pp_r_%i",i+1));
		h_proj = h_dt_vs_pp->ProjectionY();
		fitptr = h_proj->Fit("gaus","sqn");
		fitmean = fitptr->Parameters()[1];
		fitsig = fitptr->Parameters()[2];

		fitptr = h_proj->Fit("gaus","sqn","",fitmean-2*fitsig,fitmean+2*fitsig);
		fitmean = fitptr->Parameters()[1];
		fitmeanerr = fitptr->ParError(1);
		fitsig = fitptr->Parameters()[2];
		fitsigerr = fitptr->ParError(2);

		h_means->Fill(i+NMODULES+1,fitmean);
		h_means->SetBinError(i+NMODULES+1,fitmeanerr);
		h_sigmas->Fill(i+NMODULES+1,fitsig);
		h_sigmas->SetBinError(i+NMODULES+1,fitsigerr);
	}

	TCanvas *c = NULL;
	if (TVirtualPad::Pad() == NULL)
		c = new TCanvas("PSC_TW","PSC_TW",1200,600);
	else
		c = gPad->GetCanvas();
	c->Clear();
	c->Divide(2,1);
	gStyle->SetOptStat(0);

	c->cd(1);
	TPad* pad = (TPad*)c->GetPad(1);
	pad->SetLeftMargin(0.15);
	if (h_means != NULL)
		h_means->Draw("E");

	c->cd(2);
	pad = (TPad*)c->GetPad(2);
	pad->SetLeftMargin(0.15);
	if (h_sigmas != NULL)
		h_sigmas->Draw("E");
}

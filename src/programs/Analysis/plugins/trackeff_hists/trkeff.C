

void trkeff(void)
{
	gROOT->Reset();

	// If we are not already in the TRACKING directory, cd there
	TDirectory *TRACKING = (TDirectory*)gROOT->FindObject("TRACKING");
	if(TRACKING)TRACKING->cd();

	// Get pointer to trkeff tree
	TTree *trkeff = (TTree*)gROOT->FindObject("trkeff");

	// Create histograms to calculate efficiency
	TH1D *np = new TH1D("np","Numerator", 50, 0.0, 7.0);
	TH1D *dp = new TH1D("dp","Denominator", 50, 0.0, 7.0);

	trkeff->Project("np", "F.pthrown.Mag()", "F.status==0 && F.nhits_can>0");
	trkeff->Project("dp", "F.pthrown.Mag()", "F.status==0");
	
	TH1D *effp = np->Clone("effp");
	
	TCanvas *c1 = new TCanvas("c1");
	c1->SetTickx();
	c1->SetTicky();
	effp->Divide(dp);
	effp->SetMarkerStyle(2);
	effp->SetStats(0);
	effp->Draw("P");
	
	double eff = np->GetEntries()/dp->GetEntries();
	cout<<" eff = "<<eff<<endl;
}

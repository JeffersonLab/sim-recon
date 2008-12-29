

void tofHits(void)
{
	gROOT->Reset();
	gStyle->SetErrorX(0.0);

	TCanvas *c1 = new TCanvas("c1");
	c1->SetTicky();
	c1->SetTickx();
	c1->SetGridy();
	c1->Draw();
	
	TFile *f = new TFile("tof_hits.root");
	TTree *t = (TTree*)gROOT->FindObject("tofHits");
	double Nevents_thrown = (double)t->GetEntries();
	double timewindow = 1.0E-6;
	double total_integrated_time = timewindow*Nevents_thrown;
	double tagged_rate = 1.0E7;
	cout<<"Number of events thrown: "<<(int)Nevents_thrown<<endl;
	cout<<"Assuming "<<timewindow/1.0E-6<<" microsecond time window"<<endl;
	cout<<"Total beam time represented by simulation:"<<total_integrated_time<<" seconds"<<endl;
	cout<<"Assuming rates correspond to "<<tagged_rate<<" tagged photons/sec beam"<<endl;

	TH1D *hits_per_bar = new TH1D("hits_per_bar","TOF Hits per bar", 45, 0.5, 45.5);
	t->Project("hits_per_bar","bar","plane==1");
	
	TH1D *rate_per_bar = hits_per_bar->Clone("rate_per_bar");
	char str[256];
	sprintf(str, "Rate per bar at %g tags/sec (kHz)", tagged_rate);
	rate_per_bar->Scale(1.0E-3/total_integrated_time);
	ScaleError(rate_per_bar, 1.0E-3/total_integrated_time);
	rate_per_bar->SetXTitle("TOF Bar");
	rate_per_bar->SetYTitle(str);
	rate_per_bar->SetStats(0);
	rate_per_bar->SetLineColor(kRed);
	rate_per_bar->SetLineWidth(2.0);
	rate_per_bar->SetMarkerColor(kRed);
	rate_per_bar->SetMarkerStyle(20);
	rate_per_bar->SetMarkerSize(1.5);
	rate_per_bar->Draw("P");
	
	c1->SaveAs("tof_hit_rate.gif");
	c1->SaveAs("tof_hit_rate.pdf");
}

void ScaleError(TH1D *h, double scale_factor)
{
	cout<<"Scaling errors by "<<scale_factor<<endl;
	int Nbins = h->GetNbinsX();
	for(int i=1; i<=Nbins; i++){
		double val = h->GetBinError(i);
		val *= scale_factor;
		h->SetBinError(i,val);
	}
}


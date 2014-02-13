

void cdcHits(void)
{
	gROOT->Reset();
	gStyle->SetErrorX(0.0);

	TCanvas *c1 = new TCanvas("c1");
	c1->SetTicky();
	c1->SetTickx();
	c1->SetGridy();
	c1->Draw();
	
	TFile *f = new TFile("cdc_hits.root");
	TTree *t = (TTree*)gROOT->FindObject("cdcHits");
	double Nevents_thrown = (double)t->GetEntries();
	double timewindow = 1.0E-6;
	double total_integrated_time = timewindow*Nevents_thrown;
	double tagged_rate = 1.0E7;
	cout<<"Number of events thrown: "<<(int)Nevents_thrown<<endl;
	cout<<"Assuming "<<timewindow/1.0E-6<<" microsecond time window"<<endl;
	cout<<"Total beam time represented by simulation:"<<total_integrated_time<<" seconds"<<endl;
	cout<<"Assuming rates correspond to "<<tagged_rate<<" tagged photons/sec beam"<<endl;
	
	// Number of straws in each "ring"
	int n_straws[]={43,50,57,64,71,78,85,99,106,113,120,127,134,141,148,155,166,173,182,187,194,201,208,215,222};
	
	TH1D *nstraws = new TH1D("nstraws","Straws per ring", 25, 0.5, 25.5);
	for(int i=0;i<25; i++){
		nstraws->SetBinContent(i+1, n_straws[i]);
		nstraws->SetBinError(i+1, 0.0);
	}

	TH1D *hits_per_straw = new TH1D("hits_per_straw","CDC Hits per layer", 25, 0.5, 25.5);
	t->Project("hits_per_straw","ring");
	hits_per_straw->Divide(nstraws);
	
	TH1D *rate_per_straw = hits_per_straw->Clone("rate_per_straw");
	char str[256];
	sprintf(str, "Rate per straw at %g tags/sec (kHz)", tagged_rate);
	rate_per_straw->Scale(1.0E-3/total_integrated_time);
	rate_per_straw->SetXTitle("CDC layer");
	rate_per_straw->SetYTitle(str);
	rate_per_straw->SetStats(0);
	rate_per_straw->SetLineColor(kRed);
	rate_per_straw->SetLineWidth(2.0);
	rate_per_straw->SetMarkerColor(kRed);
	rate_per_straw->SetMarkerStyle(20);
	rate_per_straw->SetMarkerSize(1.5);
	rate_per_straw->Draw("P");
	
	c1->SaveAs("cdc_hit_rate.gif");
	c1->SaveAs("cdc_hit_rate.pdf");
}


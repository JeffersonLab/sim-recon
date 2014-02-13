

void fdcHits(void)
{
	gROOT->Reset();
	gStyle->SetErrorX(0.0);

	TCanvas *c1 = new TCanvas("c1");
	c1->SetTicky();
	c1->SetTickx();
	c1->SetGridy();
	c1->Draw();
	
	TFile *f = new TFile("fdc_hits.root");
	TTree *t = (TTree*)gROOT->FindObject("fdcHits");
	double Nevents_thrown = (double)t->GetEntries();
	double timewindow = 1.0E-6;
	double total_integrated_time = timewindow*Nevents_thrown;
	double tagged_rate = 1.0E7;
	cout<<"Number of events thrown: "<<(int)Nevents_thrown<<endl;
	cout<<"Assuming "<<timewindow/1.0E-6<<" microsecond time window"<<endl;
	cout<<"Total beam time represented by simulation:"<<total_integrated_time<<" seconds"<<endl;
	cout<<"Assuming rates correspond to "<<tagged_rate<<" tagged photons/sec beam"<<endl;

	// Hits per anode wire in last layer
	TH1D *hits_per_wire24 = new TH1D("hits_per_wire24","FDC Hits per anode wire", 96, 0.5, 96.5);
	t->Project("hits_per_wire24","element","layer==24 && ((plane-1)%3)==1");

	// Hits per anode wire in first layer
	TH1D *hits_per_wire1 = new TH1D("hits_per_wire1","FDC Hits per anode wire", 96, 0.5, 96.5);
	t->Project("hits_per_wire1","element","layer==1 && ((plane-1)%3)==1");
	
	// Anode Wire hit plot
	char str[256];
	sprintf(str, "Rate per wire at %g tags/sec (kHz)", tagged_rate);
	hits_per_wire24->Scale(1.0E-3/total_integrated_time);
	hits_per_wire24->SetXTitle("FDC Wire Number");
	hits_per_wire24->SetYTitle(str);
	hits_per_wire24->SetStats(0);
	hits_per_wire24->SetLineColor(kRed);
	hits_per_wire24->SetLineWidth(2.0);
	hits_per_wire24->Draw();
	
	hits_per_wire1->Scale(1.0E-3/total_integrated_time);
	hits_per_wire1->SetLineColor(kBlue);
	hits_per_wire1->SetLineWidth(2.0);
	hits_per_wire1->Draw("same");
	
	TLegend *leg = new TLegend(0.6, 0.5, 0.85, 0.75);
	leg->AddEntry(hits_per_wire1, "Upstream");
	leg->AddEntry(hits_per_wire24, "Downstream");
	leg->Draw();

	// Hits per cathode strip in last layer
	TH1D *hits_per_strip24 = new TH1D("hits_per_strip24","FDC Hits per cathode strip", 220, 0.5, 220.5);
	t->Project("hits_per_strip24","element","layer==24 && ((plane-1)%3)==0");

	// Hits per cathode strip in first layer
	TH1D *hits_per_strip1 = new TH1D("hits_per_strip1","FDC Hits per cathode strip", 220, 0.5, 220.5);
	t->Project("hits_per_strip1","element","layer==1 && ((plane-1)%3)==0");
	
	c1->SaveAs("fdc_anode_hit_rate.gif");
	c1->SaveAs("fdc_anode_hit_rate.pdf");
	
	// Cathode hit plot
	char str[256];
	sprintf(str, "Rate per strip at %g tags/sec (kHz)", tagged_rate);
	hits_per_strip24->Scale(1.0E-3/total_integrated_time);
	hits_per_strip24->SetXTitle("FDC Strip Number");
	hits_per_strip24->SetYTitle(str);
	hits_per_strip24->SetStats(0);
	hits_per_strip24->SetLineColor(kRed);
	hits_per_strip24->SetLineWidth(2.0);
	hits_per_strip24->Draw();
	
	hits_per_strip1->Scale(1.0E-3/total_integrated_time);
	hits_per_strip1->SetLineColor(kBlue);
	hits_per_strip1->SetLineWidth(2.0);
	hits_per_strip1->Draw("same");
	
	TLegend *leg = new TLegend(0.6, 0.5, 0.85, 0.75);
	leg->AddEntry(hits_per_strip1, "Upstream");
	leg->AddEntry(hits_per_strip24, "Downstream");
	leg->Draw();

	c1->SaveAs("fdc_cathode_hit_rate.gif");
	c1->SaveAs("fdc_cathode_hit_rate.pdf");
}


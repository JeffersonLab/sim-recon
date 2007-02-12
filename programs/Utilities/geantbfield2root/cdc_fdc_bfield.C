


void cdc_fdc_bfield(void)
{
	gROOT->Reset();
	
	double cdc_start = 17.0;
	double cdc_length = 175.0;
	double fdc_start = 214.9;
	double fdc_length = 16.2;
	double fdc_space = 43.8;
	
	TCanvas *c1 = new TCanvas("c1");
	c1->SetTickx();
	c1->SetTicky();
	c1->SetGridy();
	
	TH2F *axes = new TH2F("axes","GlueX Magnetic Field Map (TOSCA)", 100, -100.0, 500.0, 100, -2.5, 0.0);
	axes->SetStats(0);
	axes->SetYTitle("B_{z} (Tesla)");
	axes->SetXTitle("z-distance along beamline (cm)");
	axes->Draw();

	TFile *f = new TFile("bfield.root");
	TTree *bfield = (TTree*)gROOT->FindObject("bfield");

	bfield->SetMarkerStyle(8);
	bfield->SetMarkerSize(0.4);
	bfield->Draw("Bz:z", "r<65.0 && z<450.0","same");
	
	// CDC
	bfield->SetMarkerColor(kRed);
	char cut[256];
	sprintf(cut, "r<65.0 && z>%f && z<%f", cdc_start, cdc_start+cdc_length);
	bfield->Draw("Bz:z", cut, "same");
	TLatex *l = new TLatex(100.0, -1.0, "CDC");
	l->Draw();
	TArrow *a = new TArrow(cdc_start, -1.2, cdc_start+cdc_length, -1.2, 0.02, "|<|-|>|");
	a->SetLineColor(kRed);
	a->SetFillColor(kRed);
	a->Draw();
 
	// FDC
	bfield->SetMarkerColor(kBlue);
	TLatex *l = new TLatex(300.0, -0.3, "FDC");
	l->Draw();
	
	double fdc_ypos[]={-1.9, -1.9, -1.5, -0.5};

	for(int i=0; i<4; i++){
		double start = fdc_start + i*(fdc_length+fdc_space);
		sprintf(cut, "r<65.0 && z>%f && z<%f", start, start+fdc_length);
		bfield->Draw("Bz:z", cut, "same");

		TArrow *a = new TArrow(320.0, -0.35, start+fdc_length/2.0, fdc_ypos[i], 0.02, "-|>|");
		a->SetLineColor(kBlue);
		a->SetFillColor(kBlue);
		a->Draw();
	}
	
	// Target
	TLatex *l = new TLatex(45.0, -2.25, "Target");
	l->SetTextSize(0.03);
	l->Draw();
	TArrow *a = new TArrow(50.0, -2.15, 80.0, -2.15, 0.01, "|<|-|>|");
	a->SetLineColor(kMagenta);
	a->SetFillColor(kMagenta);
	a->Draw();
	
	// Save
	c1->SaveAs("cdc_fdc_bfield.eps");
	c1->SaveAs("cdc_fdc_bfield.gif");
}


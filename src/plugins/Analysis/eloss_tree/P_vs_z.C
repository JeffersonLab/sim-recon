

void P_vs_z(void)
{
	gROOT->Reset();
	
	TFile *f = new TFile("hd_root.root");
	TTree *geant = (TTree*)gROOT->FindObject("geant");
	TTree *dana = (TTree*)gROOT->FindObject("dana");

	TCanvas *c1 = new TCanvas();
	c1->SetTicks();
	c1->SetGrid();

	geant->SetMarkerStyle(6);
	geant->SetMarkerSize(0.1);

	geant->Draw("P:z", "P>0.95");
	dana->Draw("P:z", "", "same");
	
	TLatex *lab = new TLatex(400.0, 1.005, "single 1GeV/c proton #theta=5^{o} #phi=180^{o}");
	lab->SetTextSize(0.03);
	lab->Draw();
	
	geant->SetLineColor(geant->GetMarkerColor());
	dana->SetLineColor(dana->GetMarkerColor());
	geant->SetLineWidth(4.0);
	dana->SetLineWidth(4.0);
	
	TLegend *leg = new TLegend(0.5, 0.5, 0.7, 0.7);
	leg->SetFillColor(kWhite);
	leg->AddEntry(geant, "GEANT");
	leg->AddEntry(dana, "DANA");
	leg->Draw();

	c1->SaveAs("P_vs_z.pdf");
	c1->SaveAs("P_vs_z.gif");
}


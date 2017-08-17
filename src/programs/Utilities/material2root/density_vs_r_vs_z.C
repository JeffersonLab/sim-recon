

void density_vs_r_vs_z(void)
{
	gROOT->Reset();

	TFile *f = new TFile("material.root");
	TH2D *density_LL = (TH2D*)gROOT->FindObject("density_LL");
	TH2D *density_table = (TH2D*)gROOT->FindObject("density_table");

	TCanvas *c1 = new TCanvas("c1", "", 800, 1000);
	c1->Divide(1,2);
	
	c1->cd(1);
	gPad->SetLogz();
	gPad->SetRightMargin(0.15);
	density_LL->GetZaxis()->SetRangeUser(1.0E-3, 1.25);
	density_LL->SetStats(0);
	density_LL->SetTitle("Material Map from FindMatLL(...)");
	density_LL->GetZaxis()->SetTitle("density (g/cm^{3})");
	density_table->GetZaxis()->SetRangeUser(0.001, 1.5);
	density_LL->Draw("colz");

	c1->cd(2);
	gPad->SetLogz();
	gPad->SetRightMargin(0.15);
	density_table->GetZaxis()->SetRangeUser(1.0E-3, 1.25);
	density_table->SetStats(0);
	density_table->SetTitle("Material Map from FindMatTable(...)");
	density_table->GetZaxis()->SetTitle("density (g/cm^{3})");
	density_table->GetZaxis()->SetRangeUser(0.001, 1.5);
	density_table->Draw("colz");
	
	c1->SaveAs("density_vs_r_vs_z.png");
	c1->SaveAs("density_vs_r_vs_z.pdf");
}



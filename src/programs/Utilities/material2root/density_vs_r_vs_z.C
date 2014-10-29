

void density_vs_r_vs_z(void)
{
	gROOT->Reset();

	TFile *f = new TFile("material.root");
	TH2D *density_LL = (TH2D*)gROOT->FindObject("density_LL");
	TH2D *density_table = (TH2D*)gROOT->FindObject("density_table");

	TCanvas *c1 = new TCanvas("c1", "", 800, 1000);
	c1->Divide(1,2);
	
	c1->cd(1);
	density_LL->SetStats(0);
	density_LL->SetTitle("Material Map from FindMatLL(...)");
	density_LL->Draw("colz");

	c1->cd(2);
	density_table->SetStats(0);
	density_table->SetTitle("Material Map from FindMatTable(...)");
	density_table->Draw("colz");
	
	c1->SaveAs("density_vs_r_vs_z.gif");
	c1->SaveAs("density_vs_r_vs_z.pdf");
}



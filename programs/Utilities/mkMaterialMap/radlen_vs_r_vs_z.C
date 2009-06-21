

void radlen_vs_r_vs_z(void)
{
	gROOT->Reset();

	TFile *f = new TFile("material.root");
	TH2D *material_LL = (TH2D*)gROOT->FindObject("material_LL");
	TH2D *material_table = (TH2D*)gROOT->FindObject("material_table");

	TCanvas *c1 = new TCanvas("c1", "", 800, 1000);
	c1->Divide(1,2);
	
	c1->cd(1);
	material_LL->SetStats(0);
	material_LL->SetTitle("Material Map from FindMatLL(...)");
	material_LL->Draw("colz");

	c1->cd(2);
	material_table->SetStats(0);
	material_table->SetTitle("Material Map from FindMatTable(...)");
	material_table->Draw("colz");
	
	c1->SaveAs("radlen_vs_r_vs_z.gif");
	c1->SaveAs("radlen_vs_r_vs_z.pdf");
}



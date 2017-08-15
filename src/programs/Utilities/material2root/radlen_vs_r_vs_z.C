

void radlen_vs_r_vs_z(void)
{
	gROOT->Reset();

	TFile *f = new TFile("material.root");
	TH2D *radlen_LL = (TH2D*)gROOT->FindObject("radlen_LL");
	TH2D *radlen_table = (TH2D*)gROOT->FindObject("radlen_table");

	TCanvas *c1 = new TCanvas("c1", "", 800, 1000);
	c1->Divide(1,2);
	
	c1->cd(1);
	gPad->SetRightMargin(0.15);
	radlen_LL->GetZaxis()->SetRangeUser(0.0, 35000.0);
	radlen_LL->SetStats(0);
	radlen_LL->SetTitle("Material Map from DRootGeom::FindMatLL(...)");
	radlen_LL->Draw("colz");

	c1->cd(2);
	gPad->SetRightMargin(0.15);
	radlen_table->GetZaxis()->SetRangeUser(0.0, 35000.0);
	radlen_table->SetStats(0);
	radlen_table->SetTitle("Material Map from DGeometry::FindMat(...)");
	radlen_table->Draw("colz");
	
	c1->SaveAs("radlen_vs_r_vs_z.png");
	c1->SaveAs("radlen_vs_r_vs_z.pdf");
}



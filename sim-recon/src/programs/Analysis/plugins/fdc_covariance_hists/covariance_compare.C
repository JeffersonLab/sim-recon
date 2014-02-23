

void covariance_compare(void)
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetPadRightMargin(0.15);

	TFile *f = new TFile("hd_root.root");
	f->cd("TRACKING");
	TProfile2D *cdc_cov = (TProfile2D*)gROOT->FindObject("cdc_cov");
	TProfile2D *cdc_cov_calc = (TProfile2D*)gROOT->FindObject("cdc_cov_calc");
	
	TCanvas *c1 = new TCanvas("c1","",1300,600);
	c1->Divide(2,1);
	c1->SetLogz();
	
	c1->cd(1);
	c1->GetPad(1)->SetLogz();
	cdc_cov->GetZaxis()->SetRangeUser(0.00001, 0.115);
	cdc_cov->Draw("colz");
	
	
	c1->cd(2);
	c1->GetPad(2)->SetLogz();
	cdc_cov_calc->GetZaxis()->SetRangeUser(0.00001, 0.115);
	cdc_cov_calc->Draw("colz");
	
	c1->SaveAs("covariance_compare.gif");
	c1->SaveAs("covariance_compare.pdf");
}


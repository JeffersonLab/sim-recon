

void covariance_compare(void)
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetPadRightMargin(0.15);

	TFile *f = new TFile("hd_root.root");
	f->cd("TRACKING");
	TProfile2D *fdc_cov = (TProfile2D*)gROOT->FindObject("fdc_cov");
	TProfile2D *fdc_cov_calc = (TProfile2D*)gROOT->FindObject("fdc_cov_calc");
	
	TCanvas *c1 = new TCanvas("c1","",1300,600);
	c1->Divide(2,1);
	c1->SetLogz();
	
	c1->cd(1);
	c1->GetPad(1)->SetLogz();
	fdc_cov->GetZaxis()->SetRangeUser(0.000001, 0.115);
	fdc_cov->Draw("colz");
	
	
	c1->cd(2);
	c1->GetPad(2)->SetLogz();
	fdc_cov_calc->GetZaxis()->SetRangeUser(0.000001, 0.115);
	fdc_cov_calc->Draw("colz");
	
	c1->SaveAs("covariance_compare.gif");
	c1->SaveAs("covariance_compare.pdf");
}


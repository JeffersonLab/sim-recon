
void plotTime(string type = "pi", int p=4, int theta=4, int phi=40){

	gStyle->SetOptStat(0);

	gSystem->Unlink(Form("%s_p%d_theta%d_phi%d.gif",type.data(),p,theta,phi)); // delete old file

	//TFile *f = TFile::Open(Form("/volatile/halld/home/jrsteven/2015-3-dirc/data_single/%s_p%d_theta%d_phi%d/single_hd_root_1.root",type.data(),p,theta,phi));
	TFile *f = TFile::Open("hd_root.root");

	TH1 *hTruthHitT = (TH1*)f->Get("hTruthHitT"); 
	hTruthHitT->SetTitle("TruthHit Time; Time (ns)");

	TH3 *hTruthHit3D = (TH3*)f->Get("hTruthHitYLocWT");
	hTruthHit3D->SetTitle("TruthHit Position Y vs LocZ; Y (cm); Local Z (cm)");
	hTruthHit3D->Rebin3D(1, 1, 2);
	TH3 *hTruthHit3D_clone = (TH3*)hTruthHit3D->Clone();
	hTruthHit3D_clone->SetName("hTruthHit3D_clone");

	TCanvas *c1 = new TCanvas("cc","cc",900,1200);
	double ybound = 0.3;
	TPad *p1 = new TPad("p1", "p1", 0., ybound, 0.5, 1.); p1->Draw();
	TPad *p2 = new TPad("p1", "p1", 0.5, ybound, 1., 1.); p2->Draw();
	TPad *p3 = new TPad("p2", "p2", 0., 0., 1., ybound); p3->Draw();

	TLine *l1;

	int nTimeBins = hTruthHit3D->GetZaxis()->GetNbins();
	for(int ibin = 0; ibin<nTimeBins; ibin++){
		p1->cd();
		hTruthHit3D->GetZaxis()->SetRange(ibin+1,ibin+1);
		hTruthHit3D->Project3D("xy")->Draw("colz");
		
		p2->cd();
		hTruthHit3D_clone->GetZaxis()->SetRange(1,ibin+1);
		hTruthHit3D_clone->Project3D("xy")->Draw("colz");

		p3->cd();
		hTruthHitT->Draw();
		l1 = new TLine(hTruthHit3D->GetZaxis()->GetBinCenter(ibin+1),0,hTruthHit3D->GetZaxis()->GetBinCenter(ibin+1),hTruthHitT->GetMaximum());
		l1->Draw("same");

		c1->Print(Form("%s_p%d_theta%d_phi%d.gif+40",type.data(),p,theta,phi));
		c1->Print(Form("%s_p%d_theta%d_phi%d_%04d.ps",type.data(),p,theta,phi,ibin));
		if(ibin%10 == 0) 
			cout<<"printed time bin = "<<ibin<<" of "<<nTimeBins<<endl;
	}
	c1->Print(Form("%s_p%d_theta%d_phi%d.gif++",type.data(),p,theta,phi));

	gSystem->Exec(Form("cat %s_p%d_theta%d_phi%d*.ps | ps2pdf - %s_p%d_theta%d_phi%d.pdf",type.data(),p,theta,phi,type.data(),p,theta,phi));
	gSystem->Exec("rm *.ps");

	return;
}

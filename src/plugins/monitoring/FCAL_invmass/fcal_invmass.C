// hnamepath:  /FCAL_invmass/InvMass1
// hnamepath:  /FCAL_invmass/InvMass2
// hnamepath:  /FCAL_invmass/InvMass3
// hnamepath:  /FCAL_invmass/InvMass4
// hnamepath:  /FCAL_invmass/InvMass5
// hnamepath:  /FCAL_invmass/InvMass6
// hnamepath:  /FCAL_invmass/InvMass7
// hnamepath:  /FCAL_invmass/InvMass8
// hnamepath:  /FCAL_invmass/InvMass9


void fitHisto(TH1I* histo){

	float par_1[12]; 
	Double_t par[12] = {};  
	TF1 *g1    = new TF1("g1","gaus",0.06,0.09);
	TF1 *g2    = new TF1("g2","gaus",0.09,0.15);
	TF1 *g3    = new TF1("g3","pol2",0.15,0.28);

	TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol3(6)",0.06,0.2);
	total->SetLineColor(kRed);
	histo->Fit(g1,"nR");
	histo->Fit(g2,"nR+");
	histo->Fit(g3,"nR+");
	g1->GetParameters(&par[0]);
	g2->GetParameters(&par[3]);
	g3->GetParameters(&par[6]);
	total->SetParameters(par);

	total->SetParName(3, "Amplitude");
	total->SetParName(4, "Pi0 Mass");
	total->SetParName(5, "Pi0 Width");

	histo->Fit(total,"R+"); 

	for (int i=0; i<12; i++) {
		par_1[i] = total->GetParameter(i);
	}

	TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
	pt_300->SetFillColor(0);
	pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_1[4]*1000));
	pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_1[5]*1000));
	pt_300->AddText(Form("#sigma/M = %.3f %%",(par_1[5]/par_1[4])*100));
	pt_300->Draw();

}

void fcal_invmass(){


	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("FCAL_invmass");
	if(dir) dir->cd();

	TH1I* InvMass1 = (TH1I*)gDirectory->FindObjectAny("InvMass1");
	TH1I* InvMass2 = (TH1I*)gDirectory->FindObjectAny("InvMass2");
	TH1I* InvMass3 = (TH1I*)gDirectory->FindObjectAny("InvMass3");
	TH1I* InvMass4 = (TH1I*)gDirectory->FindObjectAny("InvMass4");
	TH1I* InvMass5 = (TH1I*)gDirectory->FindObjectAny("InvMass5");
	TH1I* InvMass6 = (TH1I*)gDirectory->FindObjectAny("InvMass6");
	TH1I* InvMass7 = (TH1I*)gDirectory->FindObjectAny("InvMass7");
	TH1I* InvMass8 = (TH1I*)gDirectory->FindObjectAny("InvMass8");
	TH1I* InvMass9 = (TH1I*)gDirectory->FindObjectAny("InvMass9");
	TH1I* qualCut_00 = (TH1I*)gDirectory->FindObjectAny("qualCut_00");
	TH1I* qualCut_03 = (TH1I*)gDirectory->FindObjectAny("qualCut_03");
	TH1I* qualCut_05 = (TH1I*)gDirectory->FindObjectAny("qualCut_05");


	if(gPad == NULL){

		TCanvas *c1 = new TCanvas( "c1", "FCAL_invmass_plot", 1200, 1200 );
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}

	if( !gPad ) return;
	TCanvas* c1 = gPad->GetCanvas();
	c1->Divide(4,3);

	if( qualCut_00 ){
		
		qualCut_00->SetStats(0);
		c1->cd(10);
		qualCut_00->Draw();
		qualCut_00->SetLineWidth(2);
		fitHisto(qualCut_00);
	}
	if( qualCut_03 ){
	
		qualCut_03->SetStats(0);
		c1->cd(11);
		qualCut_03->Draw();
		qualCut_03->SetLineWidth(2);
		fitHisto(qualCut_03);
	}


	if( qualCut_05 ){
	
		qualCut_05->SetStats(0);
		c1->cd(12);
		qualCut_05->Draw();
		qualCut_05->SetLineWidth(2);
		fitHisto(qualCut_05);
	}
	if( InvMass1 ){

		InvMass1->SetStats(0);
		c1->cd(1);
		InvMass1->Draw();
		InvMass1->SetLineWidth(2);
		fitHisto(InvMass1);
	}

	if( InvMass2 ){

		InvMass2->SetStats(0);
		c1->cd(2);
		InvMass2->Draw();
		InvMass2->SetLineWidth(2);
		fitHisto(InvMass2);
	}

	if( InvMass3 ){

		InvMass3->SetStats(0);
		c1->cd(3);
		InvMass3->Draw();
		InvMass3->SetLineWidth(2);
		fitHisto(InvMass3);
	}

	if( InvMass4 ){

		InvMass4->SetStats(0);
		c1->cd(4);
		InvMass4->Draw();
		InvMass4->SetLineWidth(2);
		fitHisto(InvMass4);
	}

	if( InvMass5 ){

		InvMass5->SetStats(0);
		c1->cd(5);
		InvMass5->Draw();
		InvMass5->SetLineWidth(2);
		fitHisto(InvMass5);
	}

	if( InvMass6 ){

		InvMass6->SetStats(0);
		c1->cd(6);
		InvMass6->Draw();
		InvMass6->SetLineWidth(2);
		fitHisto(InvMass6);
	}

	if( InvMass7 ){

		InvMass7->SetStats(0);
		c1->cd(7);
		InvMass7->Draw();
		InvMass7->SetLineWidth(2);
		fitHisto(InvMass7);
	}

	if( InvMass8 ){

		InvMass8->SetStats(0);
		c1->cd(8);
		InvMass8->Draw();
		InvMass8->SetLineWidth(2);
		fitHisto(InvMass8);
	} 

	if( InvMass9 ){

		InvMass9->SetStats(0);
		c1->cd(9);
		InvMass9->Draw();
		InvMass9->SetLineWidth(2);
		fitHisto(InvMass9);
	}
}

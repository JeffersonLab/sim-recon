// hnamepath: /highlevel/NumHighLevelObjects
// hnamepath: /highlevel/F1TDC_fADC_tdiff
//
// e-mail: davidl@jlab.org
// e-mail: pmatt@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
//

{
	TDirectory* locCurrentDir = gDirectory;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH2* NumHighLevelObjects = (TH2*)gDirectory->Get("NumHighLevelObjects");
	TH2* F1TDC_fADC_tdiff = (TH2*)gDirectory->Get("F1TDC_fADC_tdiff");

	//Get/Make Canvas
	TCanvas *c1 = NULL;
	if(TVirtualPad::Pad() == NULL)
		c1 = new TCanvas("NumHighLevelObjects", "NumHighLevelObjects", 1200, 900); //for testing
	else
		c1 = gPad->GetCanvas();

	c1->Divide(1,2);
	//c1->Draw();
	
	if(NumHighLevelObjects){
		c1->cd(1);
		gPad->SetTicks();
		gPad->SetGrid();
		gPad->SetLogz();
		gPad->SetBottomMargin(0.2);
		NumHighLevelObjects->GetXaxis()->SetLabelSize(0.05);
		NumHighLevelObjects->SetStats(0);
		NumHighLevelObjects->Draw("COLZ");
		gPad->Update();
	}

}


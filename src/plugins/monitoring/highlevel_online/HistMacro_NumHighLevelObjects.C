// hnamepath: /highlevel/NumHighLevelObjects
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
	TH2* locHist_NumHighLevel = (TH2*)gDirectory->Get("NumHighLevelObjects");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("NumHighLevelObjects", "NumHighLevelObjects", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumHighLevel != NULL)
		locHist_NumHighLevel->Draw("COLZ");
	gPad->SetLogz();
	gPad->Update();

	locCurrentDir->cd();
}


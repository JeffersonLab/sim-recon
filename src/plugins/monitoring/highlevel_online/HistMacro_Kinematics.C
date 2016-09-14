// hnamepath: /highlevel/PVsTheta_Tracks
// hnamepath: /highlevel/PhiVsTheta_Tracks

{
	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	TH2* locHist_PVsTheta_Tracks = (TH2*)gDirectory->Get("PVsTheta_Tracks");
	TH2* locHist_PhiVsTheta_Tracks = (TH2*)gDirectory->Get("PhiVsTheta_Tracks");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Kinematics", "Kinematics", 1200, 600); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 1);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_Tracks != NULL)
	{
		locHist_PVsTheta_Tracks->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Tracks->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Tracks->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Tracks->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Tracks->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PhiVsTheta_Tracks != NULL)
	{
		locHist_PhiVsTheta_Tracks->GetXaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_Tracks->GetYaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_Tracks->GetXaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_Tracks->GetYaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_Tracks->Draw("colz");
	}
}

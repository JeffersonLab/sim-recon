// hnamepath: /Independent/Hist_DetectedParticleKinematics/Beam_Photon/Momentum
// hnamepath: /Independent/Hist_EventVertexZ/2+GoodTracks/EventVertexZ
// hnamepath: /Independent/Hist_EventVertexZ/2+GoodTracks/EventVertexYVsX
// hnamepath: /Independent/Hist_EventVertexZ/2+GoodTracks/ConfidenceLevel

{
	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectedParticleKinematics");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Beam Histograms
	gDirectory->cd("Beam_Photon");
	TH1I* locHist_BeamMomentum = (TH1I*)gDirectory->Get("Momentum");

	//Goto Vertex Path
	locDirectory = (TDirectory*)locTopDirectory->FindObjectAny("Hist_EventVertex");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Beam Histograms
	gDirectory->cd("2+GoodTracks");
	TH1I* locHist_EventVertexZ = (TH1I*)gDirectory->Get("EventVertexZ");
	TH2I* locHist_EventVertexYVsX = (TH2I*)gDirectory->Get("EventVertexYVsX");
	TH1I* locHist_ConfidenceLevel = (TH1I*)gDirectory->Get("ConfidenceLevel");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("EventInfo", "EventInfo", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BeamMomentum != NULL)
	{
		locHist_BeamMomentum->GetXaxis()->SetTitleSize(0.05);
		locHist_BeamMomentum->GetYaxis()->SetTitleSize(0.05);
		locHist_BeamMomentum->GetXaxis()->SetLabelSize(0.05);
		locHist_BeamMomentum->GetYaxis()->SetLabelSize(0.05);
		locHist_BeamMomentum->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_ConfidenceLevel != NULL)
	{
		locHist_ConfidenceLevel->GetXaxis()->SetTitleSize(0.05);
		locHist_ConfidenceLevel->GetYaxis()->SetTitleSize(0.05);
		locHist_ConfidenceLevel->GetXaxis()->SetLabelSize(0.05);
		locHist_ConfidenceLevel->GetYaxis()->SetLabelSize(0.05);
		locHist_ConfidenceLevel->Draw();
		gPad->SetLogy();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_EventVertexZ != NULL)
	{
		locHist_EventVertexZ->GetXaxis()->SetTitleSize(0.05);
		locHist_EventVertexZ->GetYaxis()->SetTitleSize(0.05);
		locHist_EventVertexZ->GetXaxis()->SetLabelSize(0.05);
		locHist_EventVertexZ->GetYaxis()->SetLabelSize(0.05);
		locHist_EventVertexZ->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_EventVertexYVsX != NULL)
	{
		locHist_EventVertexYVsX->GetXaxis()->SetRangeUser(-4.0, 4.0);
		locHist_EventVertexYVsX->GetYaxis()->SetRangeUser(-4.0, 4.0);
		locHist_EventVertexYVsX->GetXaxis()->SetTitleSize(0.05);
		locHist_EventVertexYVsX->GetYaxis()->SetTitleSize(0.05);
		locHist_EventVertexYVsX->GetXaxis()->SetLabelSize(0.05);
		locHist_EventVertexYVsX->GetYaxis()->SetLabelSize(0.05);
		locHist_EventVertexYVsX->Draw("COLZ");
	}
}


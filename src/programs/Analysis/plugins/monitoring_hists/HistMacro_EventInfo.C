// hnamepath: /Independent/Hist_DetectedParticleKinematics/Beam_Photon/Momentum
// hnamepath: /Independent/Hist_DetectedParticleKinematics/EventVertexZ
// hnamepath: /Independent/Hist_DetectedParticleKinematics/EventVertexT
// hnamepath: /Independent/Hist_DetectedParticleKinematics/EventVertexYVsX

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectedParticleKinematics");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Beam_Photon");
	TH1I* locHist_BeamMomentum = (TH1I*)gDirectory->Get("Momentum");
	gDirectory->cd("..");
	TH1I* locHist_EventVertexZ = (TH1I*)gDirectory->Get("EventVertexZ");
	TH1I* locHist_EventVertexT = (TH1I*)gDirectory->Get("EventVertexT");
	TH2I* locHist_EventVertexYVsX = (TH2I*)gDirectory->Get("EventVertexYVsX");

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
	if(locHist_EventVertexZ != NULL)
	{
		locHist_EventVertexZ->GetXaxis()->SetTitleSize(0.05);
		locHist_EventVertexZ->GetYaxis()->SetTitleSize(0.05);
		locHist_EventVertexZ->GetXaxis()->SetLabelSize(0.05);
		locHist_EventVertexZ->GetYaxis()->SetLabelSize(0.05);
		locHist_EventVertexZ->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_EventVertexT != NULL)
	{
		locHist_EventVertexT->GetXaxis()->SetTitleSize(0.05);
		locHist_EventVertexT->GetYaxis()->SetTitleSize(0.05);
		locHist_EventVertexT->GetXaxis()->SetLabelSize(0.05);
		locHist_EventVertexT->GetYaxis()->SetLabelSize(0.05);
		locHist_EventVertexT->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_EventVertexYVsX != NULL)
	{
		locHist_EventVertexYVsX->GetXaxis()->SetTitleSize(0.05);
		locHist_EventVertexYVsX->GetYaxis()->SetTitleSize(0.05);
		locHist_EventVertexYVsX->GetXaxis()->SetLabelSize(0.05);
		locHist_EventVertexYVsX->GetYaxis()->SetLabelSize(0.05);
		locHist_EventVertexYVsX->Draw("COLZ");
	}
}


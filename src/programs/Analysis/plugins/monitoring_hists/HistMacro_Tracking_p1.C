// hnamepath: /Independent/Hist_Reconstruction/Tracking/NumDCHitsPerTrackVsTheta
// hnamepath: /Independent/Hist_Reconstruction/Tracking/TrackingFOM
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTrackCandidates
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumWireBasedTracks
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTimeBasedTracks

{
	//Goto Path
	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_Reconstruction");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Tracking");
	TH2I* locHist_NumDCHitsPerTrackVsTheta = (TH2I*)gDirectory->Get("NumDCHitsPerTrackVsTheta");
	TH1I* locHist_TrackingFOM = (TH1I*)gDirectory->Get("TrackingFOM");

	locDirectory = (TDirectory*)locInitDirectory->FindObjectAny("Hist_NumReconstructedObjects");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH1I* locHist_NumTrackCandidates = (TH1I*)gDirectory->Get("NumTrackCandidates");
	TH1I* locHist_NumWireBasedTracks = (TH1I*)gDirectory->Get("NumWireBasedTracks");
	TH1I* locHist_NumTimeBasedTracks = (TH1I*)gDirectory->Get("NumTimeBasedTracks");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Tracking_p1", "Tracking_p1", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumTrackCandidates != NULL)
	{
		locHist_NumTrackCandidates->GetXaxis()->SetTitleSize(0.05);
		locHist_NumTrackCandidates->GetYaxis()->SetTitleSize(0.05);
		locHist_NumTrackCandidates->GetXaxis()->SetLabelSize(0.05);
		locHist_NumTrackCandidates->GetYaxis()->SetLabelSize(0.05);
		locHist_NumTrackCandidates->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumWireBasedTracks != NULL)
	{
		locHist_NumWireBasedTracks->GetXaxis()->SetTitleSize(0.05);
		locHist_NumWireBasedTracks->GetYaxis()->SetTitleSize(0.05);
		locHist_NumWireBasedTracks->GetXaxis()->SetLabelSize(0.05);
		locHist_NumWireBasedTracks->GetYaxis()->SetLabelSize(0.05);
		locHist_NumWireBasedTracks->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumTimeBasedTracks != NULL)
	{
		locHist_NumTimeBasedTracks->GetXaxis()->SetTitleSize(0.05);
		locHist_NumTimeBasedTracks->GetYaxis()->SetTitleSize(0.05);
		locHist_NumTimeBasedTracks->GetXaxis()->SetLabelSize(0.05);
		locHist_NumTimeBasedTracks->GetYaxis()->SetLabelSize(0.05);
		locHist_NumTimeBasedTracks->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackingFOM != NULL)
	{
		locHist_TrackingFOM->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackingFOM->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackingFOM->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackingFOM->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackingFOM->Draw();
		gPad->SetLogy();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_NumDCHitsPerTrackVsTheta != NULL)
	{
		locHist_NumDCHitsPerTrackVsTheta->GetYaxis()->SetRangeUser(5.0, 30.0);
		locHist_NumDCHitsPerTrackVsTheta->GetXaxis()->SetTitleSize(0.05);
		locHist_NumDCHitsPerTrackVsTheta->GetYaxis()->SetTitleSize(0.05);
		locHist_NumDCHitsPerTrackVsTheta->GetXaxis()->SetLabelSize(0.05);
		locHist_NumDCHitsPerTrackVsTheta->GetYaxis()->SetLabelSize(0.05);
		locHist_NumDCHitsPerTrackVsTheta->Draw("COLZ");
		gPad->SetLogz();
	}
}


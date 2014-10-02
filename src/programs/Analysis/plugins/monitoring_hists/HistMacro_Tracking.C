// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/NumDCHitsPerTrackVsTheta
// hnamepath: /Independent/Hist_DetectorStudies/ReconstructedPID/Proton/Tracking_ChiSqPerDF
// hnamepath: /Independent/Hist_DetectorStudies/ReconstructedPID/Pi-/Tracking_ChiSqPerDF
// hnamepath: /Independent/Hist_DetectorStudies/ReconstructedPID/Pi+/Tracking_ChiSqPerDF

// hnamepath: /Independent/Hist_NumReconstructedObjects/NumTrackCandidates
// hnamepath: /Independent/Hist_NumReconstructedObjects/NumWireBasedTracks

{
	//Goto Path
	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Reconstruction");
	TH2I* locHist_NumDCHitsPerTrackVsTheta = (TH2I*)gDirectory->Get("NumDCHitsPerTrackVsTheta");
	gDirectory->cd("ReconstructedPID/Proton");
	TH1I* locHist_TrackingChiSqPerDF_Proton = (TH1I*)gDirectory->Get("Tracking_ChiSqPerDF");
	gDirectory->cd("../Pi+");
	TH1I* locHist_TrackingChiSqPerDF_PiPlus = (TH1I*)gDirectory->Get("Tracking_ChiSqPerDF");
	gDirectory->cd("../Pi-");
	TH1I* locHist_TrackingChiSqPerDF_PiMinus = (TH1I*)gDirectory->Get("Tracking_ChiSqPerDF");

	locDirectory = (TDirectory*)locInitDirectory->FindObjectAny("Hist_NumReconstructedObjects");
	if(!locDirectory)
		return;
	locDirectory->cd();
	TH1I* locHist_NumTrackCandidates = (TH1I*)gDirectory->Get("NumTrackCandidates");
	TH1I* locHist_NumWireBasedTracks = (TH1I*)gDirectory->Get("NumWireBasedTracks");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Kinematics"); //for testing
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
	if(locHist_NumDCHitsPerTrackVsTheta != NULL)
	{
		locHist_NumDCHitsPerTrackVsTheta->GetXaxis()->SetTitleSize(0.05);
		locHist_NumDCHitsPerTrackVsTheta->GetYaxis()->SetTitleSize(0.05);
		locHist_NumDCHitsPerTrackVsTheta->GetXaxis()->SetLabelSize(0.05);
		locHist_NumDCHitsPerTrackVsTheta->GetYaxis()->SetLabelSize(0.05);
		locHist_NumDCHitsPerTrackVsTheta->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackingChiSqPerDF_PiPlus != NULL)
	{
		locHist_TrackingChiSqPerDF_PiPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackingChiSqPerDF_PiPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackingChiSqPerDF_PiPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackingChiSqPerDF_PiPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackingChiSqPerDF_PiPlus->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackingChiSqPerDF_PiMinus != NULL)
	{
		locHist_TrackingChiSqPerDF_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackingChiSqPerDF_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackingChiSqPerDF_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackingChiSqPerDF_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackingChiSqPerDF_PiMinus->Draw();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackingChiSqPerDF_Proton != NULL)
	{
		locHist_TrackingChiSqPerDF_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackingChiSqPerDF_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackingChiSqPerDF_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackingChiSqPerDF_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackingChiSqPerDF_Proton->Draw();
	}
}


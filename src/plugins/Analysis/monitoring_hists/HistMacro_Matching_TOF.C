// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TOFTrackDistanceVsP
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TOFTrackDistanceVsTheta
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TOFTrackDeltaXVsVerticalPaddle
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TOFTrackDeltaXVsHorizontalPaddle
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TOFTrackDeltaYVsVerticalPaddle
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TOFTrackDeltaYVsHorizontalPaddle
{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatching");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("TimeBased/TOFPoint");
	TH2I* locHist_TOFTrackDistanceVsP = (TH2I*)gDirectory->Get("TOFTrackDistanceVsP");
	TH2I* locHist_TOFTrackDistanceVsTheta = (TH2I*)gDirectory->Get("TOFTrackDistanceVsTheta");
	TH2I* locHist_TrackDeltaXVsVerticalPaddle = (TH2I*)gDirectory->Get("TOFTrackDeltaXVsVerticalPaddle");
	TH2I* locHist_TrackDeltaXVsHorizontalPaddle = (TH2I*)gDirectory->Get("TOFTrackDeltaXVsHorizontalPaddle");
	TH2I* locHist_TrackDeltaYVsVerticalPaddle = (TH2I*)gDirectory->Get("TOFTrackDeltaYVsVerticalPaddle");
	TH2I* locHist_TrackDeltaYVsHorizontalPaddle = (TH2I*)gDirectory->Get("TOFTrackDeltaYVsHorizontalPaddle");

	//Get original pad margins
	double locLeftPadMargin = gStyle->GetPadLeftMargin();
	double locRightPadMargin = gStyle->GetPadRightMargin();
	double locTopPadMargin = gStyle->GetPadTopMargin();
	double locBottomPadMargin = gStyle->GetPadBottomMargin();

	//Set new pad margins
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.15);
//	gStyle->SetPadTopMargin(locTopPadMargin);
//	gStyle->SetPadBottomMargin(locBottomPadMargin);

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_TOF", "Matching_TOF", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//PADS 1 & 4: IS THE GEOMETRY & ALIGNMENT OK?
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackDeltaXVsVerticalPaddle != NULL)
	{
		locHist_TrackDeltaXVsVerticalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaXVsVerticalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaXVsVerticalPaddle->Draw("COLZ");
		locHist_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetTitleOffset(1.3);
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackDeltaYVsHorizontalPaddle != NULL)
	{
		locHist_TrackDeltaYVsHorizontalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaYVsHorizontalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaYVsHorizontalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaYVsHorizontalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaYVsHorizontalPaddle->Draw("COLZ");
		locHist_TrackDeltaYVsHorizontalPaddle->GetYaxis()->SetTitleOffset(1.3);
	}

	//PADS 2 & 5: VS KINEMATICS
	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFTrackDistanceVsP != NULL)
	{
		locHist_TOFTrackDistanceVsP->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFTrackDistanceVsP->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFTrackDistanceVsP->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFTrackDistanceVsP->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFTrackDistanceVsP->Draw("COLZ");
		locHist_TOFTrackDistanceVsP->GetYaxis()->SetTitleOffset(1.3);
		TF1* locFunc = new TF1("TOF_DCut", "exp(-1.0*[0]*x + [1]) + [2]", 0.0, 10.0);
		locFunc->SetParameters(1.1, 1.5, 6.15);
		locFunc->Draw("SAME");
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFTrackDistanceVsTheta != NULL)
	{
		locHist_TOFTrackDistanceVsTheta->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFTrackDistanceVsTheta->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFTrackDistanceVsTheta->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFTrackDistanceVsTheta->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFTrackDistanceVsTheta->Draw("COLZ");
		locHist_TOFTrackDistanceVsTheta->GetYaxis()->SetTitleOffset(1.3);
	}

	//PADS 3 & 6: IS THE CALIBRATION OK?
	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackDeltaXVsHorizontalPaddle != NULL)
	{
		locHist_TrackDeltaXVsHorizontalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaXVsHorizontalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaXVsHorizontalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaXVsHorizontalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaXVsHorizontalPaddle->Draw("COLZ");
		locHist_TrackDeltaXVsHorizontalPaddle->GetYaxis()->SetTitleOffset(1.3);
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackDeltaYVsVerticalPaddle != NULL)
	{
		locHist_TrackDeltaYVsVerticalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaYVsVerticalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaYVsVerticalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaYVsVerticalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaYVsVerticalPaddle->Draw("COLZ");
		locHist_TrackDeltaYVsVerticalPaddle->GetYaxis()->SetTitleOffset(1.3);
	}

	//Reset original pad margins
	gStyle->SetPadLeftMargin(locLeftPadMargin);
	gStyle->SetPadRightMargin(locRightPadMargin);
	gStyle->SetPadTopMargin(locTopPadMargin);
	gStyle->SetPadBottomMargin(locBottomPadMargin);
}

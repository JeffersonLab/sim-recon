// hnamepath: /Independent/Hist_DetectorMatchParams/ReconstructedPID/Proton/FCALShowerDeltaT
// hnamepath: /Independent/Hist_DetectorMatchParams/ReconstructedPID/Proton/FCALShowerDeltaTVsP
// hnamepath: /Independent/Hist_DetectorMatchParams/ReconstructedPID/Pi-/FCALShowerDeltaT
// hnamepath: /Independent/Hist_DetectorMatchParams/ReconstructedPID/Pi-/FCALShowerDeltaTVsP
// hnamepath: /Independent/Hist_DetectorPID/FCAL/BetaVsP_q-
// hnamepath: /Independent/Hist_DetectorPID/FCAL/BetaVsP_q+

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatchParams");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("ReconstructedPID/Proton");
	TH1I* locHist_FCALShowerDeltaT_Proton = (TH1I*)gDirectory->Get("FCALShowerDeltaT"); //proton
	TH2I* locHist_FCALShowerDeltaTVsP_Proton = (TH2I*)gDirectory->Get("FCALShowerDeltaTVsP"); //proton

	gDirectory->cd("../Pi-");
	TH1I* locHist_FCALShowerDeltaT_PiMinus = (TH1I*)gDirectory->Get("FCALShowerDeltaT"); //pi-
	TH2I* locHist_FCALShowerDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("FCALShowerDeltaTVsP"); //pi-

	gDirectory->cd("../../../Hist_DetectorPID/FCAL");
	TH2I* locHist_BetaVsP_QMinus = (TH2I*)gDirectory->Get("BetaVsP_q-"); //q-
	TH2I* locHist_BetaVsP_QPlus = (TH2I*)gDirectory->Get("BetaVsP_q+"); //q+

	//Beta-vs-p functions
	TF1* locBetaVsPFunc_Proton = new TF1("BetaVsPFunc_Proton", "x/sqrt(x*x + [0]*[0])", 0.0, 12.0);
	locBetaVsPFunc_Proton->SetParameter(0, 0.938272046);
	locBetaVsPFunc_Proton->SetLineWidth(2);
	locBetaVsPFunc_Proton->SetLineColor(kBlack);
	locBetaVsPFunc_Proton->SetNpx(1000);

	TF1* locBetaVsPFunc_Kaon = new TF1("BetaVsPFunc_Kaon", "x/sqrt(x*x + [0]*[0])", 0.0, 12.0);
	locBetaVsPFunc_Kaon->SetParameter(0, 0.493677);
	locBetaVsPFunc_Kaon->SetLineWidth(2);
	locBetaVsPFunc_Kaon->SetLineColor(kBlack);
	locBetaVsPFunc_Kaon->SetNpx(1000);

	TF1* locBetaVsPFunc_Pion = new TF1("BetaVsPFunc_Pion", "x/sqrt(x*x + [0]*[0])", 0.0, 12.0);
	locBetaVsPFunc_Pion->SetParameter(0, 0.13957018);
	locBetaVsPFunc_Pion->SetLineWidth(2);
	locBetaVsPFunc_Pion->SetLineColor(kBlack);
	locBetaVsPFunc_Pion->SetNpx(1000);

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("FCALReconstruction_p2", "FCALReconstruction_p2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaT_Proton != NULL)
	{
		locHist_FCALShowerDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_Proton->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaT_PiMinus != NULL)
	{
		locHist_FCALShowerDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaTVsP_Proton != NULL)
	{
		locHist_FCALShowerDeltaTVsP_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaTVsP_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaTVsP_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaTVsP_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaTVsP_Proton->Draw("COLZ");
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaTVsP_PiMinus != NULL)
	{
		locHist_FCALShowerDeltaTVsP_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaTVsP_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaTVsP_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaTVsP_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaTVsP_PiMinus->Draw("COLZ");
	}

	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QPlus != NULL)
	{
		locHist_BetaVsP_QPlus->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_BetaVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QPlus->GetListOfFunctions()->Add(locBetaVsPFunc_Proton);
		locHist_BetaVsP_QPlus->GetListOfFunctions()->Add(locBetaVsPFunc_Kaon);
		locHist_BetaVsP_QPlus->GetListOfFunctions()->Add(locBetaVsPFunc_Pion);
		locHist_BetaVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QMinus != NULL)
	{
		locHist_BetaVsP_QMinus->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_BetaVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QMinus->GetListOfFunctions()->Add(locBetaVsPFunc_Pion);
		locHist_BetaVsP_QMinus->GetListOfFunctions()->Add(locBetaVsPFunc_Kaon);
		locHist_BetaVsP_QMinus->Draw("COLZ");
	}
}


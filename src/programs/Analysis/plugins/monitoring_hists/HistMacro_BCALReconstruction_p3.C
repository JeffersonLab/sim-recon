// hnamepath: /Independent/Hist_DetectorPID/BCAL/BetaVsP_q-
// hnamepath: /Independent/Hist_DetectorPID/BCAL/BetaVsP_q0
// hnamepath: /Independent/Hist_DetectorPID/BCAL/BetaVsP_q+
// hnamepath: /Independent/Hist_DetectorPID/BCAL/EOverPVsP_q-
// hnamepath: /Independent/Hist_DetectorPID/BCAL/EOverPVsP_q+

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorPID");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("BCAL");
	TH2I* locHist_BetaVsP_QMinus = (TH2I*)gDirectory->Get("BetaVsP_q-"); //q-
	TH2I* locHist_BetaVsP_QZero = (TH2I*)gDirectory->Get("BetaVsP_q0"); //q0
	TH2I* locHist_BetaVsP_QPlus = (TH2I*)gDirectory->Get("BetaVsP_q+"); //q+
	TH2I* locHist_EOverPVsP_QPlus = (TH2I*)gDirectory->Get("EOverPVsP_q+"); //q+
	TH2I* locHist_EOverPVsP_QMinus = (TH2I*)gDirectory->Get("EOverPVsP_q-"); //q-

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
		locCanvas = new TCanvas("BCALReconstruction_p3", "BCALReconstruction_p3", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QMinus != NULL)
	{
		locHist_BetaVsP_QMinus->GetXaxis()->SetRangeUser(0.0, 2.0);
		locHist_BetaVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QMinus->GetListOfFunctions()->Add(locBetaVsPFunc_Pion);
		locHist_BetaVsP_QMinus->GetListOfFunctions()->Add(locBetaVsPFunc_Kaon);
		locHist_BetaVsP_QMinus->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QZero != NULL)
	{
		locHist_BetaVsP_QZero->GetXaxis()->SetRangeUser(0.0, 2.0);
		locHist_BetaVsP_QZero->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QZero->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QZero->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QZero->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QZero->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QPlus != NULL)
	{
		locHist_BetaVsP_QPlus->GetXaxis()->SetRangeUser(0.0, 2.0);
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
	if(locHist_EOverPVsP_QPlus != NULL)
	{
		locHist_EOverPVsP_QPlus->GetXaxis()->SetRangeUser(0.0, 2.0);
		locHist_EOverPVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_EOverPVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_EOverPVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_EOverPVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_EOverPVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_EOverPVsP_QMinus != NULL)
	{
		locHist_EOverPVsP_QMinus->GetXaxis()->SetRangeUser(0.0, 2.0);
		locHist_EOverPVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_EOverPVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_EOverPVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_EOverPVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_EOverPVsP_QMinus->Draw("COLZ");
	}
}


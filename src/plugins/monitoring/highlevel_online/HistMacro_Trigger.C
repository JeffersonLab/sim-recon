// hnamepath: /highlevel/L1GTPRate
// hnamepath: /highlevel/BCALVsFCAL_TrigBit1
// hnamepath: /highlevel/BCALVsFCAL_TrigBit6

{
	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	TH2* locHist_L1GTPRate = (TH2*)gDirectory->Get("L1GTPRate");
	TH2* locHist_BCALVsFCAL_TrigBit1 = (TH2*)gDirectory->Get("BCALVsFCAL_TrigBit1");
	TH2* locHist_BCALVsFCAL_TrigBit6 = (TH2*)gDirectory->Get("BCALVsFCAL_TrigBit6");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Kinematics", "Kinematics", 1200, 400); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 1);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_L1GTPRate != NULL)
	{
		locHist_L1GTPRate->GetXaxis()->SetTitleSize(0.05);
		locHist_L1GTPRate->GetYaxis()->SetTitleSize(0.05);
		locHist_L1GTPRate->GetXaxis()->SetLabelSize(0.05);
		locHist_L1GTPRate->GetYaxis()->SetLabelSize(0.05);
		locHist_L1GTPRate->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALVsFCAL_TrigBit1 != NULL)
	{
		locHist_BCALVsFCAL_TrigBit1->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALVsFCAL_TrigBit1->GetYaxis()->SetTitleSize(0.05);
		locHist_BCALVsFCAL_TrigBit1->Draw("colz");
		gPad->SetLogz();
		gPad->Update();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALVsFCAL_TrigBit6 != NULL)
	{
		locHist_BCALVsFCAL_TrigBit6->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALVsFCAL_TrigBit6->GetYaxis()->SetTitleSize(0.05);
		locHist_BCALVsFCAL_TrigBit6->Draw("colz");
		gPad->SetLogz();
		gPad->Update();
	}
}

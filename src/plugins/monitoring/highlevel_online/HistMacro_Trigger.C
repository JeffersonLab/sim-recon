// hnamepath: /highlevel/L1GTPRate
// hnamepath: /highlevel/BCALVsFCAL_TrigBit1
// hnamepath: /highlevel/BCALVsFCAL_TrigBit6
//
// e-mail: davidl@jlab.org
// e-mail: pmatt@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
//

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
	
	TLatex latex;
	latex.SetTextSize(0.04);
	char str[256];

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_L1GTPRate != NULL)
	{
		locHist_L1GTPRate->GetXaxis()->SetTitleSize(0.05);
		locHist_L1GTPRate->GetYaxis()->SetTitleSize(0.04);
		locHist_L1GTPRate->GetXaxis()->SetLabelSize(0.05);
		locHist_L1GTPRate->GetYaxis()->SetLabelSize(0.05);
		locHist_L1GTPRate->SetStats(0);
		locHist_L1GTPRate->Draw("colz");
		
		sprintf(str, "from %d sync events", (uint32_t)locHist_L1GTPRate->GetEntries()/8);
		latex.DrawLatex(1.0, 101.0, str);
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALVsFCAL_TrigBit1 != NULL)
	{
		locHist_BCALVsFCAL_TrigBit1->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALVsFCAL_TrigBit1->GetYaxis()->SetTitleSize(0.04);
		locHist_BCALVsFCAL_TrigBit1->SetStats(0);
		locHist_BCALVsFCAL_TrigBit1->Draw("colz");

		sprintf(str, "%d entries", (uint32_t)locHist_BCALVsFCAL_TrigBit1->GetEntries());
		latex.DrawLatex(500.0, 50000.0*1.01, str);

		gPad->SetLogz();
		gPad->Update();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALVsFCAL_TrigBit6 != NULL)
	{
		locHist_BCALVsFCAL_TrigBit6->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALVsFCAL_TrigBit6->GetYaxis()->SetTitleSize(0.04);
		locHist_BCALVsFCAL_TrigBit6->SetStats(0);
		locHist_BCALVsFCAL_TrigBit6->Draw("colz");

		sprintf(str, "%d entries", (uint32_t)locHist_BCALVsFCAL_TrigBit6->GetEntries());
		latex.DrawLatex(500.0, 50000.0*1.01, str);

		gPad->SetLogz();
		gPad->Update();
	}
}

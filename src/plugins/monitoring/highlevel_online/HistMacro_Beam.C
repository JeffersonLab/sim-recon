// hnamepath: /highlevel/RFBeamBunchPeriod
// hnamepath: /highlevel/RFBeamBunchPeriod_DFT
// hnamepath: /highlevel/BeamEnergy
//
// e-mail: davidl@jlab.org
// e-mail: pmatt@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
//

{
	// This taken from the bin contents of the BeamEnergy histogram from the
	// online monitoring root file for run 30890, an amorphous target run.
	// The macro used is in /gluonwork1/Users/davidl/2017.02.08.amorphous_norm
	Double_t amorphous_data[] = {
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,   1874227.0,   1943030.0, 
		  1755316.0,     20181.0,   1703311.0,   1606692.0,   1723045.0,   1619582.0,   1620684.0,   1600446.0,     17238.0,   1566334.0, 
		  1461708.0,   1478741.0,   1606202.0,   1500731.0,   1475408.0,     17868.0,   1404680.0,   1331043.0,   1287114.0,   1262602.0, 
		  1253865.0,   1248368.0,      9273.0,   1206261.0,     10665.0,   1211427.0,   1130631.0,   1150276.0,   1079137.0,   1061930.0, 
		  1074806.0,      8298.0,   1036731.0,   1022121.0,   1040724.0,    955157.0,         0.0,   1031693.0,      7154.0,    918333.0, 
		   955360.0,    954557.0,    925761.0,    872359.0,    861251.0,    822208.0,      5375.0,    878918.0,    828288.0,    798800.0, 
		      113.0,    770563.0,    787358.0,      5257.0,    759065.0,    749725.0,    729757.0,    682289.0,    695635.0,    694138.0, 
		   722611.0,      3874.0,    675319.0,    653009.0,    635066.0,    830412.0,    791875.0,    991135.0,      5793.0,    754065.0, 
		   975395.0,    731504.0,    754759.0,    888878.0,    899057.0,      5818.0,    858852.0,    843084.0,    825755.0,    802431.0, 
		   671132.0,      4624.0,    671460.0,    674781.0,    649224.0,    491714.0,    909804.0,    616382.0,      2849.0,    625648.0, 
		   594743.0,    552795.0,    580985.0,    983213.0,   1268615.0,   1694554.0,   1254491.0,   1398173.0,   2563137.0,   4297711.0, 
		  3220369.0,   3809785.0,   5191254.0,   2848361.0,   3186664.0,   3845752.0,   3542581.0,   3087916.0,   3581239.0,   3542077.0, 
		  2899594.0,   2127684.0,   2142266.0,   2596527.0,   3020573.0,   2498288.0,   3270618.0,   2278964.0,   1659922.0,   1178253.0, 
		   871488.0,   1039194.0,    845424.0,   1047850.0,    702899.0,   1108102.0,   1116642.0,    769841.0,   1065043.0,   1029901.0, 
		  1063093.0,   1001900.0,   1006887.0,    956399.0,    950236.0,    807695.0,    941686.0,    923643.0,   1354655.0,    888567.0, 
		   862495.0,    833182.0,    807991.0,   1339698.0,    848880.0,    904976.0,    877676.0,    844513.0,    823781.0,    787384.0, 
		  1139778.0,    738763.0,   1023720.0,    672321.0,   1006879.0,    961981.0,    902932.0,    803966.0,    773695.0,   1116473.0, 
		   687015.0,    923731.0,    654774.0,    996528.0,    588333.0,    752199.0,    845485.0,    684292.0,    565304.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
	0.0};



	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	TH1* locHist_RFBeamBunchPeriod = (TH1*)gDirectory->Get("RFBeamBunchPeriod");
	TH1* locHist_RFBeamBunchPeriod_DFT = (TH1*)gDirectory->Get("RFBeamBunchPeriod_DFT");
	TH1* locHist_BeamEnergy = (TH1*)gDirectory->Get("BeamEnergy");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Beam", "Beam", 1200, 600); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 1);

	TLatex latex;
	latex.SetTextSize(0.04);
	char str[256];

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if( (locHist_RFBeamBunchPeriod_DFT != NULL) && (locHist_RFBeamBunchPeriod != NULL))
	{
		locHist_RFBeamBunchPeriod_DFT->GetXaxis()->SetTitleSize(0.05);
		locHist_RFBeamBunchPeriod_DFT->GetXaxis()->SetLabelSize(0.05);
		locHist_RFBeamBunchPeriod_DFT->GetYaxis()->SetLabelSize(0.03);
		locHist_RFBeamBunchPeriod_DFT->SetFillStyle(3001);
		locHist_RFBeamBunchPeriod_DFT->SetFillColor(kGreen);
		locHist_RFBeamBunchPeriod_DFT->SetLineColor(kGreen-2);
		locHist_RFBeamBunchPeriod_DFT->SetLineWidth(2);
		locHist_RFBeamBunchPeriod_DFT->SetStats(0);
		locHist_RFBeamBunchPeriod_DFT->Draw();

		sprintf(str, "%d entries", (uint32_t)locHist_RFBeamBunchPeriod->GetEntries());
		latex.DrawLatex(300.0, 1.06*locHist_RFBeamBunchPeriod_DFT->GetMaximum(), str);
		
		// Plot RFBeamBunchPeriod as inset
		TPad *p = (TPad*)gDirectory->FindObjectAny("RF_DFT");
		if(!p) p = new TPad("RF_DFT", "insert", 0.1, 0.6, 0.5, 0.9);
		p->Draw();
		p->cd();

		locHist_RFBeamBunchPeriod->GetXaxis()->SetRangeUser(30.0, 90.0);
		locHist_RFBeamBunchPeriod->GetXaxis()->SetTitleSize(0.05);
		locHist_RFBeamBunchPeriod->GetXaxis()->SetLabelSize(0.05);
		locHist_RFBeamBunchPeriod->GetYaxis()->SetLabelSize(0.03);
		locHist_RFBeamBunchPeriod->SetFillStyle(3001);
		locHist_RFBeamBunchPeriod->SetFillColor(kMagenta);
		locHist_RFBeamBunchPeriod->SetStats(0);
		locHist_RFBeamBunchPeriod->Draw();
		
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BeamEnergy != NULL)
	{
		// Create normalized histogram
		TH1D* locHist_BeamEnergy_norm = (TH1D*)gDirectory->Get("BeamEnergy_norm");
		if(!locHist_BeamEnergy_norm){
			locHist_BeamEnergy_norm = new TH1D("BeamEnergy_norm", "Reconstructed photon beam energy normalized to amorphous run 30898;Beam #gamma energy (GeV)", 240, 0.0, 12.0);

			locHist_BeamEnergy_norm->GetXaxis()->SetTitleSize(0.05);
			locHist_BeamEnergy_norm->GetXaxis()->SetLabelSize(0.05);
			locHist_BeamEnergy_norm->GetYaxis()->SetLabelSize(0.03);
			locHist_BeamEnergy_norm->SetFillStyle(3001);
			locHist_BeamEnergy_norm->SetFillColor(kOrange);
			locHist_BeamEnergy_norm->SetLineColor(kRed-2);
			locHist_BeamEnergy_norm->SetLineWidth(2);
			locHist_BeamEnergy_norm->SetStats(0);
		}
		if(locHist_BeamEnergy_norm){
			for(int ibin=1; ibin<=locHist_BeamEnergy_norm->GetNbinsX(); ibin++){
				Double_t norm = amorphous_data[ibin-1];
				if( norm == 0.0) continue;

				Double_t v = (Double_t)locHist_BeamEnergy->GetBinContent(ibin);
				locHist_BeamEnergy_norm->SetBinContent(ibin, v/norm);
			}
			locHist_BeamEnergy_norm->Draw();
		}

		TPad *beamenergypad = (TPad*)gDirectory->FindObjectAny("beamenergypad");
		if(!beamenergypad) beamenergypad = new TPad("beamenergypad", "", 0.105, 0.65, 0.5, 0.895);
		beamenergypad->SetTicks();
		beamenergypad->Draw();
		beamenergypad->cd();
	
		//locHist_BeamEnergy->Rebin(5);
		locHist_BeamEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_BeamEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_BeamEnergy->GetYaxis()->SetLabelSize(0.03);
		locHist_BeamEnergy->SetFillStyle(3001);
		locHist_BeamEnergy->SetFillColor(kCyan);
		locHist_BeamEnergy->SetLineColor(kBlue);
		locHist_BeamEnergy->SetLineWidth(2);
		locHist_BeamEnergy->SetStats(0);
		locHist_BeamEnergy->Draw();

		sprintf(str, "%d entries", (uint32_t)locHist_BeamEnergy->GetEntries());
		latex.DrawLatex(2.0, 1.06*locHist_BeamEnergy->GetMaximum(), str);
	}
}


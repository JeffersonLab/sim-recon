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
	// online monitoring root file for run 30388, an amorphous target run.
	// The macro used is in /gluonwork1/Users/davidl/2017.02.08.amorphous_norm
	Double_t amorphous_data[] = {
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,
		   207112.0,    214652.0,    193329.0,    187909.0,    177206.0,    188547.0,    179962.0,    177935.0,   1785278.0,   1842415.0,
		  1669679.0,    180174.0,   1638820.0,   1539030.0,   1640370.0,   1539394.0,   1537368.0,   1510950.0,    153961.0,   1480462.0,
		  1387942.0,   1271081.0,   1511105.0,   1419369.0,   1386856.0,    141496.0,   1322292.0,   1255687.0,   1217868.0,   1198187.0,
		  1186159.0,   1184043.0,    112608.0,   1031092.0,    122066.0,   1136319.0,   1071261.0,   1088562.0,   1025939.0,   1003034.0,
		  1014784.0,     97242.0,    982312.0,    962518.0,    980330.0,    814799.0,     84245.0,    971352.0,    868050.0,     87580.0,
		   895712.0,    892575.0,    868272.0,    821409.0,    817088.0,    778239.0,     76396.0,    821478.0,    797287.0,    769617.0,
		   108929.0,    741326.0,    779274.0,     85059.0,    731231.0,    738503.0,    720861.0,    677556.0,    686818.0,    683617.0,
		   707237.0,      3443.0,    650321.0,    633095.0,    615594.0,    782139.0,    731123.0,    946509.0,     71900.0,    713521.0,
		   898371.0,    689203.0,    751601.0,    914986.0,    929665.0,    168973.0,    932448.0,    965654.0,   1266386.0,   1098035.0,
		  1190110.0,    497183.0,    928331.0,   1050795.0,    978075.0,    847976.0,   1260333.0,    976839.0,    284750.0,    826013.0,
		   867393.0,    811261.0,    842630.0,   1477289.0,   1367978.0,   1293970.0,   1202063.0,   1318935.0,   2033424.0,   3008094.0,
		  2936686.0,   3484083.0,   4766919.0,   3196669.0,   2268946.0,   3455931.0,   3252551.0,   2786139.0,   3272483.0,   3250077.0,
		  2620797.0,   1974946.0,   2031168.0,   2808866.0,   2369345.0,   2320531.0,   2981188.0,   2097775.0,   2050197.0,    674196.0,
		   875506.0,    979724.0,    846943.0,   1013793.0,    716036.0,   1093307.0,   1054607.0,   1016535.0,    741128.0,   1009375.0,
		  1022691.0,    964256.0,    927722.0,    928634.0,    912294.0,    799601.0,    897200.0,    876986.0,   1185298.0,    763340.0,
		   733589.0,    708865.0,    999313.0,    833832.0,    725050.0,    772515.0,    748708.0,    721066.0,    705554.0,    671451.0,
		   973485.0,    631075.0,    879034.0,    568478.0,    861210.0,    822981.0,    773300.0,    688667.0,    663987.0,    961948.0,
		   589918.0,    796089.0,    565022.0,    859376.0,    509337.0,    652145.0,    734002.0,    595229.0,    486093.0,         0.0,
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0
	};



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
			locHist_BeamEnergy_norm = new TH1D("BeamEnergy_norm", "Reconstructed photon beam energy normalized to amorphous run 30388;Beam #gamma energy (GeV)", 240, 0.0, 12.0);

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


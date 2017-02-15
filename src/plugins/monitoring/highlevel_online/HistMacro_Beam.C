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
	// online monitoring root file for run 30437, an amorphous target run.
	// The macro used is in /gluonwork1/Users/davidl/2017.02.08.amorphous_norm
	string amorphous_label = "Normalized to Amorphous run 30575";

	Double_t amorphous_data[] = {
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0, 
		        0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,  12499828.0,  12912764.0, 
		 11715557.0,    133235.0,  11367286.0,  10703321.0,  11490983.0,  10806950.0,  10830532.0,  10665557.0,    112860.0,  10441301.0, 
		  9760634.0,   9886114.0,  10699019.0,  10011472.0,   9843526.0,    116214.0,   9374345.0,   8877996.0,   8575177.0,   8438387.0, 
		  8362988.0,   8329770.0,     60609.0,   8045744.0,     68557.0,   8070172.0,   7546779.0,   7670028.0,   7202497.0,   7076441.0, 
		  7169419.0,     55829.0,   6914715.0,   6810160.0,   6953525.0,   6369064.0,         0.0,   6871242.0,   6130439.0,     40746.0, 
		  6370973.0,   6361449.0,   6180239.0,   5813841.0,   5752260.0,   5493445.0,     34406.0,   5858238.0,   5517410.0,   5334706.0, 
		      362.0,   5132287.0,   5248854.0,     34354.0,   5055089.0,   5005140.0,   4860634.0,   4561905.0,   4645770.0,   4636649.0, 
		  4825429.0,     24668.0,   4503646.0,   4363692.0,   4239201.0,   5545297.0,   5284855.0,   6619124.0,     37672.0,   5031279.0, 
		  6505771.0,   4924161.0,   5002071.0,   5927795.0,   6016538.0,     38722.0,   5743276.0,   5630685.0,   5512178.0,   5359186.0, 
		  4484307.0,     30038.0,   4481334.0,   4492860.0,   4324970.0,   3280944.0,   6093542.0,   4121913.0,     17526.0,   4185522.0, 
		  3974795.0,   3698084.0,   3890470.0,   8848428.0,   8423308.0,   9071154.0,   8381794.0,   9326720.0,   7098550.0,   5843469.0, 
		  5716399.0,   7166468.0,   6366792.0,   6872323.0,   5289462.0,   6540325.0,   6523947.0,   6321828.0,   6138324.0,   6111317.0, 
		  5905414.0,   6218892.0,   6205598.0,   7649013.0,   5437156.0,   4583597.0,   7403847.0,   4663690.0,   7667759.0,   5593122.0, 
		  5825496.0,   6929872.0,   5644297.0,   6997689.0,   4697671.0,   7423557.0,   7441152.0,   7176474.0,   5100355.0,   6878409.0, 
		  7100159.0,   6700049.0,   6734182.0,   6381398.0,   6349889.0,   5404426.0,   6317546.0,   6168787.0,   9036771.0,   5991137.0, 
		  5751732.0,   5560510.0,   7823666.0,   6518472.0,   5671903.0,   6052815.0,   5859889.0,   5638719.0,   5507838.0,   5254143.0, 
		  7611980.0,   4940460.0,   6867833.0,   4447819.0,   6723825.0,   6425042.0,   6036879.0,   5366979.0,   5169345.0,   7472186.0, 
		  4580756.0,   6165344.0,   4369330.0,   6627373.0,   3918911.0,   5004494.0,   5590983.0,   4536820.0,   3702407.0,         0.0, 
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
		double max = locHist_RFBeamBunchPeriod_DFT->GetMaximum()*2.0;
		locHist_RFBeamBunchPeriod_DFT->GetXaxis()->SetTitleSize(0.05);
		locHist_RFBeamBunchPeriod_DFT->GetXaxis()->SetLabelSize(0.05);
		locHist_RFBeamBunchPeriod_DFT->GetYaxis()->SetLabelSize(0.03);
		locHist_RFBeamBunchPeriod_DFT->SetFillStyle(3001);
		locHist_RFBeamBunchPeriod_DFT->SetFillColor(kGreen);
		locHist_RFBeamBunchPeriod_DFT->SetLineColor(kGreen-2);
		locHist_RFBeamBunchPeriod_DFT->SetLineWidth(2);
		locHist_RFBeamBunchPeriod_DFT->SetStats(0);
		locHist_RFBeamBunchPeriod_DFT->GetYaxis()->SetRangeUser(0.0, max);
		locHist_RFBeamBunchPeriod_DFT->Draw();

		sprintf(str, "%g entries", (double)locHist_RFBeamBunchPeriod->GetEntries());
		latex.DrawLatex(300.0, 1.02*max, str);
		
		// Plot RFBeamBunchPeriod as inset
		TPad *p = (TPad*)gDirectory->FindObjectAny("RF_DFT");
		if(!p) p = new TPad("RF_DFT", "insert", 0.16, 0.6, 0.89, 0.9);
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
			locHist_BeamEnergy_norm = new TH1D("BeamEnergy_norm", "Reconstructed Photon Beam Energy;Beam #gamma energy (GeV)", 240, 0.0, 12.0);

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
		
			// Find leftmost non-zero bin for scaling whole hist
			double scale = 0.0;
			for(int ibin=1; ibin<=locHist_BeamEnergy->GetNbinsX(); ibin++){
				if( amorphous_data[ibin-1] < 1000.0) continue;
				Double_t v = (Double_t)locHist_BeamEnergy->GetBinContent(ibin);
				if(v>0.0){
					scale = v/amorphous_data[ibin-1];
					break;
				}
			}

			// Normalize to amorphous baseline including scaling factor
			// above so that leftmost non-zero bin is always 1
			double max = 0.0;
			for(int ibin=1; ibin<=locHist_BeamEnergy_norm->GetNbinsX(); ibin++){
				Double_t norm = amorphous_data[ibin-1];
				if( norm < 1000.0) continue;
				norm *= scale;

				Double_t v = (Double_t)locHist_BeamEnergy->GetBinContent(ibin);
				locHist_BeamEnergy_norm->SetBinContent(ibin, v/norm);
				
				if(v/norm  > max) max = v/norm;
			}
			max = (max-0.75)*1.05 + 0.75;
			if(max > 3.0) max = 3.0;
			if(max < 1.4) max = 1.4;
			locHist_BeamEnergy_norm->GetYaxis()->SetRangeUser(0.75, max);
			locHist_BeamEnergy_norm->Draw();
			
			// If maximum is > 1.5 then assume this is not an amorphous run
			// and draw a label of the peak energy
			if(max > 1.5){
				double Epeak = locHist_BeamEnergy_norm->GetBinCenter(locHist_BeamEnergy_norm->GetMaximumBin());
				sprintf(str, "Epeak: %3.2f GeV", Epeak);
				latex.SetTextAlign(12);
				latex.DrawLatex(1.0, 0.5*(max-0.75)+0.75, str);
			}

			sprintf(str, "%g entries", (double)locHist_BeamEnergy->GetEntries());
			latex.SetTextAlign(22);
			latex.DrawLatex(6.0, 1.035*(max-0.75)+0.75, str);

			latex.SetTextAngle(270);
			latex.DrawLatex(12.5, 0.5*(max-0.75)+0.75, amorphous_label.c_str());
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
	}
}


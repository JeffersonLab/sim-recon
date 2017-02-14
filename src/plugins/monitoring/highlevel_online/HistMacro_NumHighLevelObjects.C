// hnamepath: /highlevel/NumHighLevelObjects
// hnamepath: /highlevel/F1TDC_fADC_tdiff
//
// e-mail: davidl@jlab.org
// e-mail: pmatt@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
//

{
	TDirectory* locCurrentDir = gDirectory;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH2* NumHighLevelObjects = (TH2*)gDirectory->Get("NumHighLevelObjects");
	TH2* F1TDC_fADC_tdiff = (TH2*)gDirectory->Get("F1TDC_fADC_tdiff");

	//Get/Make Canvas
	TCanvas *c1 = NULL;
	if(TVirtualPad::Pad() == NULL)
		c1 = new TCanvas("NumHighLevelObjects", "NumHighLevelObjects", 1200, 900); //for testing
	else
		c1 = gPad->GetCanvas();

	c1->Divide(1,2);
	//c1->Draw();
	
	if(NumHighLevelObjects){
		c1->cd(1);
		gPad->SetTicks();
		gPad->SetGrid();
		gPad->SetLogz();
		gPad->SetBottomMargin(0.2);
		NumHighLevelObjects->GetXaxis()->SetLabelSize(0.07);
		NumHighLevelObjects->SetStats(0);
		NumHighLevelObjects->SetTitle("Num. High Level Objects (all triggers)");
		NumHighLevelObjects->Draw("COLZ");
	
		// Draw mean number of hits
		TLatex latex;
		latex.SetTextColor(kWhite);
		latex.SetTextSize(0.05);
		latex.SetTextAlign(22);
		latex.SetTextAngle(90.0);
		TBox box;
		box.SetFillColor(kCyan+1);
		latex.SetTextColor(kBlack);
		TAxis *xaxis = NumHighLevelObjects->GetXaxis();
		TAxis *yaxis = NumHighLevelObjects->GetYaxis();
		for(int ibin=1; ibin<=xaxis->GetNbins(); ibin++){
			double x = xaxis->GetBinCenter(ibin);
			double y = 90.0;
			
			double sum  = 0.0;
			double sumw = 0.0;			
			for(int jbin=1; jbin<=yaxis->GetNbins(); jbin++){
 				double w = NumHighLevelObjects->GetBinContent(ibin, jbin);
 				double y = yaxis->GetBinCenter(jbin);
 				sumw += y*w;
 				sum  += w; 
			}
			
			double mean = sum==0.0 ? 0.0:sumw/sum;
			
			char str[256];
			sprintf(str, "%5.1f", mean);
			box.DrawBox(x-0.2, y-6.5, x+0.25, y+8.5);
			latex.DrawLatex(x, y, str);
		}
	}

	if(F1TDC_fADC_tdiff){
		c1->cd(2);
		gPad->SetTicks();
		gPad->SetGrid();
		gPad->SetLogz();
		gPad->SetBottomMargin(0.25);
		gPad->SetLeftMargin(0.04);
		gPad->SetRightMargin(0.08);
		F1TDC_fADC_tdiff->GetXaxis()->SetLabelSize(0.06);
		F1TDC_fADC_tdiff->GetYaxis()->SetTitleOffset(0.6);
		F1TDC_fADC_tdiff->SetStats(0);
		F1TDC_fADC_tdiff->Draw("COLZ");
	}

}


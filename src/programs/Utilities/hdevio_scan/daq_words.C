
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /DAQ/daq_words_by_type

void daq_words(void)
{

	TH1D *daq_words_by_type = (TH1D*)gROOT->FindObject("daq_words_by_type");
	if(!daq_words_by_type){
		cout << "Can't find 'daq_words_by_type' histogram! (aborting macro)" << endl;
		return;
	}
	
	TCanvas *c1 = new TCanvas("c1", "", 1600, 800);
	c1->Draw();

	if(!gPad){
		cout << "gPad==NULL!! aborting macro" << endl;
		return;
	}

	gPad->SetBottomMargin(0.3);
	gPad->SetLeftMargin(0.05);
	gPad->SetRightMargin(0.02);
	
	//c1 = new TCanvas("c1", "", 1500, 600);
	gPad->GetCanvas()->SetTicks();
	gPad->GetCanvas()->SetLogy();
	
	double max = daq_words_by_type->GetMaximum()*1.3;
	double min = 1.0;
	
	daq_words_by_type->SetStats(0);
	daq_words_by_type->SetYTitle("Total number of 32-bit words");
	daq_words_by_type->GetYaxis()->SetTitleOffset(0.6);
	daq_words_by_type->GetYaxis()->SetRangeUser(min, max);
	daq_words_by_type->SetFillColor(kBlue);
	daq_words_by_type->Draw("bar");
	
	double ylab = pow(10.0, 0.3*(log(max)-log(min)) + log(min));
	
	TLatex latex;
	latex.SetTextSize(0.030);
	latex.SetTextAlign(12);
	latex.SetTextAngle(90.0);
	
	Int_t Nbins = daq_words_by_type->GetNbinsX();
	double Nword_tot = daq_words_by_type->GetBinContent(Nbins-1);
	double sum = 0;
	for(int ibin=1; ibin<=Nbins; ibin++){
		
		double x = daq_words_by_type->GetXaxis()->GetBinCenter(ibin);
		double Nwords = daq_words_by_type->GetBinContent(ibin);
		double percent = 100.0 * Nwords / Nword_tot;
		
		if(Nwords == 0) continue;
		
		if(ibin < Nbins-2) sum += Nwords;

		char str[256];
		sprintf(str, "%3.1f%%", percent);
		if(ibin < Nbins-2) sprintf(str, "%s _{(%3.1f%% total)}", str, 100.0*sum/Nword_tot);
		latex.DrawLatex( x, ylab, str);

	}

	cout << "       sum: " << sum << endl;
	cout << " Nword_tot: " << Nword_tot << endl;
	cout << "   missing: " << (Nword_tot-sum)/Nword_tot*100.0 << "%" << endl;

	c1->SaveAs("daq_words.png");
	c1->SaveAs("daq_words.pdf");
}


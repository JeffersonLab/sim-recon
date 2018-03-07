//
// Guidance: --------------------------------------------
//
// This shows the tagger microscope (top) and hodoscope(bottom)
// occupancies.
//
// Please note that for amorphous radiator running there will be
// 3 channels in the upper plot of fADC values around colun 82
// that are nearly twice as high as its neighbors. This effect
// is not as pronounced for the diamond radiator which is shown
// in the reference plot. This is expected and should be ignored.
//
// For all other issues contact:
//
//  Alex Barnes - microscope: <see white board>
//
//  Nathan Sparks - hodoscope: <see white board>
//
//  Run Coordinator:  (757) 383-5542
//
//  David Lawrence:   x5567    (office)
//                    746-6697 (cell)
//                    595-1302 (home)
//
// End Guidance: ----------------------------------------
//
//
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/tag_num_events
// hnamepath: /occupancy/tagh_adc_occ
// hnamepath: /occupancy/tagh_tdc_occ
// hnamepath: /occupancy/tagm_adc_occ
// hnamepath: /occupancy/tagm_tdc_occ
//
// e-mail: davidl@jlab.org
// e-mail: aebarnes@jlab.org
// e-mail: somov@jlab.org
// e-mail: richard.t.jones@uconn.edu
// e-mail: tbritton@jlab.org
//
{
	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.
        TDirectory *savedir = gDirectory;
	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	double Nevents = 1.0;
	TH1I *tag_num_events = (TH1I*)gDirectory->FindObjectAny("tag_num_events");
	if(tag_num_events) Nevents = (double)tag_num_events->GetBinContent(1);

	TH2I *tagh_adc_occ = (TH2I*)gDirectory->FindObjectAny("tagh_adc_occ");
	TH2I *tagh_tdc_occ = (TH2I*)gDirectory->FindObjectAny("tagh_tdc_occ");
	TH2I *tagm_adc_occ = (TH2I*)gDirectory->FindObjectAny("tagm_adc_occ");
	TH2I *tagm_tdc_occ = (TH2I*)gDirectory->FindObjectAny("tagm_tdc_occ");

	TLegend *legend_sa = new TLegend(0.1,0.85,0.3,0.9);
	TLegend *legend_na = new TLegend(0.3,0.85,0.5,0.9);

	TLegend *legend_s = new TLegend(0.1,0.85,0.3,0.9);
	TLegend *legend_n = new TLegend(0.3,0.85,0.5,0.9);

	if(tagh_tdc_occ){
		tagh_tdc_occ->SetBarWidth(0.5);
		tagh_tdc_occ->SetBarOffset(0);
		tagh_tdc_occ->SetFillColor(kBlack);
		tagh_tdc_occ->SetStats(0);
		tagh_tdc_occ->SetTitle("TAGH column occupancy");
		tagh_tdc_occ->SetXTitle("Column number");
		tagh_tdc_occ->SetYTitle("fADC/TDC occupancy");
		tagh_tdc_occ->SetTitleSize(0.05,"X");
		tagh_tdc_occ->GetXaxis()->CenterTitle();
		tagh_tdc_occ->SetTitleSize(0.05,"Y");
		tagh_tdc_occ->GetYaxis()->CenterTitle();
	}
	
	if(tagh_adc_occ){
		tagh_adc_occ->SetBarWidth(0.5);
		tagh_adc_occ->SetBarOffset(0.5);
		tagh_adc_occ->SetFillColor(kGreen);
		tagh_adc_occ->SetStats(0);
		tagh_adc_occ->SetTitle("TAGH column occupancy");
		tagh_adc_occ->SetXTitle("Column number");
		tagh_adc_occ->SetYTitle("fADC/TDC occupancy");
		tagh_adc_occ->SetTitleSize(0.05,"X");
		tagh_adc_occ->GetXaxis()->CenterTitle();
		tagh_adc_occ->SetTitleSize(0.05,"Y");
		tagh_adc_occ->GetYaxis()->CenterTitle();
	}
	
	legend_sa->AddEntry(tagh_adc_occ,"fADC","f");
	legend_na->AddEntry(tagh_tdc_occ,"TDC","f");

	if(tagm_tdc_occ){
		tagm_tdc_occ->SetBarWidth(0.5);
		tagm_tdc_occ->SetBarOffset(0);
		tagm_tdc_occ->SetFillColor(kMagenta);
		tagm_tdc_occ->SetStats(0);
		tagm_tdc_occ->SetTitle("TAGM column occupancy");
		tagm_tdc_occ->SetXTitle("Column number");
		tagm_tdc_occ->SetYTitle("fADC/TDC occupancy");
		tagm_tdc_occ->SetTitleSize(0.05,"X");
		tagm_tdc_occ->GetXaxis()->CenterTitle();
		tagm_tdc_occ->SetTitleSize(0.05,"Y");
		tagm_tdc_occ->GetYaxis()->CenterTitle();
	}
	
	if(tagm_adc_occ){
		tagm_adc_occ->SetBarWidth(0.5);
		tagm_adc_occ->SetBarOffset(0.5);
		tagm_adc_occ->SetFillColor(kBlue);
		tagm_adc_occ->SetStats(0);
		tagm_adc_occ->SetTitle("TAGM column occupancy");
		tagm_adc_occ->SetXTitle("Column number");
		tagm_adc_occ->SetYTitle("fADC/TDC occupancy");
		tagm_adc_occ->SetTitleSize(0.05,"X");
		tagm_adc_occ->GetXaxis()->CenterTitle();
		tagm_adc_occ->SetTitleSize(0.05,"Y");
		tagm_adc_occ->GetYaxis()->CenterTitle();
	}

	legend_s->AddEntry(tagm_adc_occ,"fADC","f");
	legend_n->AddEntry(tagm_tdc_occ,"TDC","f");

	char nevents_str[256];
	sprintf(nevents_str,"%0.0f events", Nevents);
	TLatex lat;
	lat.SetTextAlign(22);
	lat.SetTextSize(0.035);

	// Just for testing
	if(gPad == NULL){
		TCanvas *c1 = new TCanvas("c1");
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}
	if(!gPad) {savedir->cd(); return;}

	TCanvas *c1 = gPad->GetCanvas();
	c1->cd(0);
	c1->Clear();

	c1->Divide(1,2);

	// ----------------- TAGM ----------------------
	TVirtualPad *pad1 = c1->cd(1);
	pad1->SetTicks();
	pad1->SetGridy();
	double ymax = 1.0;
	if(tagm_adc_occ) ymax = tagm_adc_occ->GetMaximum();
	if(tagm_tdc_occ){
		if(tagm_tdc_occ->GetMaximum() > ymax) ymax = tagm_tdc_occ->GetMaximum();
	}
	if(tagm_adc_occ){
		tagm_adc_occ->GetYaxis()->SetRangeUser(0.0, 1.05*ymax);
		tagm_adc_occ->Draw("BAR");
	}
	if(tagm_tdc_occ){
		tagm_tdc_occ->GetYaxis()->SetRangeUser(0.0, 1.05*ymax);
		tagm_tdc_occ->Draw( tagm_adc_occ==NULL ? "BAR":"BAR same" );
	}

	legend_s->Draw();
	legend_n->Draw();
	lat.DrawLatex(85.0, 1.075*ymax, nevents_str);

	// ----------------- TAGH ----------------------
	pad1 = c1->cd(2);
	pad1->SetTicks();
	pad1->SetGridy();
	ymax = 1.0;
	if(tagh_adc_occ) ymax = tagh_adc_occ->GetMaximum();
	if(tagh_tdc_occ){
		if(tagh_tdc_occ->GetMaximum() > ymax) ymax = tagh_tdc_occ->GetMaximum();
	}
	if(tagh_adc_occ){
		tagh_adc_occ->GetYaxis()->SetRangeUser(0.0, 1.05*ymax);
		tagh_adc_occ->Draw("BAR");
	}
	if(tagh_tdc_occ){
		tagh_tdc_occ->GetYaxis()->SetRangeUser(0.0, 1.05*ymax);
		tagh_tdc_occ->Draw( tagh_adc_occ==NULL ? "BAR":"BAR same" );
	}
	
	legend_sa->Draw();
	legend_na->Draw();
	lat.DrawLatex(225.0, 1.075*ymax, nevents_str);

}




// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /FDC/fdc_num_events
// hnamepath: /FDC/fdc_occ_plane_18
// hnamepath: /FDC/fdc_occ_plane_19
// hnamepath: /FDC/fdc_occ_plane_20
// hnamepath: /FDC/fdc_occ_plane_21
// hnamepath: /FDC/fdc_occ_plane_22
// hnamepath: /FDC/fdc_occ_plane_23

{
	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("FDC");
	if(dir) dir->cd();

	TH2F *fdc_occ_cell_1 = (TH2F*)gDirectory->FindObjectAny("fdc_occ_plane_18");
    TH2F *fdc_occ_cell_2 = (TH2F*)gDirectory->FindObjectAny("fdc_occ_plane_19");
    TH2F *fdc_occ_cell_3 = (TH2F*)gDirectory->FindObjectAny("fdc_occ_plane_20");
    TH2F *fdc_occ_cell_4 = (TH2F*)gDirectory->FindObjectAny("fdc_occ_plane_21");
    TH2F *fdc_occ_cell_5 = (TH2F*)gDirectory->FindObjectAny("fdc_occ_plane_22");
    TH2F *fdc_occ_cell_6 = (TH2F*)gDirectory->FindObjectAny("fdc_occ_plane_23");

	double Nevents = 1.0;
	TH1I *fdc_num_events = (TH1I*)gDirectory->FindObjectAny("fdc_num_events");
	if(fdc_num_events) Nevents = (double)fdc_num_events->GetBinContent(1);

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
	
	c1->Divide(3,2);
	
	TVirtualPad *pad1 = c1->cd(1);
	pad1->SetTicks();
    pad1->SetLogz();
	pad1->SetLeftMargin(0.15);
	pad1->SetRightMargin(0.15);
	if(fdc_occ_cell_1){
		fdc_occ_cell_1->SetStats(0);
        fdc_occ_cell_1->Scale(1./Nevents);
		fdc_occ_cell_1->Draw("colz");
	}

    TVirtualPad *pad2 = c1->cd(2);
    pad2->SetTicks();
    pad2->SetLogz();
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.15);
    if(fdc_occ_cell_2){
        fdc_occ_cell_2->SetStats(0);
        fdc_occ_cell_2->Scale(1./Nevents);
        fdc_occ_cell_2->Draw("colz");
    }

    TVirtualPad *pad3 = c1->cd(3);
    pad3->SetTicks();
    pad3->SetLogz();
    pad3->SetLeftMargin(0.15);
    pad3->SetRightMargin(0.15);
    if(fdc_occ_cell_3){
        fdc_occ_cell_3->SetStats(0);
        fdc_occ_cell_3->Scale(1./Nevents);
        fdc_occ_cell_3->Draw("colz");
    }

    TVirtualPad *pad4 = c1->cd(4);
    pad4->SetTicks();
    pad4->SetLogz();
    pad4->SetLeftMargin(0.15);
    pad4->SetRightMargin(0.15);
    if(fdc_occ_cell_4){
        fdc_occ_cell_4->SetStats(0);
        fdc_occ_cell_4->Scale(1./Nevents);
        fdc_occ_cell_4->Draw("colz");
    }

    TVirtualPad *pad5 = c1->cd(5);
    pad5->SetTicks();
    pad5->SetLogz();
    pad5->SetLeftMargin(0.15);
    pad5->SetRightMargin(0.15);
    if(fdc_occ_cell_5){
        fdc_occ_cell_5->SetStats(0);
        fdc_occ_cell_5->Scale(1./Nevents);
        fdc_occ_cell_5->Draw("colz");
    }

    TVirtualPad *pad6 = c1->cd(6);
    pad6->SetTicks();
    pad6->SetLogz();
    pad6->SetLeftMargin(0.15);
    pad6->SetRightMargin(0.15);
    if(fdc_occ_cell_6){
        fdc_occ_cell_6->SetStats(0);
        fdc_occ_cell_6->Scale(1./Nevents);
        fdc_occ_cell_6->Draw("colz");
    }

}



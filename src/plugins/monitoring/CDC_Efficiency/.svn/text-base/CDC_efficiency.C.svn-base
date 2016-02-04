
{

    //dir->cd();
    TDirectory *dir = TFile::Open("hd_root.root");
    dir = (TDirectory*)gDirectory->FindObjectAny("cdc_efficiency");
    if(!dir) return;
    dir->cd();

    TCanvas *cEfficiency = new TCanvas("cEfficiency", "cEfficiency", 900, 800);
    // Draw axes
    TH2D *axes = (TH2D *)dir->Get("axes");
    if(!axes) axes = new TH2D("axes", "CDC Efficiency", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axes->SetStats(0);
    axes->Fill(100,100); // without this, the color ramp is not drawn
    axes->GetZaxis()->SetRangeUser(minScale, maxScale);
    axes->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiency->SaveAs("cEfficiency.png");

    TCanvas *cEfficiencyDOCA0 = new TCanvas("cEfficiencyDOCA0", "cEfficiencyDOCA0", 900, 800);
    // Draw axes
    TH2D *axesDOCA0 = (TH2D *)dir->Get("axesDOCA0");
    if(!axesDOCA0) axesDOCA0 = new TH2D("axesDOCA0", "CDC Efficiency for DOCA #in [0.0, 0.1 cm]", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axesDOCA0->SetStats(0);
    axesDOCA0->Fill(100,100); // without this, the color ramp is not drawn
    axesDOCA0->GetZaxis()->SetRangeUser(minScale, maxScale);
    axesDOCA0->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]DOCA0", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]DOCA0", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiencyDOCA0->SaveAs("cEfficiencyDOCA0.png");

    TCanvas *cEfficiencyDOCA1 = new TCanvas("cEfficiencyDOCA1", "cEfficiencyDOCA1", 900, 800);
    // Draw axes
    TH2D *axesDOCA1 = (TH2D *)dir->Get("axesDOCA1");
    if(!axesDOCA1) axesDOCA1 = new TH2D("axesDOCA1", "CDC Efficiency for DOCA #in [0.1, 0.2 cm]", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axesDOCA1->SetStats(0);
    axesDOCA1->Fill(100,100); // without this, the color ramp is not drawn
    axesDOCA1->GetZaxis()->SetRangeUser(minScale, maxScale);
    axesDOCA1->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]DOCA1", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]DOCA1", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiencyDOCA1->SaveAs("cEfficiencyDOCA1.png");

    TCanvas *cEfficiencyDOCA2 = new TCanvas("cEfficiencyDOCA2", "cEfficiencyDOCA2", 900, 800);
    // Draw axes
    TH2D *axesDOCA2 = (TH2D *)dir->Get("axesDOCA2");
    if(!axesDOCA2) axesDOCA2 = new TH2D("axesDOCA2", "CDC Efficiency for DOCA #in [0.2, 0.3 cm]", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axesDOCA2->SetStats(0);
    axesDOCA2->Fill(100,100); // without this, the color ramp is not drawn
    axesDOCA2->GetZaxis()->SetRangeUser(minScale, maxScale);
    axesDOCA2->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]DOCA2", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]DOCA2", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiencyDOCA2->SaveAs("cEfficiencyDOCA2.png");

    TCanvas *cEfficiencyDOCA3 = new TCanvas("cEfficiencyDOCA3", "cEfficiencyDOCA3", 900, 800);
    // Draw axes
    TH2D *axesDOCA3 = (TH2D *)dir->Get("axesDOCA3");
    if(!axesDOCA3) axesDOCA3 = new TH2D("axesDOCA3", "CDC Efficiency for DOCA #in [0.3, 0.4 cm]", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axesDOCA3->SetStats(0);
    axesDOCA3->Fill(100,100); // without this, the color ramp is not drawn
    axesDOCA3->GetZaxis()->SetRangeUser(minScale, maxScale);
    axesDOCA3->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]DOCA3", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]DOCA3", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiencyDOCA3->SaveAs("cEfficiencyDOCA3.png");

    TCanvas *cEfficiencyDOCA4 = new TCanvas("cEfficiencyDOCA4", "cEfficiencyDOCA4", 900, 800);
    // Draw axes
    TH2D *axesDOCA4 = (TH2D *)dir->Get("axesDOCA4");
    if(!axesDOCA4) axesDOCA4 = new TH2D("axesDOCA4", "CDC Efficiency for DOCA #in [0.4, 0.5 cm]", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axesDOCA4->SetStats(0);
    axesDOCA4->Fill(100,100); // without this, the color ramp is not drawn
    axesDOCA4->GetZaxis()->SetRangeUser(minScale, maxScale);
    axesDOCA4->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]DOCA4", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]DOCA4", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiencyDOCA4->SaveAs("cEfficiencyDOCA4.png");

    TCanvas *cEfficiencyDOCA5 = new TCanvas("cEfficiencyDOCA5", "cEfficiencyDOCA5", 900, 800);
    // Draw axes
    TH2D *axesDOCA5 = (TH2D *)dir->Get("axesDOCA5");
    if(!axesDOCA5) axesDOCA5 = new TH2D("axesDOCA5", "CDC Efficiency for DOCA #in [0.5, 0.6 cm]", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axesDOCA5->SetStats(0);
    axesDOCA5->Fill(100,100); // without this, the color ramp is not drawn
    axesDOCA5->GetZaxis()->SetRangeUser(minScale, maxScale);
    axesDOCA5->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]DOCA5", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]DOCA5", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiencyDOCA5->SaveAs("cEfficiencyDOCA5.png");

    TCanvas *cEfficiencyDOCA6 = new TCanvas("cEfficiencyDOCA6", "cEfficiencyDOCA6", 900, 800);
    // Draw axes
    TH2D *axesDOCA6 = (TH2D *)dir->Get("axesDOCA6");
    if(!axesDOCA6) axesDOCA6 = new TH2D("axesDOCA6", "CDC Efficiency for DOCA #in [0.6, 0.7 cm]", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axesDOCA6->SetStats(0);
    axesDOCA6->Fill(100,100); // without this, the color ramp is not drawn
    axesDOCA6->GetZaxis()->SetRangeUser(minScale, maxScale);
    axesDOCA6->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]DOCA6", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]DOCA6", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiencyDOCA6->SaveAs("cEfficiencyDOCA6.png");

    TCanvas *cEfficiencyDOCA7 = new TCanvas("cEfficiencyDOCA7", "cEfficiencyDOCA7", 900, 800);
    // Draw axes
    TH2D *axesDOCA7 = (TH2D *)dir->Get("axesDOCA7");
    if(!axesDOCA7) axesDOCA7 = new TH2D("axesDOCA7", "CDC Efficiency for DOCA #in [0.7, 0.78 cm]", 100, -65.0, 65.0, 100, -65.0, 65.0);

    Float_t minScale = 0.5; Float_t maxScale = 1.0;
    axesDOCA7->SetStats(0);
    axesDOCA7->Fill(100,100); // without this, the color ramp is not drawn
    axesDOCA7->GetZaxis()->SetRangeUser(minScale, maxScale);
    axesDOCA7->Draw("colz");

    for(unsigned int iring=1; iring<=28; iring++){
        char hname[256];
        sprintf(hname, "cdc_measured_ring[%d]DOCA7", iring);
        TH1 *h = (TH1*)(dir->Get(hname));
        char hname2[256];
        sprintf(hname2, "cdc_expected_ring[%d]DOCA7", iring);
        TH1 *h2 = (TH1*)(dir->Get(hname2));

        if(h && h2){
            h->Divide(h2);
            //sprintf(hname, "cdc_measured_ring[%d]", iring);
            //TH1 *hh = (TH1*)h->Clone(hname);
            //hh->Scale(1.0/Nevents);
            h->SetMinimum(minScale);
            h->SetMaximum(maxScale);
            h->SetStats(0);
            h->Draw("same col pol");  // draw remaining histos without overwriting color palette
        }
    }
    cEfficiencyDOCA7->SaveAs("cEfficiencyDOCA7.png");

    /*
       Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs Path Length",
       dx,
       "Measured Hits",
       100, 0 , 4.0);
       Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs DOCA",
       distanceToWire,
       "Measured Hits",
       100, 0 , 0.78);
       Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs Tracking FOM",
       TMath::Prob((*trackIter)->chisq, (*trackIter)->Ndof),
       "Measured Hits",
       100, 0 , 1.0);
       Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs theta",
       (*trackIter)->momentum().Theta()*TMath::RadToDeg(),
       "Measured Hits",
       100, 0, 180);
       Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs p",
       (*trackIter)->pmag(),
       "Measured Hits",
       100, 0 , 4.0);
       Fill2DHistogram("cdc_efficiency", "", "Measured hits p Vs Theta",
       (*trackIter)->momentum().Theta()*TMath::RadToDeg(), (*trackIter)->pmag(),
       "Measured Hits",
       100, 0, 180, 100, 0 , 4.0);
       */
    // Draw the other efficiency plots
    TCanvas *cPathLength = new TCanvas("cPathLength", "cPathLength", 800, 600);
    cPathLength->SetGridx();
    cPathLength->SetGridy();
    TH1I *MeasPathLength = (TH1I*)(dir->Get("Measured Hits Vs Path Length"));
    TH1I *ExpPathLength = (TH1I*)(dir->Get("Expected Hits Vs Path Length"));
    if(MeasPathLength && ExpPathLength){
        //EffPathLength->Draw();
        TGraphAsymmErrors *EffPathLength = new  TGraphAsymmErrors(MeasPathLength, ExpPathLength, "ac");
        EffPathLength->Draw("ap");
        EffPathLength->SetMinimum(0.0);
        EffPathLength->SetMaximum(1.0);
        EffPathLength->SetTitle("CDC Per Straw Efficiency Vs. dX");
        EffPathLength->GetXaxis()->SetTitle("Path Length Through Straw [cm]");
        EffPathLength->GetYaxis()->SetTitle("Efficiency");
        cPathLength->SaveAs("cPathLength.png");
    }

    TCanvas *cDOCA = new TCanvas("cDOCA", "cDOCA", 800, 600);
    cDOCA->SetGridx();
    cDOCA->SetGridy();
    TH1I *MeasDOCA = (TH1I*)(dir->Get("Measured Hits Vs DOCA"));
    TH1I *ExpDOCA = (TH1I*)(dir->Get("Expected Hits Vs DOCA"));
    if(MeasDOCA && ExpDOCA){
        //EffDOCA->Draw();
        TGraphAsymmErrors *EffDOCA = new  TGraphAsymmErrors(MeasDOCA, ExpDOCA, "ac");
        EffDOCA->Draw("ap");
        EffDOCA->SetMinimum(0.0);
        EffDOCA->SetMaximum(1.0);
        EffDOCA->SetTitle("CDC Per Straw Efficiency Vs. DOCA");
        EffDOCA->GetXaxis()->SetTitle("Closest distance between track and wire [cm]");
        EffDOCA->GetYaxis()->SetTitle("Efficiency");
        cDOCA->SaveAs("cDOCA.png");
    }

    TCanvas *cTrackingFOM = new TCanvas("cTrackingFOM", "cTrackingFOM", 800, 600);
    cTrackingFOM->SetGridx();
    cTrackingFOM->SetGridy();
    TH1I *MeasTrackingFOM = (TH1I*)(dir->Get("Measured Hits Vs Tracking FOM"));
    TH1I *ExpTrackingFOM = (TH1I*)(dir->Get("Expected Hits Vs Tracking FOM"));
    if(MeasTrackingFOM && ExpTrackingFOM){
        //EffTrackingFOM->Draw();
        TGraphAsymmErrors *EffTrackingFOM = new  TGraphAsymmErrors(MeasTrackingFOM, ExpTrackingFOM, "ac");
        EffTrackingFOM->Draw("ap");
        EffTrackingFOM->SetMinimum(0.0);
        EffTrackingFOM->SetMaximum(1.0);
        EffTrackingFOM->SetTitle("CDC Per Straw Efficiency Vs. Tracking FOM");
        EffTrackingFOM->GetXaxis()->SetTitle("Tracking FOM");
        EffTrackingFOM->GetYaxis()->SetTitle("Efficiency");
        cTrackingFOM->SaveAs("cTrackingFOM.png");
    }

    TCanvas *ctheta = new TCanvas("ctheta", "ctheta", 800, 600);
    ctheta->SetGridx();
    ctheta->SetGridy();
    TH1I *Meastheta = (TH1I*)(dir->Get("Measured Hits Vs theta"));
    TH1I *Exptheta = (TH1I*)(dir->Get("Expected Hits Vs theta"));
    if(Meastheta && Exptheta){
        //Efftheta->Draw();
        TGraphAsymmErrors *Efftheta = new  TGraphAsymmErrors(Meastheta, Exptheta, "ac");
        Efftheta->Draw("ap");
        Efftheta->SetMinimum(0.0);
        Efftheta->SetMaximum(1.0);
        Efftheta->SetTitle("CDC Per Straw Efficiency Vs. #theta");
        Efftheta->GetXaxis()->SetTitle("Track #theta [deg.]");
        Efftheta->GetYaxis()->SetTitle("Efficiency");
        ctheta->SaveAs("ctheta.png");
    }

    TCanvas *cp = new TCanvas("cp", "cp", 800, 600);
    cp->SetGridx();
    cp->SetGridy();
    TH1I *Measp = (TH1I*)(dir->Get("Measured Hits Vs p"));
    TH1I *Expp = (TH1I*)(dir->Get("Expected Hits Vs p"));
    if(Measp && Expp){
        //Effp->Draw();
        TGraphAsymmErrors *Effp = new  TGraphAsymmErrors(Measp, Expp, "ac");
        Effp->Draw("ap");
        Effp->SetMinimum(0.0);
        Effp->SetMaximum(1.0);
        Effp->SetTitle("CDC Per Straw Efficiency Vs. p");
        Effp->GetXaxis()->SetTitle("Track Momentum [GeV/c]");
        Effp->GetYaxis()->SetTitle("Efficiency");
        cp->SaveAs("cp.png");
    }


}

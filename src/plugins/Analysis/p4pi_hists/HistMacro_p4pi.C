{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locReactionDirectory;
	if((TDirectory*)locTopDirectory->FindObjectAny("p4pi_hists") != 0)
	  locReactionDirectory = (TDirectory*)locTopDirectory->FindObjectAny("p4pi_hists");
	else
	  return;

	TDirectory *locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Hist_MissingMassSquared");
	if(!locDirectory)
	  return;
	locDirectory->cd();
	TH1I* locHist_MM2 = (TH1I*)gROOT->FindObject("MissingMassSquared");
	
 	locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Hist_InvariantMass_FourPi");
	if(!locDirectory)
	  return;
	locDirectory->cd();
	TH1I* locHist_4Pi = (TH1I*)gROOT->FindObject("InvariantMass");

 	locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Hist_2DInvariantMass_TwoPi_vs_TwoPi");
	if(!locDirectory)
	  return;
	locDirectory->cd();
	TH2I* locHist_Dalitz = (TH2I*)gROOT->FindObject("2DInvariantMass");
	
 	locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Hist_InvariantMass_ProtonPip");
	if(!locDirectory)
	  return;
	locDirectory->cd();
	TH1I* locHist_ProtonPip = (TH1I*)gROOT->FindObject("InvariantMass");
	
 	locDirectory = (TDirectory*)locReactionDirectory->FindObjectAny("Hist_InvariantMass_ProtonPim");
	if(!locDirectory)
	  return;
	locDirectory->cd();
	TH1I* locHist_ProtonPim = (TH1I*)gROOT->FindObject("InvariantMass");
	
 	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
	  locCanvas = new TCanvas("p4pi_hists", "p4pi_hists", 1200, 800); //for testing
	else
	  locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);
	
	//Draw
	locCanvas->cd(1);
	if(locHist_MM2 != NULL) {	
	  locHist_MM2->GetXaxis()->SetRangeUser(-0.04,0.04);
	  locHist_MM2->Draw();
	}

	locCanvas->cd(2);
	if(locHist_4Pi != NULL) {	
	  locHist_4Pi->Draw();
	}
	
	locCanvas->cd(3);
	if(locHist_Dalitz != NULL) {	
	  locHist_Dalitz->Draw("colz");
	}
	
	locCanvas->cd(4);
	if(locHist_ProtonPip != NULL && locHist_ProtonPim != NULL) {	
	  locHist_ProtonPip->Draw();
	  locHist_ProtonPim->SetLineColor(kRed);
	  locHist_ProtonPim->Draw("same");
	}
	
}


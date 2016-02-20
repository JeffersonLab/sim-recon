// hnamepath: /TAGGER_online/TaggerEnergy_DeltaTSC

{
	TDirectory *locTopDirectory = gDirectory;
	TDirectory *locDirectory = (TDirectory*)locTopDirectory->FindObjectAny("TAGGER");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get diamond data for display
	TH2D* locHist_TaggerEnergy_DeltaTSC = (TH2D*)gDirectory->Get("TaggerEnergy_DeltaTSC");
	int lowTime = locHist_TaggerEnergy_DeltaTSC->GetXaxis()->FindBin(5.);
	int highTime = locHist_TaggerEnergy_DeltaTSC->GetXaxis()->FindBin(20.);
	int lowTimeAcc = locHist_TaggerEnergy_DeltaTSC->GetXaxis()->FindBin(-20.);
	int highTimeAcc = locHist_TaggerEnergy_DeltaTSC->GetXaxis()->FindBin(-5.);

	//Get amorphous data from reference file
	TFile *f = TFile::Open("/group/halld/Users/jrsteven/hist_amorph_reference.root");
	if(!f) return;
	TH2D* locHist_TaggerEnergy_DeltaTSC_Amorph = (TH2D*)f->Get("TAGGER/TaggerEnergy_DeltaTSC");
	if(!locHist_TaggerEnergy_DeltaTSC_Amorph) return;

  	TH1F* locHist_TaggerEnergy = (TH1F*)locHist_TaggerEnergy_DeltaTSC->ProjectionY("TaggerEnergy", lowTime, highTime);
	locHist_TaggerEnergy->SetTitle("Diamond");
	TH1F* locHist_TaggerEnergy_Amorph = (TH1F*)locHist_TaggerEnergy_DeltaTSC_Amorph->ProjectionY("TaggerEnergy_Amorph", lowTime, highTime);
	locHist_TaggerEnergy_Amorph->SetTitle("Amorphous Reference");
  	TH1F* locHist_TaggerEnergyAcc = (TH1F*)locHist_TaggerEnergy_DeltaTSC->ProjectionY("TaggerEnergyAcc", lowTimeAcc, highTimeAcc);
  	TH1F* locHist_TaggerEnergyAcc_Amorph = (TH1F*)locHist_TaggerEnergy_DeltaTSC_Amorph->ProjectionY("TaggerEnergyAcc_Amorph", lowTimeAcc, highTimeAcc);
 	locHist_TaggerEnergy->Rebin(); locHist_TaggerEnergy_Amorph->Rebin();
	locHist_TaggerEnergyAcc->Rebin(); locHist_TaggerEnergyAcc_Amorph->Rebin();  
	
	locHist_TaggerEnergy->Add(locHist_TaggerEnergyAcc, -1.);
	locHist_TaggerEnergy_Amorph->Add(locHist_TaggerEnergyAcc_Amorph, -1.);

	// make enhancment and scale to lowest energy bin
	TH1F* locHist_TaggerEnhancement = (TH1F*)locHist_TaggerEnergy->Clone();
	locHist_TaggerEnhancement->Divide(locHist_TaggerEnergy_Amorph); 
	Double_t scaleFactor = 1.;
	for(int i=0; i<locHist_TaggerEnhancement->GetXaxis()->GetNbins(); i++){
		if(locHist_TaggerEnhancement->GetBinContent(i)) {
			scaleFactor = locHist_TaggerEnhancement->GetBinContent(i);
			break;
		}
	}
	locHist_TaggerEnhancement->Scale(1./scaleFactor);
	locHist_TaggerEnhancement->SetTitle("Enhancement: Diamond/Amorphous");


	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TaggerEnhancement", "TaggerEnhancement", 1200, 500); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 1);

	//Draw
	locCanvas->cd(1);
	if(locHist_TaggerEnergy != NULL) {	
		locHist_TaggerEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_TaggerEnergy->GetYaxis()->SetTitleSize(0.05);
		locHist_TaggerEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_TaggerEnergy->GetYaxis()->SetLabelSize(0.05);
		locHist_TaggerEnergy->Draw();
	}

	locCanvas->cd(2);
	if(locHist_TaggerEnergy_Amorph != NULL) {         
                locHist_TaggerEnergy_Amorph->GetXaxis()->SetTitleSize(0.05);
                locHist_TaggerEnergy_Amorph->GetYaxis()->SetTitleSize(0.05);
                locHist_TaggerEnergy_Amorph->GetXaxis()->SetLabelSize(0.05);
                locHist_TaggerEnergy_Amorph->GetYaxis()->SetLabelSize(0.05);
                locHist_TaggerEnergy_Amorph->Draw();
        }

	locCanvas->cd(3);
	locHist_TaggerEnhancement->GetXaxis()->SetTitleSize(0.05);
        locHist_TaggerEnhancement->GetYaxis()->SetTitleSize(0.05);
        locHist_TaggerEnhancement->GetXaxis()->SetLabelSize(0.05);
        locHist_TaggerEnhancement->GetYaxis()->SetLabelSize(0.05);
	locHist_TaggerEnhancement->Draw();
}


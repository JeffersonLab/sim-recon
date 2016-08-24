

void plot_module_histos(int module, char filename[255])
{
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);

	TFile *_file0 = TFile::Open(filename);

	//	TDirectory *main = gDirectory;  // save current directory
    //Goto Path
    TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("bcalgainratio");
    if(!locDirectory)
        return;
    locDirectory->cd("channels");

	char name[255];
	int pad=0;

	TCanvas *histos = new TCanvas("histos","Fits",1000,1000);
	histos->Divide(4,4,0.001,0.001);
	TH2I *histograms[16];
	for (int layer=1; layer<=4; layer++) {
		for (int sector=1; sector<=4; sector++) {
			pad++;
			histos->cd(pad);
			sprintf(name,"logintratiovsZ_%02i%i%i",module,layer,sector);
			histograms[pad-1] = (TH2I*)gDirectory->Get(name);
			//			histograms[pad-1]->
			histograms[pad-1]->Draw();
		}
	}

	pad=0;
	TCanvas *Ehistos = new TCanvas("Ehistos","Fits",1000,1000);
	Ehistos->Divide(4,4,0.001,0.001);
	TH2I *Ehistograms[16];
	for (int layer=1; layer<=4; layer++) {
		for (int sector=1; sector<=4; sector++) {
			pad++;
			Ehistos->cd(pad);
			sprintf(name,"EvsZ_%02i%i%i",module,layer,sector);
			Ehistograms[pad-1] = (TH2I*)gDirectory->Get(name);
			Ehistograms[pad-1]->Rebin2D(4,4);
			//			histograms[pad-1]->
			Ehistograms[pad-1]->Draw("colz");
		}
	}

}


 


void plot_results(char filename[255], int module=1)
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
    locDirectory->cd();

    // load directories
    TH1I *hist_attenlength  = (TH1I *)gDirectory->Get("hist_attenlength");
    TH1I *hist_gainratio  = (TH1I *)gDirectory->Get("hist_gainratio");
    TH1I *hist_attenlength_err  = (TH1I *)gDirectory->Get("hist_attenlength_err");
    TH1I *hist_gainratio_err  = (TH1I *)gDirectory->Get("hist_gainratio_err");
    TH1I *hist_attenlength_relerr  = (TH1I *)gDirectory->Get("hist_attenlength_relerr");
    TH1I *hist_gainratio_relerr  = (TH1I *)gDirectory->Get("hist_gainratio_relerr");

    TH2I *hist2D_attenlength  = (TH2I *)gDirectory->Get("hist2D_attenlength");
    TH2I *hist2D_gainratio  = (TH2I *)gDirectory->Get("hist2D_gainratio");
    TH2I *EvsZ_layer1  = (TH2I *)gDirectory->Get("EvsZ_layer1");
    TH2I *EvsZ_layer2  = (TH2I *)gDirectory->Get("EvsZ_layer2");
    TH2I *EvsZ_layer3  = (TH2I *)gDirectory->Get("EvsZ_layer3");
    TH2I *EvsZ_layer4  = (TH2I *)gDirectory->Get("EvsZ_layer4");
    

	TCanvas *results = new TCanvas("results","Results of fit",800,800);
	results->Divide(2,2,0.001,0.001);

	results->cd(1);
	hist_attenlength->Draw();
	results->cd(2);
	hist_gainratio->Draw();
	results->cd(3);
	hist2D_attenlength->Draw("colz");
	results->cd(4);
	hist2D_gainratio->Draw("colz");	

	TCanvas *results_err = new TCanvas("results_err","Error of fit",800,800);
	results_err->Divide(2,2,0.001,0.001);

	results_err->cd(1);
	hist_attenlength_err->Draw();
	results_err->cd(2);
	hist_gainratio_err->Draw();
	results_err->cd(3);
	hist_attenlength_relerr->Draw();
	results_err->cd(4);
	hist_gainratio_relerr->Draw();


	TCanvas *layerE = new TCanvas("layerE","E vs Z",800,800);
	layerE->Divide(2,2,0.001,0.001);

	layerE->cd(1);
	EvsZ_layer1->GetYaxis()->SetRangeUser(0,0.02);
	EvsZ_layer1->Draw("colz");
	layerE->cd(2);
	EvsZ_layer2->GetYaxis()->SetRangeUser(0,0.02);
	EvsZ_layer2->Draw("colz");
	layerE->cd(3);
	EvsZ_layer3->GetYaxis()->SetRangeUser(0,0.02);
	EvsZ_layer3->Draw("colz");
	layerE->cd(4);
	EvsZ_layer4->GetYaxis()->SetRangeUser(0,0.02);
	EvsZ_layer4->Draw("colz");


	TCanvas *layerE_prof = new TCanvas("layerE_prof","E vs Z",800,800);
	layerE_prof->Divide(2,2,0.001,0.001);

	layerE_prof->cd(1);
	EvsZ_layer1->ProfileX()->Draw();
	layerE_prof->cd(2);
	EvsZ_layer2->ProfileX()->Draw();
	layerE_prof->cd(3);
	EvsZ_layer3->ProfileX()->Draw();
	layerE_prof->cd(4);
	EvsZ_layer4->ProfileX()->Draw();






}


 

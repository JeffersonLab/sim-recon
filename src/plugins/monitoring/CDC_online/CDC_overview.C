
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /CDC/cdc_raw_intpp
// hnamepath: /CDC/cdc_raw_t
// hnamepath: /CDC/cdc_ped
// hnamepath: /CDC/cdc_windata_ped
// hnamepath: /CDC/cdc_raw_intpp_vs_n
// hnamepath: /CDC/cdc_raw_t_vs_n
// hnamepath: /CDC/cdc_ped_vs_n
// hnamepath: /CDC/cdc_windata_ped_vs_n

{

  TDirectory *CDCdir = (TDirectory*)gDirectory->FindObjectAny("CDC");

  if(!CDCdir) return;
	
  CDCdir->cd();

  TH1D *cdc_raw_intpp   = (TH1D*)gDirectory->FindObjectAny("cdc_raw_intpp");
  TH1D *cdc_raw_t   = (TH1D*)gDirectory->FindObjectAny("cdc_raw_t");
  TH1I *cdc_ped = (TH1I*)gDirectory->FindObjectAny("cdc_ped");
  TH1I *cdc_windata_ped = (TH1I*)gDirectory->FindObjectAny("cdc_windata_ped");

  TH2D *cdc_raw_intpp_vs_n   = (TH2D*)gDirectory->FindObjectAny("cdc_raw_intpp_vs_n");
  TH2D *cdc_raw_t_vs_n   = (TH2D*)gDirectory->FindObjectAny("cdc_raw_t_vs_n");
  TH2I *cdc_ped_vs_n = (TH2I*)gDirectory->FindObjectAny("cdc_ped_vs_n");
  TH2I *cdc_windata_ped_vs_n = (TH2I*)gDirectory->FindObjectAny("cdc_windata_ped_vs_n");



	// Just for testing
	if(gPad == NULL){
		TCanvas *c1 = new TCanvas("c1");
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}
	
	if(!gPad) return;



	TCanvas *c1 = gPad->GetCanvas();

	gStyle->SetOptStat(1000000011);

	TGaxis::SetMaxDigits(4);

        c1->SetTitle("CDC overview");
	c1->Divide(4, 2, 0.001, 0.001);
	

	if (cdc_ped) {
  	  c1->cd(1);
          cdc_ped->SetLabelSize(0.025,"xyz");
          cdc_ped->SetTitleSize(0.03,"xy");
          cdc_ped->Draw();
	  gPad->Update();

          TPaveStats *st = (TPaveStats*)cdc_ped->FindObject("stats");
          st->SetOptStat(1000001111);
          st->SetY1NDC(0.85);
          gPad->GetCanvas()->Modified();
        }


	if(cdc_windata_ped) {
  	  c1->cd(2);
          cdc_windata_ped->SetLabelSize(0.025,"xyz");
          cdc_windata_ped->SetTitleSize(0.03,"xy");
	  cdc_windata_ped->Draw();
          gPad->Update();

          TPaveStats *st1 = (TPaveStats*)cdc_windata_ped->FindObject("stats");
          st1->SetOptStat(1000001111);
          st1->SetY1NDC(0.85);
          gPad->GetCanvas()->Modified();
        }


	if(cdc_raw_t) {
  	  c1->cd(3);
          cdc_raw_t->SetLabelSize(0.025,"xyz");
          cdc_raw_t->SetTitleSize(0.03,"xy");
          cdc_raw_t->Draw();
	}



	if(cdc_raw_intpp) {
  	  c1->cd(4);
          cdc_raw_intpp->SetLabelSize(0.025,"xyz");
          cdc_raw_intpp->SetTitleSize(0.03,"xy");
          cdc_raw_intpp->Draw();
        }



	if(cdc_ped_vs_n) {
  	  c1->cd(5);
          gPad->SetLogz();
          cdc_ped_vs_n->SetLabelSize(0.025,"xyz");
          cdc_ped_vs_n->SetTitleSize(0.03,"xy");
          cdc_ped_vs_n->SetTitleOffset(1.3,"y");
          cdc_ped_vs_n->Draw("colz");
        }


	if(cdc_windata_ped_vs_n) {
  	  c1->cd(6);
          gPad->SetLogz();
          cdc_windata_ped_vs_n->SetLabelSize(0.025,"xyz");
          cdc_windata_ped_vs_n->SetTitleSize(0.03,"xy");
          cdc_windata_ped_vs_n->SetTitleOffset(1.3,"y");
          cdc_windata_ped_vs_n->Draw("colz");
        }


	if(cdc_raw_t_vs_n) {
  	  c1->cd(7);
          gPad->SetLogz();
          cdc_raw_t_vs_n->SetLabelSize(0.025,"xyz");
          cdc_raw_t_vs_n->SetTitleSize(0.03,"xy");
          cdc_raw_t_vs_n->SetTitleOffset(1.3,"y");
          cdc_raw_t_vs_n->Draw("colz");
        }


	if(cdc_raw_intpp_vs_n) {
  	  c1->cd(8);
          gPad->SetLogz();
          cdc_raw_intpp_vs_n->SetLabelSize(0.025,"xyz");
          cdc_raw_intpp_vs_n->SetTitleSize(0.03,"xy");
          cdc_raw_intpp_vs_n->SetTitleOffset(1.3,"y");
	  cdc_raw_intpp_vs_n->Draw("colz");
        }


}

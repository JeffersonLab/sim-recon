
// File: bcal_hist_eff.C
// Created: 03/12/15
// Creator: Elton Smith
// Purpose: Display bcal efficiency plots
// Uses output of BCAL_Eff/DEventProcessor_BCAL_Eff.cc

// hnamepath: /bcal_eff/h1eff_eff
// hnamepath: /bcal_eff/h1eff_cellideff
// hnamepath: /bcal_eff/h1eff2_eff2
// hnamepath: /bcal_eff/h1eff2_cellideff2

// hnamepath: /bcal_eff/h1eff_layer
// hnamepath: /bcal_eff/h1eff_layertot
// hnamepath: /bcal_eff/h1eff2_layer
// hnamepath: /bcal_eff/h1eff2_layertot
// hnamepath: /bcal_eff/h1eff_cellid
// hnamepath: /bcal_eff/h1eff_cellidtot
// hnamepath: /bcal_eff/h1eff2_cellid
// hnamepath: /bcal_eff/h1eff2_cellidtot

{
	gStyle->SetPalette(1,0);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);

	// get histograms from bcal_eff subdirectory

   TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal_eff");
   if(!dir){
		cout << "Can't find bcal_eff TDirectory!" << endl;
		return;
	}

   TH1F *h1eff_eff         = (TH1F*)dir->Get("h1eff_eff"         );
   TH1F *h1eff_cellideff   = (TH1F*)dir->Get("h1eff_cellideff"   );
   TH1F *h1eff2_eff2       = (TH1F*)dir->Get("h1eff2_eff2"       );
   TH1F *h1eff2_cellideff2 = (TH1F*)dir->Get("h1eff2_cellideff2" );

   TH1F *h1eff_layer       = (TH1F*)dir->Get("h1eff_layer"       );
   TH1F *h1eff_layertot    = (TH1F*)dir->Get("h1eff_layertot"    );
   TH1F *h1eff2_layer      = (TH1F*)dir->Get("h1eff2_layer"      );
   TH1F *h1eff2_layertot   = (TH1F*)dir->Get("h1eff2_layertot"   );
   TH1F *h1eff_cellid      = (TH1F*)dir->Get("h1eff_cellid"      );
   TH1F *h1eff_cellidtot   = (TH1F*)dir->Get("h1eff_cellidtot"   );
   TH1F *h1eff2_cellid     = (TH1F*)dir->Get("h1eff2_cellid"     );
   TH1F *h1eff2_cellidtot  = (TH1F*)dir->Get("h1eff2_cellidtot"  );
	
	if( !h1eff_eff         ) { cout << "Can't find h1eff_eff!"         << endl; return; }
	if( !h1eff_cellideff   ) { cout << "Can't find h1eff_cellideff!"   << endl; return; }
	if( !h1eff2_eff2       ) { cout << "Can't find h1eff2_eff2!"       << endl; return; }
	if( !h1eff2_cellideff2 ) { cout << "Can't find h1eff2_cellideff2!" << endl; return; }

	if( !h1eff_layer       ) { cout << "Can't find h1eff_layer!"       << endl; return; }
	if( !h1eff_layertot    ) { cout << "Can't find h1eff_layertot!"    << endl; return; }
	if( !h1eff2_layer      ) { cout << "Can't find h1eff2_layer!"      << endl; return; }
	if( !h1eff2_layertot   ) { cout << "Can't find h1eff2_layertot!"   << endl; return; }
	if( !h1eff_cellid      ) { cout << "Can't find h1eff_cellid!"      << endl; return; }
	if( !h1eff_cellidtot   ) { cout << "Can't find h1eff_cellidtot!"   << endl; return; }
	if( !h1eff2_cellid     ) { cout << "Can't find h1eff2_cellid!"     << endl; return; }
	if( !h1eff2_cellidtot  ) { cout << "Can't find h1eff2_cellidtot!"  << endl; return; }

	// The efficiency hist definitions are grabbed from the remote
	// process, but we calculate the ratio here using the cumulative
	// histos. Note that the BCAL_Eff plugin does this in the fini
	// method, but not during filling so the efficiency histos we grab
	// for display in ROOTSpy are empty anyway.
	h1eff_eff->Divide(h1eff_layer,h1eff_layertot,1,1,"B");
	h1eff2_eff2->Divide(h1eff2_layer,h1eff2_layertot,1,1,"B");
	h1eff_cellideff->Divide(h1eff_cellid,h1eff_cellidtot,1,1,"B");
	h1eff2_cellideff2->Divide(h1eff2_cellid,h1eff2_cellidtot,1,1,"B");

   //
   if(gPad == NULL){
	  TCanvas *c1 = new TCanvas("c1","c1 bcal_hist_eff ",200,10,700,700);
	  c1->cd(0);
	  c1->Draw();
	  c1->Update();
	}
	
   if(!gPad) return;
   TCanvas *c1 = gPad->GetCanvas();

   c1->SetBorderMode(0);
   c1->SetFillColor(0);

   Double_t ymin=0.5;
   // Double_t ymin=0.;
   Double_t ymax=1.1;

   c1->SetGridx();
   c1->SetGridy();
   c1->SetBorderMode(0);
   c1->SetFillColor(0);

   c1->Divide(2,2);
   c1->cd(1);

   h1eff_eff->SetTitle("");
   //h1eff_eff->GetXaxis()->SetRangeUser(xmin,xmax);
   h1eff_eff->GetYaxis()->SetRangeUser(ymin,ymax);
   h1eff_eff->GetXaxis()->SetTitleSize(0.05);
   h1eff_eff->GetYaxis()->SetTitleSize(0.05);
   //h1eff_eff->GetXaxis()->SetTitle("Layer");
   h1eff_eff->GetYaxis()->SetTitle("Efficiency Clusters");
   h1eff_eff->SetLineColor(2);
   h1eff_eff->Draw("");


   c1->cd(2);

   h1eff_cellideff->SetTitle("");
   // h1eff_cellideff->GetXaxis()->SetRangeUser(xmin,xmax);
   h1eff_cellideff->GetYaxis()->SetRangeUser(ymin,ymax);
   h1eff_cellideff->GetXaxis()->SetTitleSize(0.05);
   h1eff_cellideff->GetYaxis()->SetTitleSize(0.05);
   //h1eff_cellideff->GetXaxis()->SetTitle("Layer");
   h1eff_cellideff->GetYaxis()->SetTitle("Efficiency Clusters");
   h1eff_cellideff->SetLineColor(2);
   h1eff_cellideff->Draw("");


   c1->cd(3);

   h1eff2_eff2->SetTitle("");
   // h1eff2_eff2->GetXaxis()->SetRangeUser(xmin,xmax);
   h1eff2_eff2->GetYaxis()->SetRangeUser(ymin,ymax);
   h1eff2_eff2->GetXaxis()->SetTitleSize(0.05);
   h1eff2_eff2->GetYaxis()->SetTitleSize(0.05);
   //h1eff2_eff2->GetXaxis()->SetTitle("Layer");
   h1eff2_eff2->GetYaxis()->SetTitle("Efficiency Enhanced Hits");
   h1eff2_eff2->SetLineColor(2);
   h1eff2_eff2->Draw("");


   c1->cd(4);

   h1eff2_cellideff2->SetTitle("");
   // h1eff2_cellideff2->GetXaxis()->SetRangeUser(xmin,xmax);
   h1eff2_cellideff2->GetYaxis()->SetRangeUser(ymin,ymax);
   h1eff2_cellideff2->GetXaxis()->SetTitleSize(0.05);
   h1eff2_cellideff2->GetYaxis()->SetTitleSize(0.05);
   //h1eff2_cellideff2->GetXaxis()->SetTitle("Layer");
   h1eff2_cellideff2->GetYaxis()->SetTitle("Efficiency Enhanced Hits");
   h1eff2_cellideff2->SetLineColor(2);
   h1eff2_cellideff2->Draw("");

}



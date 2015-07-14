void bcal_hist_eff_summary(void)
{
// read in histograms plotting BCAL efficiencies and print history to file.
//

#include <TRandom.h>
#include <fstream.h>

gROOT->Reset();
//TTree *Bfield = (TTree *) gROOT->FindObject("Bfield");
gStyle->SetPalette(1,0);
gStyle->SetOptStat(kFALSE);
// gStyle->SetOptStat(11111111);
gStyle->SetPadRightMargin(0.15);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
//

   char string[256];
   Int_t j,jj;
   Int_t ndx;
   // Int_t npts=7;   // run 1506
   // Int_t npts=42;   // run 2400
   // Int_t npts=423;    // run 2439
   // Int_t npts=205;    // run 2439
   Int_t npts=18;    // run 2179

   // define history histogram

   Double_t xmin=0;
   Double_t xmax=npts;
   //Double_t ymin=0.5;
   // Double_t ymax=1.1;
   Double_t ymin=0.;
   Double_t ymax=1.2;
  TH1F *hsummary_layer1 = new TH1F("hsummary_layer1","Summary layer 1",npts,xmin,xmax);
  TH1F *hsummary_layer2 = new TH1F("hsummary_layer2","Summary layer 2",npts,xmin,xmax);
  TH1F *hsummary_layer3 = new TH1F("hsummary_layer3","Summary layer 3",npts,xmin,xmax);
  TH1F *hsummary_layer4 = new TH1F("hsummary_layer4","Summary layer 4",npts,xmin,xmax);

  TH1F *hsummary2_layer1 = new TH1F("hsummary2_layer1","Summary layer 1",npts,xmin,xmax);
  TH1F *hsummary2_layer2 = new TH1F("hsummary2_layer2","Summary layer 2",npts,xmin,xmax);
  TH1F *hsummary2_layer3 = new TH1F("hsummary2_layer3","Summary layer 3",npts,xmin,xmax);
  TH1F *hsummary2_layer4 = new TH1F("hsummary2_layer4","Summary layer 4",npts,xmin,xmax);

   // read histograms from file

   for (ndx=0; ndx < npts; ndx++) {

   sprintf(string,"hd_rawdata_002179_%03d_.root",ndx);


   // if (ndx == 1) continue;  // run 1506
   /*if (ndx == 2) continue;
   if (ndx == 6) continue;
   if (ndx == 12) continue;;
   if (ndx == 39) continue;*/ 

   /*if (ndx == 10) continue;     // file #3 is empty
   if (ndx == 92) continue;     // file #3 is empty
   if (ndx == 95) continue;     // file #3 is empty
   if (ndx == 141) continue;     // file #3 is empty
   if (ndx == 176) continue;     // file #3 is empty
   if (ndx == 225) continue;     // file #3 is empty
   if (ndx == 226) continue;     // file #3 is empty
   if (ndx == 227) continue;     // file #3 is empty
   if (ndx >= 230 && ndx<=254 ) continue;     // file #3 is empty
   if (ndx == 257) continue;     // file #3 is empty
   if (ndx == 258) continue;     // file #3 is empty
   if (ndx == 260) continue;     // file #3 is empty
   if (ndx == 262) continue;     // file #3 is empty
   if (ndx == 264) continue;     // file #3 is empty
   if (ndx == 278) continue;     // file #3 is empty
   if (ndx == 281) continue;     // file #3 is empty
   if (ndx == 284) continue;     // file #3 is empty
   if (ndx == 285) continue;     // file #3 is empty
   if (ndx == 286) continue;     // file #3 is empty
   if (ndx == 287) continue;     // file #3 is empty
   if (ndx == 289) continue;     // file #3 is empty
   if (ndx == 290) continue;     // file #3 is empty
   if (ndx == 291) continue;     // file #3 is empty
   if (ndx == 293) continue;     // file #3 is empty
   if (ndx == 294) continue;     // file #3 is empty
   if (ndx == 295) continue;     // file #3 is empty
   if (ndx == 298) continue;     // file #3 is empty
   if (ndx == 299) continue;     // file #3 is empty
   if (ndx == 304) continue;     // file #3 is empty
   if (ndx == 306) continue;     // file #3 is empty*/

   TFile *in=NULL;
   in = new TFile(string,"read");
   if (!in) {
     printf ("Histogram input filename=%s does not open\n",string);
     continue;
   }
   else if {
     printf ("Histogram input filename=%s OK\n",string);
   }

   // cd into bcal_eff directory
   TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal_eff");
   if(dir) dir->cd();

   TH1F *h1eff_eff= NULL;
   h1eff_eff = (TH1F*)gDirectory->FindObjectAny("h1eff_eff");
   TH1F *h1eff_cellideff = (TH1F*))gDirectory->FindObjectAny("h1eff_cellideff");
   TH1F *h1eff2_eff2 = (TH1F*))gDirectory->FindObjectAny("h1eff2_eff2");
   TH1F *h1eff2_cellideff2 = (TH1F*))gDirectory->FindObjectAny("h1eff2_cellideff2");

   if (!h1eff_eff == NULL !! h1eff_eff->GetEntries() <=0) continue;

   // retrieving information from the histogram

   Int_t ndim = h1eff_eff->GetNbinsX();
   Double_t xlo = h1eff_eff->GetBinLowEdge(1);
   Double_t width = h1eff_eff->GetBinWidth(1);
   Double_t xhi = xlo + width*ndim;
   printf ("\nbcal_hist_eff_summary: ndx=%d, ndim=%d, xlo=%f, width=%f, xhi=%d\n\n",ndx,ndim,xlo,width,xhi);

      for(j=0;j<ndim;j++) {
      Double_t xbin = xlo + j*width;
      Double_t content = h1eff_eff->GetBinContent(j+1);
      Double_t error = h1eff_eff->GetBinError(j+1);
      printf ("j=%d, xbin=%f, content=%f\n",j,xbin,content);

     // fill summary histogram

          if (j == 1) {
	hsummary_layer1->SetBinContent(ndx+1,content);
                    hsummary_layer1->SetBinError(ndx+1,error);
	}
         if (j == 2) {
	hsummary_layer2->SetBinContent(ndx+1,content);
                    hsummary_layer2->SetBinError(ndx+1,error);
	}
         if (j == 3) {
	hsummary_layer3->SetBinContent(ndx+1,content);
                    hsummary_layer3->SetBinError(ndx+1,error);
	}
         if (j == 4) {
	hsummary_layer4->SetBinContent(ndx+1,content);
                    hsummary_layer4->SetBinError(ndx+1,error);
	}
    
   }



   // retrieving information from the histogram

   Int_t ndim = h1eff2_eff2->GetNbinsX();
   Double_t xlo = h1eff2_eff2->GetBinLowEdge(1);
   Double_t width = h1eff2_eff2->GetBinWidth(1);
   Double_t xhi = xlo + width*ndim;
   printf ("\nbcal_hist_eff2_summary: ndx=%d, ndim=%d, xlo=%f, width=%f, xhi=%d\n\n",ndx,ndim,xlo,width,xhi);

      for(j=0;j<ndim;j++) {
      Double_t xbin = xlo + j*width;
      Double_t content = h1eff2_eff2->GetBinContent(j+1);
      Double_t error = h1eff2_eff2->GetBinError(j+1);
      printf ("j=%d, xbin=%f, content=%f\n",j,xbin,content);

     // fill summary histogram

          if (j == 1) {
	hsummary2_layer1->SetBinContent(ndx+1,content);
                    hsummary2_layer1->SetBinError(ndx+1,error);
	}
         if (j == 2) {
	hsummary2_layer2->SetBinContent(ndx+1,content);
                    hsummary2_layer2->SetBinError(ndx+1,error);
	}
         if (j == 3) {
	hsummary2_layer3->SetBinContent(ndx+1,content);
                    hsummary2_layer3->SetBinError(ndx+1,error);
	}
         if (j == 4) {
	hsummary2_layer4->SetBinContent(ndx+1,content);
                    hsummary2_layer4->SetBinError(ndx+1,error);
	}
    
   }

      in->Close();
   
   }   // end loop over histogram files

   //
   TCanvas *c1 = new TCanvas("c1","c1 bcal_hist_eff_summary ",200,10,700,700);
   c1->SetBorderMode(0);
   c1->SetFillColor(0);

   c1->SetGridx();
   c1->SetGridy();
   c1->SetBorderMode(0);
   c1->SetFillColor(0);

  TLegend *leg = new TLegend(0.6,0.8,0.85,0.95);
  leg->AddEntry(hsummary_layer1,"Layer 1","p");
  leg->AddEntry(hsummary_layer2,"Layer 2","p");
  leg->AddEntry(hsummary_layer3,"Layer 3","p");
  leg->AddEntry(hsummary_layer4,"Layer 4","p");

   hsummary_layer1->SetTitle("");
   //hsummary_layer1->GetXaxis()->SetRangeUser(xmin,xmax);
   hsummary_layer1->GetYaxis()->SetRangeUser(ymin,ymax);
   hsummary_layer1->GetXaxis()->SetTitleSize(0.05);
   hsummary_layer1->GetYaxis()->SetTitleSize(0.05);
   hsummary_layer1->GetXaxis()->SetTitle("File Number");
   hsummary_layer1->GetYaxis()->SetTitle("Efficiency Clusters");
   hsummary_layer1->SetLineColor(2);
   hsummary_layer1->SetMarkerColor(2);
   hsummary_layer1->SetMarkerStyle(20);
   hsummary_layer1->Draw("p");

   hsummary_layer2->SetLineColor(4);
   hsummary_layer2->SetMarkerColor(4);
   hsummary_layer2->SetMarkerStyle(20);
   hsummary_layer2->Draw("samep");

   hsummary_layer3->SetLineColor(1);
   hsummary_layer3->SetMarkerColor(1);
   hsummary_layer3->SetMarkerStyle(20);
   hsummary_layer3->Draw("samep");

   hsummary_layer4->SetLineColor(3);
   hsummary_layer4->SetMarkerColor(3);
   hsummary_layer4->SetMarkerStyle(20);
   hsummary_layer4->Draw("samep");

   sprintf(string,"hd_rawdata_002179");
   printf ("Histogram input filename=%s\n",string);
   t1 = new TLatex(0.15,0.92,string);
   t1->SetNDC();
   t1->SetTextSize(0.03);
   t1->Draw();

  leg->Draw();

   //
   TCanvas *c2 = new TCanvas("c2","c2 bcal_hist_eff_summary ",200,10,700,700);
   c2->SetBorderMode(0);
   c2->SetFillColor(0);

   c2->SetGridx();
   c2->SetGridy();
   c2->SetBorderMode(0);
   c2->SetFillColor(0);

  TLegend *leg = new TLegend(0.6,0.80,0.85,0.95);
  leg->AddEntry(hsummary2_layer1,"Layer 1","p");
  leg->AddEntry(hsummary2_layer2,"Layer 2","p");
  leg->AddEntry(hsummary2_layer3,"Layer 3","p");
  leg->AddEntry(hsummary2_layer4,"Layer 4","p");

   hsummary2_layer1->SetTitle("");
   //hsummary2_layer1->GetXaxis()->SetRangeUser(xmin,xmax);
   hsummary2_layer1->GetYaxis()->SetRangeUser(ymin,ymax);
   hsummary2_layer1->GetXaxis()->SetTitleSize(0.05);
   hsummary2_layer1->GetYaxis()->SetTitleSize(0.05);
   hsummary2_layer1->GetXaxis()->SetTitle("File Number");
   hsummary2_layer1->GetYaxis()->SetTitle("Efficiency Enhanced Hits");
   hsummary2_layer1->SetLineColor(2);
   hsummary2_layer1->SetMarkerColor(2);
   hsummary2_layer1->SetMarkerStyle(20);
   hsummary2_layer1->Draw("p");

   hsummary2_layer2->SetLineColor(4);
   hsummary2_layer2->SetMarkerColor(4);
   hsummary2_layer2->SetMarkerStyle(20);
   hsummary2_layer2->Draw("samep");

   hsummary2_layer3->SetLineColor(1);
   hsummary2_layer3->SetMarkerColor(1);
   hsummary2_layer3->SetMarkerStyle(20);
   hsummary2_layer3->Draw("samep");

   hsummary2_layer4->SetLineColor(3);
   hsummary2_layer4->SetMarkerColor(3);
   hsummary2_layer4->SetMarkerStyle(20);
   hsummary2_layer4->Draw("samep");

   sprintf(string,"hd_rawdata_002179");
   printf ("Histogram input filename=%s\n",string);
   t1 = new TLatex(0.15,0.92,string);
   t1->SetNDC();
   t1->SetTextSize(0.03);
   t1->Draw();

  leg->Draw();



   sprintf(string,"bcal_hist_eff_summary_002179.pdf(");
   c1->SaveAs(string);

   sprintf(string,"bcal_hist_eff_summary_002179.pdf)");
   c2->SaveAs(string);

}


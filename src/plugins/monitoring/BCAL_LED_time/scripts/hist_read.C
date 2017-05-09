void hist_read(TString filename)
{
// Read histogram file and output offset constants to file
//

#include <TRandom.h>

gROOT->Reset();
//TTree *Bfield = (TTree *) gROOT->FindObject("Bfield");
gStyle->SetPalette(1,0);
gStyle->SetOptStat(kFALSE);
gStyle->SetOptStat(11111111);
gStyle->SetPadRightMargin(0.15);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
//

   char string[256];
   Int_t const nids=767;
   Double_t const rms0=1;
   Double_t const offset_down=23+10;
   Double_t const offset_up=-23+10;
   Int_t fail0=0;
   Int_t fail1=0;

   // define histograms
   TH1F *sigma = new TH1F ("sigma","Sigma",100,0,5);

   // read histograms from file

   // TString filename = "hd_rawdata_030494_fix";
   // TString filename = "hd_rawdata_030964_005";
   TString title = filename;

  // open file for output

   TString outfile = "dat/"+filename+".dat";
   cout << "Opening text output file: " << outfile.Data() << endl;
   TString outfile2 = "dat/"+filename+".dat2";
   cout << "Opening constant output file: " << outfile.Data() << endl;

   ofstream filedat;
   filedat.open (outfile.Data(), ios::out);
   ofstream filedat2;
   filedat2.open (outfile2.Data(), ios::out);

   TFile *in = new TFile(filename+".root","read");
        
   TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("BCAL_LED_time");
   if(dir) {
     dir->cd();
   }
   else {
     cout << "*** Cannot find directory" << endl;
   }

   TH1F *low_bias_down_time_vchannel  = (TH1F*)gDirectory->FindObjectAny("low_bias_down_time_vchannel");
   TH1F *low_bias_up_time_vchannel  = (TH1F*)gDirectory->FindObjectAny("low_bias_up_time_vchannel");
   TH1F *high_bias_down_time_vchannel  = (TH1F*)gDirectory->FindObjectAny("high_bias_down_time_vchannel");
   TH1F *high_bias_up_time_vchannel  = (TH1F*)gDirectory->FindObjectAny("high_bias_up_time_vchannel");

   // retrieving information from the histogram

   Int_t ndim = low_bias_down_time_vchannel->GetNbinsX();
   Double_t xlo = low_bias_down_time_vchannel->GetBinLowEdge(1);
   Double_t width = low_bias_down_time_vchannel->GetBinWidth(1);
   printf ("\nhist_read: ndim=%d, xlo=%f, width=%f\n\n",ndim,xlo,width);

   for(Int_t j=0;j<ndim;j++) {

      Int_t chid = low_bias_down_time_vchannel->GetBinCenter(j+1);

      // compute module, layer, sector

      Int_t module = chid/16 + 1;
      Int_t layer = (chid - (module-1)*16)/4 + 1;
      Int_t sector = (chid - (module-1)*16 - (layer-1)*4) +1; 
      // cout << " chid=" << chid << " module=" << module << " layer=" << layer << " sector=" << sector << endl;
      Double_t tdiff_low_down = low_bias_down_time_vchannel->GetBinContent(j+1);
      Double_t tdiff_low_up = low_bias_up_time_vchannel->GetBinContent(j+1);
      Double_t tdiff_high_down = high_bias_down_time_vchannel->GetBinContent(j+1);
      Double_t tdiff_high_up = high_bias_up_time_vchannel->GetBinContent(j+1);

      // subtract offset only if non-zero data are present

      tdiff_low_down = tdiff_low_down != 0? tdiff_low_down + offset_down : 0;
      tdiff_low_up =  tdiff_low_up != 0? tdiff_low_up + offset_up : 0;
      tdiff_high_down = tdiff_high_down != 0? tdiff_high_down  + offset_down : 0;
      tdiff_high_up = tdiff_high_up != 0? tdiff_high_up  + offset_up : 0;

      // compute rms
      // Double_t sum = abs(tdiff_low_down) + abs(tdiff_low_up) + abs(tdiff_high_down) + abs(tdiff_high_up);
      Double_t sum = tdiff_low_down + tdiff_low_up + tdiff_high_down + tdiff_high_up;
      Double_t sum2 = tdiff_low_down*tdiff_low_down + tdiff_low_up*tdiff_low_up + tdiff_high_down*tdiff_high_down + tdiff_high_up*tdiff_high_up;
      Double_t sum1 = 0;
      sum1 = abs(tdiff_low_down)>0? sum1+1 : sum1;
      sum1 = abs(tdiff_low_up)>0? sum1+1 : sum1;
      sum1 = abs(tdiff_high_down)>0? sum1+1 : sum1;
      sum1 = abs(tdiff_high_up)>0? sum1+1 : sum1;
      Double_t rms = sum1 > 0? sqrt(sum2/sum1 - sum*sum/sum1/sum1): 0;
      sigma->Fill(rms);

      sprintf (string,"chid low_down  low_up high_down  high_up rms %4d  %4d  %4d  %4d %10.2f %10.2f %10.2f %10.2f %10.4f",chid,module,layer,sector,tdiff_low_down,tdiff_low_up,tdiff_high_down,tdiff_high_up,rms);
      if (rms > rms0 | rms ==0) {
	// cout << "***hist_read: chid=" << chid << " rms=" << rms << endl;
	printf ("%s   ***\n",string);
	if (rms > rms0) fail1++;
	if (rms == 0) fail0++;
	filedat << string << "   ***" <<  endl;
      }
      else {
	filedat << string <<  endl;
      }
      // output average value to constant file. Add one line for each end:
      Double_t average = sum1 > 0? sum/sum1 : 0;
      filedat2 << chid*2 << "\t" << 0 << endl;
      filedat2 << chid*2 + 1 << "\t" << average << endl;
      
   }

   cout << "*** Fails: N(rms==0):" << fail0 << ", N(rms > rms0):" << fail1 << endl;
   filedat << "*** Fails: N(rms==0):" << fail0 << ", N(rms > rms0):" << fail1 << endl;
     
   //
   TCanvas *c0 = new TCanvas("c0","c0 hist_read",200,10,1200,700);
   c0->Divide(3,2);

   c0->cd(1);

   low_bias_down_time_vchannel->SetTitle(title);
   // low_bias_down_time_vchannel->GetXaxis()->SetRangeUser(xmin,xmax);
   // low_bias_down_time_vchannel->GetYaxis()->SetRangeUser(ymin,ymax);
   low_bias_down_time_vchannel->GetXaxis()->SetTitleSize(0.05);
   low_bias_down_time_vchannel->GetYaxis()->SetTitleSize(0.05);
   low_bias_down_time_vchannel->GetXaxis()->SetTitle("Low Down - chid");
   low_bias_down_time_vchannel->GetYaxis()->SetTitle("Time Difference (Down-Up) (ns)");
   low_bias_down_time_vchannel->SetMarkerColor(2);
   low_bias_down_time_vchannel->SetMarkerStyle(20);
   low_bias_down_time_vchannel->SetMarkerSize(0.5);
   low_bias_down_time_vchannel->Draw("");

   /*sprintf(string,"CMS Phys Lett B 716 (2012) 30");
   printf("string=%s\n",string);
   t1 = new TLatex(0.15,0.92,string);
   t1->SetNDC();
   t1->SetTextSize(0.03);
   t1->Draw();*/

   c0->cd(2);

   low_bias_up_time_vchannel->SetTitle("");
   // low_bias_up_time_vchannel->GetXaxis()->SetRangeUser(xmin,xmax);
   // low_bias_up_time_vchannel->GetYaxis()->SetRangeUser(ymin,ymax);
   low_bias_up_time_vchannel->GetXaxis()->SetTitleSize(0.05);
   low_bias_up_time_vchannel->GetYaxis()->SetTitleSize(0.05);
   low_bias_up_time_vchannel->GetXaxis()->SetTitle("Low Up - chid");
   low_bias_up_time_vchannel->GetYaxis()->SetTitle("Time Difference (Down-Up) (ns)");
   low_bias_up_time_vchannel->SetMarkerColor(2);
   low_bias_up_time_vchannel->SetMarkerStyle(20);
   low_bias_up_time_vchannel->SetMarkerSize(0.5);
   low_bias_up_time_vchannel->Draw("");

   c0->cd(3);

   high_bias_down_time_vchannel->SetTitle("");
   // high_bias_down_time_vchannel->GetXaxis()->SetRangeUser(xmin,xmax);
   // high_bias_down_time_vchannel->GetYaxis()->SetRangeUser(ymin,ymax);
   high_bias_down_time_vchannel->GetXaxis()->SetTitleSize(0.05);
   high_bias_down_time_vchannel->GetYaxis()->SetTitleSize(0.05);
   high_bias_down_time_vchannel->GetXaxis()->SetTitle("High Down - chid");
   high_bias_down_time_vchannel->GetYaxis()->SetTitle("Time Difference (Down-Up) (ns)");
   high_bias_down_time_vchannel->SetMarkerColor(2);
   high_bias_down_time_vchannel->SetMarkerStyle(20);
   high_bias_down_time_vchannel->SetMarkerSize(0.5);
   high_bias_down_time_vchannel->Draw("");

   c0->cd(4);

   high_bias_up_time_vchannel->SetTitle("");
   // high_bias_up_time_vchannel->GetXaxis()->SetRangeUser(xmin,xmax);
   // high_bias_up_time_vchannel->GetYaxis()->SetRangeUser(ymin,ymax);
   high_bias_up_time_vchannel->GetXaxis()->SetTitleSize(0.05);
   high_bias_up_time_vchannel->GetYaxis()->SetTitleSize(0.05);
   high_bias_up_time_vchannel->GetXaxis()->SetTitle("High Up - chid");
   high_bias_up_time_vchannel->GetYaxis()->SetTitle("Time Difference (Down-Up) (ns)");
   high_bias_up_time_vchannel->SetMarkerColor(2);
   high_bias_up_time_vchannel->SetMarkerStyle(20);
   high_bias_up_time_vchannel->SetMarkerSize(0.5);
   high_bias_up_time_vchannel->Draw("");

   c0->cd(5);
   gPad->SetLogy();

   sigma->SetTitle("");
   // sigma->GetXaxis()->SetRangeUser(xmin,xmax);
   // sigma->GetYaxis()->SetRangeUser(ymin,ymax);
   sigma->GetXaxis()->SetTitleSize(0.05);
   sigma->GetYaxis()->SetTitleSize(0.05);
   sigma->GetYaxis()->SetTitle("Entries");
   sigma->GetXaxis()->SetTitle("Sigma between configurations");
   sigma->SetLineColor(2);
   sigma->Draw("");


   c0->SaveAs("plots/"+filename+"_plots.pdf");

   filedat.close();   // close text file
   filedat2.close();   // close text file


   // in->Close();

}


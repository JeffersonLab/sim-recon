// File: ST_Monitoring_Expert_waveform.C
// Created: 05/18/2015
// Creator: Mahmoud Kamel, mkame006@fiu.edu
// Purpose: Displaying waveforms histograms for online monitoring purposes.  Designed for ST Experts.

{
 // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_lowlevel/waveforms");
   if(dir) dir->cd();

  // Declare the number of channels 
  const int NCHANNELS = 30;

  // Define arrays of dynamic pointers to 1D histograms
  TH1I** h_amp_vs_sampl_chan = new TH1I*[NCHANNELS];
  TH1I** h_amp_vs_sampl_chan150 = new TH1I*[NCHANNELS];
  TH1I** h_amp_vs_sampl_chan1000 = new TH1I*[NCHANNELS];
  TH1I** h_amp_vs_sampl_chan2000 = new TH1I*[NCHANNELS];
  TH1I** h_amp_vs_sampl_chan3000 = new TH1I*[NCHANNELS];
  TH1I** h_amp_vs_sampl_chan4000 = new TH1I*[NCHANNELS];
  for(unsigned int i = 0; i < NCHANNELS; i++)
    {
        // Grab 1D histograms for root file 
        h_amp_vs_sampl_chan[i]       = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan_%i", i+1));
        h_amp_vs_sampl_chan150[i]    = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan150_%i", i+1));
        h_amp_vs_sampl_chan1000[i]   = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan1000_%i", i+1));
        h_amp_vs_sampl_chan2000[i]   = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan2000_%i", i+1));
        h_amp_vs_sampl_chan3000[i]   = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan3000_%i", i+1));
        h_amp_vs_sampl_chan4000[i]   = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan4000_%i", i+1));
    }


  // Create the canvas c1
  if(gPad == NULL)
    {
      TCanvas *c1 = new TCanvas("c1","Start Counter Expert Waveform Histograms( 100 < Pulse height <= 150)", 200, 10, 600, 480);
      c1->cd(0);
      c1->Draw();
      c1->Update();
    }
  
  if(!gPad) return;
  TCanvas *c1 = gPad->GetCanvas();
  c1->Divide(3,2);

      // f1ADC250 multiplicity histograms
      c1->cd(1);
      gPad->SetTicks();
      gStyle->SetOptStat(10);
	        
      //set margins
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.03);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.1);
      //set axis label size
      h_amp_vs_sampl_chan150[3]->SetLabelSize(0.04,"X");
      h_amp_vs_sampl_chan150[3]->SetLabelSize(0.04,"Y");
      h_amp_vs_sampl_chan150[3]->GetXaxis()->CenterTitle();
      h_amp_vs_sampl_chan150[3]->GetYaxis()->CenterTitle();
      h_amp_vs_sampl_chan150[3]->SetTitleSize(0.06,"X");
      h_amp_vs_sampl_chan150[3]->SetTitleSize(0.06,"Y");
      
      //set marker,&  marker size and color
      h_amp_vs_sampl_chan150[3]->SetMarkerStyle(8);
      h_amp_vs_sampl_chan150[3]->SetMarkerColor(4);
      h_amp_vs_sampl_chan150[3]->SetMarkerSize(0.5);
	    
      if(h_amp_vs_sampl_chan150[3]) h_amp_vs_sampl_chan150[3]->Draw("PC");
      h_amp_vs_sampl_chan150[3]->SetMinimum(50);
      h_amp_vs_sampl_chan150[3]->SetMaximum(150);
  
      // f1ADC250 multiplicity histograms
      c1->cd(2);
      gPad->SetTicks();
      gStyle->SetOptStat(10);
	        
      //set margins
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.03);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.1);
      //set axis label size
      h_amp_vs_sampl_chan[3]->SetLabelSize(0.04,"X");
      h_amp_vs_sampl_chan[3]->SetLabelSize(0.04,"Y");
      h_amp_vs_sampl_chan[3]->GetXaxis()->CenterTitle();
      h_amp_vs_sampl_chan[3]->GetYaxis()->CenterTitle();
      h_amp_vs_sampl_chan[3]->SetTitleSize(0.06,"X");
      h_amp_vs_sampl_chan[3]->SetTitleSize(0.06,"Y");
      
      //set marker,&  marker size and color
      h_amp_vs_sampl_chan[3]->SetMarkerStyle(8);
      h_amp_vs_sampl_chan[3]->SetMarkerColor(4);
      h_amp_vs_sampl_chan[3]->SetMarkerSize(0.5);
	    
      if(h_amp_vs_sampl_chan[3]) h_amp_vs_sampl_chan[3]->Draw("PC");
      h_amp_vs_sampl_chan[3]->SetMinimum(50);
      h_amp_vs_sampl_chan[3]->SetMaximum(1000);
  
  
	  c1->cd(3);
	  gPad->SetTicks();
	  gStyle->SetOptStat(10);
	  //set margins
	  gPad->SetLeftMargin(0.12);
	  gPad->SetRightMargin(0.03);
	  gPad->SetBottomMargin(0.14);
	  gPad->SetTopMargin(0.1);
	  //set axis label size
	  h_amp_vs_sampl_chan1000[3]->SetLabelSize(0.04,"X");
	  h_amp_vs_sampl_chan1000[3]->SetLabelSize(0.04,"Y");
	  h_amp_vs_sampl_chan1000[3]->GetXaxis()->CenterTitle();
	  h_amp_vs_sampl_chan1000[3]->GetYaxis()->CenterTitle();
	  h_amp_vs_sampl_chan1000[3]->SetTitleSize(0.06,"X");
	  h_amp_vs_sampl_chan1000[3]->SetTitleSize(0.06,"Y");
	  //set marker,&  marker size and color
	  h_amp_vs_sampl_chan1000[3]->SetMarkerStyle(8);
	  h_amp_vs_sampl_chan1000[3]->SetMarkerColor(4);
	  h_amp_vs_sampl_chan1000[3]->SetMarkerSize(0.5);
	  if(h_amp_vs_sampl_chan1000[3]) h_amp_vs_sampl_chan1000[3]->Draw("PC");
	  h_amp_vs_sampl_chan1000[3]->SetMinimum(50);
	  h_amp_vs_sampl_chan1000[3]->SetMaximum(2000);
	

	  c1->cd(4);
	  gPad->SetTicks();
	  gStyle->SetOptStat(10);
	  //set margins
	  gPad->SetLeftMargin(0.12);
	  gPad->SetRightMargin(0.03);
	  gPad->SetBottomMargin(0.14);
	  gPad->SetTopMargin(0.1);
	  //set axis label size
	  h_amp_vs_sampl_chan2000[3]->SetLabelSize(0.04,"X");
	  h_amp_vs_sampl_chan2000[3]->SetLabelSize(0.04,"Y");
	  h_amp_vs_sampl_chan2000[3]->GetXaxis()->CenterTitle();
	  h_amp_vs_sampl_chan2000[3]->GetYaxis()->CenterTitle();
	  h_amp_vs_sampl_chan2000[3]->SetTitleSize(0.06,"X");
	  h_amp_vs_sampl_chan2000[3]->SetTitleSize(0.06,"Y");
	  //set marker,&  marker size and color
	  h_amp_vs_sampl_chan2000[3]->SetMarkerStyle(8);
	  h_amp_vs_sampl_chan2000[3]->SetMarkerColor(4);
	  h_amp_vs_sampl_chan2000[3]->SetMarkerSize(0.5);
	  if(h_amp_vs_sampl_chan2000[3]) h_amp_vs_sampl_chan2000[3]->Draw("PC");
	  h_amp_vs_sampl_chan2000[3]->SetMinimum(50);
	  h_amp_vs_sampl_chan2000[3]->SetMaximum(3000);


	  c1->cd(5);
	  gPad->SetTicks();
	  gStyle->SetOptStat(10);
	  //set margins
	  gPad->SetLeftMargin(0.12);
	  gPad->SetRightMargin(0.03); 
	  gPad->SetBottomMargin(0.14);
	  gPad->SetTopMargin(0.1);
	  //set axis label size
	  h_amp_vs_sampl_chan3000[3]->SetLabelSize(0.04,"X");
	  h_amp_vs_sampl_chan3000[3]->SetLabelSize(0.04,"Y");
	  h_amp_vs_sampl_chan3000[3]->GetXaxis()->CenterTitle();
	  h_amp_vs_sampl_chan3000[3]->GetYaxis()->CenterTitle();
	  h_amp_vs_sampl_chan3000[3]->SetTitleSize(0.06,"X");
	  h_amp_vs_sampl_chan3000[3]->SetTitleSize(0.06,"Y");
	  //set marker,&  marker size and color
	  h_amp_vs_sampl_chan3000[3]->SetMarkerStyle(8);
	  h_amp_vs_sampl_chan3000[3]->SetMarkerColor(4);
	  h_amp_vs_sampl_chan3000[3]->SetMarkerSize(0.5);
	  if(h_amp_vs_sampl_chan3000[3]) h_amp_vs_sampl_chan3000[3]->Draw("PC");
	  h_amp_vs_sampl_chan3000[3]->SetMinimum(50);
	  h_amp_vs_sampl_chan3000[3]->SetMaximum(4000);
       

	  c1->cd(6);
	  gPad->SetTicks();
	  gStyle->SetOptStat(10);
	  //set margins
	  gPad->SetLeftMargin(0.12);
	  gPad->SetRightMargin(0.03);
	  gPad->SetBottomMargin(0.14);
	  gPad->SetTopMargin(0.1);
	  //set axis label size
	  h_amp_vs_sampl_chan4000[3]->SetLabelSize(0.04,"X");
	  h_amp_vs_sampl_chan4000[3]->SetLabelSize(0.04,"Y");
	  h_amp_vs_sampl_chan4000[3]->GetXaxis()->CenterTitle();
	  h_amp_vs_sampl_chan4000[3]->GetYaxis()->CenterTitle();
	  h_amp_vs_sampl_chan4000[3]->SetTitleSize(0.06,"X");
	  h_amp_vs_sampl_chan4000[3]->SetTitleSize(0.06,"Y");
	  //set marker,&  marker size and color
	  h_amp_vs_sampl_chan4000[3]->SetMarkerStyle(8);
	  h_amp_vs_sampl_chan4000[3]->SetMarkerColor(4);
	  h_amp_vs_sampl_chan4000[3]->SetMarkerSize(0.5);
	  if(h_amp_vs_sampl_chan4000[3]) h_amp_vs_sampl_chan4000[3]->Draw("PC");
	  h_amp_vs_sampl_chan4000[3]->SetMinimum(50);
	  h_amp_vs_sampl_chan4000[3]->SetMaximum(9000);
	

      
}

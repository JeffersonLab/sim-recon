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
      TH1I *h_amp_vs_sampl_chan[i]       = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan_%i", i+1));
      TH1I *h_amp_vs_sampl_chan150[i]    = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan150_%i", i+1));
      TH1I *h_amp_vs_sampl_chan1000[i]   = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan1000_%i", i+1));
      TH1I *h_amp_vs_sampl_chan2000[i]   = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan2000_%i", i+1));
      TH1I *h_amp_vs_sampl_chan3000[i]   = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan3000_%i", i+1));
      TH1I *h_amp_vs_sampl_chan4000[i]   = (TH1I*)gDirectory->FindObjectAny(Form("amp_vs_sampl_chan4000_%i", i+1));
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
  c1->Divide(5,6);

  for(unsigned int i = 0; i < NCHANNELS; i++)
    {
      // f1ADC250 multiplicity histograms
      c1->cd(i+1);
      gPad->SetTicks();
      gStyle->SetOptStat(10);
	        
      //set margins
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.02);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.1);
      //set axis label size
      h_amp_vs_sampl_chan150[i]->SetLabelSize(0.06,"X");
      h_amp_vs_sampl_chan150[i]->SetLabelSize(0.06,"Y");
      h_amp_vs_sampl_chan150[i]->GetXaxis()->CenterTitle();
      h_amp_vs_sampl_chan150[i]->GetYaxis()->CenterTitle();
      h_amp_vs_sampl_chan150[i]->SetTitleSize(0.06,"X");
      h_amp_vs_sampl_chan150[i]->SetTitleSize(0.06,"Y");
      
      //set marker,&  marker size and color
      h_amp_vs_sampl_chan150[i]->SetMarkerStyle(8);
      h_amp_vs_sampl_chan150[i]->SetMarkerColor(4);
      h_amp_vs_sampl_chan150[i]->SetMarkerSize(0.5);
	    
      if(h_amp_vs_sampl_chan150[i]) h_amp_vs_sampl_chan150[i]->Draw("PC");
      h_amp_vs_sampl_chan150[i]->SetMinimum(50);
      h_amp_vs_sampl_chan150[i]->SetMaximum(150);
    }
  // Create the canvas c2
  TCanvas *c2 = new TCanvas("c2","Start Counter Expert Waveform Histograms( 150 < Pulse height <= 1000)", 200, 10, 600, 480);
  c2->cd(0);
  c2->Draw();
  c2->Update();  
  if(!gPad) return;
  TCanvas *c2 = gPad->GetCanvas();
  c2->Divide(5,6);

  for(unsigned int i = 0; i < NCHANNELS; i++)
    {
      // f1ADC250 multiplicity histograms
      c2->cd(i+1);
      gPad->SetTicks();
      gStyle->SetOptStat(10);
	        
      //set margins
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.02);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.1);
      //set axis label size
      h_amp_vs_sampl_chan[i]->SetLabelSize(0.06,"X");
      h_amp_vs_sampl_chan[i]->SetLabelSize(0.06,"Y");
      h_amp_vs_sampl_chan[i]->GetXaxis()->CenterTitle();
      h_amp_vs_sampl_chan[i]->GetYaxis()->CenterTitle();
      h_amp_vs_sampl_chan[i]->SetTitleSize(0.06,"X");
      h_amp_vs_sampl_chan[i]->SetTitleSize(0.06,"Y");
      
      //set marker,&  marker size and color
      h_amp_vs_sampl_chan[i]->SetMarkerStyle(8);
      h_amp_vs_sampl_chan[i]->SetMarkerColor(4);
      h_amp_vs_sampl_chan[i]->SetMarkerSize(0.5);
	    
      if(h_amp_vs_sampl_chan[i]) h_amp_vs_sampl_chan[i]->Draw("PC");
      h_amp_vs_sampl_chan[i]->SetMinimum(50);
      h_amp_vs_sampl_chan[i]->SetMaximum(1000);
    }
  // Create the canvas c3
  TCanvas *c3 = new TCanvas("c3","Start Counter Expert Waveform Histograms(1000 < Pulse height <= 2000)", 200, 10, 600, 480);
      c3->cd(0);
      c3->Draw();
      c3->Update();
      if(!gPad) return;
      TCanvas *c3 = gPad->GetCanvas();
      c3->Divide(5,6);
      
      for(unsigned int i = 0; i < NCHANNELS; i++)
	{
	  // f1ADC250 multiplicity histograms
	  c3->cd(i+1);
	  gPad->SetTicks();
	  gStyle->SetOptStat(10);
	  //set margins
	  gPad->SetLeftMargin(0.12);
	  gPad->SetRightMargin(0.02);
	  gPad->SetBottomMargin(0.14);
	  gPad->SetTopMargin(0.1);
	  //set axis label size
	  h_amp_vs_sampl_chan1000[i]->SetLabelSize(0.06,"X");
	  h_amp_vs_sampl_chan1000[i]->SetLabelSize(0.06,"Y");
	  h_amp_vs_sampl_chan1000[i]->GetXaxis()->CenterTitle();
	  h_amp_vs_sampl_chan1000[i]->GetYaxis()->CenterTitle();
	  h_amp_vs_sampl_chan1000[i]->SetTitleSize(0.06,"X");
	  h_amp_vs_sampl_chan1000[i]->SetTitleSize(0.06,"Y");
	  //set marker,&  marker size and color
	  h_amp_vs_sampl_chan1000[i]->SetMarkerStyle(8);
	  h_amp_vs_sampl_chan1000[i]->SetMarkerColor(4);
	  h_amp_vs_sampl_chan1000[i]->SetMarkerSize(0.5);
	  if(h_amp_vs_sampl_chan1000[i]) h_amp_vs_sampl_chan1000[i]->Draw("PC");
	  h_amp_vs_sampl_chan1000[i]->SetMinimum(50);
	  h_amp_vs_sampl_chan1000[i]->SetMaximum(2000);
	}
// Create the canvas c4
  TCanvas *c4 = new TCanvas("c4","Start Counter Expert Waveform Histograms (2000 < Pulse height <= 3000)", 200, 10, 600, 480);
      c4->cd(0);
      c4->Draw();
      c4->Update();
      if(!gPad) return;
      TCanvas *c4 = gPad->GetCanvas();
      c4->Divide(5,6);
      
      for(unsigned int i = 0; i < NCHANNELS; i++)
	{
	  // f1ADC250 multiplicity histograms
	  c4->cd(i+1);
	  gPad->SetTicks();
	  gStyle->SetOptStat(10);
	  //set margins
	  gPad->SetLeftMargin(0.12);
	  gPad->SetRightMargin(0.02);
	  gPad->SetBottomMargin(0.14);
	  gPad->SetTopMargin(0.1);
	  //set axis label size
	  h_amp_vs_sampl_chan2000[i]->SetLabelSize(0.06,"X");
	  h_amp_vs_sampl_chan2000[i]->SetLabelSize(0.06,"Y");
	  h_amp_vs_sampl_chan2000[i]->GetXaxis()->CenterTitle();
	  h_amp_vs_sampl_chan2000[i]->GetYaxis()->CenterTitle();
	  h_amp_vs_sampl_chan2000[i]->SetTitleSize(0.06,"X");
	  h_amp_vs_sampl_chan2000[i]->SetTitleSize(0.06,"Y");
	  //set marker,&  marker size and color
	  h_amp_vs_sampl_chan2000[i]->SetMarkerStyle(8);
	  h_amp_vs_sampl_chan2000[i]->SetMarkerColor(4);
	  h_amp_vs_sampl_chan2000[i]->SetMarkerSize(0.5);
	  if(h_amp_vs_sampl_chan2000[i]) h_amp_vs_sampl_chan2000[i]->Draw("PC");
	  h_amp_vs_sampl_chan2000[i]->SetMinimum(50);
	  h_amp_vs_sampl_chan2000[i]->SetMaximum(3000);
	}
// Create the canvas c5
  TCanvas *c5 = new TCanvas("c5","Start Counter Expert Waveform Histograms(3000 < Pulse height <= 4000)", 200, 10, 600, 480);
      c5->cd(0);
      c5->Draw();
      c5->Update();
      if(!gPad) return;
      TCanvas *c5 = gPad->GetCanvas();
      c5->Divide(5,6);
      
      for(unsigned int i = 0; i < NCHANNELS; i++)
	{
	  // f1ADC250 multiplicity histograms
	  c5->cd(i+1);
	  gPad->SetTicks();
	  gStyle->SetOptStat(10);
	  //set margins
	  gPad->SetLeftMargin(0.12);
	  gPad->SetRightMargin(0.02);
	  gPad->SetBottomMargin(0.14);
	  gPad->SetTopMargin(0.1);
	  //set axis label size
	  h_amp_vs_sampl_chan3000[i]->SetLabelSize(0.06,"X");
	  h_amp_vs_sampl_chan3000[i]->SetLabelSize(0.06,"Y");
	  h_amp_vs_sampl_chan3000[i]->GetXaxis()->CenterTitle();
	  h_amp_vs_sampl_chan3000[i]->GetYaxis()->CenterTitle();
	  h_amp_vs_sampl_chan3000[i]->SetTitleSize(0.06,"X");
	  h_amp_vs_sampl_chan3000[i]->SetTitleSize(0.06,"Y");
	  //set marker,&  marker size and color
	  h_amp_vs_sampl_chan3000[i]->SetMarkerStyle(8);
	  h_amp_vs_sampl_chan3000[i]->SetMarkerColor(4);
	  h_amp_vs_sampl_chan3000[i]->SetMarkerSize(0.5);
	  if(h_amp_vs_sampl_chan3000[i]) h_amp_vs_sampl_chan3000[i]->Draw("PC");
	  h_amp_vs_sampl_chan3000[i]->SetMinimum(50);
	  h_amp_vs_sampl_chan3000[i]->SetMaximum(4000);
	}
// Create the canvas c6
  TCanvas *c6 = new TCanvas("c6","Start Counter Expert Waveform Histograms(4000 < Pulse height)", 200, 10, 600, 480);
      c6->cd(0);
      c6->Draw();
      c6->Update();
      if(!gPad) return;
      TCanvas *c6 = gPad->GetCanvas();
      c6->Divide(5,6);
      
      for(unsigned int i = 0; i < NCHANNELS; i++)
	{
	  // f1ADC250 multiplicity histograms
	  c6->cd(i+1);
	  gPad->SetTicks();
	  gStyle->SetOptStat(10);
	  //set margins
	  gPad->SetLeftMargin(0.12);
	  gPad->SetRightMargin(0.02);
	  gPad->SetBottomMargin(0.14);
	  gPad->SetTopMargin(0.1);
	  //set axis label size
	  h_amp_vs_sampl_chan4000[i]->SetLabelSize(0.06,"X");
	  h_amp_vs_sampl_chan4000[i]->SetLabelSize(0.06,"Y");
	  h_amp_vs_sampl_chan4000[i]->GetXaxis()->CenterTitle();
	  h_amp_vs_sampl_chan4000[i]->GetYaxis()->CenterTitle();
	  h_amp_vs_sampl_chan4000[i]->SetTitleSize(0.06,"X");
	  h_amp_vs_sampl_chan4000[i]->SetTitleSize(0.06,"Y");
	  //set marker,&  marker size and color
	  h_amp_vs_sampl_chan4000[i]->SetMarkerStyle(8);
	  h_amp_vs_sampl_chan4000[i]->SetMarkerColor(4);
	  h_amp_vs_sampl_chan4000[i]->SetMarkerSize(0.5);
	  if(h_amp_vs_sampl_chan4000[i]) h_amp_vs_sampl_chan4000[i]->Draw("PC");
	  h_amp_vs_sampl_chan4000[i]->SetMinimum(50);
	  h_amp_vs_sampl_chan4000[i]->SetMaximum(9000);
	}

      
}

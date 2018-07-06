
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/digihits_trig1
// hnamepath: /occupancy/digihits_trig3
// hnamepath: /occupancy/digihits_trig4
// hnamepath: /occupancy/digihits_scale_factors
//
// e-mail: davidl@jlab.org
// e-mail: tbritton@jlab.org
//

{

// The following are empty versions of routines defined in RootSpy
// compiled executables. These are defined here for when this
// macro is run outside that context.
#ifndef ROOTSPY_MACROS
#define rs_SetFlag(A) cout<<"rs_SetFlag ignored outside of RootSpy context"<<endl
#define rs_GetFlag(A) 0
#define rs_ResetHisto(A) cout<<"rs_ResetHisto ignored outside of RootSpy context"<<endl
#define rs_RestoreHisto(A) cout<<"rs_RestoreHisto ignored outside of RootSpy context"<<endl
#define InsertSeriesData(A) cout<<"InsertSeriesData ignored outside of RootSpy context"<<endl
#define InsertSeriesMassFit(A,B,C,D,E,F) cout<<"InsertSeriesMassFit ignored outside of RootSpy context"<<endl
#endif

	// RootSpy saves the current directory and style before
	// calling the macro and restores it after so it is OK to
	// change them and not change them back.
	TDirectory *savedir = gDirectory;

	// First get EventInfo which is used to get unix_time for time series DB
	double unix_time =  0;
	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(dir) dir->cd();
	TH1* EventInfo           = (TH1*)gDirectory->Get("EventInfo");
	if(EventInfo){
		Double_t Nunix = EventInfo->GetBinContent(1);
		if(Nunix>0.0){
			Double_t sum_unix_time = EventInfo->GetBinContent(2);
			unix_time = (sum_unix_time/Nunix);
			time_t t = (time_t)unix_time;
			cout << ctime(&t);
		}
	}

	// Get all other histos from occupancy
	savedir->cd();
	dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	TH2I *digihits_trig1 = (TH2I*)gDirectory->FindObjectAny("digihits_trig1");
	TH2I *digihits_trig3 = (TH2I*)gDirectory->FindObjectAny("digihits_trig3");
	TH2I *digihits_trig4 = (TH2I*)gDirectory->FindObjectAny("digihits_trig4");
	TH1F *digihits_scale_factors = (TH1F*)gDirectory->FindObjectAny("digihits_scale_factors");

	// Just for testing
	if(gPad == NULL){
		TCanvas *c1 = new TCanvas("c1", "", 1200, 900);
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}
	if(!gPad) {savedir->cd(); return;}
	
	// Some detectors are scaled in order to fit them on the same
	// graph. The scale factors are communicated via the digihits_scale_factors
	// histogram. Since RootSpy sums this histogram we must divide the bin contents
	// by the last bin which keeps track of how many hists were added.
	map<string,double> digihitbinmap; // bin number
	map<string,double> digihitsclmap; // Scale number of hits by this (0 means don't scale)
	if(digihits_scale_factors){
		double Nhists_summed = digihits_scale_factors->GetBinContent(digihits_scale_factors->GetNbinsX());
		TAxis *xaxis = digihits_scale_factors->GetXaxis();
		for(int ibin=1; ibin<xaxis->GetNbins(); ibin++){
			string lab = xaxis->GetBinLabel(ibin);
			digihitbinmap[lab] = (double)ibin;
			digihitsclmap[lab] = digihits_scale_factors->GetBinContent(ibin)/Nhists_summed;
		}
	}

	TCanvas *c1 = gPad->GetCanvas();

	c1->cd(0);	
	if(digihits_trig4){
	
		TPad *pad1 = (TPad*)gDirectory->FindObjectAny("digihitspad2");
		if(!pad1) pad1 = new TPad("digihitspad2", "", 0.0, 0.0, 1.0, 0.5);
		pad1->Draw();
		pad1->cd();
		pad1->SetTicks();
		pad1->SetRightMargin(0.10);
		pad1->SetLeftMargin(0.05);
		pad1->SetBottomMargin(0.15);
		pad1->SetTopMargin(0.05);
		pad1->SetLogz();
		digihits_trig4->SetStats(0);
		digihits_trig4->SetTitle(" ");
		digihits_trig4->GetXaxis()->SetLabelSize(0.06);
		digihits_trig4->GetYaxis()->SetTitleSize(0.06);
		digihits_trig4->GetYaxis()->SetTitleOffset(0.35);
		digihits_trig4->Draw("colz");
		
		// Draw trigger label
		TLatex latex;
		latex.SetTextSize(0.07);
		latex.SetTextAlign(12);
		latex.DrawLatex(12.0, 115.0, "Trig 4: PS");
		
		// Draw any non-zero and non-one scale factors
		latex.SetTextSize(0.04);
		latex.SetTextAngle(90.0);
		latex.SetTextAlign(22);
		latex.SetTextColor(kWhite);
		TBox box;
		box.SetFillColor(kGray+1);
		TAxis *xaxis = digihits_trig4->GetXaxis();
		for(auto p : digihitsclmap){
			if(p.second == 0.0) continue;
			double x = xaxis->GetBinCenter((int)digihitbinmap[p.first]);
			double y = 70.0;
			char str[256];
			if( p.second>1.0 ){
				sprintf(str, "#divide%2.0f", p.second);
			}else{
				sprintf(str, "#times%2.0f", 1.0/p.second);
			}
			box.DrawBox(x-0.3, y-9.0, x+0.33, y+9.0);
			latex.DrawLatex(x, y, str);
		}
		
		// Draw mean number of hits
		box.SetFillColor(kGray);
		latex.SetTextColor(kBlack);
		TAxis *yaxis = digihits_trig4->GetYaxis();
		stringstream ss;
		ss << "digihits,trig=4 ";
		for(auto p : digihitbinmap){
			int ibin = (int)digihitbinmap[p.first];
			double x = xaxis->GetBinCenter(ibin);
			double y = 139.0;
			
			double sum  = 0.0;
			double sumw = 0.0;			
			for(int jbin=1; jbin<=yaxis->GetNbins(); jbin++){
 				double w = digihits_trig4->GetBinContent(ibin, jbin);
 				double y = yaxis->GetBinLowEdge(jbin);
 				sumw += y*w;
 				sum  += w; 
			}
			
			double mean = sum==0.0 ? 0.0:sumw/sum;
			double scale = digihitsclmap[p.first];
			if(scale != 0.0) mean *= scale;
			
			char str[256];
			if( mean > 0.5 ){
				sprintf(str, "%5.1f", mean);
			}else{
				sprintf(str, "%4.3f", mean);
			}
			box.DrawBox(x-0.3, y-10.0, x+0.33, y+10.0);
			latex.DrawLatex(x, y, str);
			
			// Build time series string
			if(digihits_scale_factors){
				string lab = digihits_scale_factors->GetXaxis()->GetBinLabel(ibin);
				if( ibin>1 ) ss << ",";
				ss << lab << "=" << mean;
			}
		}
		
		// Only write out once every 100000 objects
		if(digihits_trig1->Integral()>100000){
			if(unix_time!=0.0) ss<<" "<<(uint64_t)(unix_time*1.0E9);  // time is in units of ns
			InsertSeriesData( ss.str() );
		
			// Optionally reset the histogram so next fit is independent of this one
			if(rs_GetFlag("RESET_AFTER_FIT")) rs_ResetHisto("/occupancy/digihits_trig4");
		}
		
	}	

	c1->cd(0);	
	if(digihits_trig1){
	
		TPad *pad1 = (TPad*)gDirectory->FindObjectAny("digihitspad1");
		if(!pad1) pad1 = new TPad("digihitspad1", "", 0.0, 0.5, 1.0, 1.0);
		pad1->Draw();
		pad1->cd();
		pad1->SetTicks();
		pad1->SetRightMargin(0.10);
		pad1->SetLeftMargin(0.05);
		pad1->SetBottomMargin(0.15);
		pad1->SetTopMargin(0.05);
		pad1->SetLogz();
		digihits_trig1->SetStats(0);
		digihits_trig1->SetTitle(" ");
		digihits_trig1->GetXaxis()->SetLabelSize(0.06);
		digihits_trig1->GetYaxis()->SetTitleSize(0.06);
		digihits_trig1->GetYaxis()->SetTitleOffset(0.35);
		digihits_trig1->Draw("colz");
		
		// Draw trigger label
		TLatex latex;
		latex.SetTextSize(0.07);
		latex.SetTextAlign(12);
		latex.DrawLatex(12.0, 115.0, "Trig 1: FCAL+BCAL");
		
		// Draw any non-zero and non-one scale factors
		latex.SetTextSize(0.04);
		latex.SetTextAngle(90.0);
		latex.SetTextAlign(22);
		latex.SetTextColor(kWhite);
		TBox box;
		box.SetFillColor(kGray+1);
		TAxis *xaxis = digihits_trig1->GetXaxis();
		for(auto p : digihitsclmap){
			if(p.second == 0.0) continue;
			double x = xaxis->GetBinCenter((int)digihitbinmap[p.first]);
			double y = 70.0;
			char str[256];
			if( p.second>1.0 ){
				sprintf(str, "#divide%2.0f", p.second);
			}else{
				sprintf(str, "#times%2.0f", 1.0/p.second);
			}
			box.DrawBox(x-0.3, y-9.0, x+0.33, y+9.0);
			latex.DrawLatex(x, y, str);
		}
		
		// Draw mean number of hits
		box.SetFillColor(kGray);
		latex.SetTextColor(kBlack);
		TAxis *yaxis = digihits_trig1->GetYaxis();
		stringstream ss;
		ss << "digihits,trig=1 ";
		for(auto p : digihitbinmap){
			int ibin = (int)digihitbinmap[p.first];
			double x = xaxis->GetBinCenter(ibin);
			double y = 139.0;
			
			double sum  = 0.0;
			double sumw = 0.0;
			for(int jbin=1; jbin<=yaxis->GetNbins(); jbin++){
				double w = digihits_trig1->GetBinContent(ibin, jbin);
				double y = yaxis->GetBinLowEdge(jbin);
				sumw += y*w;
				sum  += w; 
			}
			
			double mean = sum==0.0 ? 0.0:sumw/sum;
			double scale = digihitsclmap[p.first];
			if(scale>0.0) mean *= scale;
			
			char str[256];
			if( mean > 0.5 ){
				sprintf(str, "%5.1f", mean);
			}else{
				sprintf(str, "%4.3f", mean);
			}
			box.DrawBox(x-0.3, y-10.0, x+0.33, y+10.0);
			latex.DrawLatex(x, y, str);
			
			// Build time series string
			if(digihits_scale_factors){
				string lab = digihits_scale_factors->GetXaxis()->GetBinLabel(ibin);
				if( ibin>1 ) ss << ",";
				ss  << lab << "=" << mean;
			}
		}
		
		// Only write out once every 100000 objects
		if(digihits_trig1->Integral()>100000){
			if(unix_time!=0.0) ss<<" "<<(uint64_t)(unix_time*1.0E9);  // time is in units of ns
			InsertSeriesData( ss.str() );
		
			// Optionally reset the histogram so next fit is independent of this one
			if(rs_GetFlag("RESET_AFTER_FIT")) rs_ResetHisto("/occupancy/digihits_trig1");
		}
	}
}



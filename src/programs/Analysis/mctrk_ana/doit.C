

void doit(void)
{
	gROOT->Reset();

	int Ny = dp_over_p_vs_theta->GetNbinsY();
	dp_over_p_vs_theta->ProjectionX();
	dp_over_p_vs_theta->FitSlicesY(0,1,Ny,0);
	
	TAxis *xaxis = dp_over_p_vs_theta->GetXaxis();
	int Nx = xaxis->GetNbins();
	float xmin = xaxis->GetXmin();
	float xmax = xaxis->GetXmax();

	TH1F *perr = new TH1F("perr","perr",Nx, xmin, xmax);
	TH1F *perr_deg = new TH1F("perr","perr",Nx, xmin*57.3, xmax*57.3);

	int i;
	for(i=1;i<=Nx;i++){

		float sig1;
		sig1 = dp_over_p_vs_theta_2->GetBinContent(i);

		// Make sure there's enough events
		int n;
		n = dp_over_p_vs_theta_px->GetBinContent(i);
		if(n<50)sig1 = 0.0;

		perr->SetBinContent(i,sig1);
		perr_deg->SetBinContent(i,sig1);
		cout<<" "<<i<<" : "<<sig1<<"  n="<<n<<endl;
	}
	
	perr->SetMarkerStyle(20);
	perr->SetMarkerColor(kRed);
	perr_deg->SetMarkerStyle(20);
	perr_deg->SetMarkerColor(kRed);
	//dp_over_p_vs_theta->SetMarkerStyle(22);
	dp_over_p_vs_theta->SetMarkerColor(kBlack);
	
	//TH2F *axes = new TH2F("axes",";#theta(radians);#Deltap/p",10,0.0,160.0/57.3,10,0.0,0.1);
	//axes->Draw();
	//dp_over_p_vs_theta->Draw("Psame");
	//perr->Draw("Psame");
	
	TH2F *axes = new TH2F("axes",";#theta(degrees);#Deltap/p",10,0.0,160.0,10,0.0,0.1);
	axes->Draw();
	perr_deg->Draw("Psame");
	
	// Show FDC and CDC limits
	drawalllimits();
}

void drawalllimits(void)
{
	drawlimits(1.15, 20.1, kBlue);
	drawlimits(0.88, 15.3, kBlue);
	drawlimits(0.71, 12.4, kBlue);
	drawlimits(0.59, 10.3, kBlue);

	drawlimits(6.0, 108.0, kRed);
}

void drawlimits(double lo, double hi, Color_t lcolor)
{
	TLine *fdcloA = new TLine(lo, 0.0, lo, 0.1);
	TLine *fdchiA = new TLine(hi, 0.0, hi, 0.1);
	fdcloA->SetLineColor(lcolor);
	fdchiA->SetLineColor(lcolor);
	fdcloA->Draw();
	fdchiA->Draw();
}

void doitp(void)
{
	gROOT->Reset();

	int Ny = dp_over_p_vs_p->GetNbinsY();
	dp_over_p_vs_p->ProjectionX();
	dp_over_p_vs_p->FitSlicesY(0,1,Ny,0);
	
	TAxis *xaxis = dp_over_p_vs_p->GetXaxis();
	int Nx = xaxis->GetNbins();
	float xmin = xaxis->GetXmin();
	float xmax = xaxis->GetXmax();

	TH1F *perr = new TH1F("perr","perr",Nx, xmin, xmax);

	int i;
	for(i=1;i<=Nx;i++){

		float sig1;
		sig1 = dp_over_p_vs_p_2->GetBinContent(i);

		// Make sure there's enough events
		int n;
		n = dp_over_p_vs_p_px->GetBinContent(i);
		if(n<50)sig1 = 0.0;

		perr->SetBinContent(i,sig1);
		cout<<" "<<i<<" : "<<sig1<<"  n="<<n<<endl;
	}
	
	perr->SetMarkerStyle(20);
	perr->SetMarkerColor(kRed);
	//dp_over_p_vs_theta->SetMarkerStyle(22);
	dp_over_p_vs_p->SetMarkerColor(kBlack);
	
	TH2F *axes = new TH2F("axes",";#p(GeV/c);#Deltap/p",10,0.0,10.0,10,0.0,0.1);
	axes->Draw();
	dp_over_p_vs_p->Draw("Psame");
	perr->Draw("Psame");
}


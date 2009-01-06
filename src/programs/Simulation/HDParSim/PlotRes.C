

//----------------
// PlotRes
//----------------
void PlotRes(void)
{
	gROOT->Reset();
	gStyle->SetPalette(1,NULL);
	
	TCanvas *c1 = new TCanvas("c1");
	c1->SetTickx();
	c1->SetTicky();
	
	TFile *f = new TFile("hd_res_charged.root");

	PlotEfficiency();
	PlotNevents();
	PlotPtotRes();
	PlotPtransRes();
	PlotPtotResVsTheta();
	PlotPtransResVsTheta();
}

//----------------
// AddStandardLabels
//----------------
void AddStandardLabels(TH2F *axes=NULL)
{
	// This will draw a label or two on the
	// current plot using the NDC coordinates.
	// It is put here to make sure all plots have
	// a consistent labeling.
	
	// Date, Author
	TLatex *lab = new TLatex(0.7, 0.7, "Feb. 22, 2008 DL");
	ConvertFromNDC(lab, axes);
	lab->SetTextSize(0.03);
	lab->SetTextAlign(33);
	lab->Draw();

	// SVN Revision
	lab = new TLatex(0.7, 0.645, "svn revision: 3355");
	ConvertFromNDC(lab, axes);
	lab->SetTextSize(0.02);
	lab->SetTextAlign(31);
	lab->Draw();
	
	// Event type
	lab = new TLatex(0.575, 0.59, "Single #pi^{+} with MULS and LOSS on");
	ConvertFromNDC(lab, axes);
	lab->SetTextSize(0.03);
	lab->SetTextAlign(31);
	lab->Draw();

	// Candidate type
	//lab = new TLatex(-0.2, 0.580, "Candidates from THROWN values");
	lab = new TLatex(0.52, 0.65, "Candidates from THROWN values");
	ConvertFromNDC(lab, axes);
	lab->SetTextSize(0.02);
	lab->SetTextAlign(32);
	lab->Draw();
}

//----------------
// PlotPtransRes
//----------------
void PlotPtransRes(void)
{
	TH2D* dpt_over_pt_vs_p_vs_theta_2 = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_2");
	//dpt_over_pt_vs_p_vs_theta_2->GetZaxis()->SetRangeUser(0.0, 8.0);		
	dpt_over_pt_vs_p_vs_theta_2->SetStats(0);
	dpt_over_pt_vs_p_vs_theta_2->SetTitle("#sigma_{#Deltap_{tran}/p_{tran}} vs. p_{tot} vs. #theta");
	dpt_over_pt_vs_p_vs_theta_2->SetXTitle("#theta Angle (degrees)");
	dpt_over_pt_vs_p_vs_theta_2->SetYTitle("Total Momentum (GeV/c)");
	dpt_over_pt_vs_p_vs_theta_2->Draw("cont4z");
	//PlotBinGrid(dpt_over_pt_vs_p_vs_theta_2, true);
	
	AddStandardLabels();
	
	c1->SaveAs("dpt_over_pt_vs_p_vs_theta.gif");
	c1->SaveAs("dpt_over_pt_vs_p_vs_theta.pdf");
}

//----------------
// PlotPtotRes
//----------------
void PlotPtotRes(void)
{
	TH2D* dp_over_p_vs_p_vs_theta_2 = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_2");
	//dp_over_p_vs_p_vs_theta_2->GetZaxis()->SetRangeUser(0.0, 8.0);		
	dp_over_p_vs_p_vs_theta_2->SetStats(0);
	dp_over_p_vs_p_vs_theta_2->SetTitle("#sigma_{#Deltap/p} vs. p_{tot} vs. #theta");
	dp_over_p_vs_p_vs_theta_2->SetXTitle("#theta Angle (degrees)");
	dp_over_p_vs_p_vs_theta_2->SetYTitle("Total Momentum (GeV/c)");
	dp_over_p_vs_p_vs_theta_2->Draw("cont4z");
	//PlotBinGrid(dp_over_p_vs_p_vs_theta_2, true);
	
	AddStandardLabels();
	
	c1->SaveAs("dp_over_p_vs_p_vs_theta.gif");
	c1->SaveAs("dp_over_p_vs_p_vs_theta.pdf");
}

//----------------
// PlotNevents
//----------------
void PlotNevents(void)
{
	TH2D* eff_vs_p_vs_theta_numerator = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta_numerator");
	//eff_vs_p_vs_theta_numerator->GetZaxis()->SetRangeUser(0.8, 1.0);		
	eff_vs_p_vs_theta_numerator->SetStats(0);
	eff_vs_p_vs_theta_numerator->SetTitle("Number of Reconstructed Tracks per bin");
	eff_vs_p_vs_theta_numerator->SetXTitle("#theta Angle (degrees)");
	eff_vs_p_vs_theta_numerator->SetYTitle("Total Momentum (GeV/c)");
	eff_vs_p_vs_theta_numerator->Draw("cont4z");
	//PlotBinGrid(eff_vs_p_vs_theta_numerator, true);
	
	AddStandardLabels();
	
	c1->SaveAs("nrecon_event_vs_p_vs_theta.gif");
	c1->SaveAs("nrecon_event_vs_p_vs_theta.pdf");
}

//----------------
// PlotEfficiency
//----------------
void PlotEfficiency(void)
{
	TH2D* eff_vs_p_vs_theta = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta");
	eff_vs_p_vs_theta->GetZaxis()->SetRangeUser(0.8, 1.0);		
	eff_vs_p_vs_theta->SetStats(0);
	eff_vs_p_vs_theta->SetXTitle("#theta Angle (degrees)");
	eff_vs_p_vs_theta->SetYTitle("Total Momentum (GeV/c)");
	eff_vs_p_vs_theta->Draw("cont4z");
	//PlotBinGrid(eff_vs_p_vs_theta, true);
	
	AddStandardLabels();
	
	c1->SaveAs("eff_vs_p_vs_theta.gif");
	c1->SaveAs("eff_vs_p_vs_theta.pdf");
}

//----------------
// PlotPtotResVsTheta
//----------------
void PlotPtotResVsTheta(void)
{
	TCanvas *c1 = new TCanvas("c1");
	c1->SetTickx();
	c1->SetTicky();

	//----------- Draw emtpy axes with correct limits
	TH2F *axes = new TH2F("axes","#sigma_{#Deltap/p} vs. #theta", 100, 0.0, 150.0, 100, 0.0, 7.0);
	axes->SetStats(0);
	axes->SetXTitle("#theta Angle (degrees)");
	axes->SetYTitle("Total Momentum Resolution (%)");
	c1->SetGridy();
	axes->Draw();

	AddStandardLabels(axes);

	//----------- Create Inset
	TPad *pad = new TPad("inset","",0.20, 0.50, 0.60, 0.89);
	pad->SetTickx();
	pad->SetTicky();
	pad->SetGridy();
	pad->Draw();
	pad->cd();

	// Draw emtpy axes with right limits
	TH2F *axes = new TH2F("axes_inset","#sigma_{#Deltap/p} vs. #theta", 100, 0.0, 40.0, 100, 0.0, 6.0);
	axes->SetStats(0);
	axes->SetXTitle("#theta Angle (degrees)");
	axes->SetYTitle("Total Momentum Resolution (%)");
	c1->SetGridy();
	axes->Draw();

	// Get 3D histogram
	TH3D* dp_over_p_vs_p_vs_theta = (TH3D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta");

	TAxis *thetaaxis = dp_over_p_vs_p_vs_theta->GetXaxis();
	int Nbinstheta = thetaaxis->GetNbins();
	double thetamin = thetaaxis->GetBinLowEdge(1);
	double thetamax = thetaaxis->GetBinLowEdge(Nbinstheta+1);

	// Legend
	TLegend *leg = new TLegend(0.60, 0.6, 0.78, 0.85);

	// Loop over momenta
	int color[]={kBlue, kRed, kGreen, kMagenta, kBlack};
	int Ncolors=4;
	int i=0;
	for(double p=1.0; p<=7.0; p+=2.0, i++){
		TAxis *paxis = dp_over_p_vs_p_vs_theta->GetYaxis();
		int ipmin = paxis->FindBin(p-0.15);
		int ipmax = paxis->FindBin(p+0.15);

		// Create a 2D histo to hold results
		char hname[256];
		sprintf(hname, "dp_over_p_vs_theta%1.1f", p);
		TH1D *dp_over_p_vs_theta = new TH1D(hname, "", Nbinstheta, thetamin, thetamax);
		
		// Loop over theta bins and project out p-range onto z-axis
		for(int bin=1; bin<=Nbinstheta; bin++){
			TH1D *h = dp_over_p_vs_p_vs_theta->ProjectionZ("_pz", bin, bin, ipmin, ipmax);
			h->Fit("gaus", "0Q");
			TF1 *fun = h->GetFunction("gaus");
			double sigma = fun->GetParameter(2);
			dp_over_p_vs_theta->SetBinContent(bin, sigma);
		}

		
		// Overlay histo
		c1->cd();
		dp_over_p_vs_theta->SetTitle("#sigma_{#Deltap/p} vs. #theta");
		dp_over_p_vs_theta->SetMarkerStyle(8);
		dp_over_p_vs_theta->SetMarkerSize(0.75);
		dp_over_p_vs_theta->SetMarkerColor(color[i%Ncolors]);
		dp_over_p_vs_theta->SetLineColor(color[i%Ncolors]);
		dp_over_p_vs_theta->Draw("PLsame");

		char str[256];
		sprintf(str, "p=%1.1f GeV/c", p);
		leg->AddEntry(dp_over_p_vs_theta, str);

		// Overlay on inset
		pad->cd();
		char hname2[256];
		sprintf(hname2, "%s_inset", hname);
		TH1D *h2 = dp_over_p_vs_theta->Clone("hname2");
		h2->SetMarkerSize(0.5);
		h2->Draw("PLsame");
		
	}
	
	c1->cd();
	leg->Draw();
	
	char fname[256];
	sprintf(fname,"dp_over_p_vs_theta.gif");
	c1->SaveAs(fname);
	sprintf(fname,"dp_over_p_vs_theta.pdf");
	c1->SaveAs(fname);
}

//----------------
// PlotPtransResVsTheta
//----------------
void PlotPtransResVsTheta(void)
{
	TCanvas *c1 = new TCanvas("c1");
	c1->SetTickx();
	c1->SetTicky();

	//----------- Draw emtpy axes with correct limits
	TH2F *axes = new TH2F("axes","#sigma_{#Deltap_{trans}/p_{trans}} vs. #theta", 100, 0.0, 150.0, 100, 0.0, 7.0);
	axes->SetStats(0);
	axes->SetXTitle("#theta Angle (degrees)");
	axes->SetYTitle("Total Momentum Resolution (%)");
	c1->SetGridy();
	axes->Draw();

	AddStandardLabels(axes);

	//----------- Create Inset
	TPad *pad = new TPad("inset","",0.20, 0.50, 0.60, 0.89);
	pad->SetTickx();
	pad->SetTicky();
	pad->SetGridy();
	pad->Draw();
	pad->cd();

	// Draw emtpy axes with right limits
	TH2F *axes = new TH2F("axes_inset","#sigma_{#Deltap_{trans}/p_{trans}} vs. #theta", 100, 0.0, 40.0, 100, 0.0, 6.0);
	axes->SetStats(0);
	axes->SetXTitle("#theta Angle (degrees)");
	axes->SetYTitle("Total Momentum Resolution (%)");
	c1->SetGridy();
	axes->Draw();

	// Get 3D histogram
	TH3D* dpt_over_pt_vs_p_vs_theta = (TH3D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta");

	TAxis *thetaaxis = dpt_over_pt_vs_p_vs_theta->GetXaxis();
	int Nbinstheta = thetaaxis->GetNbins();
	double thetamin = thetaaxis->GetBinLowEdge(1);
	double thetamax = thetaaxis->GetBinLowEdge(Nbinstheta+1);

	// Legend
	TLegend *leg = new TLegend(0.60, 0.6, 0.78, 0.85);

	// Loop over momenta
	int color[]={kBlue, kRed, kGreen, kMagenta, kBlack};
	int Ncolors=4;
	int i=0;
	for(double p=1.0; p<=7.0; p+=2.0, i++){
		TAxis *paxis = dpt_over_pt_vs_p_vs_theta->GetYaxis();
		int ipmin = paxis->FindBin(p-0.15);
		int ipmax = paxis->FindBin(p+0.15);

		// Create a 2D histo to hold results
		char hname[256];
		sprintf(hname, "dpt_over_pt_vs_theta%1.1f", p);
		TH1D *dpt_over_pt_vs_theta = new TH1D(hname, "", Nbinstheta, thetamin, thetamax);
		
		// Loop over theta bins and project out p-range onto z-axis
		for(int bin=1; bin<=Nbinstheta; bin++){
			TH1D *h = dpt_over_pt_vs_p_vs_theta->ProjectionZ("_pz", bin, bin, ipmin, ipmax);
			h->Fit("gaus", "0Q");
			TF1 *fun = h->GetFunction("gaus");
			double sigma = fun->GetParameter(2);
			dpt_over_pt_vs_theta->SetBinContent(bin, sigma);
		}

		
		// Overlay histo
		c1->cd();
		dpt_over_pt_vs_theta->SetTitle("#sigma_{#Deltap_{trans}/p_{trans}} vs. #theta");
		dpt_over_pt_vs_theta->SetMarkerStyle(8);
		dpt_over_pt_vs_theta->SetMarkerSize(0.75);
		dpt_over_pt_vs_theta->SetMarkerColor(color[i%Ncolors]);
		dpt_over_pt_vs_theta->SetLineColor(color[i%Ncolors]);
		dpt_over_pt_vs_theta->Draw("PLsame");

		char str[256];
		sprintf(str, "p=%1.1f GeV/c", p);
		leg->AddEntry(dpt_over_pt_vs_theta, str);

		// Overlay on inset
		pad->cd();
		char hname2[256];
		sprintf(hname2, "%s_inset", hname);
		TH1D *h2 = dpt_over_pt_vs_theta->Clone("hname2");
		h2->SetMarkerSize(0.5);
		h2->Draw("PLsame");
		
	}
	
	c1->cd();
	leg->Draw();
	
	char fname[256];
	sprintf(fname,"dpt_over_pt_vs_theta.gif");
	c1->SaveAs(fname);
	sprintf(fname,"dpt_over_pt_vs_theta.pdf");
	c1->SaveAs(fname);
}

//----------------
// PlotBinGrid
//----------------
void PlotBinGrid(TH2D *h, bool useNDC=false)
{
	// Plot a grid showing the bins of the 2D histogram
	// on top of the current plot. If the value of "useNDC"
	// is true, then the coordinates are converted to a coordinate
	// system with 0,0 at the center and the edges at +/-1.

	TAxis *xaxis = h->GetXaxis();
	int Nbinsx = xaxis->GetNbins();
	double xmin = xaxis->GetBinLowEdge(1);
	double xmax = xmin + xaxis->GetBinLowEdge(Nbinsx);

	TAxis *yaxis = h->GetYaxis();
	int Nbinsy = yaxis->GetNbins();
	double ymin = yaxis->GetBinLowEdge(1);
	double ymax = ymin + yaxis->GetBinLowEdge(Nbinsy);

	// Lines spanning Y
	for(int bin=1; bin<=Nbinsx; bin++){
		double x = xaxis->GetBinLowEdge(bin);
		double x1 = x;
		double x2 = x;
		double y1 = ymin;
		double y2 = ymax;
		if(useNDC)Convert(xmin, ymin, xmax, ymax, x1, y1);
		if(useNDC)Convert(xmin, ymin, xmax, ymax, x2, y2);
		TLine *l = new TLine(x1, y1, x2, y2);
		l->SetLineColor(15);
		l->Draw();
	}

	// Lines spanning X
	for(int bin=1; bin<=Nbinsy; bin++){
		double y = yaxis->GetBinLowEdge(bin);
		double x1 = xmin;
		double x2 = xmax;
		double y1 = y;
		double y2 = y;
		if(useNDC)Convert(xmin, ymin, xmax, ymax, x1, y1);
		if(useNDC)Convert(xmin, ymin, xmax, ymax, x2, y2);
		TLine *l = new TLine(x1, y1, x2, y2);
		l->SetLineColor(15);
		l->Draw();
	}


}

//----------------
// Convert
//----------------
void Convert(double xmin, double ymin, double xmax, double ymax, double &x, double &y)
{
	// Convert from a NDC (I think) to the coordinates system specified
	// by the xmin, ... values.
	// This is because when the 2D histo from the file is drawn, it seems
	// to leave the default coordinate system for lines to be one such that
	// the center of the pad is 0,0 and the edges something like +/-1.0.
	// The values here (0.5 and 1.15) were obtained empirically.
	double xx = ((x-xmin)/(xmax-xmin)-0.5)*1.15;
	double yy = ((y-ymin)/(ymax-ymin)-0.5)*1.15;
	
	x = xx;
	y = yy;
}

//----------------
// ConvertFromNDC
//----------------
void ConvertFromNDC(TLatex *obj, TH2F *h=NULL)
{
	// Bugs in ROOT make it hard to plot labels consistently.
	// For 1D plots, the histogram axes define the coordinate
	// system. For 2D plots, we seem to be forced to use the
	// NDC. There does not seem to be an obvious way to tell
	// which we're using so we pass the information in in the
	// form of the "axes" histogram. If it is not NULL, then
	// we use it to define the limits. Otherwise, we do nothing.

	if(h==NULL)return;

	TAxis *xaxis = h->GetXaxis();
	int Nbinsx = xaxis->GetNbins();
	double xmin = xaxis->GetBinLowEdge(1);
	double xmax = xmin + xaxis->GetBinLowEdge(Nbinsx);

	TAxis *yaxis = h->GetYaxis();
	int Nbinsy = yaxis->GetNbins();
	double ymin = yaxis->GetBinLowEdge(1);
	double ymax = ymin + yaxis->GetBinLowEdge(Nbinsy);
	
	
	double x = obj->GetX();
	double y = obj->GetY();
	cout<<" in:  x="<<x<<"  y="<<y<<endl;
	
	x = xmin + (xmax-xmin)*(0.5+x/1.15);
	y = ymin + (ymax-ymin)*(0.5+y/1.15);
	
	obj->SetX(x);
	obj->SetY(y);

	cout<<"out:  x="<<x<<"  y="<<y<<endl;
}


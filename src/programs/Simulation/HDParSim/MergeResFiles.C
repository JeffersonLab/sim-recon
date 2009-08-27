
// Read the resolution and efficiency histograms from fileA and fileB
// and merge them into a new file. The bins of the new histograms will
// be filled with the values from the histograms in fileA for
// theta<theta_border and from the histograms in fileB for 
// theta>=theta_border.
//
// This is intended to merge results from 2 different simulation
// data sets.

void MergeResFiles(const char *fileA="hd_res_charged_thrown_fdc.root", const char *fileB="hd_res_charged_thrown_cdc.root", double theta_border=38.0)
{
	// Open fileA and get pointers to histograms
	TFile *f = new TFile(fileA);
	TH3D *dp_over_p_vs_p_vs_thetaA = (TH3D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta");
	TH2D *dp_over_p_vs_p_vs_theta_0A = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_0");
	TH2D *dp_over_p_vs_p_vs_theta_1A = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_1");
	TH2D *dp_over_p_vs_p_vs_theta_2A = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_2");
	TH2D *dp_over_p_vs_p_vs_theta_chi2A = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_chi2");
	TH3D *dpt_over_pt_vs_p_vs_thetaA = (TH3D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta");
	TH2D *dpt_over_pt_vs_p_vs_theta_0A = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_0");
	TH2D *dpt_over_pt_vs_p_vs_theta_1A = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_1");
	TH2D *dpt_over_pt_vs_p_vs_theta_2A = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_2");
	TH2D *dpt_over_pt_vs_p_vs_theta_chi2A = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_chi2");
	TH3D *dtheta_vs_p_vs_thetaA = (TH3D*)gROOT->FindObject("dtheta_vs_p_vs_theta");
	TH2D *dtheta_vs_p_vs_theta_0A = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta_0");
	TH2D *dtheta_vs_p_vs_theta_1A = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta_1");
	TH2D *dtheta_vs_p_vs_theta_2A = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta_2");
	TH2D *dtheta_vs_p_vs_theta_chi2A = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta_chi2");
	TH3D *dphi_vs_p_vs_thetaA = (TH3D*)gROOT->FindObject("dphi_vs_p_vs_theta");
	TH2D *dphi_vs_p_vs_theta_0A = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_0");
	TH2D *dphi_vs_p_vs_theta_1A = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_1");
	TH2D *dphi_vs_p_vs_theta_2A = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_2");
	TH2D *dphi_vs_p_vs_theta_chi2A = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_chi2");
	TH2D *eff_vs_p_vs_thetaA = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta");
	TH2D *eff_vs_p_vs_theta_numeratorA = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta_numerator");
	TH2D *eff_vs_p_vs_theta_denominatorA = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta_denominator");

	// Open fileB and get pointers to histograms
	TFile *f = new TFile(fileB);
	TH3D *dp_over_p_vs_p_vs_thetaB = (TH3D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta");
	TH2D *dp_over_p_vs_p_vs_theta_0B = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_0");
	TH2D *dp_over_p_vs_p_vs_theta_1B = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_1");
	TH2D *dp_over_p_vs_p_vs_theta_2B = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_2");
	TH2D *dp_over_p_vs_p_vs_theta_chi2B = (TH2D*)gROOT->FindObject("dp_over_p_vs_p_vs_theta_chi2");
	TH3D *dpt_over_pt_vs_p_vs_thetaB = (TH3D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta");
	TH2D *dpt_over_pt_vs_p_vs_theta_0B = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_0");
	TH2D *dpt_over_pt_vs_p_vs_theta_1B = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_1");
	TH2D *dpt_over_pt_vs_p_vs_theta_2B = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_2");
	TH2D *dpt_over_pt_vs_p_vs_theta_chi2B = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_chi2");
	TH3D *dtheta_vs_p_vs_thetaB = (TH3D*)gROOT->FindObject("dtheta_vs_p_vs_theta");
	TH2D *dtheta_vs_p_vs_theta_0B = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta_0");
	TH2D *dtheta_vs_p_vs_theta_1B = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta_1");
	TH2D *dtheta_vs_p_vs_theta_2B = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta_2");
	TH2D *dtheta_vs_p_vs_theta_chi2B = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta_chi2");
	TH3D *dphi_vs_p_vs_thetaB = (TH3D*)gROOT->FindObject("dphi_vs_p_vs_theta");
	TH2D *dphi_vs_p_vs_theta_0B = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_0");
	TH2D *dphi_vs_p_vs_theta_1B = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_1");
	TH2D *dphi_vs_p_vs_theta_2B = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_2");
	TH2D *dphi_vs_p_vs_theta_chi2B = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_chi2");
	TH2D *eff_vs_p_vs_thetaB = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta");
	TH2D *eff_vs_p_vs_theta_numeratorB = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta_numerator");
	TH2D *eff_vs_p_vs_theta_denominatorB = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta_denominator");

	// Open output file
	TFile *fout = new TFile("hd_res_charged.root","RECREATE","Created by MergeResFiles.C");
	
	// Merge histos into a new histo
	MergeHistograms(dp_over_p_vs_p_vs_thetaA, dp_over_p_vs_p_vs_thetaB, theta_border);	
	MergeHistograms(dp_over_p_vs_p_vs_theta_0A, dp_over_p_vs_p_vs_theta_0B, theta_border);	
	MergeHistograms(dp_over_p_vs_p_vs_theta_1A, dp_over_p_vs_p_vs_theta_1B, theta_border);	
	MergeHistograms(dp_over_p_vs_p_vs_theta_2A, dp_over_p_vs_p_vs_theta_2B, theta_border);	
	MergeHistograms(dp_over_p_vs_p_vs_theta_chi2A, dp_over_p_vs_p_vs_theta_chi2B, theta_border);	
	MergeHistograms(dpt_over_pt_vs_p_vs_thetaA, dpt_over_pt_vs_p_vs_thetaB, theta_border);	
	MergeHistograms(dpt_over_pt_vs_p_vs_theta_0A, dpt_over_pt_vs_p_vs_theta_0B, theta_border);	
	MergeHistograms(dpt_over_pt_vs_p_vs_theta_1A, dpt_over_pt_vs_p_vs_theta_1B, theta_border);	
	MergeHistograms(dpt_over_pt_vs_p_vs_theta_2A, dpt_over_pt_vs_p_vs_theta_2B, theta_border);	
	MergeHistograms(dpt_over_pt_vs_p_vs_theta_chi2A, dpt_over_pt_vs_p_vs_theta_chi2B, theta_border);	
	MergeHistograms(dtheta_vs_p_vs_thetaA, dtheta_vs_p_vs_thetaB, theta_border);	
	MergeHistograms(dtheta_vs_p_vs_theta_0A, dtheta_vs_p_vs_theta_0B, theta_border);	
	MergeHistograms(dtheta_vs_p_vs_theta_1A, dtheta_vs_p_vs_theta_1B, theta_border);	
	MergeHistograms(dtheta_vs_p_vs_theta_2A, dtheta_vs_p_vs_theta_2B, theta_border);	
	MergeHistograms(dtheta_vs_p_vs_theta_chi2A, dtheta_vs_p_vs_theta_chi2B, theta_border);	
	MergeHistograms(dphi_vs_p_vs_thetaA, dphi_vs_p_vs_thetaB, theta_border);	
	MergeHistograms(dphi_vs_p_vs_theta_0A, dphi_vs_p_vs_theta_0B, theta_border);	
	MergeHistograms(dphi_vs_p_vs_theta_1A, dphi_vs_p_vs_theta_1B, theta_border);	
	MergeHistograms(dphi_vs_p_vs_theta_2A, dphi_vs_p_vs_theta_2B, theta_border);	
	MergeHistograms(dphi_vs_p_vs_theta_chi2A, dphi_vs_p_vs_theta_chi2B, theta_border);	
	MergeHistograms(eff_vs_p_vs_thetaA, eff_vs_p_vs_thetaB, theta_border);	
	MergeHistograms(eff_vs_p_vs_theta_numeratorA, eff_vs_p_vs_theta_numeratorB, theta_border);	
	MergeHistograms(eff_vs_p_vs_theta_denominatorA, eff_vs_p_vs_theta_denominatorB, theta_border);	
	
	// Close output file
	fout->Write();
	delete fout;
}

TH3D* MergeHistograms(TH3D *hA, TH3D *hB, double theta_border)
{
	// Make clone of histogram B
	TH3D *h = hB->Clone(hA->GetName());
	
	TAxis *xaxis = h->GetXaxis();
	double binwidth = xaxis->GetBinWidth(1);
	
	// Loop over X (theta) bins
	for(int xbin=1; xbin<xaxis->GetNbins(); xbin++){
		if(xaxis->GetBinLowEdge(xbin)+binwidth > theta_border)break;
		
		// Loop over Y bins
		TAxis *yaxis = h->GetYaxis();
		for(int ybin=1; ybin<yaxis->GetNbins(); ybin++){
			// Loop over Z bins
			TAxis *zaxis = h->GetZaxis();
			for(int zbin=1; zbin<zaxis->GetNbins(); zbin++){
				h->SetBinContent(xbin, ybin, zbin, hA->GetBinContent(xbin, ybin,zbin));
				h->SetBinError(xbin, ybin, zbin, hA->GetBinError(xbin, ybin,zbin));
			}		
		}
	}

}

TH2D* MergeHistograms(TH2D *hA, TH2D *hB, double theta_border)
{
	TH2D *h = hB->Clone(hA->GetName());

	TAxis *xaxis = h->GetXaxis();
	double binwidth = xaxis->GetBinWidth(1);
	
	// Loop over X (theta) bins
	for(int xbin=1; xbin<xaxis->GetNbins(); xbin++){
		if(xaxis->GetBinLowEdge(xbin)+binwidth > theta_border)break;
		
		// Loop over Y bins
		TAxis *yaxis = h->GetYaxis();
		for(int ybin=1; ybin<yaxis->GetNbins(); ybin++){
			h->SetBinContent(xbin, ybin, hA->GetBinContent(xbin, ybin));
			h->SetBinError(xbin, ybin, hA->GetBinError(xbin, ybin));
		}
	}
}


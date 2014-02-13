void trkeff(int pid)
{
	gROOT->Reset();

	Double_t Red[5]   = { 0.5,1.00,,0,0, 0 };
	Double_t Green[5] = { 0.1, 0.3, 1., 0., 0.};
	Double_t Blue[5]  = { 0.1, 0.35, 0.65,0.8,1.0 };	
	Double_t Length[5] = { 0.00, 0.3, 0.5, 0.85, 1.00 };
	Int_t nb=50;
	TColor::CreateGradientColorTable(5,Length,Red,Green,Blue,nb);

	// If we are not already in the TRACKING directory, cd there
	TDirectory *TRACKING = (TDirectory*)gROOT->FindObject("TRACKING");
	if(TRACKING)TRACKING->cd();

	// Get pointer to trkeff tree
	TTree *trkeff = (TTree*)gROOT->FindObject("trkeff");

	// Create histograms to calculate efficiency
	TH2D *np = new TH2D("np","Numerator", 70,0,140,70,-0.05,6.95);
	TH2D *dp = new TH2D("dp","Denominator",70,0,140,70,-0.05,6.95);	
	TH2D *npcut = new TH2D("npcut","numerator",70,0,140,70,-0.05,6.95);

	trkeff->Project("np", "F.pthrown.Mag():180./3.14159*F.pthrown.Theta()", "F.trktb.Nfdc>0 || F.trktb.Ncdc>0");
	trkeff->Project("dp", "F.pthrown.Mag():180/3.14159*F.pthrown.Theta()", "");
	trkeff->Project("npcut", "F.pthrown.Mag():180./3.14159*F.pthrown.Theta()", "(F.trktb.Nfdc>0 || F.trktb.Ncdc>0) && TMath::Prob(F.trktb.trk_chisq,F.trktb.trk_Ndof)>0.01");
	
	TH2D *effp = np->Clone("effp");
	
	TCanvas *c1 = new TCanvas("c1");
	c1->SetTickx();
	c1->SetTicky();
	effp->Divide(dp);
	effp->SetMarkerStyle(2);
	effp->SetStats(0);
	//gStyle->SetPalette(1,0);
	effp->SetContour(nb);
	if (pid==9)
	  effp->SetTitle("Pion reconstruction efficiency");
	else{
	  effp->SetTitle("Proton reconstruction efficiency");
	  effp->SetAxisRange(0,88,"X");
	  effp->SetAxisRange(0,2,"Y");
	}
	effp->SetXTitle("#theta [degrees]");
	effp->SetYTitle("p [GeV/c]");

	// Find the region where the efficiency drops through 50% and fit 
	// the resulting ridge with a polynomial.
	TH2D *ridge = new TH2D("ridge","ridge", 70,0,140,70,-0.05,6.95);
	for (unsigned int i=2;i<65;i++){
	  for (unsigned int j=2;j<70;j++){
	    double frac=effp->GetBinContent(i,j);	  
	    if (frac>0.0){
	      double min_weight=exp(-(0.8-0.5)*(0.8-0.5)/1);
	      double weight=exp(-(frac-0.5)*(frac-0.5)/1);
	      if (weight>min_weight){
		ridge->SetBinContent(i,j,1000.*(weight-min_weight));
	   
	      }
	    }
	  }
	}
	ridge->Fit("pol3");
	effp->Draw("colz");
	double par[4];
	TF1 *f1=ridge->GetFunction("pol3");
	f1->GetParameters(par);
	f1->SetRange(0.95,130.5);
	f1->SetLineColor(5);
	f1->Draw("same");
	
	double pmax=6.95;
	if (pid==14) pmax=2.05;
	double p1=par[0]+par[1]+par[2]+par[3];
	double p2=par[0]+par[1]*130.+par[2]*130.*130.+par[3]*130.*130.*130.;
	TLine *line1 = new TLine(1,p1,1,pmax);
	line1->SetLineColor(5);
	line1->SetLineWidth(2);
	line1->Draw();
	if (pid==9){
	  TLine *line2 = new TLine(130,p2,130,pmax);
	  line2->SetLineColor(5);
	  line2->SetLineWidth(2);
	  line2->Draw();
	}
	if (pid==9) c1->SaveAs("pion_efficiency.png");
	else c1->SaveAs("proton_efficiency.png");


	TH2D *effpcut = npcut->Clone("effpcut");
	
	TCanvas *c2 = new TCanvas("c2");
	c2->SetTickx();
	c2->SetTicky();
	effpcut->Divide(dp);
	effpcut->SetMarkerStyle(2);
	effpcut->SetStats(0);


	effpcut->SetContour(nb);
	if (pid==9)
	  effpcut->SetTitle("Pion reconstruction efficiency, P(#chi^2)>0.01");
	else{
	  effpcut->SetTitle("Proton reconstruction efficiency, P(#chi^2)>0.01");
	  effpcut->SetAxisRange(0,88,"X");
	  effpcut->SetAxisRange(0,2,"Y");
	}
	effpcut->SetXTitle("#theta [degrees]");
	effpcut->SetYTitle("p [GeV/c]");
	effpcut->Draw("colz");
	f1->Draw("same");
	line1->Draw();
	if (pid==9) line2->Draw();

	if (pid==9) c2->SaveAs("pion_efficiency_with_chi2_cut.png");
	else c2->SaveAs("proton_efficiency_with_chi2_cut.png");

	if (pid==9){
	  ofstream ofile("pion.txt");
	  ofile << "<html>" <<endl;
	  ofile<<"<h1>Single track reconstruction</h1>" <<endl;
	  ofile<<"<p> Yellow lines indicate recommended fiducial cuts </p>" 
	    <<endl;
	  ofile<<"<h2>Pion events</h6>" <<endl;
	}
	else{
	  ofstream ofile("proton.txt");
	  ofile<<"<h2>Proton events</h6>" <<endl;
	}
	double eff = np->GetEntries()/dp->GetEntries();
	ofile <<"<p>Overall efficiency = "<<eff<<"<br/>"<<endl;
	double effcut = npcut->GetEntries()/dp->GetEntries();
	ofile <<"  Efficiency (Prob>0.01) = "<<effcut<<"<br/>" << endl;
	ofile << "Momentum cut: p=" << par[0] << "&ensp;+&ensp;" 
	      << par[1] << "&thinsp;&theta;&ensp;+&ensp;" 
	      << par[2] << "&thinsp;&theta;<sup>2</sup>&ensp;+&ensp;" << par[3]
	      <<"&thinsp;&theta;<sup>3</sup></p>" <<endl;
	if (pid==9){
	  ofile <<"<img src=\"pion_efficiency.png\">" <<endl;
	  ofile <<"<img src=\"pion_efficiency_with_chi2_cut.png\">" <<endl;
	}
	else{
	  ofile <<"<img src=\"proton_efficiency.png\">" <<endl;
	  ofile <<"<img src=\"proton_efficiency_with_chi2_cut.png\">" <<endl;
	  ofile << "</html>" <<endl;
	}
	
	ofile.close();
}

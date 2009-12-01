


void gen_covariance_code(const char *fname="hd_root.root")
{
	gROOT->Reset();
	
	TFile *f = new TFile(fname);
	f->cd("TRACKING");
	TH2D *fdc_cov = (TH2D*)gROOT->FindObject("fdc_cov");
	TH2D *fdc_cath_cov = (TH2D*)gROOT->FindObject("fdc_cath_cov");
	
	// Print C++ routine for returning covariance of drift-time
	// residuals for FDC due to MULS as a function of two FDC
	// layers.
	ostream &mout = cout; // Maybe we want this to go straight to a file?
	
	mout<<endl;
	mout<<endl;
	mout<<"// The following was auto-generated from the gen_covariance_code.C"<<endl;
	mout<<"// macro using a ROOT file generated with the fdc_covariance_tree"<<endl;
	mout<<"// plugin."<<endl;
	mout<<endl;
	mout<<endl;
	mout<<"double GetFDCCovariance(int layer1, int layer2);"<<endl;
	mout<<"double GetFDCCathodeCovariance(int layer1, int layer2);"<<endl;
	mout<<endl;
	mout<<endl;
	mout<<"//-------------------------"<<endl;
	mout<<"// GetFDCCovariance"<<endl;
	mout<<"//-------------------------"<<endl;
	mout<<"double GetFDCCovariance(int layer1, int layer2)"<<endl;
	mout<<"{"<<endl;
	mout<<"	if(layer1<1 || layer2>24 || layer2<1 || layer2>24)return 0.0;"<<endl;
	mout<<"	if(layer2<layer1){"<<endl;
	mout<<"		int tmp = layer1;"<<endl;
	mout<<"		layer1 = layer2;"<<endl;
	mout<<"		layer2 = tmp;"<<endl;
	mout<<"	}"<<endl;
	mout<<endl;
	mout<<"	switch(layer1){"<<endl;
	for(int layer1=1; layer1<=24; layer1++){
		mout<<"		case "<<layer1<<":"<<endl;
		for(int layer2=layer1; layer2<=24; layer2++){
			mout<<"			if(layer2=="<<layer2<<")return "<<fdc_cov->GetBinContent(layer1, layer2)<<";"<<endl;
		}
		mout<<"			break; // layer "<<layer1<<endl;
	}
	mout<<"	} // switch for layer1 "<<layer1<<endl;
	mout<<endl;
	mout<<"	return 0.0;"<<endl;
	mout<<"}"<<endl;

	mout<<endl;
	mout<<endl;
	mout<<"// NOTE: At this point, the following is incorrect. It does not"<<endl;
	mout<<"// Properly account for the Lorentz deflections that are in the"<<endl;
	mout<<"// simulated data file. It's better to use the above for the"<<endl;
	mout<<"// MULS error along the wire as well for the time being."<<endl;
	mout<<"//    April 20, 2009  DL"<<endl;
	mout<<endl;
	mout<<endl;
	mout<<"//-------------------------"<<endl;
	mout<<"// GetFDCCathodeCovariance"<<endl;
	mout<<"//-------------------------"<<endl;
	mout<<"double GetFDCCathodeCovariance(int layer1, int layer2)"<<endl;
	mout<<"{"<<endl;
	mout<<"	if(layer1<1 || layer2>24 || layer2<1 || layer2>24)return 0.0;"<<endl;
	mout<<"	if(layer2<layer1){"<<endl;
	mout<<"		int tmp = layer1;"<<endl;
	mout<<"		layer1 = layer2;"<<endl;
	mout<<"		layer2 = tmp;"<<endl;
	mout<<"	}"<<endl;
	mout<<endl;
	mout<<"	switch(layer1){"<<endl;
	for(int layer1=1; layer1<=24; layer1++){
		mout<<"		case "<<layer1<<":"<<endl;
		for(int layer2=layer1; layer2<=24; layer2++){
			mout<<"			if(layer2=="<<layer2<<")return "<<fdc_cath_cov->GetBinContent(layer1, layer2)<<";"<<endl;
		}
		mout<<"			break; // layer "<<layer1<<endl;
	}
	mout<<"	} // switch for layer1 "<<layer1<<endl;
	mout<<endl;
	mout<<"	return 0.0;"<<endl;
	mout<<"}"<<endl;
}




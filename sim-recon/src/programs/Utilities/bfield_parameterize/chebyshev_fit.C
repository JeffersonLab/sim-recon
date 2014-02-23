
// Chebyshev polynomial functions
Double_t T0(Double_t x){return 1;}
Double_t T1(Double_t x){return x;}
Double_t T2(Double_t x){return 2*pow(x,2)-1;}
Double_t T3(Double_t x){return 4*pow(x,3)-3*x;}
Double_t T4(Double_t x){return 8*pow(x,4)-8*pow(x,2)+1;}
Double_t T5(Double_t x){return 16*pow(x,5)-20*pow(x,3)+5*x;}
Double_t T6(Double_t x){return 32*pow(x,6)-48*pow(x,4)+18*pow(x,2)-1;}
Double_t T7(Double_t x){return 64*pow(x,7)-112*pow(x,5)+56*pow(x,3)-7*x;}
Double_t T8(Double_t x){return 128*pow(x,8)-256*pow(x,6)+160*pow(x,4)-32*pow(x,2)+1;}
Double_t T9(Double_t x){return 256*pow(x,9)-576*pow(x,7)+432*pow(x,5)-120*pow(x,3)+9*x;}


//---------------------
// chebyshev_fit
//---------------------
void chebyshev_fit(void)
{
	gROOT->Reset();
	TColor::CreateColorWheel();
	
	ChebyshevTestFit();
}

//---------------------
// ChebyshevTestFit
//---------------------
void ChebyshevTestFit()
{
	TCanvas *c1 = new TCanvas("c1");
	c1->SetGrid();
	c1->SetTicks();
	
	TH1D *h = new TH1D("h","A test histogram", 1000, -1.0, 1.0);
	h->FillRandom("gaus", 1000000);
	
	TF1 *f = chebyshev_FindBestFunction(h);

	h->Draw();
	f->SetLineColor(kRed);
	f->Draw("same");
	
	UInt_t order = f->GetNpar() - 1;
	cout<<"Chebyshev polynomial fit to data of order: "<<order<<endl;
}

//---------------------
// chebyshev_Fit
//---------------------
TF1* chebyshev_Fit(TH1D *h, UInt_t order=9)
{
	double xmin = h->GetXaxis()->GetXmin();
	double xmax = h->GetXaxis()->GetXmax();
	TF1 *f = chebyshev_MakeTF1(xmin, xmax, order);

	h->Fit(f, "0Q");
	return f;
}

//---------------------
// chebyshev_FindBestFunction
//---------------------
TF1* chebyshev_FindBestFunction(TH1D *h, UInt_t max_order=9)
{
	// Verify histogram limits
	if(h->GetXaxis()->GetXmin()!=-1.0 || h->GetXaxis()->GetXmax()!=+1.0){
		cerr<<"-- ERROR! The limits of the histo passed to chebyshev_FindBestFunction"<<endl;
		cerr<<"-- ERROR! MUST be -1.0 to +1.0 ("<<h->GetXaxis()->GetXmin()<<" and "<<h->GetXaxis()->GetXmax()<<" passed)"<<endl;
	}

	// Loop over all orders keeping track of the one with the best chisq
	UInt_t best_order = 0;
	Double_t best_chisq_per_dof = 1.0E6;
	for(UInt_t order=0; order<=max_order; order++){
		TF1 *f = chebyshev_MakeTF1(order);
		h->Fit(f, "0Q");
		Double_t chisq_per_dof = f->GetChisquare()/(double)f->GetNDF();
		cout<<"order "<<order<<": chisq="<<f->GetChisquare()<<"  NDF="<<f->GetNDF()<<"  chisq/NDF="<<chisq_per_dof<<endl;
		if(chisq_per_dof<best_chisq_per_dof){
			best_chisq_per_dof = chisq_per_dof;
			best_order = order;
		}
	}
	
	// Regenerate and fit with best function
	TF1 *best_func = chebyshev_MakeTF1(best_order);
	h->Fit(best_func, "0Q");
	
	return best_func;
}

//---------------------
// chebyshev_MakeTF1
//---------------------
TF1* chebyshev_MakeTF1(Double_t xmin, Double_t xmax, UInt_t order)
{
	// ROOT gives a "Too many operators !" error when we try and define
	// past order 9. Truncate it with a warning if too high of an order
	// is specified.
	if(order>9){
		cerr<<"--WARNING! The maximum order that can be used is 9"<<endl;
		cerr<<"--WARNING! order="<<order<<" will be reduced to 9"<<endl;
		order = 9; 
	}

	// Build function string of given order
	stringstream func_str;
	for(UInt_t i=0; i<=order; i++){
		if(i!=0)func_str<<"+";
		func_str<<"["<<i+2<<"]*T"<<i<<"((x-[0])/[1])";
	}

	// Create Function
	stringstream funcname;
	funcname<<"polT"<<order;
	TF1 *myfunc = new TF1(funcname.str().c_str(), func_str.str().c_str(), xmin, xmax);

	// Parmeters used to convert x to -1 to 1 coordinates
	double xmid  = (xmax+xmin)/2.0;
	double xnorm = (xmax-xmin)/2.0;

	// Set initial parameter values and parameter names
	myfunc->SetParName(0, "xmid");
	myfunc->SetParName(1, "xnorm");
	myfunc->FixParameter(0, xmid);
	myfunc->FixParameter(1, xnorm);
	for(UInt_t i=0; i<=order; i++){
		stringstream fname;
		fname<<"T"<<i;
		myfunc->SetParName(i+2, fname.str().c_str());
		myfunc->SetParameter(i+2, 1.0);
	}

	return myfunc;
}





#include <TFile.h>
#include <TH1.h>

// FORTRAN routines
extern "C"{
void cobrems_(float* Emax, float* Epeak, float* Dist);
float dntdx_(float* x);
};

// Wrapper function for total
double dNtdx(double x)
{
	float xx = x;
	return (double)dntdx_(&xx);
}

int main(int narg, char* argv[])
{
	// Initialize coherent brem table
	float Emax = 12.0;
	float Epeak = 9.0;
	float Dist = 76.0;
	cobrems_(&Emax, &Epeak, &Dist);

	// Open ROOT file
	TFile *f = new TFile("cobrems.root","RECREATE");
	
	// Create histogram
	TH1D *cobrem_vs_E = new TH1D("cobrem_vs_E", "Coherent Bremstrahlung vs. E_{#gamma}", 100, 0.0, 12.0);
	
	// Fill histogram
	for(int i=1; i<=cobrem_vs_E->GetNbinsX(); i++){
		double x = cobrem_vs_E->GetBinCenter(i)/Emax;
		double y = dNtdx(x);
		cobrem_vs_E->SetBinContent(i, y);
	}
	
	// Close out ROOT file
	f->Write();
	delete f;

	return 0;
}



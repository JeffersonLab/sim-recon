
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TString.h"
#include "TF1.h"

void fit_ADC_4ns(TString filename, int runnumber) {

    bool writeout=0;

    TFile *infile = TFile::Open(filename.Data());
    if (!infile) {
        printf("Could not find %s\n",filename.Data());
        return;
    }
        
    char name[255];
    TH2F* ZvsDeltat = NULL;
    
    TF1* pol1f1 = new TF1("pol1f1","212 + [0]*(x - [1])",-30,30);

    char outfilename[255];
    sprintf(outfilename,"BCAL_ADC_4ns_correction_%i.txt",runnumber);
    FILE *outfile = fopen(outfilename, "w");

    for(int module = 0; module < 48; ++module){
        for(int layer = 0; layer < 4 ; ++layer){
            for(int sector = 0; sector < 4; ++sector){

                // load 2D histogram
                sprintf(name,"BCAL_ADC_Deltat/ZvsDeltat/M%02dL%dS%d", module+1, layer+1, sector+1);
                ZvsDeltat = (TH2F*)infile->Get(name);
                int correction = 0;
                if (!ZvsDeltat) {
                    printf("M%02dL%dS%d  missing\n",module+1, layer+1, sector+1);
                } else {
                    ZvsDeltat->Fit(pol1f1,"Q");

                    float p0 = pol1f1->GetParameter(0);
                    float dp0 = pol1f1->GetParError(0);
                    float p1 = pol1f1->GetParameter(1);
                    float dp1 = pol1f1->GetParError(1);
                    printf("M%02dL%dS%d  %6.2f +-%6.2f %6.2f +-%6.2f\n", 
                           module+1, layer+1, sector+1, p0, dp0, p1, dp1);
                    if (p1>2) correction = 4;
                }
                fprintf(outfile,"0\n%i\n",correction);
            } 
        }
    }
    fclose(outfile);
}

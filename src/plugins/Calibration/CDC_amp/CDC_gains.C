void CDC_gains(int EXIT_EARLY=0) {

  // set EXIT_EARLY to 1 to fit the sum histograms only, 0 for the individual straw gains

  // Fits histograms from CDC_amp plugin to estimate CDC gains
  // writes digi_scales/ascale to cdc_new_ascale.txt
  // writes wire_gains to cdc_new_wiregains.txt
  // writes fitted histograms to cdc_amphistos.root
  //
  // wire_gains are gain correction constants for individual straws, with respect to sum histogram
  // fit landau to histogram for sum of all straws and histograms for each individual straw
  // set wire_gains to scale individual fit mpvs to the mpv for the sum of all straws histogram.

  // digi_scales/ascale is the overall gain correction factor for the chamber
  // Reference values of IDEALMPV and ASCALE are chosen from a high gain (low pressure) run 11621  
  // The new ascale is ASCALE*IDEALMPV/sum fit mpv
  //
  // (For the first round of gain calibrations, ASCALE was chosen to match the simulation charge histogram.
  // For the next round, ASCALE was chosen to match the dedx histograms from the first round).

  // CDC_amp provides 3 2D histograms of pulse amplitude vs straw
  // The histograms are projected here into groups of 1D histograms, one per straw
  //
  // asum:   igroup 0 : all hits
  // atsum:  igroup 1 : tracked hits, use these if not enough stats to use group 2
  // attsum: igroup 2 : tracked hits for z=50-80 and theta 85-95.

  // attsum provide the most Landau-like fits

  // Need 80+ run files to get good fits for attsum.  2 run files are enough for atsum

  // NSJ 18 Nov 2016



  // reference values

  const float IDEALMPV=32.385; //mpv for tracked hits with restricted z,theta from low pressure run 11621
  const float ASCALE=0.176;   //ascale for tracked ztheta hits from 011621

  const int USEGROUP=2; // use group 2 (attsum) to calc gain consts

  // fit limits
  const int AL0=26; //amp fit lower limit asum (untracked hits histo)
  const int AL1=20; //amp fit lower limit atsum
  const int AL2=20; //amp fit lower limit attsum
  const int AH=400; //amp fit upper limit landau
  const int AH0=60; //amp fit upper limit gaussian

  const int MINCOUNTS=5000; //counts required to fit histogram
  

  gStyle->SetOptFit(1);
  gStyle->SetFuncWidth(1);

  TF1 *f = new TF1("f","landau");
  f->SetLineColor(6);

  TF1 *g = new TF1("g","gaus");
  g->SetLineColor(434);
  g->SetRange(AL0,AH0);

  int a_fitstat,q_fitstat;

  float newascale = 0;

  float thismpv = 0;

  // get histograms
   
  TDirectory *fmain = (TDirectory*)gDirectory->FindObjectAny("CDC_amp");
  if (!fmain) printf("Cannot find directory CDC_amp\n"); 
  if (!fmain) return;
  fmain->cd();


  TH1I *asum = (TH1I*)fmain->Get("asum");
  TH1I *atsum = (TH1I*)fmain->Get("atsum");
  TH1I *attsum = (TH1I*)fmain->Get("attsum");


  new TCanvas;
  printf("\nasum fit, all hits:\n");
  a_fitstat = asum->Fit(g,"R");
  if (a_fitstat) printf("Bad fit to asum\n\n");

  new TCanvas;
  printf("\natsum fit, tracked hits:\n");
  f->SetRange(AL1,AH);
  a_fitstat = atsum->Fit(f,"R");
  if (a_fitstat) printf("Bad fit to asum\n\n");

  if (USEGROUP==1 && !a_fitstat) {
    thismpv = f->GetParameter(1);
    newascale = ASCALE*IDEALMPV/thismpv;
    printf("\nnew digi_scales/ascale should be %.3f\n\n",newascale);
  } 

  new TCanvas;
  printf("\nattsum fit, tracked hits, restricted z and theta:\n");
  f->SetRange(AL2,AH);
  a_fitstat = attsum->Fit(f,"R");

  if (USEGROUP==2 && !a_fitstat) {
    thismpv = f->GetParameter(1);
    newascale = ASCALE*IDEALMPV/thismpv;
    printf("\nnew digi_scales/ascale should be %.3f\n\n",newascale);
  }

  if (a_fitstat) printf("\nSum histogram fit error, exiting script\n");
  if (a_fitstat) return;


  FILE *outfile = fopen("cdc_new_ascale.txt","w");
  fprintf(outfile,"%.3f 0.8\n",newascale); //0.8 is the time scale, for in the 2nd column of ccdb's CDC/digi_scales table
  fclose(outfile);


  TH1D *qsum = (TH1D*)fmain->Get("qsum");
  TH1D *qtsum = (TH1D*)fmain->Get("qtsum");
  TH1D *qttsum = (TH1D*)fmain->Get("qttsum");


  TFile *hfile = new TFile("cdc_amphistos.root","RECREATE");

  asum->Write();
  atsum->Write();
  attsum->Write();

  qsum->Write();
  qtsum->Write();
  qttsum->Write();

  if (EXIT_EARLY==1) hfile->Write();
  if (EXIT_EARLY==1) hfile->Close();
  if (EXIT_EARLY==1) return;


  TTree *atstats = new TTree("atstats","fit stats for tracked hits");
  TTree *attstats = new TTree("attstats","fit stats for tracked hits with restricted z and theta");

  int n,ring;  //straw number (1-3522), ring number (1-28)

  int a_n;  //count
  double a_mean;
  double a_c;  //const
  double a_mpv;
  double a_mpverr;
  double a_sig;
  double a_chisq;


  atstats->Branch("n",&n,"n/I");
  atstats->Branch("hits",&a_n,"hits/I");
  atstats->Branch("mean",&a_mean,"mean/D");
  atstats->Branch("c",&a_c,"c/D");
  atstats->Branch("mpv",&a_mpv,"mpv/D");
  atstats->Branch("mpverr",&a_mpverr,"mpverr/D");
  atstats->Branch("sig",&a_sig,"sig/D");
  atstats->Branch("chisq",&a_chisq,"chisq/D");
  atstats->Branch("fitstat",&a_fitstat,"fitstat/I");

  attstats->Branch("n",&n,"n/I");
  attstats->Branch("hits",&a_n,"hits/I");
  attstats->Branch("mean",&a_mean,"mean/D");
  attstats->Branch("c",&a_c,"c/D");
  attstats->Branch("mpv",&a_mpv,"mpv/D");
  attstats->Branch("mpverr",&a_mpverr,"mpverr/D");
  attstats->Branch("sig",&a_sig,"sig/D");
  attstats->Branch("chisq",&a_chisq,"chisq/D");
  attstats->Branch("fitstat",&a_fitstat,"fitstat/I");

  int i,igroup;  

  TH1D *ahisto;
  TH2I *anhisto;
  char htitle[300];

  int fitlowerlimit;
  int badfit, nofit; //counters
  int bincont, lastbin, nbins; // to flag low-gain preamp channels

  float wiregain;



  hfile->cd();

  TDirectory *hmain = gDirectory;

  outfile = fopen("cdc_new_wiregains.txt","w");

  new TCanvas;
  
  for (igroup=1; igroup<3; igroup++) {   //untracked individual straws are not fittable for inner rings

    printf("\n-----------  igroup %i -----------\n\n",igroup);

    hmain->cd();

    //if (igroup==0) gDirectory->mkdir("all")->cd();
    if (igroup==1) gDirectory->mkdir("tracked")->cd();
    if (igroup==2) gDirectory->mkdir("theta")->cd();

    TDirectory *hsub = gDirectory;

    fmain->cd();

    //if (igroup==0) anhisto = (TH2I*)fmain->Get("an");
    if (igroup==1) anhisto = (TH2I*)fmain->Get("atn");
    if (igroup==2) anhisto = (TH2I*)fmain->Get("attn");

    //if (igroup==0) fitlowerlimit = AL0;
    if (igroup==1) fitlowerlimit = AL1;
    if (igroup==2) fitlowerlimit = AL2;
    
    f->SetRange(fitlowerlimit,AH);

    badfit = 0;
    nofit = 0;


    for (i=1; i<3523; i++) {   

      a_c = 0;
      a_mpv = 0;
      a_mpverr = 0;
      a_sig = 0;
      a_chisq = 0;
      a_fitstat = 0;

      n = i;

      ahisto = anhisto->ProjectionY(Form("a[%i]",i),i+1,i+1);  //straw 1 is in bin 2

      //if (igroup==0) sprintf(htitle,"Amplitude, straw %i; amplitude-pedestal",i);
      if (igroup==1) sprintf(htitle,"Amplitude, tracked hits, straw %i; amplitude-pedestal",i);
      if (igroup==2) sprintf(htitle,"Amplitude, tracked hits, z=50 to 80 cm, theta=85 to 95 degrees, straw %i; amplitude-pedestal",i);

      ahisto->SetTitle(htitle);

      a_n = ahisto->GetEntries();
      a_mean = ahisto->GetMean();
 
      a_fitstat = -1;
      if (a_n > MINCOUNTS) {
        a_fitstat = ahisto->Fit(f,"RQ");  //landau   LL is for low stats
      }

      if (a_fitstat==0) { 

        a_c      = f->GetParameter(0);
        a_mpv    = f->GetParameter(1);
        a_mpverr = f->GetParError(1);
        a_sig    = f->GetParameter(2);

        if (f->GetNDF()>0) a_chisq  = f->GetChisquare()/f->GetNDF();

      } else if (a_fitstat>0){
        printf("straw %i fitstatus %i fit MPV %f\n",i,a_fitstat,f->GetParameter(1));
      } else if (a_fitstat<0){
        printf("straw %i not enough counts \n",i);
        nofit++;
      } else if (a_fitstat==4){
        printf("straw %i unconverged fit \n",i);
      }

      if (a_fitstat>0) badfit++;
    
      // flag low-gain preamp channels
      if (a_n > MINCOUNTS) {    // sum last 10 bins before amp histo usually ends 
        nbins = ahisto->GetNbinsX();
        bincont = 0;
        for (int j=nbins-50; j < nbins; j++) bincont += ahisto->GetBinContent(j);
        if (bincont==0) {
          lastbin = nbins-1;
          while (ahisto->GetBinContent(lastbin)==0) lastbin--;
        }
      }
    
      if (a_mpv<0) printf("straw %i fit MPV %f\n",i,a_mpv);
      if (a_mpv>0 && a_mpv<fitlowerlimit) printf("straw %i MPV %f below fit lower limit\n",i,a_mpv);
      if (a_mpv>0 && a_mean>3.3*a_mpv) printf("straw %i mean %f MPV %f\n",i,a_mean,a_mpv);
      if (bincont==0 && a_n>MINCOUNTS) printf("straw %i no amplitudes above %i\n",i,lastbin*(int)ahisto->GetBinWidth(0));


      //if (igroup==0) astats->Fill();
      if (igroup==1) atstats->Fill();
      if (igroup==2) attstats->Fill();

      if (igroup==USEGROUP) {
        wiregain=0;
        if (a_mpv>0) wiregain = thismpv/a_mpv;
        fprintf(outfile,"%.3f\n",wiregain);
      }

      hsub->cd();

      ahisto->Write();

    }

    printf("\nfitstatus 4 count: %i \n",badfit);
    printf("no fit count: %i \n",nofit);

  }

  fclose(outfile);

  hfile->Write();

}

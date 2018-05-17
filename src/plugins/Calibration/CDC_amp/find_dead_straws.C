

void find_dead_straws(void) {

  // run this on CDC_amp output root file containing histo attn

  // scans 2d histogram of tracked hits from the target attn 
  // identifies straws with very few hits
  // writes text file straw_eff.txt containing 1 for live straws, 0 for dead straws

 
  TDirectory *fmain = (TDirectory*)gDirectory->FindObjectAny("CDC_amp");
  if (!fmain) printf("Cannot find directory CDC_amp\n"); 
  if (!fmain) return;
  fmain->cd();
 
  TH2I *anhisto = (TH2I*)fmain->Get("attn");
  if (!anhisto) printf("Cannot find histo attn\n");
  if (!anhisto) return;


  FILE *outfile = fopen("straw_eff.txt","w");


  int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
 
  //  int Ntotal[28] = {42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313,3522};


  TH1D *ahisto;

  int Nfirst = 0;


  for (int i=0; i<28; i++) {
    
   ahisto = anhisto->ProjectionY(Form("a[%i]",i),Nfirst+2,Nfirst+Nstraws[i]+1);  //straw 1 is in bin 2
   //   cout << "ring " << i+1 << " straws " << Nfirst+1 << " to " << Nfirst+Nstraws[i] << endl;

   int nring = ahisto->GetEntries();
   
   int mean_n = (int)(nring/(float)Nstraws[i]);

   //   cout << "Ring " << i << " total counts " << nring << " mean per straw " << mean_n << endl;

   for (int j=Nfirst; j<Nfirst+Nstraws[i]; j++) {   
     // if (j==Nfirst) cout << "Start with straw " << j+1 << "End with straw " << Nfirst+Nstraws[i] << endl;

      ahisto = anhisto->ProjectionY(Form("a[%i]",i), j+2, j+2);  //straw 1 is in bin 2

      double a_n = ahisto->GetEntries();

      int eff = 1;

      if (a_n < 0.25 * mean_n) eff = 0;

      if (!eff) printf("n=%i counts %.0f mean for ring is %i\n",j+1,a_n,mean_n);

      fprintf(outfile,"%i\n",eff);

    }

    Nfirst += Nstraws[i];

  }

  fclose(outfile);

}

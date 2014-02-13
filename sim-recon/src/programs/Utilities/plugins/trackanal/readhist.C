{

  char fnam[128] = "histograms_trackanal.root";
  TFile f1(fnam);

  TList * hist_list = (TList*)f1.Get("hist_list");
  TList * hist2d_list = (TList*)f1.Get("hist2d_list");

  Int_t Idx = hist_list->GetSize();
  TH1F * hist[Idx+1];
  cout<<"1D histograms: "<<Idx<<endl;
  for (Int_t k=0;k<Idx+1;k++)
    hist[k] = (TH1F *)hist_list->At(k);

  Idx = hist2d_list->GetSize();
  TH1F * hist2d[Idx+1];
  cout<<"2D histograms: "<<Idx<<endl;
  for (Int_t k=0;k<Idx+1;k++)
    hist2d[k] = (TH1F *)hist2d_list->At(k);

}

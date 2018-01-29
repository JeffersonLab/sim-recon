struct devilTreePT_t {
  //Int_t nGen;
  Double_t eGamma;
  Double_t weight;
  Double_t recoilE;
  Double_t recoilPx;
  Double_t recoilPy;
  Double_t recoilPz;
  Double_t electronE;
  Double_t electronPx;
  Double_t electronPy;
  Double_t electronPz;
  Double_t positronE;
  Double_t positronPx;
  Double_t positronPy;
  Double_t positronPz;
};

void setBranchesT1(TTree *t1, devilTreePT_t *devilTree){
  t1->Branch("devilTree.eGamma",&devilTree->eGamma,"eGamma/D");
  t1->Branch("devilTree.weight",&devilTree->weight,"weight/D");

  t1->Branch("devilTree.recoilE",&devilTree->recoilE,"recoilE/D");
  t1->Branch("devilTree.recoilPx",&devilTree->recoilPx,"recoilPx/D");
  t1->Branch("devilTree.recoilPy",&devilTree->recoilPy,"recoilPy/D");
  t1->Branch("devilTree.recoilPz",&devilTree->recoilPz,"recoilPz/D");

  t1->Branch("devilTree.electronE",&devilTree->electronE,"electronE/D");
  t1->Branch("devilTree.electronPx",&devilTree->electronPx,"electronPx/D");
  t1->Branch("devilTree.electronPy",&devilTree->electronPy,"electronPy/D");
  t1->Branch("devilTree.electronPz",&devilTree->electronPz,"electronPz/D");

  t1->Branch("devilTree.positronE",&devilTree->positronE,"positronE/D");
  t1->Branch("devilTree.positronPx",&devilTree->positronPx,"positronPx/D");
  t1->Branch("devilTree.positronPy",&devilTree->positronPy,"positronPy/D");
  t1->Branch("devilTree.positronPz",&devilTree->positronPz,"positronPz/D");
}

void getBranchesT1(TTree *t1,devilTreePT_t *devilTree){
  t1->SetBranchAddress("devilTree.eGamma",&devilTree->eGamma);
  t1->SetBranchAddress("devilTree.vsWeight",&devilTree->weight);

  t1->SetBranchAddress("devilTree.recoilE",&devilTree->recoilE);
  t1->SetBranchAddress("devilTree.recoilPx",&devilTree->recoilPx);
  t1->SetBranchAddress("devilTree.recoilPy",&devilTree->recoilPy);
  t1->SetBranchAddress("devilTree.recoilPz",&devilTree->recoilPz);

  t1->SetBranchAddress("devilTree.electronE",&devilTree->electronE);
  t1->SetBranchAddress("devilTree.electronPx",&devilTree->electronPx);
  t1->SetBranchAddress("devilTree.electronPy",&devilTree->electronPy);
  t1->SetBranchAddress("devilTree.electronPz",&devilTree->electronPz);

  t1->SetBranchAddress("devilTree.positronE",&devilTree->positronE);
  t1->SetBranchAddress("devilTree.positronPx",&devilTree->positronPx);
  t1->SetBranchAddress("devilTree.positronPy",&devilTree->positronPy);
  t1->SetBranchAddress("devilTree.positronPz",&devilTree->positronPz);
}



#include <string>
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Riostream.h"
#include "TH1.h"
#include "TH2.h" 
#include "TRandom3.h"
#include "qDevilLib.h"
#include "devilTreePT.h"
#include "HddmOut.h"

//STRUCTURE TO KEEP THE CONFIGURATION SETTINGS
struct genSettings_t {
  int beamType;         //Type of beam (1 -> single energy; 2-> bremstrahlung spectrum)
  int polDir;           //beam polarization direction (0-> unpolarized; 1-> pol in x; 2->pol in y)
  double eGammaInit;    //incident photon energy
  double eLower;        //incident photon energy min (only used for beamType = 1)
  double eUpper;        //incident photon energy max (only used for beamType = 1)
  int corrYes;          //do screening and radiative corrections if corrYes = 1
  int nToGen;           //number of events to generate
  int prescale;         //number of events between printing to terminal  
  int rSeed;            //seed for random number generator
  int tOut;             //type of output file (1->ROOT; 2->HDDM)
  int reaction;         //reaction (2->pair; 3->triplet)
  char outFile[80];     //name of output file
  char inFileBrem[80];  //name of input root file that contains spectra histogram cobrem_vs_E
};

//FUNCTION PROTOTYPES
double ampSqPT(int type, int polDir, TLorentzVector target, TLorentzVector beam,
                 TLorentzVector recoil,TLorentzVector q1,TLorentzVector q2);
void printUsage(genSettings_t genSettings, int goYes);

int main(int argc, char **argv){

  char *argptr;
  //SET THE DEFAULT CONFIGURATION SETTINGS
  genSettings_t genSettings;
  genSettings.reaction     = 2;
  genSettings.tOut         = 1;
  genSettings.polDir       = 2;
  genSettings.beamType     = 1;
  genSettings.eGammaInit   = 9.0;
  genSettings.eLower       = 8.0;
  genSettings.eUpper       = 9.0;
  genSettings.corrYes      = 1;
  genSettings.nToGen       = 100000;
  genSettings.prescale     = 1000;
  genSettings.rSeed        = 103;
  sprintf(genSettings.outFile,"genOut.root");
  sprintf(genSettings.inFileBrem,"cobrems.root");

  char rootFile[80];
  sprintf(rootFile,"genOut.root");

  char hddmFile[80];
  sprintf(hddmFile,"genOut.hddm");

  int outFileSet = 0;
  //COMMAND LINE PARSING
  for (int i=1; i<argc; i++) {
    argptr = argv[i];
    if (*argptr == '-') {
      argptr++;
      switch (*argptr) {
      case 'r':
        genSettings.rSeed = atoi(++argptr);
        break;
      case 'R':
        genSettings.reaction = atoi(++argptr);
        break;
      case 'n':
        genSettings.nToGen = atoi(++argptr);
        break;
      case 'p':
        genSettings.prescale = atoi(++argptr);
        break;
      case 'P':
        genSettings.polDir = atoi(++argptr);
        break;
      case 'b':
        genSettings.beamType = atoi(++argptr);
        break;
      case 't':
        genSettings.tOut = atoi(++argptr);
        break;
      case 'e':
        genSettings.eGammaInit = atof(++argptr);
        break;
      case 'l':
        genSettings.eLower = atof(++argptr);
        break;
      case 'u':
        genSettings.eUpper = atof(++argptr);
        break;
      case 'o':
	outFileSet = 1;
        strcpy(genSettings.outFile,++argptr);
        break;
      case 's':
        strcpy(genSettings.inFileBrem,++argptr);
        break;
	case 'h':
	// Do nothing right now
        break;
      default:
        fprintf(stderr,"Unrecognized argument: [-%s]\n\n",argptr);
        printUsage(genSettings,0);
        break;
      }
    }
  }


  if (genSettings.tOut == 1){
    if (outFileSet == 1) sprintf(rootFile,genSettings.outFile);
    sprintf(hddmFile,"/dev/null");
  }
  if (genSettings.tOut == 2){
    if (outFileSet == 1) sprintf(hddmFile,genSettings.outFile);
    sprintf(rootFile,"/dev/null");
  }

    
  for (int i=1; i<argc; i++) {
    argptr = argv[i];
    if (*argptr == '-') {
      argptr++;
      switch (*argptr) {
      case 'h':
        printUsage(genSettings,0);
        break;
      }
    }
  }

  double eGamma = genSettings.eGammaInit;

  //GET THE HISTOGRAM FOR COHERENT BREMSTRAHLUNG SPECTRUM
  TH1D* hGvsE;
  TFile *inCoBrem=new TFile(genSettings.inFileBrem); //Using spectrum from Richard
  hGvsE=(TH1D*)inCoBrem->Get("cobrem_vs_E");
  TH1D* hGvsEout = (TH1D*)hGvsE->Clone("hGvsEout");
  hGvsEout->Reset();
  hGvsEout->Rebin(30);
  int eBinLow = hGvsE->GetXaxis()->FindBin(genSettings.eLower);
  int eBinHigh = hGvsE->GetXaxis()->FindBin(genSettings.eUpper);
  hGvsE->GetXaxis()->SetRange(eBinLow,eBinHigh);
  double gMax = hGvsE->GetMaximum();

  //GET THE TRIPLET TO PAIR FRACTION HISTOGRAM NEEDED FOR RADIATIVE CORRECTIONS
  TH1D* hcsFraction;
  TFile *inCSfrac=new TFile("csFraction.root"); 
  hcsFraction=(TH1D*)inCSfrac->Get("hcsFraction");

  //DEFINE OUTPUT FILE

  HddmOut hddmGo(hddmFile);
  int evtNumber = 0;
  
  TFile *fout = new TFile(rootFile,"RECREATE");  

  // DEFINE TREE TO STORE THE DATA (SEE qDevilLib.h)
  TTree *t1 = new TTree("t1","genDevilPairs");
  devilTreePT_t devilTree; 
  setBranchesT1(t1, &devilTree); 

  //PRINT OUT THE SETTINGS
  printUsage(genSettings,1);

  double PIval =2*atan2(1,0);
  double alphaQED = 1.0/137.036;
  double hbarcSqr = 389.37966;
  double crossSection = 0.0;
  double fullWeight;
  double mElectron = 0.51099907e-3;
  double mProton = 0.938;
  //SETTING THE CUT PARAMETERS IS A 
  //BALANCING ACT BETWEEN SPEED OF CONVERGENCE
  //AND MAKING SURE THAT THE PHASE SPACE IS COMPLETE
  double Mcut= 5.0e-3; 
  if (genSettings.reaction == 3) Mcut= 20.0e-3; 
  if (genSettings.reaction == 2) Mcut= 1.0; 
  double qRcut = 1.0e-3; 
  if (genSettings.reaction == 2) qRcut= 2.0;

  //DEFINE SOME FOUR VECTORS
  TLorentzVector beam(0.0,0.0,eGamma,eGamma);
  TLorentzVector target;
  if (genSettings.reaction == 2) target.SetPxPyPzE(0.0,0.0,0.0,mProton);
  if (genSettings.reaction == 3) target.SetPxPyPzE(0.0,0.0,0.0,mElectron);
  TLorentzVector wVec = beam + target;
  TLorentzVector q1; //electron from pair
  TLorentzVector q2; //positron from pair
  TLorentzVector recoil;
  TLorentzVector q12;
  TLorentzVector q23;
  TLorentzVector moTransfer;

  if (genSettings.reaction == 2) recoil.SetPxPyPzE(0.0,0.0,0.0,mProton);
  if (genSettings.reaction == 3) recoil.SetPxPyPzE(0.0,0.0,0.0,mElectron);

  //DEFINE SOME HISTGRAMS
  double wMax = sqrt(target.Mag2() + 2*eGamma*target.Mag());
  if (genSettings.beamType == 2) {
    wMax = sqrt(target.Mag2() + 2*genSettings.eUpper*target.Mag());
  }
  double m12MinSq = pow(q1.Mag() + q2.Mag(),2);
  double m23MinSq = pow(q2.Mag() + recoil.Mag(),2);
  double m12MaxSq = pow(wMax - recoil.Mag(),2);
  double m23MaxSq = pow(wMax - q1.Mag(),2);
  int nBinXY = sqrt(genSettings.nToGen/100);
  TH2D* hPhaseSpaceR = new TH2D("hPhaseSpaceR","",nBinXY,m12MinSq,m12MaxSq,nBinXY,m23MinSq,m23MaxSq);


  TRandom3 *random = new TRandom3;
  random->SetSeed(genSettings.rSeed);
  double sigmaVal;
  PhasePT phaseGenR; //phase space object (based off Richards calculations)
  phaseGenR.SetRCut(qRcut);
  phaseGenR.SetM12Cut(Mcut);
  phaseGenR.SetBeam(beam);
  phaseGenR.SetTarget(target);
  double phaseSpaceWeight = 0.0;
  double sum=0.0;
  double sum2=0.0;

  double sumTest1=0.0;
  double sumTest2=0.0;
  double sumTest3=0.0;
  double sumTest4=0.0;
  double sumTest5=0.0;

  int genVal;
  int nTest = 0;
  double yMax,yVal,testValY;
  int eBin;
  int nGen = 0;
  int nSkip = 0;
  //nSkip = 79992;//ASDF
  for (int nGenTmp = 1; nGenTmp <= genSettings.nToGen; nGenTmp++) {
    genVal = -1;
    //MONTE CARLO THE COHERENT BREMSTRAHLUNG SPECTRUM
    if (genSettings.beamType == 2){ 
      yMax = gMax*1.02;
      yVal = 0.0;
      testValY = yMax + 10.0;
      while(testValY>yVal){//Monte Carlo the event to get brem spectrum
	eGamma = random->Uniform(genSettings.eLower,genSettings.eUpper); //Grab a photon energy
	testValY = random->Uniform(0.0,yMax); //Grab a test value                  
	eBin = hGvsE->GetXaxis()->FindBin(eGamma);
	yVal = hGvsE->GetBinContent(eBin);
      }
      hGvsEout->Fill(eGamma);
      beam.SetPxPyPzE(0,0,eGamma,eGamma);
      wVec = beam + target;
      phaseGenR.SetBeam(beam);
    }
    while (genVal < 0) {
      if (nGenTmp >= nSkip) nTest++;
      //GENERATE THE PHASE SPACE EVENT
      genVal = phaseGenR.Gen(random);
      if (genVal == -1) sumTest1 += 1.0;
      if (genVal == -2) sumTest2 += 1.0;
      if (genVal == -3) sumTest3 += 1.0;
      if (genVal == -4) sumTest4 += 1.0;
      if (genVal == -5) sumTest5 += 1.0;
      //NOTE: if (genVal < 0) then phase space not physical
    }
    if (nGenTmp < nSkip) continue;
    nGen++;

    //GET PHASE SPACE WEIGHT
    phaseSpaceWeight = phaseGenR.GetWeight();

    //GET THE FOUR-VECTORS
    beam       = phaseGenR.GetBeam();
    target     = phaseGenR.GetTarget();
    q1         = phaseGenR.GetQ1(); 
    q2         = phaseGenR.GetQ2();
    q12        = phaseGenR.GetQ12();
    q23        = phaseGenR.GetQ23();
    recoil     = phaseGenR.GetRecoil();
    moTransfer = recoil - target;

    //FILL PHASE SPACE HISTOGRAM
    hPhaseSpaceR->Fill(q12.Mag2(),q23.Mag2(),phaseSpaceWeight);
    
    //NOTE: Richard's calculation for the phase space gives small recoil momentum
    //very high probability. To account for this, the phase space weight is large
    //when the recoil is large. The problem is that these rare events have a large
    //weight. The rare events should not be a problem, except that one needs to 
    //generate a huge number of events to get the proper value. Instead, in the case
    //of triplet production, we can notice that every diagram has a corresponding
    //diagram that has the recoil electon momentum switched with the momentum of the 
    //produced electron. This means that the cross section for the phase space where 
    // q1.P > recoil.P will give identical results to the cross section for the
    //phase space where q1.P < recoil.P. Because of this symmetry between
    //q1.P and recoil.P, we can calculate the cross sections for the case
    //q1.P > recoil.P and just make sure that we account for the "lost" phase space.
    //Since we already have to account for the unphysical parts of the phase space
    //generation, accounting for the case where a1.P < recoil.P does not require
    //any additional work at this point in the code. For pair production, luckily
    //the rare events (high momentum protons) are so rare that I have not seen
    //any instability in the cross section results due to the rare large-weight 
    //events.
    if (genSettings.reaction == 3 && recoil.P() > q1.P()) continue;
    //If you comment out the line above and generate triplets, you should
    //notice that everything is fine with the cross sections except that
    //at some point (usualy after a very large number of events) there is
    //a spike in the cross section. 
    
    //ADDITIONAL FACTORS TO GET CROSS SECTION (micro barns)
    double fluxFactor = 4*beam.E()*(target.P() + target.E());
    double rhoFactor = 1.0/(8*recoil.E()*q12.P());
    double piFactor = pow(2*PIval,4-9)*pow(4*PIval,3);

    //SCREENING AND RADIATIVE CORRECTIONS
    double bohrRadius = (1.0/alphaQED)*(1.0/mElectron);
    double fH = 1.0/pow(1 + pow(bohrRadius*moTransfer.P()/2.0,2),2);

    double sHfactor = 1.0;
    double sHfactorPair = pow(1.0 - fH,2); //Screening for pair production
    double sHfactorTrip = 1.0 - pow(fH,2); //Screening for triplet production

    //RADIATIVE CORRECTIONS FOR PAIRS IS SIMPLY A COMMON FACTOR :)
    double radFactor = 1.0;
    double radFactorDelta = 0.0093;
    double radFactorPair =  1.0+radFactorDelta;  
    //THE RADIATIVE CORRECTION FOR TRIPLETS IS THE SAME MAG AS FOR PAIRS 
    //THIS MEANS WE HAVE TO GET THE FRACTION OF TRIPLETS TO PAIRS AND
    //SCALE THE radFactorDelta USING THE TRIPLET TO PAIR FRACTION
    int    eFracBin = hcsFraction->GetXaxis()->FindBin(beam.E());
    double csFrac = hcsFraction->GetBinContent(eFracBin);
    double radFactorTrip = 1.0 + radFactorDelta/csFrac;
    
    //GET THE CROSS SECTION
    if (genSettings.reaction == 2) { //PAIR PRODUCTION
      radFactor = radFactorPair;
      sHfactor = sHfactorPair;
      crossSection = ampSqPT(2,genSettings.polDir,target,beam,recoil,q1,q2)*hbarcSqr*pow(alphaQED,3)
	/ fluxFactor * rhoFactor * piFactor;
    }
    if (genSettings.reaction == 3) { //TRIPLET PRODUCTION
      radFactor = radFactorTrip;
      sHfactor = sHfactorTrip;
      crossSection = ampSqPT(3,genSettings.polDir,target,beam,recoil,q1,q2)*hbarcSqr*pow(alphaQED,3)
	/ fluxFactor * rhoFactor * piFactor;
    }

    if (genSettings.corrYes == 1) {
      crossSection *= sHfactor;
      crossSection *= radFactor;  
    }

    //GET INITIAL fullWeight
    fullWeight = crossSection*phaseSpaceWeight;

    //CHECK FOR SOMETHING BAD
    if (fullWeight >= 0 || fullWeight <=0){
      //do nothing
    } else {
      cout<<"!!!!Something Bad!!!!"<<endl;
    }
    //CORRECT fullWeight FOR UNPHYSICAL PHASE SPACE EVENTS THAT ARE REMOVED
    Double_t phaseSpaceCorrection = nGen*1.0/(1.0*nTest);
    fullWeight *= phaseSpaceCorrection;
    
    //FILL THE TREE
    devilTree.eGamma   = eGamma;
    devilTree.weight   = fullWeight;
    devilTree.recoilE  = recoil.E();
    devilTree.recoilPx = recoil.Px();
    devilTree.recoilPy = recoil.Py();
    devilTree.recoilPz = recoil.Pz();
    devilTree.electronE  = q1.E();
    devilTree.electronPx = q1.Px();
    devilTree.electronPy = q1.Py();
    devilTree.electronPz = q1.Pz();
    devilTree.positronE  = q2.E();
    devilTree.positronPx = q2.Px();
    devilTree.positronPy = q2.Py();
    devilTree.positronPz = q2.Pz();
    t1->Fill(); 

    //CALCULATING CROSS SECTION ERROR (SAME WAY THAT RICHARD DOES) AND PRINT RESULT
    sum += fullWeight;
    sum2 += pow(fullWeight,2);

    if (nGen/genSettings.prescale*genSettings.prescale == nGen) {
      cout <<"Integrated cross section after " << nGen << " events : "
           << sum/nGen << " +/- " <<sqrt(sum2-pow(sum,2)/nGen)/nGen<< sigmaVal <<" micro barns"<< endl;
    }
    //HDDM STUFF
    tmpEvt_t tmpEvt;
    tmpEvt.beam = beam;
    tmpEvt.target = target;
    tmpEvt.q1 = q1;
    tmpEvt.q2 = q2;
    tmpEvt.recoil = recoil;
    tmpEvt.nGen = 3;
    tmpEvt.rxn = genSettings.reaction;
    tmpEvt.weight = fullWeight;
    evtNumber++;
    if (genSettings.tOut == 2) hddmGo.write(tmpEvt,evtNumber);
  }
  //WRITE THE TREE AND HISTOGRAMS TO FILE
  if (genSettings.tOut == 1){
    fout->cd();
    t1->Write();       
    hPhaseSpaceR->Write();
    hGvsEout->Write(); 
  }
  fout->Close();     
  cout<<"All done. Bye"<<endl;
}

void printUsage(genSettings_t genSettings, int goYes){
  if (goYes == 0){
    fprintf(stderr,"\nSWITCHES:\n");
    fprintf(stderr,"-h\tPrint this message\n");
    fprintf(stderr,"-R<arg>\tReaction:\n");
    fprintf(stderr,"\t\t-R2 = Pair production off of proton\n");
    fprintf(stderr,"\t\t-R3 = Triplet production\n");
    fprintf(stderr,"-n<arg>\tNumber of events to generate\n");
    fprintf(stderr,"-r<arg>\tUser defined random number seed\n");

    fprintf(stderr,"-t<arg>\tType of output file\n");
    fprintf(stderr,"\t\t-t1 = ROOT file\n");
    fprintf(stderr,"\t\t-t2 = HDDM file\n");

    fprintf(stderr,"-p<arg>\tPrescale factor (number of events between printing to terminal)\n");
    fprintf(stderr,"-P<arg>\tPhoton beam polarization direction:\n");
    fprintf(stderr,"\t\t-P0 = Unpolarized\n");
    fprintf(stderr,"\t\t-P2 = Polarized in x-direction (100 percent)\n");
    fprintf(stderr,"\t\t-P3 = Polarized in y-direction (100 percent)\n");
    fprintf(stderr,"-b<arg>\tBeam type:\n");
    fprintf(stderr,"\t\t-b1 = Single photon energy\n");
    fprintf(stderr,"\t\t-b2 = Bremstrahlung spectra\n");
    fprintf(stderr,"-e<arg>\tPhoton energy in GeV.                  ONLY USED IF -b1\n");
    fprintf(stderr,"-l<arg>\tMinimum incident photon energy in GeV. ONLY USED IF -b2\n");
    fprintf(stderr,"-u<arg>\tMaximum incident photon energy in GeV. ONLY USED IF -b2\n");
    fprintf(stderr,"-s<arg>\tfile with histogram cobrem_vs_E.       ONLY USED IF -b2\n");
    fprintf(stderr,"-o<arg>\tOutFile name\n");

    cout<<""<<endl;
    cout<<"The above switches overide the default setting."<<endl;
    cout<<""<<endl;
    
    if (genSettings.beamType == 1){
      cout<<"The current operation is equivalent to the command:"<<endl;
      cout<<"genDevilPT"
	  <<" -n"<<genSettings.nToGen
	  <<" -R"<<genSettings.reaction
#ifdef WITHHDDM
	  <<" -t"<<genSettings.tOut
#endif
	  <<" -r"<<genSettings.rSeed
	  <<" -p"<<genSettings.prescale
	  <<" -b"<<genSettings.beamType
	  <<" -P"<<genSettings.polDir
	  <<" -e"<<genSettings.eGammaInit
	//<<" -l"<<genSettings.eLower
	//<<" -u"<<genSettings.eUpper
	  <<" -o"<<genSettings.outFile
	//<<" -s"<<genSettings.inFileBrem
	  <<"\n"<<endl;
    } else {
      cout<<"The current operation is equivalent to the command:"<<endl;
      cout<<"genDevilPT"
	  <<" -n"<<genSettings.nToGen
	  <<" -R"<<genSettings.reaction
	  <<" -t"<<genSettings.tOut
	  <<" -r"<<genSettings.rSeed
	  <<" -p"<<genSettings.prescale
	  <<" -b"<<genSettings.beamType
	  <<" -P"<<genSettings.polDir
	//<<" -e"<<genSettings.eGammaInit
	  <<" -l"<<genSettings.eLower
	  <<" -u"<<genSettings.eUpper
	  <<" -o"<<genSettings.outFile
	  <<" -s"<<genSettings.inFileBrem
	  <<"\n"<<endl;
    }
    cout<<""<<endl;
    cout<<"NOTE 1-> If providing custom photon spectrum histogram:"<<endl; 
    cout<<"         * The file containing the histogram is provided"<<endl; 
    cout<<"           as an argument to the -s switch"<<endl;
    cout<<"         * The histogram MUST be of type TH1D"<<endl;
    cout<<"         * The histogram MUST have energy units of GeV for x-axis"<<endl;
    cout<<"         * The histogram MUST be named cobrem_vs_E"<<endl;
    cout<<"         * The energy bin widths in the histogram are not important"<<endl;
    cout<<""<<endl;
//
//    cout<<"NOTE 2-> To be able to output in HDDM mode you must compile like:"<<endl; 
//    cout<<"         make HDDM=yes"<<endl;
//    cout<<""<<endl;

  }

  if (goYes == 1){
    if(genSettings.reaction == 2) cout<<"\nReaction: Pair production off of proton"<<endl;
    if(genSettings.reaction == 3) cout<<"\nReaction: Triplet production"<<endl;
    cout<<"Generating "<<genSettings.nToGen<<" events"
    <<" with output file "<<genSettings.outFile<<endl;
    if(genSettings.beamType == 1){
      cout<<"Beam type: Single photon energy = "<<genSettings.eGammaInit<<" GeV"<<endl;
    }
    if(genSettings.beamType == 2){
      cout<<"Beam type: Bremstrahlung spectrum from histogram "<<"cobrem_vs_E"<<endl;
      cout<<"           within file "<<genSettings.inFileBrem<<endl;
      cout<<"           Using energy range from "
	  <<genSettings.eLower<<" GeV to "<<genSettings.eUpper<<" GeV"<<endl;
    }
    cout<<"Random number seed = "<<genSettings.rSeed<<endl;
    cout<<""<<endl;
  }

  if (goYes == 0) exit(0);
}

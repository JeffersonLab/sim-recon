#include<TGClient.h>
#include<TCanvas.h>
#include<TF1.h>
#include<TRandom.h>
#include<TGButton.h>
#include<TGFrame.h>
#include<TRootEmbeddedCanvas.h>
#include<RQ_OBJECT.h>

class MyMainFrame{

  RQ_OBJECT("MyMainFrame")

  private: 
    TGMainFrame     *fMain;
    TRootEmbeddedCanvas *fEcanvas;
    TView *fView;
    Int_t fThetaView = -60;
    Int_t fPhiView = 0;
    Int_t fPsiView = 0;

    int colorPad = 41;
    int colorBCAL = 20;
    int colorCDC = 4;

  public:
    MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~MyMainFrame();
    void DrawSITE();
    void DrawCDC();
    void DrawBCAL();
    void ThetaPlus();
    void PhiPlus();
    void ThetaMinus();
    void PhiMinus();
    void SideView();
};

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h)
{
  //Create a main frame
  fMain = new TGMainFrame(p, w, h);

  TGHorizontalFrame *mFrame = new TGHorizontalFrame(fMain,600,600);

  TGVerticalFrame *vframe1 = new TGVerticalFrame(mFrame, 40, 40);
  TGVerticalFrame *vframe2 = new TGVerticalFrame(mFrame, 40, 40);


  TGTextButton *exit2 = new TGTextButton(vframe1, "&Exit2", "gApplication->Terminate(0)");
  vframe1->AddFrame(exit2, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  //Create canvas widget
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas", mFrame, 400, 400);

  //Create a horizontal frame widget with buttons

  TGTextButton *drawSITE = new TGTextButton(vframe2, "&SITE");
  drawSITE->Connect("Clicked()", "MyMainFrame", this, "DrawSITE()");
  vframe2->AddFrame(drawSITE, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawCDC = new TGTextButton(vframe2, "&CDC");
  drawCDC->Connect("Clicked()", "MyMainFrame", this, "DrawCDC()");
  vframe2->AddFrame(drawCDC, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawBCAL = new TGTextButton(vframe2, "&BCAL");
  drawBCAL->Connect("Clicked()", "MyMainFrame", this, "DrawBCAL()");
  vframe2->AddFrame(drawBCAL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *phiM = new TGTextButton(vframe2, "&phi-");
  phiM->Connect("Clicked()", "MyMainFrame", this, "PhiMinus()");
  vframe2->AddFrame(phiM, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *phiP = new TGTextButton(vframe2, "&phi+");
  phiP->Connect("Clicked()", "MyMainFrame", this, "PhiPlus()");
  vframe2->AddFrame(phiP, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *thetaM = new TGTextButton(vframe2, "&theta-");
  thetaM->Connect("Clicked()", "MyMainFrame", this, "ThetaMinus()");
  vframe2->AddFrame(thetaM, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *thetaP = new TGTextButton(vframe2, "&theta+");
  thetaP->Connect("Clicked()", "MyMainFrame", this, "ThetaPlus()");
  vframe2->AddFrame(thetaP, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *sideV = new TGTextButton(vframe2, "&Side View");
  sideV->Connect("Clicked()", "MyMainFrame", this, "SideView()");
  vframe2->AddFrame(sideV, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *exit = new TGTextButton(vframe2, "&Exit", "gApplication->Terminate(0)");
  vframe2->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  mFrame->AddFrame(vframe1, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
  mFrame->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 1));
  mFrame->AddFrame(vframe2, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));

  fMain->AddFrame(mFrame, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 1));

  cerr << "Here.....0" << endl;

  // Set a name to the main frame
  fMain->SetWindowName("Simple Example");

  // Map all subwindows of main frame
  fMain->MapSubwindows();


  // Initialize the layout algorithm
  fMain->Resize(fMain->GetDefaultSize());

  // Map main frame
  fMain->MapWindow();

  //gROOT->LoadMacro("hdgeant.C");
  gROOT->LoadMacro("hddsroot.C");
  float xmin[3] = {-100, -100, -100};
  float xmax[3] = {100, 100, 100};
  fView = new TView(xmin, xmax, 1);

  //hdgeant();
  hddsroot();

  this->DrawSITE();


}

void MyMainFrame::DrawSITE()
{
  // Draws function graphics in randomly choosen interval
  gGeoManager->SetVisLevel(3);
  gGeoManager->GetMasterVolume()->Draw();

  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    

  TCanvas *fCanvas = fEcanvas->GetCanvas();
  fCanvas->cd();
  fCanvas->Update();

}

void MyMainFrame::DrawBCAL()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *BCAL = gGeoManager->GetVolume("BCAL");
  gGeoManager->SetTopVolume(BCAL);
  gGeoManager->SetVisLevel(1);
  gGeoManager->GetTopVolume()->Draw();

  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    

  TCanvas *fCanvas = fEcanvas->GetCanvas();
  fCanvas->cd();
  fCanvas->Update();

}

void MyMainFrame::DrawCDC()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *CDC = gGeoManager->GetVolume("CDC");
  gGeoManager->SetTopVolume(CDC);
  gGeoManager->SetVisLevel(1);
  gGeoManager->GetTopVolume()->Draw();

  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    

  TCanvas *fCanvas = fEcanvas->GetCanvas();
  fCanvas->cd();
  fCanvas->Update();

}

void MyMainFrame::PhiPlus()
{
  fPhiView+=10;
  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

void MyMainFrame::ThetaPlus()
{
  fThetaView+=10;
  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

void MyMainFrame::PhiMinus()
{
  fPhiView-=10;
  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

void MyMainFrame::ThetaMinus()
{
  fThetaView-=10;
  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

void MyMainFrame::SideView()
{
    fView->Side();
    this->fPhiView = 0;
    this->fThetaView = 90;

  //this->fPhiView = 0;
  //this->fThetaView = 90;
  //this->fPsiView = 90;
  //fView->RotateView(this->fPhiView,this->fThetaView);
  //fView->SetPsi(this->fPsiView);      
  //fView->Zoom();    
}

MyMainFrame::~MyMainFrame()
{
  // Clean up used widgets: frames, buttons, layouthints
  fMain->Cleanup();
  delete fMain;
}

void hdexample()
{
  // Popup the GUI.....
  new MyMainFrame(gClient->GetRoot(), 800, 800);
}

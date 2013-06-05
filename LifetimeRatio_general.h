#include <iostream>
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TExec.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>

using namespace std;

TColor* ApplyZAxisColor(TColor *pal) {
  //palette settings - completely independent to histograms
  TColor::InitializeColors();
  const Int_t nRGBs = 3;
  const Int_t NCont = 20;

  Double_t stops[nRGBs] = { 0.00, 0.50, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.50, 1.00 };
  Double_t green[nRGBs] = { 0.00, 0.50, 1.00 };
  Double_t blue[nRGBs]  = { 0.50, 0.50, 0.00 };
  pal = new TColor();
  pal->CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
  return pal;
}

void ApplyHistStyle(TH1* h, int i, int j) {
  int colorArr[] = {kBlack, kRed+1, kOrange+7, kGreen+2, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
  int markerFullArr[6] = {kFullCircle, kFullSquare, kFullStar, kFullTriangleUp, kFullTriangleDown, 33};
//  int markerOpenArr[6] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenTriangleUp, kOpenTriangleDown, kOpenDiamond};
  
  h->SetLineColor(colorArr[j]);
  h->SetLineWidth(1.800);
  h->SetMarkerColor(colorArr[j]);
  h->SetMarkerStyle(markerFullArr[i]);
  h->SetMarkerSize(1.200);
  if (i < 2) h->SetMarkerSize(1.000);
  if (i == 2) h->SetMarkerSize(1.700);
  if (i == 5) h->SetMarkerSize(1.500);
  
  h->SetTitleSize(0.048,"XYZ");
  h->SetLabelSize(0.048,"XYZ");
}

void ApplyHistStyle(TH2* h) {
  h->SetTitleSize(0.048,"XYZ");
  h->SetLabelSize(0.048,"XYZ");
}

void ApplyLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.2);
}

void Draw2DHistsComm1(TH1 *h, TExec *pal) {
  gPad->SetRightMargin(0.125);
  h->Draw("col");
  pal->Draw();
  h->Draw("colz same");
}

void Draw2DHistsComm2(TH1 *h, TExec *pal, bool range=false) {
  gPad->SetRightMargin(0.150);
  if (range) h->GetZaxis()->SetRangeUser(0.5,1.5);
  h->Draw("col");
  pal->Draw();
  h->Draw("colz same");
}

void Draw1DHistsComm1(TH1 *h) {
  h->Reset("ICES");
  h->GetYaxis()->SetRangeUser(0.5,1.5);
  ApplyHistStyle(h,0,0);
}


class LifetimeCheck {
  private:
    TTree *myTree;
    
    vector<double> rapbins;
    vector<UInt_t> runbins;

    int nBins;
    double xmin, xmax;
    double xArray[20];
    int nBinsErr;
    double xminErr, xmaxErr;
    double xErrArray[20];
    int nBinsEta;
    double xminEta, xmaxEta;
    double xEtaArray[20];
    int nBinsPhi;
    double xminPhi, xmaxPhi;
    double xPhiArray[20];
    int nBinsPT;
    double xminPT, xmaxPT;
    double xPTArray[20];
    int nBinsY;
    double xminY, xmaxY;
    double xYArray[20];
    int nBinsDxy;
    double xminDxy, xmaxDxy;
    double xDxyArray[20];
    int nBinsDz;
    double xminDz, xmaxDz;
    double xDzArray[20];
    int nBinsChi2Trk;
    double xminChi2Trk, xmaxChi2Trk;
    double xChi2TrkArray[20];
    int nBinsChi2Glb;
    double xminChi2Glb, xmaxChi2Glb;
    double xChi2GlbArray[20];
    int nBinsNTrkHits;
    double xminNTrkHits, xmaxNTrkHits;
    double xNTrkHitsArray[20];
    int nBinsNPixWMea;
    double xminNPixWMea, xmaxNPixWMea;
    double xNPixWMeaArray[20];

    double mmin, mmax, ptmin, ptmax;
    static const int nHists=10;
    unsigned int nRuns;
    bool FillCowboy;

    // J/psi kinematic variables
    TH1D *ctauP[nHists], *ctauM[nHists], *ratio[nHists];
    TH1D *ctauErrP[nHists], *ctauErrM[nHists], *ratioErr[nHists];
    TH1D *pTP[nHists], *pTM[nHists], *ratioPT[nHists];
    TH1D *yPA, *yMA, *ratioYA;
    TH1D *phiP[nHists], *phiM[nHists], *ratioPhi[nHists];
    
    TH2D *ctauCtauErrP[nHists], *ctauCtauErrM[nHists], *ratioCtauCtauErr[nHists];
    TH2D *ctauPTP[nHists], *ctauPTM[nHists], *ratioCtauPT[nHists];
    TH2D *ctauYPA, *ctauYMA, *ratioCtauYA;
    TH2D *ctauPhiP[nHists], *ctauPhiM[nHists], *ratioCtauPhi[nHists];
    
    TH2D *yPhiP, *yPhiM, *ratioYPhi;

    // Single muon kinematic variables
    TH1D *etaP[nHists], *etaM[nHists], *ratioEta[nHists];
    TH1D *phiMuP[nHists], *phiMuM[nHists], *ratioPhiMu[nHists];
    TH1D *dxyP[nHists], *dxyM[nHists], *ratioDxy[nHists];
    TH1D *dzP[nHists], *dzM[nHists], *ratioDz[nHists];
    TH1D *chi2TrkP[nHists], *chi2TrkM[nHists], *ratioChi2Trk[nHists];
    TH1D *chi2GlbP[nHists], *chi2GlbM[nHists], *ratioChi2Glb[nHists];
    TH1D *nTrkHitsP[nHists], *nTrkHitsM[nHists], *ratioNTrkHits[nHists];
    TH1D *nPixWMeaP[nHists], *nPixWMeaM[nHists], *ratioNPixWMea[nHists];
    
    TH2D *ctauEtaP[nHists], *ctauEtaM[nHists], *ratioCtauEta[nHists];
    TH2D *ctauPhiMuP[nHists], *ctauPhiMuM[nHists], *ratioCtauPhiMu[nHists];
    TH2D *ctauDxyP[nHists], *ctauDxyM[nHists], *ratioCtauDxy[nHists];
    TH2D *ctauDzP[nHists], *ctauDzM[nHists], *ratioCtauDz[nHists];
    TH2D *ctauChi2TrkP[nHists], *ctauChi2TrkM[nHists], *ratioCtauChi2Trk[nHists];
    TH2D *ctauChi2GlbP[nHists], *ctauChi2GlbM[nHists], *ratioCtauChi2Glb[nHists];
    TH2D *ctauNTrkHitsP[nHists], *ctauNTrkHitsM[nHists], *ratioCtauNTrkHits[nHists];
    TH2D *ctauNPixWMeaP[nHists], *ctauNPixWMeaM[nHists], *ratioCtauNPixWMea[nHists];

    TH2D *etaPhiMuP[nHists], *etaPhiMuM[nHists], *ratioEtaPhiMu[nHists];
    
    // Run by run checking
    TH1D *runPhiP[nHists], *runPhiM[nHists], *runRatioPhi[nHists];
    TH1D *runPhiMuP[nHists], *runPhiMuM[nHists], *runRatioPhiMu[nHists];
    TH1D *runNTrkHitsP[nHists], *runNTrkHitsM[nHists], *runRatioNTrkHits[nHists];
  
    // Titles of histograms
    string titleRaw, titleRawErr, titleEta, titlePhi, titlePhiMu, titlePT, titleY, titleDxy, titleDz, titleChi2Trk, titleChi2Glb, titleNTrkHits, titleNPixWMea;
    string titleRatio, titleRatioErr, titleRatioEta, titleRatioPhi, titleRatioPhiMu, titleRatioPT, titleRatioY, titleRatioDxy, titleRatioDz, titleRatioChi2Trk, titleRatioChi2Glb, titleRatioNTrkHits, titleRatioNPixWMea;
    string title2DCtauErr, title2DEta, title2DPhi, title2DPhiMu, title2DPT, title2DY, title2DDxy, title2DDz, title2DChi2Trk, title2DChi2Glb, title2DYPhi, title2DEtaPhiMu, title2DNTrkHits, title2DNPixWMea;
    
  public:
    LifetimeCheck(TTree *tree, bool fillCowboy);
    LifetimeCheck(TFile *input, const char *name, bool fillCowboy);
    void Init(bool fillCowboy);
    void CreateHistos(const char *name);
    void ReadTree();
    void SetHistStyle(int a, int b);
    void WriteToFile(TFile *output);
    void Draw1DHists(TCanvas *canv[], bool drawAxis);
    void Draw2DHists(TCanvas *canv[]);
    void DrawRunCompHists(TCanvas *canv[]);
    void LegendOn2DHists(TCanvas *canv[], string strname[]);
    bool IsCowboy(int charge1, double phi1, int charge2, double phi2);
    void GetLimits(int t_nBins[], double t_xmin[], double t_xmax[], vector<double> *t_rapbins, vector<UInt_t> *t_runbins);

};

bool LifetimeCheck::IsCowboy(int charge1, double phi1, int charge2, double phi2) {
  double dPhi2Mu = phi1-phi2;
  while (dPhi2Mu > TMath::Pi()) dPhi2Mu -= 2*TMath::Pi();
  while (dPhi2Mu <= -1*TMath::Pi()) dPhi2Mu += 2*TMath::Pi();

  if (charge1*dPhi2Mu > 0) return true;
  return false;
}

void LifetimeCheck::Init(bool fillCowboy) {
  nBins=50, xmin=-1, xmax=2;
  nBinsErr=30;
  xminErr=0, xmaxErr=0.3;
  nBinsEta=12;
  xminEta=0.0, xmaxEta=2.4;
  nBinsPhi=15;
  xminPhi=-1*TMath::Pi(), xmaxPhi=TMath::Pi();
  nBinsPT=15;
  xminPT=0, xmaxPT=30;
  nBinsY=12;
  xminY=0, xmaxY=2.4;
  nBinsDxy=10;
  xminDxy=0, xmaxDxy=0.1;
  nBinsDz=40;
  xminDz=0, xmaxDz=0.4;
  nBinsChi2Trk=40;
  xminChi2Trk=0, xmaxChi2Trk=4;
  nBinsChi2Glb=200;
  xminChi2Glb=0, xmaxChi2Glb=20;
  nBinsNTrkHits=32;
  xminNTrkHits=0, xmaxNTrkHits=32;
  nBinsNPixWMea=6;
  xminNPixWMea=0, xmaxNPixWMea=6;

//  mmin = 2.95;
//  mmax = 3.25;
  mmin = 2.6;
  mmax = 3.5;

  ptmin=6.5;
  ptmax=30.0;

  FillCowboy = fillCowboy;
  
  const double xarray[]={xmin,-0.6,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.9,xmax};
  for (unsigned int i=0; i<sizeof(xarray)/sizeof(double); i++) {
    xArray[i] = xarray[i];
  }
  nBins = sizeof(xarray)/sizeof(double)-1;
  
  const double xerrarray[]={xminErr,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,xmaxErr};
  for (unsigned int i=0; i<sizeof(xerrarray)/sizeof(double); i++) {
    xErrArray[i] = xerrarray[i];
  }
  nBinsErr = sizeof(xerrarray)/sizeof(double)-1;

  const double xptarray[]={xminPT,2,4,5,6,7,8,9,10,11,12,14,16,20,xmaxPT};
  for (unsigned int i=0; i<sizeof(xptarray)/sizeof(double); i++) {
    xPTArray[i] = xptarray[i];
  }
  nBinsPT = sizeof(xptarray)/sizeof(double)-1;

  const double xdxyarray[]={xminDxy,0.005,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.07,xmaxDxy};
  for (unsigned int i=0; i<sizeof(xdxyarray)/sizeof(double); i++) {
    xDxyArray[i] = xdxyarray[i];
  }
  nBinsDxy = sizeof(xdxyarray)/sizeof(double)-1;

  const double xdzarray[]={xminDz,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.10,0.15,0.20,xmaxDz};
  for (unsigned int i=0; i<sizeof(xdzarray)/sizeof(double); i++) {
    xDzArray[i] = xdzarray[i];
  }
  nBinsDz = sizeof(xdzarray)/sizeof(double)-1;

  const double xchi2glbarray[]={xminChi2Glb,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10,xmaxChi2Glb};
  for (unsigned int i=0; i<sizeof(xchi2glbarray)/sizeof(double); i++) {
    xChi2GlbArray[i] = xchi2glbarray[i];
  }
  nBinsChi2Glb = sizeof(xchi2glbarray)/sizeof(double)-1;

  const double xchi2trkarray[]={xminChi2Trk,0.3,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,2.0,2.5,xmaxChi2Trk};
  for (unsigned int i=0; i<sizeof(xchi2trkarray)/sizeof(double); i++) {
    xChi2TrkArray[i] = xchi2trkarray[i];
  }
  nBinsChi2Trk = sizeof(xchi2trkarray)/sizeof(double)-1;

  const double xntrkhitsarray[]={xminNTrkHits,10,11,12,13,14,15,16,17,18,19,20,21,22,24,xmaxNTrkHits};
  for (unsigned int i=0; i<sizeof(xntrkhitsarray)/sizeof(double); i++) {
    xNTrkHitsArray[i] = xntrkhitsarray[i];
  }
  nBinsNTrkHits = sizeof(xntrkhitsarray)/sizeof(double)-1;

  const double xnpixwmeaarray[]={xminNPixWMea,1,2,3,4,5,xmaxNPixWMea};
  for (unsigned int i=0; i<sizeof(xnpixwmeaarray)/sizeof(double); i++) {
    xNPixWMeaArray[i] = xnpixwmeaarray[i];
  }
  nBinsNPixWMea = sizeof(xnpixwmeaarray)/sizeof(double)-1;

  titleRaw =";#font[12]{l}_{J/#psi} (mm);Counts";
  titleRawErr =";#font[12]{l}_{J/#psi} error (mm);Counts";
  titleEta =";#mu |#eta|;Counts";
  titlePhi =";J/#psi #phi (rad);Counts";
  titlePhiMu =";#mu #phi (rad);Counts";
  titlePT =";J/#psi p_{T} (GeV/c);Counts";
  titleY =";J/#psi |y|;Counts";
  titleDxy = ";#mu D_{xy};Counts";
  titleDz = ";#mu D_{z};Counts";
  titleChi2Trk = ";Tracker #mu #chi^{2}/ndof;Counts";
  titleChi2Glb = ";Global #mu #chi^{2}/ndof;Counts";
  titleNTrkHits = ";Number of tracker hits;Counts";
  titleNPixWMea = ";Number of pixel layers with measurement;Counts";

  titleRatio =";#font[12]{l}_{J/#psi} (mm);Counts (y>0)/(y<0)";
  titleRatioErr =";#font[12]{l}_{J/#psi} error (mm);Counts (y>0)/(y<0)";
  titleRatioEta =";#mu |#eta|;Counts (y>0)/(y<0)";
  titleRatioPhi =";J/#psi #phi (rad);Counts (y>0)/(y<0)";
  titleRatioPhiMu =";#mu #phi (rad);Counts (y>0)/(y<0)";
  titleRatioPT =";J/#psi p_{T} (GeV/c);Counts (y>0)/(y<0)";
  titleRatioY =";J/#psi |y|;Counts (y>0)/(y<0)";
  titleRatioDxy = ";#mu D_{xy};Counts (y>0)/(y<0)";
  titleRatioDz = ";#mu D_{z};Counts (y>0)/(y<0)";
  titleRatioChi2Trk = ";Tracker #mu #chi^{2}/ndof;Counts (y>0)/(y<0)";
  titleRatioChi2Glb = ";Global #mu #chi^{2}/ndof;Counts (y>0)/(y<0)";
  titleRatioNTrkHits = ";Number of tracker hits;Counts (y>0)/(y<0)";
  titleRatioNPixWMea = ";Number of pixel layers with measurement;Counts (y>0)/(y<0)";
  
  title2DCtauErr =";#font[12]{l}_{J/#psi} (mm);#font[12]{l}_{J/#psi} error (mm)";
  title2DEta =";#font[12]{l}_{J/#psi} (mm);#mu |#eta|";
  title2DPhi =";#font[12]{l}_{J/#psi} (mm);J/#psi #phi (rad)";
  title2DPhiMu =";#font[12]{l}_{J/#psi} (mm);#mu #phi (rad)";
  title2DPT =";#font[12]{l}_{J/#psi} (mm);J/#psi p_{T} (GeV/c)";
  title2DY =";#font[12]{l}_{J/#psi} (mm);J/#psi |y|";
  title2DDxy = ";#font[12]{l}_{J/#psi} (mm);#mu D_{xy}";
  title2DDz = ";#font[12]{l}_{J/#psi} (mm);#mu D_{z}";
  title2DChi2Trk = ";#font[12]{l}_{J/#psi} (mm);Tracker #mu #chi^{2}/ndof";
  title2DChi2Glb = ";#font[12]{l}_{J/#psi} (mm);Global #mu #chi^{2}/ndof";
  title2DNTrkHits = ";#font[12]{l}_{J/#psi} (mm);Number of tracker hits";
  title2DNPixWMea = ";#font[12]{l}_{J/#psi} (mm);Number of pixel layers with measurement";

  title2DYPhi = ";J/#psi |y|;J/#psi #phi (rad)";
  title2DEtaPhiMu = ";#mu |#eta|;#mu #phi (rad)";

  const double raparray[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  for (unsigned int i=0; i<sizeof(raparray)/sizeof(double); i++) {
    rapbins.push_back(raparray[i]);
  }

  const UInt_t runarray[] = {181611, 182133, 182324, 182664, 182960, 183127};
  for (unsigned int i=0; i<sizeof(runarray)/sizeof(UInt_t); i++) {
    runbins.push_back(runarray[i]);
  }
  nRuns = runbins.size()-1;
}

LifetimeCheck::LifetimeCheck(TTree *tree, bool fillCowboy) {
  myTree = tree;

  Init(fillCowboy);
}


LifetimeCheck::LifetimeCheck(TFile *input, const char *name, bool fillCowboy) {
  Init(fillCowboy);

  for (unsigned int i=0; i<rapbins.size(); i++) {
    ratio[i] = (TH1D*)input->Get(Form("%s_ctauRatio_%d",name,i));
    ratioErr[i] = (TH1D*)input->Get(Form("%s_ctauErrRatio_%d",name,i));
    ratioEta[i] = (TH1D*)input->Get(Form("%s_etaRatio_%d",name,i));
    ratioPhi[i] = (TH1D*)input->Get(Form("%s_phiRatio_%d",name,i));
    ratioPhiMu[i] = (TH1D*)input->Get(Form("%s_phiMuRatio_%d",name,i));
    ratioPT[i] = (TH1D*)input->Get(Form("%s_pTRatio_%d",name,i));
    ratioDxy[i] = (TH1D*)input->Get(Form("%s_dxyRatio_%d",name,i));
    ratioDz[i] = (TH1D*)input->Get(Form("%s_dzRatio_%d",name,i));
    ratioChi2Trk[i] = (TH1D*)input->Get(Form("%s_chi2TrkRatio_%d",name,i));
    ratioChi2Glb[i] = (TH1D*)input->Get(Form("%s_chi2GlbRatio_%d",name,i));
    ratioNTrkHits[i] = (TH1D*)input->Get(Form("%s_nTrkHitsRatio_%d",name,i));
    ratioNPixWMea[i] = (TH1D*)input->Get(Form("%s_nPixWMeaRatio_%d",name,i));
    
    ctauP[i] = (TH1D*)input->Get(Form("%s_ctauPlus_%d",name,i));
    ctauM[i] = (TH1D*)input->Get(Form("%s_ctauMinus_%d",name,i));
    ctauErrP[i] = (TH1D*)input->Get(Form("%s_ctauErrPlus_%d",name,i));
    ctauErrM[i] = (TH1D*)input->Get(Form("%s_ctauErrMinus_%d",name,i));
    etaP[i] = (TH1D*)input->Get(Form("%s_etaPlus_%d",name,i));
    etaM[i] = (TH1D*)input->Get(Form("%s_etaMinus_%d",name,i));
    phiP[i] = (TH1D*)input->Get(Form("%s_phiPlus_%d",name,i));
    phiM[i] = (TH1D*)input->Get(Form("%s_phiMinus_%d",name,i));
    phiMuP[i] = (TH1D*)input->Get(Form("%s_phiMuPlus_%d",name,i));
    phiMuM[i] = (TH1D*)input->Get(Form("%s_phiMuMinus_%d",name,i));
    pTP[i] = (TH1D*)input->Get(Form("%s_pTPlus_%d",name,i));
    pTM[i] = (TH1D*)input->Get(Form("%s_pTMinus_%d",name,i));
    dxyP[i] = (TH1D*)input->Get(Form("%s_dxyPlus_%d",name,i));
    dxyM[i] = (TH1D*)input->Get(Form("%s_dxyMinus_%d",name,i));
    dzP[i] = (TH1D*)input->Get(Form("%s_dzPlus_%d",name,i));
    dzM[i] = (TH1D*)input->Get(Form("%s_dzMinus_%d",name,i));
    chi2TrkP[i] = (TH1D*)input->Get(Form("%s_chi2TrkPlus_%d",name,i));
    chi2TrkM[i] = (TH1D*)input->Get(Form("%s_chi2TrkMinus_%d",name,i));
    chi2GlbP[i] = (TH1D*)input->Get(Form("%s_chi2GlbPlus_%d",name,i));
    chi2GlbM[i] = (TH1D*)input->Get(Form("%s_chi2GlbMinus_%d",name,i));
    nTrkHitsP[i] = (TH1D*)input->Get(Form("%s_nTrkHitsPlus_%d",name,i));
    nTrkHitsM[i] = (TH1D*)input->Get(Form("%s_nTrkHitsMinus_%d",name,i));
    nPixWMeaP[i] = (TH1D*)input->Get(Form("%s_nPixWMeaPlus_%d",name,i));
    nPixWMeaM[i] = (TH1D*)input->Get(Form("%s_nPixWMeaMinus_%d",name,i));

    ctauCtauErrP[i] = (TH2D*)input->Get(Form("%s_ctauCtauErrPlus_%d",name,i));
    ctauCtauErrM[i] = (TH2D*)input->Get(Form("%s_ctauCtauErrMinus_%d",name,i));
    ratioCtauCtauErr[i] = (TH2D*)input->Get(Form("%s_ctauCtauErrRatio_%d",name,i));

    ctauPhiP[i] = (TH2D*)input->Get(Form("%s_ctauPhiPlus_%d",name,i));
    ctauPhiM[i] = (TH2D*)input->Get(Form("%s_ctauPhiMinus_%d",name,i));
    ratioCtauPhi[i] = (TH2D*)input->Get(Form("%s_ctauPhiRatio_%d",name,i));

    ctauPhiMuP[i] = (TH2D*)input->Get(Form("%s_ctauPhiMuPlus_%d",name,i));
    ctauPhiMuM[i] = (TH2D*)input->Get(Form("%s_ctauPhiMuMinus_%d",name,i));
    ratioCtauPhiMu[i] = (TH2D*)input->Get(Form("%s_ctauPhiMuRatio_%d",name,i));

    ctauEtaP[i] = (TH2D*)input->Get(Form("%s_ctauEtaPlus_%d",name,i));
    ctauEtaM[i] = (TH2D*)input->Get(Form("%s_ctauEtaMinus_%d",name,i));
    ratioCtauEta[i] = (TH2D*)input->Get(Form("%s_ctauEtaRatio_%d",name,i));

    ctauPTP[i] = (TH2D*)input->Get(Form("%s_ctauPTPlus_%d",name,i));
    ctauPTM[i] = (TH2D*)input->Get(Form("%s_ctauPTMinus_%d",name,i));
    ratioCtauPT[i] = (TH2D*)input->Get(Form("%s_ctauPTRatio_%d",name,i));

    ctauDxyP[i] = (TH2D*)input->Get(Form("%s_ctauDxyPlus_%d",name,i));
    ctauDxyM[i] = (TH2D*)input->Get(Form("%s_ctauDxyMinus_%d",name,i));
    ratioCtauDxy[i] = (TH2D*)input->Get(Form("%s_ctauDxyRatio_%d",name,i));

    ctauDzP[i] = (TH2D*)input->Get(Form("%s_ctauDzPlus_%d",name,i));
    ctauDzM[i] = (TH2D*)input->Get(Form("%s_ctauDzMinus_%d",name,i));
    ratioCtauDz[i] = (TH2D*)input->Get(Form("%s_ctauDzRatio_%d",name,i));

    ctauChi2TrkP[i] = (TH2D*)input->Get(Form("%s_ctauChi2TrkPlus_%d",name,i));
    ctauChi2TrkM[i] = (TH2D*)input->Get(Form("%s_ctauChi2TrkMinus_%d",name,i));
    ratioCtauChi2Trk[i] = (TH2D*)input->Get(Form("%s_ctauChi2TrkRatio_%d",name,i));

    ctauChi2GlbP[i] = (TH2D*)input->Get(Form("%s_ctauChi2GlbPlus_%d",name,i));
    ctauChi2GlbM[i] = (TH2D*)input->Get(Form("%s_ctauChi2GlbMinus_%d",name,i));
    ratioCtauChi2Glb[i] = (TH2D*)input->Get(Form("%s_ctauChi2GlbRatio_%d",name,i));

    ctauNTrkHitsP[i] = (TH2D*)input->Get(Form("%s_ctauNTrkHitsPlus_%d",name,i));
    ctauNTrkHitsM[i] = (TH2D*)input->Get(Form("%s_ctauNTrkHitsMinus_%d",name,i));
    ratioCtauNTrkHits[i] = (TH2D*)input->Get(Form("%s_ctauNTrkHitsRatio_%d",name,i));

    ctauNPixWMeaP[i] = (TH2D*)input->Get(Form("%s_ctauNPixWMeaPlus_%d",name,i));
    ctauNPixWMeaM[i] = (TH2D*)input->Get(Form("%s_ctauNPixWMeaMinus_%d",name,i));
    ratioCtauNPixWMea[i] = (TH2D*)input->Get(Form("%s_ctauNPixWMeaRatio_%d",name,i));

    etaPhiMuP[i] = (TH2D*)input->Get(Form("%s_etaPhiMuPlus_%d",name,i));
    etaPhiMuM[i] = (TH2D*)input->Get(Form("%s_etaPhiMuMinus_%d",name,i));
    ratioEtaPhiMu[i] = (TH2D*)input->Get(Form("%s_etaPhiMuRatio_%d",name,i));
  }

  yPA = (TH1D*)input->Get(Form("%s_yPlusAll",name));
  yMA = (TH1D*)input->Get(Form("%s_yMinusAll",name));
  ratioYA = (TH1D*)input->Get(Form("%s_yRatioAll",name));

  ctauYPA = (TH2D*)input->Get(Form("%s_ctauYPlusAll",name));
  ctauYMA = (TH2D*)input->Get(Form("%s_ctauYMinusAll",name));
  ratioCtauYA = (TH2D*)input->Get(Form("%s_ctauYAllRatio",name));

  yPhiP = (TH2D*)input->Get(Form("%s_yPhiPlus",name));
  yPhiM = (TH2D*)input->Get(Form("%s_yPhiMinus",name));
  ratioYPhi = (TH2D*)input->Get(Form("%s_yPhiRatio",name));

  for (unsigned int i=0; i<nRuns; i++) {
    runPhiP[i] = (TH1D*)input->Get(Form("%s_runPhiP_%d",name,i));
    runPhiM[i] = (TH1D*)input->Get(Form("%s_runPhiM_%d",name,i));
    runRatioPhi[i] = (TH1D*)input->Get(Form("%s_runRatioPhi_%d",name,i));
    runPhiMuP[i] = (TH1D*)input->Get(Form("%s_runPhiMuP_%d",name,i));
    runPhiMuM[i] = (TH1D*)input->Get(Form("%s_runPhiMuM_%d",name,i));
    runRatioPhiMu[i] = (TH1D*)input->Get(Form("%s_runRatioPhiMu_%d",name,i));
    runNTrkHitsP[i] = (TH1D*)input->Get(Form("%s_runNTrkHitsP_%d",name,i));
    runNTrkHitsM[i] = (TH1D*)input->Get(Form("%s_runNTrkHitsM_%d",name,i));
    runRatioNTrkHits[i] = (TH1D*)input->Get(Form("%s_runRatioTrkHits_%d",name,i));
  }

}

void LifetimeCheck::CreateHistos(const char *name) {
  for (unsigned int i=0; i<rapbins.size(); i++) {
//    ratio[i] = new TH1D(Form("%s_ctauRatio_%d",name,i),titleRatio.c_str(),nBins,xmin,xmax); ratio[i]->Sumw2();
//    ratioErr[i] = new TH1D(Form("%s_ctauErrRatio_%d",name,i),titleRatioErr.c_str(),nBinsErr,xminErr,xmaxErr); ratioErr[i]->Sumw2();
    ratio[i] = new TH1D(Form("%s_ctauRatio_%d",name,i),titleRatio.c_str(),nBins,xArray); ratio[i]->Sumw2();
    ratioErr[i] = new TH1D(Form("%s_ctauErrRatio_%d",name,i),titleRatioErr.c_str(),nBinsErr,xErrArray); ratioErr[i]->Sumw2();
    ratioEta[i] = new TH1D(Form("%s_etaRatio_%d",name,i),titleRatioEta.c_str(),nBinsEta,xminEta,xmaxEta); ratioEta[i]->Sumw2();
    ratioPhi[i] = new TH1D(Form("%s_phiRatio_%d",name,i),titleRatioPhi.c_str(),nBinsPhi,xminPhi,xmaxPhi); ratioPhi[i]->Sumw2();
    ratioPhiMu[i] = new TH1D(Form("%s_phiMuRatio_%d",name,i),titleRatioPhiMu.c_str(),nBinsPhi,xminPhi,xmaxPhi); ratioPhiMu[i]->Sumw2();
//    ratioPT[i] = new TH1D(Form("%s_pTRatio_%d",name,i),titleRatioPT.c_str(),nBinsPT,xminPT,xmaxPT); ratioPT[i]->Sumw2();
    ratioPT[i] = new TH1D(Form("%s_pTRatio_%d",name,i),titleRatioPT.c_str(),nBinsPT,xPTArray); ratioPT[i]->Sumw2();
    ratioDxy[i] = new TH1D(Form("%s_dxyRatio_%d",name,i),titleRatioDxy.c_str(),nBinsDxy,xDxyArray); ratioDxy[i]->Sumw2();
    ratioDz[i] = new TH1D(Form("%s_dzRatio_%d",name,i),titleRatioDz.c_str(),nBinsDz,xDzArray); ratioDz[i]->Sumw2();
    ratioChi2Trk[i] = new TH1D(Form("%s_chi2TrkRatio_%d",name,i),titleRatioChi2Trk.c_str(),nBinsChi2Trk,xChi2TrkArray); ratioChi2Trk[i]->Sumw2();
    ratioChi2Glb[i] = new TH1D(Form("%s_chi2GlbRatio_%d",name,i),titleRatioChi2Glb.c_str(),nBinsChi2Glb,xChi2GlbArray); ratioChi2Glb[i]->Sumw2();
    ratioNTrkHits[i] = new TH1D(Form("%s_nTrkHitsRatio_%d",name,i),titleRatioNTrkHits.c_str(),nBinsNTrkHits,xNTrkHitsArray); ratioNTrkHits[i]->Sumw2();
    ratioNPixWMea[i] = new TH1D(Form("%s_nPixWMeaRatio_%d",name,i),titleRatioNPixWMea.c_str(),nBinsNPixWMea,xNPixWMeaArray); ratioNPixWMea[i]->Sumw2();
    
//    ctauP[i] = new TH1D(Form("%s_ctauPlus_%d",name,i),titleRaw.c_str(),nBins,xmin,xmax); ctauP[i]->Sumw2();
//    ctauM[i] = new TH1D(Form("%s_ctauMinus_%d",name,i),titleRaw.c_str(),nBins,xmin,xmax); ctauM[i]->Sumw2();
//    ctauErrP[i] = new TH1D(Form("%s_ctauErrPlus_%d",name,i),titleRawErr.c_str(),nBinsErr,xminErr,xmaxErr); ctauErrP[i]->Sumw2();
//    ctauErrM[i] = new TH1D(Form("%s_ctauErrMinus_%d",name,i),titleRawErr.c_str(),nBinsErr,xminErr,xmaxErr); ctauErrM[i]->Sumw2();
    ctauErrP[i] = new TH1D(Form("%s_ctauErrPlus_%d",name,i),titleRawErr.c_str(),nBinsErr,xErrArray); ctauErrP[i]->Sumw2();
    ctauErrM[i] = new TH1D(Form("%s_ctauErrMinus_%d",name,i),titleRawErr.c_str(),nBinsErr,xErrArray); ctauErrM[i]->Sumw2();
    ctauP[i] = new TH1D(Form("%s_ctauPlus_%d",name,i),titleRaw.c_str(),nBins,xArray); ctauP[i]->Sumw2();
    ctauM[i] = new TH1D(Form("%s_ctauMinus_%d",name,i),titleRaw.c_str(),nBins,xArray); ctauM[i]->Sumw2();
    etaP[i] = new TH1D(Form("%s_etaPlus_%d",name,i),titleEta.c_str(),nBinsEta,xminEta,xmaxEta); etaP[i]->Sumw2();
    etaM[i] = new TH1D(Form("%s_etaMinus_%d",name,i),titleEta.c_str(),nBinsEta,xminEta,xmaxEta); etaM[i]->Sumw2();
    phiP[i] = new TH1D(Form("%s_phiPlus_%d",name,i),titlePhi.c_str(),nBinsPhi,xminPhi,xmaxPhi); phiP[i]->Sumw2();
    phiM[i] = new TH1D(Form("%s_phiMinus_%d",name,i),titlePhi.c_str(),nBinsPhi,xminPhi,xmaxPhi); phiM[i]->Sumw2();
    phiMuP[i] = new TH1D(Form("%s_phiMuPlus_%d",name,i),titlePhiMu.c_str(),nBinsPhi,xminPhi,xmaxPhi); phiMuP[i]->Sumw2();
    phiMuM[i] = new TH1D(Form("%s_phiMuMinus_%d",name,i),titlePhiMu.c_str(),nBinsPhi,xminPhi,xmaxPhi); phiMuM[i]->Sumw2();
//    pTP[i] = new TH1D(Form("%s_pTPlus_%d",name,i),titlePT.c_str(),nBinsPT,xminPT,xmaxPT); pTP[i]->Sumw2();
//    pTM[i] = new TH1D(Form("%s_pTMinus_%d",name,i),titlePT.c_str(),nBinsPT,xminPT,xmaxPT); pTM[i]->Sumw2();
    pTP[i] = new TH1D(Form("%s_pTPlus_%d",name,i),titlePT.c_str(),nBinsPT,xPTArray); pTP[i]->Sumw2();
    pTM[i] = new TH1D(Form("%s_pTMinus_%d",name,i),titlePT.c_str(),nBinsPT,xPTArray); pTM[i]->Sumw2();
    dxyP[i] = new TH1D(Form("%s_dxyPlus_%d",name,i),titleDxy.c_str(),nBinsDxy,xDxyArray); dxyP[i]->Sumw2();
    dxyM[i] = new TH1D(Form("%s_dxyMinus_%d",name,i),titleDxy.c_str(),nBinsDxy,xDxyArray); dxyM[i]->Sumw2();
    dzP[i] = new TH1D(Form("%s_dzPlus_%d",name,i),titleDz.c_str(),nBinsDz,xDzArray); dzP[i]->Sumw2();
    dzM[i] = new TH1D(Form("%s_dzMinus_%d",name,i),titleDz.c_str(),nBinsDz,xDzArray); dzM[i]->Sumw2();
    chi2TrkP[i] = new TH1D(Form("%s_chi2TrkPlus_%d",name,i),titleChi2Trk.c_str(),nBinsChi2Trk,xChi2TrkArray); chi2TrkP[i]->Sumw2();
    chi2TrkM[i] = new TH1D(Form("%s_chi2TrkMinus_%d",name,i),titleChi2Trk.c_str(),nBinsChi2Trk,xChi2TrkArray); chi2TrkM[i]->Sumw2();
    chi2GlbP[i] = new TH1D(Form("%s_chi2GlbPlus_%d",name,i),titleChi2Glb.c_str(),nBinsChi2Glb,xChi2GlbArray); chi2GlbP[i]->Sumw2();
    chi2GlbM[i] = new TH1D(Form("%s_chi2GlbMinus_%d",name,i),titleChi2Glb.c_str(),nBinsChi2Glb,xChi2GlbArray); chi2GlbM[i]->Sumw2();
    nTrkHitsP[i] = new TH1D(Form("%s_nTrkHitsPlus_%d",name,i),titleNTrkHits.c_str(),nBinsNTrkHits,xNTrkHitsArray); nTrkHitsP[i]->Sumw2();
    nTrkHitsM[i] = new TH1D(Form("%s_nTrkHitsMinus_%d",name,i),titleNTrkHits.c_str(),nBinsNTrkHits,xNTrkHitsArray); nTrkHitsM[i]->Sumw2();
    nPixWMeaP[i] = new TH1D(Form("%s_nPixWMeaPlus_%d",name,i),titleNPixWMea.c_str(),nBinsNPixWMea,xNPixWMeaArray); nPixWMeaP[i]->Sumw2();
    nPixWMeaM[i] = new TH1D(Form("%s_nPixWMeaMinus_%d",name,i),titleNPixWMea.c_str(),nBinsNPixWMea,xNPixWMeaArray); nPixWMeaM[i]->Sumw2();

//    ctauCtauErrP[i] = new TH2D(Form("%s_ctauCtauErrPlus_%d",name,i),title2DCtauErr.c_str(),nBins,xmin,xmax,nBinsErr,xminErr,xmaxErr); ctauCtauErrP[i]->Sumw2();
//    ctauCtauErrM[i] = new TH2D(Form("%s_ctauCtauErrMinus_%d",name,i),title2DCtauErr.c_str(),nBins,xmin,xmax,nBinsErr,xminErr,xmaxErr); ctauCtauErrM[i]->Sumw2();
//    ratioCtauCtauErr[i] = new TH2D(Form("%s_ctauCtauErrRatio_%d",name,i),title2DCtauErr.c_str(),nBins,xmin,xmax,nBinsErr,xminErr,xmaxErr); ratioCtauCtauErr[i]->Sumw2();
    ctauCtauErrP[i] = new TH2D(Form("%s_ctauCtauErrPlus_%d",name,i),title2DCtauErr.c_str(),nBins,xmin,xmax,nBinsErr,xErrArray); ctauCtauErrP[i]->Sumw2();
    ctauCtauErrM[i] = new TH2D(Form("%s_ctauCtauErrMinus_%d",name,i),title2DCtauErr.c_str(),nBins,xmin,xmax,nBinsErr,xErrArray); ctauCtauErrM[i]->Sumw2();
    ratioCtauCtauErr[i] = new TH2D(Form("%s_ctauCtauErrRatio_%d",name,i),title2DCtauErr.c_str(),nBins,xmin,xmax,nBinsErr,xErrArray); ratioCtauCtauErr[i]->Sumw2();

    ctauPhiP[i] = new TH2D(Form("%s_ctauPhiPlus_%d",name,i),title2DPhi.c_str(),nBins,xmin,xmax,nBinsPhi,xminPhi,xmaxPhi); ctauPhiP[i]->Sumw2();
    ctauPhiM[i] = new TH2D(Form("%s_ctauPhiMinus_%d",name,i),title2DPhi.c_str(),nBins,xmin,xmax,nBinsPhi,xminPhi,xmaxPhi); ctauPhiM[i]->Sumw2();
    ratioCtauPhi[i] = new TH2D(Form("%s_ctauPhiRatio_%d",name,i),title2DPhi.c_str(),nBins,xmin,xmax,nBinsPhi,xminPhi,xmaxPhi); ratioCtauPhi[i]->Sumw2();

    ctauPhiMuP[i] = new TH2D(Form("%s_ctauPhiMuPlus_%d",name,i),title2DPhiMu.c_str(),nBins,xmin,xmax,nBinsPhi,xminPhi,xmaxPhi); ctauPhiMuP[i]->Sumw2();
    ctauPhiMuM[i] = new TH2D(Form("%s_ctauPhiMuMinus_%d",name,i),title2DPhiMu.c_str(),nBins,xmin,xmax,nBinsPhi,xminPhi,xmaxPhi); ctauPhiMuM[i]->Sumw2();
    ratioCtauPhiMu[i] = new TH2D(Form("%s_ctauPhiMuRatio_%d",name,i),title2DPhiMu.c_str(),nBins,xmin,xmax,nBinsPhi,xminPhi,xmaxPhi); ratioCtauPhiMu[i]->Sumw2();

    ctauEtaP[i] = new TH2D(Form("%s_ctauEtaPlus_%d",name,i),title2DEta.c_str(),nBins,xmin,xmax,nBinsEta,xminEta,xmaxEta); ctauEtaP[i]->Sumw2();
    ctauEtaM[i] = new TH2D(Form("%s_ctauEtaMinus_%d",name,i),title2DEta.c_str(),nBins,xmin,xmax,nBinsEta,xminEta,xmaxEta); ctauEtaM[i]->Sumw2();
    ratioCtauEta[i] = new TH2D(Form("%s_ctauEtaRatio_%d",name,i),title2DEta.c_str(),nBins,xmin,xmax,nBinsEta,xminEta,xmaxEta); ratioCtauEta[i]->Sumw2();

//    ctauPTP[i] = new TH2D(Form("%s_ctauPTPlus_%d",name,i),title2DPT.c_str(),nBins,xmin,xmax,nBinsPT,xminPT,xmaxPT); ctauPTP[i]->Sumw2();
//    ctauPTM[i] = new TH2D(Form("%s_ctauPTMinus_%d",name,i),title2DPT.c_str(),nBins,xmin,xmax,nBinsPT,xminPT,xmaxPT); ctauPTM[i]->Sumw2();
//    ratioCtauPT[i] = new TH2D(Form("%s_ctauPTRatio_%d",name,i),title2DPT.c_str(),nBins,xmin,xmax,nBinsPT,xminPT,xmaxPT); ratioCtauPT[i]->Sumw2();
    ctauPTP[i] = new TH2D(Form("%s_ctauPTPlus_%d",name,i),title2DPT.c_str(),nBins,xmin,xmax,nBinsPT,xPTArray); ctauPTP[i]->Sumw2();
    ctauPTM[i] = new TH2D(Form("%s_ctauPTMinus_%d",name,i),title2DPT.c_str(),nBins,xmin,xmax,nBinsPT,xPTArray); ctauPTM[i]->Sumw2();
    ratioCtauPT[i] = new TH2D(Form("%s_ctauPTRatio_%d",name,i),title2DPT.c_str(),nBins,xmin,xmax,nBinsPT,xPTArray); ratioCtauPT[i]->Sumw2();

    ctauDxyP[i] = new TH2D(Form("%s_ctauDxyPlus_%d",name,i),title2DDxy.c_str(),nBins,xmin,xmax,nBinsDxy,xDxyArray); ctauDxyP[i]->Sumw2();
    ctauDxyM[i] = new TH2D(Form("%s_ctauDxyMinus_%d",name,i),title2DDxy.c_str(),nBins,xmin,xmax,nBinsDxy,xDxyArray); ctauDxyM[i]->Sumw2();
    ratioCtauDxy[i] = new TH2D(Form("%s_ctauDxyRatio_%d",name,i),title2DDxy.c_str(),nBins,xmin,xmax,nBinsDxy,xDxyArray); ratioCtauDxy[i]->Sumw2();

    ctauDzP[i] = new TH2D(Form("%s_ctauDzPlus_%d",name,i),title2DDz.c_str(),nBins,xmin,xmax,nBinsDz,xDzArray); ctauDzP[i]->Sumw2();
    ctauDzM[i] = new TH2D(Form("%s_ctauDzMinus_%d",name,i),title2DDz.c_str(),nBins,xmin,xmax,nBinsDz,xDzArray); ctauDzM[i]->Sumw2();
    ratioCtauDz[i] = new TH2D(Form("%s_ctauDzRatio_%d",name,i),title2DDz.c_str(),nBins,xmin,xmax,nBinsDz,xDzArray); ratioCtauDz[i]->Sumw2();

    ctauChi2TrkP[i] = new TH2D(Form("%s_ctauChi2TrkPlus_%d",name,i),title2DChi2Trk.c_str(),nBins,xmin,xmax,nBinsChi2Trk,xChi2TrkArray); ctauChi2TrkP[i]->Sumw2();
    ctauChi2TrkM[i] = new TH2D(Form("%s_ctauChi2TrkMinus_%d",name,i),title2DChi2Trk.c_str(),nBins,xmin,xmax,nBinsChi2Trk,xChi2TrkArray); ctauChi2TrkM[i]->Sumw2();
    ratioCtauChi2Trk[i] = new TH2D(Form("%s_ctauChi2TrkRatio_%d",name,i),title2DChi2Trk.c_str(),nBins,xmin,xmax,nBinsChi2Trk,xChi2TrkArray); ratioCtauChi2Trk[i]->Sumw2();

    ctauChi2GlbP[i] = new TH2D(Form("%s_ctauChi2GlbPlus_%d",name,i),title2DChi2Glb.c_str(),nBins,xmin,xmax,nBinsChi2Glb,xChi2GlbArray); ctauChi2GlbP[i]->Sumw2();
    ctauChi2GlbM[i] = new TH2D(Form("%s_ctauChi2GlbMinus_%d",name,i),title2DChi2Glb.c_str(),nBins,xmin,xmax,nBinsChi2Glb,xChi2GlbArray); ctauChi2GlbM[i]->Sumw2();
    ratioCtauChi2Glb[i] = new TH2D(Form("%s_ctauChi2GlbRatio_%d",name,i),title2DChi2Glb.c_str(),nBins,xmin,xmax,nBinsChi2Glb,xChi2GlbArray); ratioCtauChi2Glb[i]->Sumw2();

    ctauNTrkHitsP[i] = new TH2D(Form("%s_ctauNTrkHitsPlus_%d",name,i),title2DNTrkHits.c_str(),nBins,xmin,xmax,nBinsNTrkHits,xNTrkHitsArray); ctauNTrkHitsP[i]->Sumw2();
    ctauNTrkHitsM[i] = new TH2D(Form("%s_ctauNTrkHitsMinus_%d",name,i),title2DNTrkHits.c_str(),nBins,xmin,xmax,nBinsNTrkHits,xNTrkHitsArray); ctauNTrkHitsM[i]->Sumw2();
    ratioCtauNTrkHits[i] = new TH2D(Form("%s_ctauNTrkHitsRatio_%d",name,i),title2DNTrkHits.c_str(),nBins,xmin,xmax,nBinsNTrkHits,xNTrkHitsArray); ratioCtauNTrkHits[i]->Sumw2();

    ctauNPixWMeaP[i] = new TH2D(Form("%s_ctauNPixWMeaPlus_%d",name,i),title2DNPixWMea.c_str(),nBins,xmin,xmax,nBinsNPixWMea,xNPixWMeaArray); ctauNPixWMeaP[i]->Sumw2();
    ctauNPixWMeaM[i] = new TH2D(Form("%s_ctauNPixWMeaMinus_%d",name,i),title2DNPixWMea.c_str(),nBins,xmin,xmax,nBinsNPixWMea,xNPixWMeaArray); ctauNPixWMeaM[i]->Sumw2();
    ratioCtauNPixWMea[i] = new TH2D(Form("%s_ctauNPixWMeaRatio_%d",name,i),title2DNPixWMea.c_str(),nBins,xmin,xmax,nBinsNPixWMea,xNPixWMeaArray); ratioCtauNPixWMea[i]->Sumw2();

    etaPhiMuP[i] = new TH2D(Form("%s_etaPhiMuPlus_%d",name,i),title2DEtaPhiMu.c_str(),nBinsEta,xminEta,xmaxEta,nBinsPhi,xminPhi,xmaxPhi); etaPhiMuP[i]->Sumw2();
    etaPhiMuM[i] = new TH2D(Form("%s_etaPhiMuMinus_%d",name,i),title2DEtaPhiMu.c_str(),nBinsEta,xminEta,xmaxEta,nBinsPhi,xminPhi,xmaxPhi); etaPhiMuM[i]->Sumw2();
    ratioEtaPhiMu[i] = new TH2D(Form("%s_etaPhiMuRatio_%d",name,i),title2DEtaPhiMu.c_str(),nBinsEta,xminEta,xmaxEta,nBinsPhi,xminPhi,xmaxPhi); ratioEtaPhiMu[i]->Sumw2();

  }

  yPA = new TH1D(Form("%s_yPlusAll",name),titleY.c_str(),nBinsY,xminY,xmaxY); yPA->Sumw2();
  yMA = new TH1D(Form("%s_yMinusAll",name),titleY.c_str(),nBinsY,xminY,xmaxY); yMA->Sumw2();
  ratioYA = new TH1D(Form("%s_yRatioAll",name),titleRatioY.c_str(),nBinsY,xminY,xmaxY); ratioYA->Sumw2();
  
  ctauYPA = new TH2D(Form("%s_ctauYPlusAll",name),title2DY.c_str(),nBins,xmin,xmax,nBinsY,xminY,xmaxY); ctauYPA->Sumw2();
  ctauYMA = new TH2D(Form("%s_ctauYMinusAll",name),title2DY.c_str(),nBins,xmin,xmax,nBinsY,xminY,xmaxY); ctauYMA->Sumw2();
  ratioCtauYA = new TH2D(Form("%s_ctauYAllRatio",name),title2DY.c_str(),nBins,xmin,xmax,nBinsY,xminY,xmaxY); ratioCtauYA->Sumw2();
  
  yPhiP = new TH2D(Form("%s_yPhiPlus",name),title2DYPhi.c_str(),nBinsY,xminY,xmaxY,nBinsPhi,xminPhi,xmaxPhi); yPhiP->Sumw2();
  yPhiM = new TH2D(Form("%s_yPhiMinus",name),title2DYPhi.c_str(),nBinsY,xminY,xmaxY,nBinsPhi,xminPhi,xmaxPhi); yPhiM->Sumw2();
  ratioYPhi = new TH2D(Form("%s_yPhiRatio",name),title2DYPhi.c_str(),nBinsY,xminY,xmaxY,nBinsPhi,xminPhi,xmaxPhi); ratioYPhi->Sumw2();

  for (unsigned int i=0; i<nRuns; i++) {
    runPhiP[i] = new TH1D(Form("%s_runPhiP_%d",name,i),titlePhi.c_str(),nBinsPhi,xminPhi,xmaxPhi); runPhiP[i]->Sumw2();
    runPhiM[i] = new TH1D(Form("%s_runPhiM_%d",name,i),titlePhi.c_str(),nBinsPhi,xminPhi,xmaxPhi); runPhiM[i]->Sumw2();
    runRatioPhi[i] = new TH1D(Form("%s_runRatioPhi_%d",name,i),titleRatioPhi.c_str(),nBinsPhi,xminPhi,xmaxPhi); runRatioPhi[i]->Sumw2();
    runPhiMuP[i] = new TH1D(Form("%s_runPhiMuP_%d",name,i),titlePhiMu.c_str(),nBinsPhi,xminPhi,xmaxPhi); runPhiMuP[i]->Sumw2();
    runPhiMuM[i] = new TH1D(Form("%s_runPhiMuM_%d",name,i),titlePhiMu.c_str(),nBinsPhi,xminPhi,xmaxPhi); runPhiMuM[i]->Sumw2();
    runRatioPhiMu[i] = new TH1D(Form("%s_runRatioPhiMu_%d",name,i),titleRatioPhiMu.c_str(),nBinsPhi,xminPhi,xmaxPhi); runRatioPhiMu[i]->Sumw2();
    runNTrkHitsP[i] = new TH1D(Form("%s_runNTrkHitsP_%d",name,i),titleNTrkHits.c_str(),nBinsNTrkHits,xNTrkHitsArray); runNTrkHitsP[i]->Sumw2();
    runNTrkHitsM[i] = new TH1D(Form("%s_runNTrkHitsM_%d",name,i),titleNTrkHits.c_str(),nBinsNTrkHits,xNTrkHitsArray); runNTrkHitsM[i]->Sumw2();
    runRatioNTrkHits[i] = new TH1D(Form("%s_runRatioTrkHits_%d",name,i),titleRatioNTrkHits.c_str(),nBinsNTrkHits,xNTrkHitsArray); runRatioNTrkHits[i]->Sumw2();
  }
  
  cout << "Memory test: " << ctauEtaM[0] << endl;

}

void LifetimeCheck::ReadTree() {
  float Reco_QQ_ctau[10000], Reco_QQ_ctauErr[10000];
  float Reco_QQ_mupl_dxy[10000], Reco_QQ_mumi_dxy[10000];
  float Reco_QQ_mupl_dz[10000], Reco_QQ_mumi_dz[10000];
  float Reco_QQ_mupl_norChi2_inner[10000], Reco_QQ_mumi_norChi2_inner[10000];
  float Reco_QQ_mupl_norChi2_global[10000], Reco_QQ_mumi_norChi2_global[10000];
  int Reco_QQ_mupl_nTrkHits[10000], Reco_QQ_mumi_nTrkHits[10000];
  int Reco_QQ_mupl_nPixWMea[10000], Reco_QQ_mumi_nPixWMea[10000];
  
  UInt_t runNb;
  int Reco_QQ_trig[10000], HLTriggers, Reco_QQ_size;
  TClonesArray *Reco_QQ_4mom, *Reco_QQ_mupl_4mom, *Reco_QQ_mumi_4mom;
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;

  myTree->SetBranchAddress("runNb",&runNb);
  myTree->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4mom);
  myTree->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4mom);
  myTree->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4mom);
  myTree->SetBranchAddress("Reco_QQ_mupl_dxy",&Reco_QQ_mupl_dxy);
  myTree->SetBranchAddress("Reco_QQ_mumi_dxy",&Reco_QQ_mumi_dxy);
  myTree->SetBranchAddress("Reco_QQ_mupl_dz",&Reco_QQ_mupl_dz);
  myTree->SetBranchAddress("Reco_QQ_mumi_dz",&Reco_QQ_mumi_dz);
  myTree->SetBranchAddress("Reco_QQ_mupl_norChi2_inner",&Reco_QQ_mupl_norChi2_inner);
  myTree->SetBranchAddress("Reco_QQ_mumi_norChi2_inner",&Reco_QQ_mumi_norChi2_inner);
  myTree->SetBranchAddress("Reco_QQ_mupl_norChi2_global",&Reco_QQ_mupl_norChi2_global);
  myTree->SetBranchAddress("Reco_QQ_mumi_norChi2_global",&Reco_QQ_mumi_norChi2_global);
  myTree->SetBranchAddress("Reco_QQ_mupl_nTrkHits",&Reco_QQ_mupl_nTrkHits);
  myTree->SetBranchAddress("Reco_QQ_mumi_nTrkHits",&Reco_QQ_mumi_nTrkHits);
  myTree->SetBranchAddress("Reco_QQ_mupl_nPixWMea",&Reco_QQ_mupl_nPixWMea);
  myTree->SetBranchAddress("Reco_QQ_mumi_nPixWMea",&Reco_QQ_mumi_nPixWMea);
  myTree->SetBranchAddress("Reco_QQ_ctau",&Reco_QQ_ctau);
  myTree->SetBranchAddress("Reco_QQ_ctauErr",&Reco_QQ_ctauErr);
  myTree->SetBranchAddress("Reco_QQ_size",&Reco_QQ_size);
  myTree->SetBranchAddress("Reco_QQ_trig",&Reco_QQ_trig);
  myTree->SetBranchAddress("HLTriggers",&HLTriggers);

//  for (int ev=0; ev<100000; ev++) {
  for (int ev=0; ev<myTree->GetEntries(); ev++) {
    if (ev%100000==0) cout << ">>>>> EVENT " << ev << " / " << myTree->GetEntries() <<  endl;
    myTree->GetEntry(ev);

    for (int idx=0; idx<Reco_QQ_size; idx++) {
      TLorentzVector *Jpsi = (TLorentzVector*)Reco_QQ_4mom->At(idx);
      TLorentzVector *muP = (TLorentzVector*)Reco_QQ_mupl_4mom->At(idx);
      TLorentzVector *muM = (TLorentzVector*)Reco_QQ_mumi_4mom->At(idx);

      if (IsCowboy(+1, muP->Phi(), -1, muM->Phi()) && !FillCowboy) continue;

      if ( (Reco_QQ_trig[idx]&1)==1 && (HLTriggers&1)==1 && Reco_QQ_ctau[idx]>xmin && Reco_QQ_ctau[idx]<xmax && Jpsi->Rapidity()>-2.4 && Jpsi->Rapidity()<2.4 && Jpsi->M()>mmin && Jpsi->M()<mmax && Jpsi->Pt()>=ptmin && Jpsi->Pt()<ptmax) {
        for (unsigned int i=0; i<rapbins.size(); i++) {
          if ( ((i<rapbins.size()-1) && Jpsi->Rapidity()<-1*rapbins[i] && Jpsi->Rapidity()>=-1*rapbins[i+1]) ||
               ((i==rapbins.size()-1) && Jpsi->Rapidity()<0 && Jpsi->Rapidity()>=-1*xmaxEta)
             ) {
            ctauM[i]->Fill(Reco_QQ_ctau[idx]);
            ctauErrM[i]->Fill(Reco_QQ_ctauErr[idx]);
            etaM[i]->Fill(TMath::Abs(muP->Eta()));
            etaM[i]->Fill(TMath::Abs(muM->Eta()));
            phiM[i]->Fill(Jpsi->Phi());
            phiMuM[i]->Fill(muP->Phi());
            phiMuM[i]->Fill(muM->Phi());
            pTM[i]->Fill(Jpsi->Pt());
            dxyM[i]->Fill(Reco_QQ_mupl_dxy[idx]);
            dxyM[i]->Fill(Reco_QQ_mumi_dxy[idx]);
            dzM[i]->Fill(Reco_QQ_mupl_dz[idx]);
            dzM[i]->Fill(Reco_QQ_mumi_dz[idx]);
            chi2TrkM[i]->Fill(Reco_QQ_mupl_norChi2_inner[idx]);
            chi2TrkM[i]->Fill(Reco_QQ_mumi_norChi2_inner[idx]);
            chi2GlbM[i]->Fill(Reco_QQ_mupl_norChi2_global[idx]);
            chi2GlbM[i]->Fill(Reco_QQ_mumi_norChi2_global[idx]);
            nTrkHitsM[i]->Fill(Reco_QQ_mupl_nTrkHits[idx]);
            nTrkHitsM[i]->Fill(Reco_QQ_mumi_nTrkHits[idx]);
            nPixWMeaM[i]->Fill(Reco_QQ_mupl_nPixWMea[idx]);
            nPixWMeaM[i]->Fill(Reco_QQ_mumi_nPixWMea[idx]);
            
            ctauCtauErrM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_ctauErr[idx]);
            ctauPhiM[i]->Fill(Reco_QQ_ctau[idx],Jpsi->Phi());
            ctauPhiMuM[i]->Fill(Reco_QQ_ctau[idx],muP->Phi());
            ctauPhiMuM[i]->Fill(Reco_QQ_ctau[idx],muM->Phi());
            ctauEtaM[i]->Fill(Reco_QQ_ctau[idx],TMath::Abs(muP->Eta()));
            ctauEtaM[i]->Fill(Reco_QQ_ctau[idx],TMath::Abs(muM->Eta()));
            ctauPTM[i]->Fill(Reco_QQ_ctau[idx],Jpsi->Pt());
            ctauDxyM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_dxy[idx]);
            ctauDxyM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_dxy[idx]);
            ctauDzM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_dz[idx]);
            ctauDzM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_dz[idx]);
            ctauChi2TrkM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_norChi2_inner[idx]);
            ctauChi2TrkM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_norChi2_inner[idx]);
            ctauChi2GlbM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_norChi2_global[idx]);
            ctauChi2GlbM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_norChi2_global[idx]);
            ctauNTrkHitsM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_nTrkHits[idx]);
            ctauNTrkHitsM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_nTrkHits[idx]);
            ctauNPixWMeaM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_nPixWMea[idx]);
            ctauNPixWMeaM[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_nPixWMea[idx]);

            etaPhiMuM[i]->Fill(TMath::Abs(muP->Eta()),muP->Phi());
            etaPhiMuM[i]->Fill(TMath::Abs(muM->Eta()),muM->Phi());
//            etaPhiMuM[i]->Fill(muP->Eta(),muP->Phi());
//            etaPhiMuM[i]->Fill(muM->Eta(),muM->Phi());
          } else if ( ((i<rapbins.size()-1) && Jpsi->Rapidity()>=rapbins[i] && Jpsi->Rapidity()<rapbins[i+1]) ||
                      ((i==rapbins.size()-1) && Jpsi->Rapidity()>=0 && Jpsi->Rapidity()<xmaxEta)
                    ) {
            ctauP[i]->Fill(Reco_QQ_ctau[idx]);
            ctauErrP[i]->Fill(Reco_QQ_ctauErr[idx]);
            etaP[i]->Fill(TMath::Abs(muP->Eta()));
            etaP[i]->Fill(TMath::Abs(muM->Eta()));
            phiP[i]->Fill(Jpsi->Phi());
            phiMuP[i]->Fill(muP->Phi());
            phiMuP[i]->Fill(muM->Phi());
            pTP[i]->Fill(Jpsi->Pt());
            dxyP[i]->Fill(Reco_QQ_mupl_dxy[idx]);
            dxyP[i]->Fill(Reco_QQ_mumi_dxy[idx]);
            dzP[i]->Fill(Reco_QQ_mupl_dz[idx]);
            dzP[i]->Fill(Reco_QQ_mumi_dz[idx]);
            chi2TrkP[i]->Fill(Reco_QQ_mupl_norChi2_inner[idx]);
            chi2TrkP[i]->Fill(Reco_QQ_mumi_norChi2_inner[idx]);
            chi2GlbP[i]->Fill(Reco_QQ_mupl_norChi2_global[idx]);
            chi2GlbP[i]->Fill(Reco_QQ_mumi_norChi2_global[idx]);
            nTrkHitsP[i]->Fill(Reco_QQ_mupl_nTrkHits[idx]);
            nTrkHitsP[i]->Fill(Reco_QQ_mumi_nTrkHits[idx]);
            nPixWMeaP[i]->Fill(Reco_QQ_mupl_nPixWMea[idx]);
            nPixWMeaP[i]->Fill(Reco_QQ_mumi_nPixWMea[idx]);
            
            ctauCtauErrP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_ctauErr[idx]);
            ctauPhiP[i]->Fill(Reco_QQ_ctau[idx],Jpsi->Phi());
            ctauPhiMuP[i]->Fill(Reco_QQ_ctau[idx],muP->Phi());
            ctauPhiMuP[i]->Fill(Reco_QQ_ctau[idx],muM->Phi());
            ctauEtaP[i]->Fill(Reco_QQ_ctau[idx],TMath::Abs(muP->Eta()));
            ctauEtaP[i]->Fill(Reco_QQ_ctau[idx],TMath::Abs(muM->Eta()));
            ctauPTP[i]->Fill(Reco_QQ_ctau[idx],Jpsi->Pt());
            ctauDxyP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_dxy[idx]);
            ctauDxyP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_dxy[idx]);
            ctauDzP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_dz[idx]);
            ctauDzP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_dz[idx]);
            ctauChi2TrkP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_norChi2_inner[idx]);
            ctauChi2TrkP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_norChi2_inner[idx]);
            ctauChi2GlbP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_norChi2_global[idx]);
            ctauChi2GlbP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_norChi2_global[idx]);
            ctauNTrkHitsP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_nTrkHits[idx]);
            ctauNTrkHitsP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_nTrkHits[idx]);
            ctauNPixWMeaP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mupl_nPixWMea[idx]);
            ctauNPixWMeaP[i]->Fill(Reco_QQ_ctau[idx],Reco_QQ_mumi_nPixWMea[idx]);

            etaPhiMuP[i]->Fill(TMath::Abs(muP->Eta()),muP->Phi());
            etaPhiMuP[i]->Fill(TMath::Abs(muM->Eta()),muM->Phi());
//            etaPhiMuP[i]->Fill(muP->Eta(),muP->Phi());
//            etaPhiMuP[i]->Fill(muM->Eta(),muM->Phi());
          } 
        } //end: for loop of nRaps
       
        if (Jpsi->Rapidity()>=-1*xmaxEta && Jpsi->Rapidity()<0) {
          yMA->Fill(TMath::Abs(Jpsi->Rapidity()));
          ctauYMA->Fill(Reco_QQ_ctau[idx],TMath::Abs(Jpsi->Rapidity()));
          yPhiM->Fill(TMath::Abs(Jpsi->Rapidity()),Jpsi->Phi());

          for (unsigned int i=0; i<nRuns; i++) {
            if ( runNb>=runbins[i] && runNb<runbins[i+1] ) {
              runPhiM[i]->Fill(Jpsi->Phi());
              runPhiMuM[i]->Fill(muP->Phi());
              runPhiMuM[i]->Fill(muM->Phi());
              runNTrkHitsM[i]->Fill(Reco_QQ_mupl_nTrkHits[idx]);
              runNTrkHitsM[i]->Fill(Reco_QQ_mumi_nTrkHits[idx]);
            }
          }
        } else if (Jpsi->Rapidity()>=0 && Jpsi->Rapidity()<xmaxEta) {
          yPA->Fill(TMath::Abs(Jpsi->Rapidity()));
          ctauYPA->Fill(Reco_QQ_ctau[idx],TMath::Abs(Jpsi->Rapidity()));
          yPhiP->Fill(TMath::Abs(Jpsi->Rapidity()),Jpsi->Phi());

          for (unsigned int i=0; i<nRuns; i++) {
            if ( runNb>=runbins[i] && runNb<runbins[i+1] ) {
              runPhiP[i]->Fill(Jpsi->Phi());
              runPhiMuP[i]->Fill(muP->Phi());
              runPhiMuP[i]->Fill(muM->Phi());
              runNTrkHitsP[i]->Fill(Reco_QQ_mupl_nTrkHits[idx]);
              runNTrkHitsP[i]->Fill(Reco_QQ_mumi_nTrkHits[idx]);
            }
          }
        }

      } //end: selection condition test
    } //end: Reco_QQ_size
  } //end: myTree loop

  for (unsigned int i=0; i<rapbins.size(); i++) {
    ratio[i]->Divide(ctauP[i],ctauM[i]);
    ratioErr[i]->Divide(ctauErrP[i],ctauErrM[i]);
    ratioEta[i]->Divide(etaP[i],etaM[i]);
    ratioPhi[i]->Divide(phiP[i],phiM[i]);
    ratioPhiMu[i]->Divide(phiMuP[i],phiMuM[i]);
    ratioPT[i]->Divide(pTP[i],pTM[i]);
    ratioDxy[i]->Divide(dxyP[i],dxyM[i]);
    ratioDz[i]->Divide(dzP[i],dzM[i]);
    ratioChi2Trk[i]->Divide(chi2TrkP[i],chi2TrkM[i]);
    ratioChi2Glb[i]->Divide(chi2GlbP[i],chi2GlbM[i]);
    ratioNTrkHits[i]->Divide(nTrkHitsP[i],nTrkHitsM[i]);
    ratioNPixWMea[i]->Divide(nPixWMeaP[i],nPixWMeaM[i]);
    
    ratioCtauCtauErr[i]->Divide(ctauCtauErrP[i],ctauCtauErrM[i]);
    ratioCtauPhi[i]->Divide(ctauPhiP[i],ctauPhiM[i]);
    ratioCtauPhiMu[i]->Divide(ctauPhiMuP[i],ctauPhiMuM[i]);
    ratioCtauEta[i]->Divide(ctauEtaP[i],ctauEtaM[i]);
    ratioCtauPT[i]->Divide(ctauPTP[i],ctauPTM[i]);
    ratioCtauDxy[i]->Divide(ctauDxyP[i],ctauDxyM[i]);
    ratioCtauDz[i]->Divide(ctauDzP[i],ctauDzM[i]);
    ratioCtauChi2Trk[i]->Divide(ctauChi2TrkP[i],ctauChi2TrkM[i]);
    ratioCtauChi2Glb[i]->Divide(ctauChi2GlbP[i],ctauChi2GlbM[i]);
    ratioCtauNTrkHits[i]->Divide(ctauNTrkHitsP[i],ctauNTrkHitsM[i]);
    ratioCtauNPixWMea[i]->Divide(ctauNPixWMeaP[i],ctauNPixWMeaM[i]);
    
    ratioEtaPhiMu[i]->Divide(etaPhiMuP[i],etaPhiMuM[i]);
  }
  ratioYA->Divide(yPA,yMA);
  ratioCtauYA->Divide(ctauYPA,ctauYMA);
  ratioYPhi->Divide(yPhiP,yPhiM);

  for (unsigned int i=0; i<nRuns; i++) {
    runRatioPhi[i]->Divide(runPhiP[i],runPhiM[i]);
    runRatioPhiMu[i]->Divide(runPhiMuP[i],runPhiMuM[i]);
    runRatioNTrkHits[i]->Divide(runNTrkHitsP[i],runNTrkHitsM[i]);
  }

}

void LifetimeCheck::GetLimits(int t_nBins[], double t_xmin[], double t_xmax[], vector<double> *t_rapbins, vector<UInt_t> *t_runbins){
    t_nBins[0] = nBins;
      t_xmin[0] = xmin;
      t_xmax[0] = xmax;
    t_nBins[1] = nBinsErr;
      t_xmin[1] = xminErr;
      t_xmax[1] = xmaxErr;
    t_nBins[2] = nBinsEta;
      t_xmin[2] = xminEta;
      t_xmax[2] = xmaxEta;
    t_nBins[3] = nBinsPhi;
      t_xmin[3] = xminPhi;
      t_xmax[3] = xmaxPhi;
    t_nBins[4] = nBinsPhi;
      t_xmin[4] = xminPhi;
      t_xmax[4] = xmaxPhi;
    t_nBins[5] = nBinsPT;
      t_xmin[5] = xminPT;
      t_xmax[5] = xmaxPT;
    t_nBins[6] = nBinsDxy;
      t_xmin[6] = xminDxy;
      t_xmax[6] = xmaxDxy;
    t_nBins[7] = nBinsDz;
      t_xmin[7] = xminDz;
      t_xmax[7] = xmaxDz;
    t_nBins[8] = nBinsChi2Trk;
      t_xmin[8] = xminChi2Trk;
      t_xmax[8] = xmaxChi2Trk;
    t_nBins[9] = nBinsChi2Glb;
      t_xmin[9] = xminChi2Glb;
      t_xmax[9] = xmaxChi2Glb;
    t_nBins[10] = nBinsNTrkHits;
      t_xmin[10] = xminNTrkHits;
      t_xmax[10] = xmaxNTrkHits;
    t_nBins[11] = nBinsNPixWMea;
      t_xmin[11] = xminNPixWMea;
      t_xmax[11] = xmaxNPixWMea;
    t_nBins[12] = nBinsY;
      t_xmin[12] = xminY;
      t_xmax[12] = xmaxY;

    t_xmin[13] = mmin;
    t_xmax[13] = mmax;
    
    *t_rapbins = rapbins;
    *t_runbins = runbins;
}

void LifetimeCheck::SetHistStyle(int a, int b) {
  for (unsigned int i=0; i<rapbins.size(); i++) {
    ApplyHistStyle(ctauP[i],a,b);
    ApplyHistStyle(ctauErrP[i],a,b);
    ApplyHistStyle(etaP[i],a,b);
    ApplyHistStyle(phiP[i],a,b);
    ApplyHistStyle(phiMuP[i],a,b);
    ApplyHistStyle(pTP[i],a,b);
    ApplyHistStyle(dxyP[i],a,b);
    ApplyHistStyle(dzP[i],a,b);
    ApplyHistStyle(chi2TrkP[i],a,b);
    ApplyHistStyle(chi2GlbP[i],a,b);
    ApplyHistStyle(nTrkHitsP[i],a,b);
    ApplyHistStyle(nPixWMeaP[i],a,b);
    ApplyHistStyle(ctauM[i],a,b);
    ApplyHistStyle(ctauErrM[i],a,b);
    ApplyHistStyle(etaM[i],a,b);
    ApplyHistStyle(phiM[i],a,b);
    ApplyHistStyle(phiMuM[i],a,b);
    ApplyHistStyle(pTM[i],a,b);
    ApplyHistStyle(dxyM[i],a,b);
    ApplyHistStyle(dzM[i],a,b);
    ApplyHistStyle(chi2TrkM[i],a,b);
    ApplyHistStyle(chi2GlbM[i],a,b);
    ApplyHistStyle(nTrkHitsM[i],a,b);
    ApplyHistStyle(nPixWMeaM[i],a,b);

    ApplyHistStyle(ratio[i],a,b);
    ApplyHistStyle(ratioErr[i],a,b);
    ApplyHistStyle(ratioEta[i],a,b);
    ApplyHistStyle(ratioPhi[i],a,b);
    ApplyHistStyle(ratioPhiMu[i],a,b);
    ApplyHistStyle(ratioPT[i],a,b);
    ApplyHistStyle(ratioDxy[i],a,b);
    ApplyHistStyle(ratioDz[i],a,b);
    ApplyHistStyle(ratioChi2Trk[i],a,b);
    ApplyHistStyle(ratioChi2Glb[i],a,b);
    ApplyHistStyle(ratioNTrkHits[i],a,b);
    ApplyHistStyle(ratioNPixWMea[i],a,b);

    ApplyHistStyle(ctauCtauErrP[i]);
    ApplyHistStyle(ctauCtauErrM[i]);
    ApplyHistStyle(ctauEtaP[i]);
    ApplyHistStyle(ctauEtaM[i]);
    ApplyHistStyle(ctauPhiP[i]);
    ApplyHistStyle(ctauPhiM[i]);
    ApplyHistStyle(ctauPhiMuP[i]);
    ApplyHistStyle(ctauPhiMuM[i]);
    ApplyHistStyle(ctauPTP[i]);
    ApplyHistStyle(ctauPTM[i]);
    ApplyHistStyle(ctauDxyP[i]);
    ApplyHistStyle(ctauDxyM[i]);
    ApplyHistStyle(ctauDzP[i]);
    ApplyHistStyle(ctauDzM[i]);
    ApplyHistStyle(ctauChi2TrkP[i]);
    ApplyHistStyle(ctauChi2TrkM[i]);
    ApplyHistStyle(ctauChi2GlbP[i]);
    ApplyHistStyle(ctauChi2GlbM[i]);
    ApplyHistStyle(ctauNTrkHitsP[i]);
    ApplyHistStyle(ctauNTrkHitsM[i]);
    ApplyHistStyle(ctauNPixWMeaP[i]);
    ApplyHistStyle(ctauNPixWMeaM[i]);
    
    ApplyHistStyle(etaPhiMuP[i]);
    ApplyHistStyle(etaPhiMuM[i]);

    ApplyHistStyle(ratioCtauCtauErr[i]);
    ApplyHistStyle(ratioCtauEta[i]);
    ApplyHistStyle(ratioCtauPhi[i]);
    ApplyHistStyle(ratioCtauPhiMu[i]);
    ApplyHistStyle(ratioCtauPT[i]);
    ApplyHistStyle(ratioCtauDxy[i]);
    ApplyHistStyle(ratioCtauDz[i]);
    ApplyHistStyle(ratioCtauChi2Trk[i]);
    ApplyHistStyle(ratioCtauChi2Glb[i]);
    ApplyHistStyle(ratioCtauNTrkHits[i]);
    ApplyHistStyle(ratioCtauNPixWMea[i]);
    
    ApplyHistStyle(ratioEtaPhiMu[i]);
  }

  ApplyHistStyle(yPA,a,b);
  ApplyHistStyle(yMA,a,b);
  ApplyHistStyle(ratioYA,a,b);
  ApplyHistStyle(ctauYPA);
  ApplyHistStyle(ctauYMA);
  ApplyHistStyle(ratioCtauYA);
  ApplyHistStyle(yPhiP);
  ApplyHistStyle(yPhiM);
  ApplyHistStyle(ratioYPhi);

  for (unsigned int i=0; i<nRuns; i++) {
    ApplyHistStyle(runRatioPhi[i],1,i+1);
    ApplyHistStyle(runRatioPhiMu[i],1,i+1);
    ApplyHistStyle(runRatioNTrkHits[i],1,i+1);
  }

}

void LifetimeCheck::LegendOn2DHists(TCanvas *canv2[], string strname[]) {
  TLatex *lat = new TLatex();
  lat->SetTextSize(0.065);
  
  double xleg[] = {0.400, 0.050};
  double yleg[13] = {
    0.250, 2.000, 1.500, 1.500,
    23.00, 0.070, 0.300, 3.000,
    16.00, 27.000, 5.100, 1.500,
    0.310
  };
  // 0th element(ctauY) isn't presented on yleg[] but should be presented here
  string legstr[] = {
    "#font[12]{l}_{J/#psi} vs J/#psi |y|",
    "#font[12]{l}_{J/#psi} vs #font[12]{l}_{J/#psi} error",
    "#font[12]{l}_{J/#psi} vs #mu #eta",
    "#font[12]{l}_{J/#psi} vs J/#psi #phi",
    "#font[12]{l}_{J/#psi} vs #mu #phi",
    "#font[12]{l}_{J/#psi} vs J/#psi p_{T}",
    "#font[12]{l}_{J/#psi} vs #mu D_{xy}",
    "#font[12]{l}_{J/#psi} vs #mu D_{z}",
    "#font[12]{l}_{J/#psi} vs Tracker #mu #chi^{2}/dof",
    "#font[12]{l}_{J/#psi} vs Global #mu #chi^{2}/dof",
    "#font[12]{l}_{J/#psi} vs # tracker hits",
    "#font[12]{l}_{J/#psi} vs # pixel layers",
    "#mu #eta vs #mu #phi",
    "J/#psi |y| vs J/#psi #phi"
  };

  for (unsigned int i=0; i<24; i++) {
    for (unsigned int j=0; j<rapbins.size()-1; j++) {
      canv2[i]->cd(j+1);
      char name[512];
      sprintf(name,"%.1f<|y|<%.1f",rapbins[j],rapbins[j+1]);
      lat->DrawLatex(xleg[0],yleg[i/2],name);
    } // end: rapbins loop
  }
 
 // only for ctauY 
 for (unsigned int j=0; j<3; j++) {
  char name[512];
  if (j==0) sprintf(name,"0.0<y<2.4");
  else if (j==1) sprintf(name,"-2.4<y<0.0");
  else if (j==2) sprintf(name,"#left(#frac{Count y>0}{Count y<0}#right)");
  
  canv2[24]->cd(j+1);
  lat->DrawLatex(xleg[1],2.000,name);
 }
 
 // all other integrated plots and ratio plot
  for (unsigned int i=25; i<38; i++) {
    for (unsigned int j=0; j<3; j++) {
      char name[512];
      if (j==0) sprintf(name,"0.0<y<2.4");
      else if (j==1) sprintf(name,"-2.4<y<0.0");
      else if (j==2) sprintf(name,"#left(#frac{Count y>0}{Count y<0}#right)");
      
      canv2[i]->cd(j+1);
      lat->DrawLatex(xleg[1],yleg[i-25],name);
    }
  }  // end: 2 last canvases loop

  TLatex *lat2 = new TLatex();
  lat2->SetNDC();
  
  for (unsigned int i=24; i<38; i++) {
    canv2[i]->cd(4);
    lat2->SetTextColor(kBlack);
    lat2->SetTextSize(0.120);
    lat2->DrawLatex(xleg[1],0.82,strname[0].c_str());
    lat2->DrawLatex(xleg[1],0.70,strname[1].c_str());
    lat2->SetTextSize(0.090);
    lat2->DrawLatex(xleg[1],0.57,strname[2].c_str());
    lat2->DrawLatex(xleg[1],0.44,strname[3].c_str());
    
    lat2->SetTextSize(0.100);
    lat2->SetTextColor(kBlue);
    lat2->DrawLatex(xleg[1],yleg[12],legstr[i-24].c_str()); 
  }
}

void LifetimeCheck::Draw2DHists(TCanvas *canv[]) {
  TExec *palette1 = new TExec("palette1","gStyle->SetPalette(1)");
//  TExec *palette54 = new TExec("palette54","gStyle->SetPalette(53)");

  TH2D **hArr[24]={
    ctauCtauErrP, ctauCtauErrM, ctauEtaP, ctauEtaM,
    ctauPhiP, ctauPhiM, ctauPhiMuP, ctauPhiMuM,
    ctauPTP, ctauPTM, ctauDxyP, ctauDxyM,
    ctauDzP, ctauDzM, ctauChi2TrkP, ctauChi2TrkM,
    ctauChi2GlbP, ctauChi2GlbM, ctauNTrkHitsP, ctauNTrkHitsM,
    ctauNPixWMeaP, ctauNPixWMeaM, etaPhiMuP, etaPhiMuM
  };

  for (unsigned int i=0; i<rapbins.size()-1; i++) {
    for (unsigned int j=0; j<24; j++) {
      canv[j]->cd(i+1);
      Draw2DHistsComm1(hArr[j][i], palette1);
    }
  }
  canv[24]->cd(1);
  Draw2DHistsComm2(ctauYPA, palette1);
  canv[24]->cd(2);
  Draw2DHistsComm2(ctauYMA, palette1);
  canv[24]->cd(3);
  Draw2DHistsComm2(ratioCtauYA, palette1, true);

  TH2D **hArrR[12]={
    ratioCtauCtauErr, ratioCtauEta,
    ratioCtauPhi, ratioCtauPhiMu,
    ratioCtauPT, ratioCtauDxy,
    ratioCtauDz, ratioCtauChi2Trk,
    ratioCtauChi2Glb, ratioCtauNTrkHits,
    ratioCtauNPixWMea, ratioEtaPhiMu,
  };

  unsigned int ind=rapbins.size()-1;
  for (unsigned int i=0; i<12; i++) {
    for (unsigned int j=0; j<2; j++) {
      canv[25+i]->cd(j+1);
      Draw2DHistsComm2(hArr[i*2+j][ind], palette1);
    }
    canv[25+i]->cd(3);
    Draw2DHistsComm2(hArrR[i][ind], palette1, true);
  }

  canv[37]->cd(1);
  Draw2DHistsComm2(yPhiP, palette1);
  canv[37]->cd(2);
  Draw2DHistsComm2(yPhiM, palette1);
  canv[37]->cd(3);
  Draw2DHistsComm2(ratioYPhi, palette1, true);

}

void LifetimeCheck::Draw1DHists(TCanvas *canv[], bool drawAxis) {
  TH1D **ratioH[12] = {
    ratio, ratioErr, ratioEta, ratioPhi,
    ratioPhiMu, ratioPT, ratioDxy, ratioDz,
    ratioChi2Trk, ratioChi2Glb, ratioNTrkHits, ratioNPixWMea
  };
  TH1D **numeH[12] = {
    ctauP, ctauErrP, etaP, phiP,
    phiMuP, pTP, dxyP, dzP,
    chi2TrkP, chi2GlbP, nTrkHitsP, nPixWMeaP
  };
  TH1D **denoH[12] = {
    ctauM, ctauErrM, etaM, phiM,
    phiMuM, pTM, dxyM, dzM,
    chi2TrkM, chi2GlbM, nTrkHitsM, nPixWMeaM
  };
  const unsigned int nCanv = 12;

  // Ghost graph is created for axis
  TH1D *ratioG[12];
  for (unsigned int i=0; i<nCanv; i++) {
    ratioG[i] = (TH1D*)ratioH[i][0]->Clone();
    Draw1DHistsComm1(ratioG[i]);
  }
 
  // Ratio plots for detailed rapidity bins are drawn
  for (unsigned int j=0; j<nCanv; j++) {
    for (unsigned int i=0; i<rapbins.size()-1; i++) {
      canv[j]->cd(i+1);
      if (drawAxis) ratioG[j]->Draw();
      ratioH[j][i]->Draw("pe, same");
    }
  }

  // Integrated over all rap ranges, ratio plots
  canv[12]->cd(1);
  gPad->SetLogy(true);
  if (drawAxis) yPA->DrawNormalized("pe");
  else yPA->DrawNormalized("pe, same");
  canv[12]->cd(2);
  gPad->SetLogy(true);
  if (drawAxis) yMA->DrawNormalized("pe");
  else yMA->DrawNormalized("pe, same");
  canv[12]->cd(3);
  gPad->SetLogy(false);
  ratioYA->GetYaxis()->SetRangeUser(0,2);
  if (drawAxis) ratioYA->Draw();
  else ratioYA->Draw("pe, same");

  for (unsigned int j=0; j<nCanv; j++) {
    canv[13+j]->cd(1);
    gPad->SetLogy(true);
    if (drawAxis) numeH[j][rapbins.size()-1]->DrawNormalized("pe");
    numeH[j][rapbins.size()-1]->DrawNormalized("pe, same");
    canv[13+j]->cd(2);
    gPad->SetLogy(true);
    if (drawAxis) denoH[j][rapbins.size()-1]->DrawNormalized("pe");
    denoH[j][rapbins.size()-1]->DrawNormalized("pe, same");
    canv[13+j]->cd(3);
    gPad->SetLogy(false);
    if (drawAxis) ratioG[j]->Draw();
    ratioH[j][rapbins.size()-1]->Draw("pe, same");
  }
  gPad->SetLogy(false);
}

void LifetimeCheck::DrawRunCompHists(TCanvas *canv[]) {
  ApplyHistStyle(ratioPhi[rapbins.size()-1],0,0);
  ratioPhi[rapbins.size()-1]->SetMarkerStyle(kOpenCircle);
  ratioPhi[rapbins.size()-1]->SetMarkerSize(2.5);
  ratioPhi[rapbins.size()-1]->GetYaxis()->SetRangeUser(0.5,1.5);
  ApplyHistStyle(ratioPhiMu[rapbins.size()-1],0,0);
  ratioPhiMu[rapbins.size()-1]->SetMarkerStyle(kOpenCircle);
  ratioPhiMu[rapbins.size()-1]->SetMarkerSize(2.5);
  ratioPhiMu[rapbins.size()-1]->GetYaxis()->SetRangeUser(0.5,1.5);
  ApplyHistStyle(ratioNTrkHits[rapbins.size()-1],0,0);
  ratioNTrkHits[rapbins.size()-1]->SetMarkerStyle(kOpenCircle);
  ratioNTrkHits[rapbins.size()-1]->SetMarkerSize(2.5);
  ratioNTrkHits[rapbins.size()-1]->GetYaxis()->SetRangeUser(0.5,1.5);

  TLatex *lat = new TLatex();
  lat->SetTextSize(0.065);  lat->SetTextColor(kBlack);

  int color[] = {kBlack, kRed, kOrange, kGreen, kCyan, kBlue};
  for (unsigned int i=0; i<nRuns; i++) {
    canv[0]->cd(i+1);
    ratioPhi[rapbins.size()-1]->Draw("pe");
    runRatioPhi[i]->SetMarkerSize(1.7);
    runRatioPhi[i]->Draw("c,l,same");
    lat->DrawLatex(-1.400,1.350,Form("%u-%u",runbins[i],runbins[i+1]-1));
    canv[1]->cd(i+1);
    ratioPhiMu[rapbins.size()-1]->Draw("pe");
    runRatioPhiMu[i]->SetMarkerSize(1.7);
    runRatioPhiMu[i]->Draw("c,l,same");
    lat->DrawLatex(-1.400,1.350,Form("%u-%u",runbins[i],runbins[i+1]-1));
    canv[2]->cd(i+1);
    ratioNTrkHits[rapbins.size()-1]->Draw("pe");
    runRatioNTrkHits[i]->SetMarkerSize(1.7);
    runRatioNTrkHits[i]->Draw("l,c,same");
    lat->DrawLatex(8.000,1.350,Form("%u-%u",runbins[i],runbins[i+1]-1));
  }

}

void LifetimeCheck::WriteToFile(TFile *output) {
  output->cd();
  for (unsigned int i=0; i<rapbins.size(); i++) {
    ctauP[i]->Write(); ctauM[i]->Write(); ratio[i]->Write();
    ctauErrP[i]->Write(); ctauErrM[i]->Write(); ratioErr[i]->Write();
    etaP[i]->Write(); etaM[i]->Write(); ratioEta[i]->Write();
    phiP[i]->Write(); phiM[i]->Write(); ratioPhi[i]->Write();
    phiMuP[i]->Write(); phiMuM[i]->Write(); ratioPhiMu[i]->Write();
    pTP[i]->Write(); pTM[i]->Write(); ratioPT[i]->Write();
    dxyP[i]->Write(); dxyM[i]->Write(); ratioDxy[i]->Write();
    dzP[i]->Write(); dzM[i]->Write(); ratioDz[i]->Write();
    chi2TrkP[i]->Write(); chi2TrkM[i]->Write(); ratioChi2Trk[i]->Write();
    chi2GlbP[i]->Write(); chi2GlbM[i]->Write(); ratioChi2Glb[i]->Write();
    nTrkHitsP[i]->Write(); nTrkHitsM[i]->Write(); ratioNTrkHits[i]->Write();
    nPixWMeaP[i]->Write(); nPixWMeaM[i]->Write(); ratioNPixWMea[i]->Write();
    
    ctauCtauErrP[i]->Write(); ctauCtauErrM[i]->Write(); ratioCtauCtauErr[i]->Write();
    ctauEtaP[i]->Write(); ctauEtaM[i]->Write(); ratioCtauEta[i]->Write();
    ctauPhiP[i]->Write(); ctauPhiM[i]->Write(); ratioCtauPhi[i]->Write();
    ctauPhiMuP[i]->Write(); ctauPhiMuM[i]->Write(); ratioCtauPhiMu[i]->Write();
    ctauPTP[i]->Write(); ctauPTM[i]->Write(); ratioCtauPT[i]->Write();
    ctauDxyP[i]->Write(); ctauDxyM[i]->Write(); ratioCtauDxy[i]->Write();
    ctauDzP[i]->Write(); ctauDzM[i]->Write(); ratioCtauDz[i]->Write();
    ctauChi2TrkP[i]->Write(); ctauChi2TrkM[i]->Write(); ratioCtauChi2Trk[i]->Write();
    ctauChi2GlbP[i]->Write(); ctauChi2GlbM[i]->Write(); ratioCtauChi2Glb[i]->Write();
    ctauNTrkHitsP[i]->Write(); ctauNTrkHitsM[i]->Write(); ratioCtauNTrkHits[i]->Write();
    ctauNPixWMeaP[i]->Write(); ctauNPixWMeaM[i]->Write(); ratioCtauNPixWMea[i]->Write();
    
    etaPhiMuP[i]->Write(); etaPhiMuM[i]->Write(); ratioEtaPhiMu[i]->Write();
  }
  yPA->Write(); yMA->Write(); ratioYA->Write();
  ctauYPA->Write(); ctauYMA->Write(); ratioCtauYA->Write();
  yPhiP->Write(); yPhiM->Write(); ratioYPhi->Write();

  for (unsigned int i=0; i<nRuns; i++) {
    runRatioPhi[i]->Write();  runPhiP[i]->Write();  runPhiM[i]->Write();
    runRatioPhiMu[i]->Write();  runPhiMuP[i]->Write();  runPhiMuM[i]->Write();
    runRatioNTrkHits[i]->Write();  runNTrkHitsP[i]->Write();  runNTrkHitsM[i]->Write();
  }
}

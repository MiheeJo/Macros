#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>

#include <RooFit.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

using namespace std;
using namespace RooFit;

int main(void) {
  TFile finput("tnp_Ana_Z0_Signal_total_0120_Glb.root");
  if (finput.IsZombie()) return 1;

  finput.cd("MuonID/Glb_eta");
  TDirectory *dir = gDirectory;
  RooDataSet *ds = (RooDataSet*)dir->Get("fit_eff");
  if (ds != NULL) {
    const RooArgSet *argset = ds->get();
    argset->Print("v"); cout << endl;
    RooRealVar *eta = (RooRealVar*)argset->find("eta");
    RooRealVar *eff = (RooRealVar*)argset->find("efficiency");

    // Fill TGraphAsymmErrors
    const int nbins = eta->getBinning().numBins();
    const double *x = eta->getBinning().array();
    double ty[nbins], tyhi[nbins], tylo[nbins];

    for (int i=0; i<nbins; i++) {
      ds->get(i);
      ty[i] = eff->getVal();
      tyhi[i] = eff->getErrorHi();
      tylo[i] = eff->getErrorLo(); 
    }

    const double *y = ty; 
    const double *yhi = tyhi;
    const double *ylo = tylo;

    TGraphAsymmErrors *graph = new TGraphAsymmErrors(nbins,x,y,0,0,ylo,yhi);
    graph->Draw("apz");

    // Fill histogram
    TH1D *htmp = new TH1D("final histogram","",nbins,x);
    for (int i=0; i<ds->numEntries(); i++) {
      ds->get(i);
      htmp->SetBinContent(htmp->FindBin(eta->getVal()),eff->getVal());
    }

    /*TCanvas *c = new TCanvas("test","test",800,600);
      htmp->Draw();
      c->SaveAs("test.png");*/
  }
  else
    cout << "Error: read Dataset!" << endl;

  finput.Close();
  return 0;
} 

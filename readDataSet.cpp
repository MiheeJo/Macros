#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TH1D.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
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
    RooRealVar *x = (RooRealVar*)argset->find("eta");
    RooRealVar *y = (RooRealVar*)argset->find("efficiency");
    TH1D *htmp = new TH1D("final histogram","",10,-2.4,2.4);
    for (int i=0; i<ds->numEntries(); i++) {
      ds->get(i);
      htmp->SetBinContent(htmp->FindBin(x->getVal()),y->getVal());
    }

    TCanvas *c = new TCanvas("test","test",800,600);
    htmp->Draw();
    c->SaveAs("test.png");
  }
  else
    cout << "Error: read Dataset!" << endl;

  finput.Close();
  return 0;
} 

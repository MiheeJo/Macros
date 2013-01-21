//root -b -l datasets/RAA/default_bit1/Data2011_cent0-100_dPhi0.000-1.571.root datasets_RegIt/default_bit1/Data2011_cent0-100_dPhi0.000-1.571.root

#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooPlot.h>

using namespace RooFit;

void drawYComp() {
  RooWorkspace *work = new RooWorkspace("workspace");

  RooDataSet *prompt = (RooDataSet*)_file0->Get("dataJpsi");
  RooDataSet *regit = (RooDataSet*)_file1->Get("dataJpsi");
  RooDataSet *Prompt = prompt->reduce("Jpsi_Pt>6.5 && Jpsi_Pt<30 && Jpsi_Y>-2.4 && Jpsi_Y<2.4 && Jpsi_Ct>-1 && Jpsi_Ct<2");
  RooDataSet *Regit = regit->reduce("Jpsi_Pt>6.5 && Jpsi_Pt<30 && Jpsi_Y>-2.4 && Jpsi_Y<2.4 && Jpsi_Ct>-1 && Jpsi_Ct<2");
  Prompt->SetName("prompt");  work->import(*Prompt);
  Regit->SetName("regit");   work->import(*Regit);

  RooDataSet *yPlusPrompt = Prompt->reduce("Jpsi_Y>0 && Jpsi_Y<2.4");
  RooDataSet *yMinusPrompt = Prompt->reduce("Jpsi_Y<0 && Jpsi_Y>-2.4");
  RooDataSet *yPlusRegit = Regit->reduce("Jpsi_Y>0 && Jpsi_Y<2.4");
  RooDataSet *yMinusRegit = Regit->reduce("Jpsi_Y<0 && Jpsi_Y>-2.4");

  TCanvas *c0 = new TCanvas("test","test",600,600);
  c0->SetLogy(1); c0->Draw();

  RooPlot *cframe = work->var("Jpsi_Ct")->frame();
  prompt->plotOn(cframe,DataError(RooAbsData::SumW2));
  regit->plotOn(cframe,DataError(RooAbsData::SumW2),MarkerColor(kRed));
  cframe->Draw();  cframe->SetMinimum(0.5);
  c0->SaveAs("prompt_regit_ctau.pdf");
  delete cframe; c0->Clear();

  cframe = work->var("Jpsi_Ct")->frame();
  yPlusPrompt->plotOn(cframe,DataError(RooAbsData::SumW2));
  yMinusPrompt->plotOn(cframe,DataError(RooAbsData::SumW2),MarkerColor(kRed));
  c0->cd(); cframe->Draw();  cframe->SetMinimum(0.5);
  c0->SaveAs("prompt_yPyM_ctau.pdf");
  delete cframe; c0->Clear();

  cframe = work->var("Jpsi_Ct")->frame();
  yPlusRegit->plotOn(cframe,DataError(RooAbsData::SumW2));
  yMinusRegit->plotOn(cframe,DataError(RooAbsData::SumW2),MarkerColor(kRed));
  c0->cd(); cframe->Draw();  cframe->SetMinimum(0.5);
  c0->SaveAs("regit_yPyM_ctau.pdf");
  delete cframe;
  
  const char title[3][] = {"y<0","y>0","ratio"};
  TH1D *hPromptReco[3] *hRegIt[3];
  for (int i=0; i<3; i++) {
    char str[512];
    sprintf(str,"PromptReco_%s",title[i]);
    hPromptReco[i] = new TH1D(str,";#font[12]{l}_{J/#psi}[mm];",60,-1,2);
    sprintf(str,"RegIt_%s",title[i]);
    hRegIt[i] = new TH1D(str,";#font[12]{l}_{J/#psi}[mm];",60,-1,2);
  }

}


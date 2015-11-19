#include <TSystem.h>
#include <TROOT.h>
#include <TTree.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <iostream>

using namespace std;

//MC : ZEEMM, Data : MinBiasUPC
int SuImpByDirectory(string sampleName="MinBiasUPC", string verName="7_6_0", string refName="7_6_0_pre7")
{
  gROOT->Macro("~/JpsiStyle.C");
  gStyle->SetOptStat(1);
  gSystem->mkdir(Form("./%s__%s_vs_%s",sampleName.c_str(),verName.c_str(),refName.c_str()), 1);
  
  //TFile *froot = new TFile("DQM_V0001_R000000001__RelValZEEMM_13_HI__CMSSW_7_6_0-76X_mcRun2_HeavyIon_v11-v1__DQMIO.root");
  TFile *froot = new TFile("DQM_V0001_R000182124__HIMinBiasUPC__CMSSW_7_6_0-76X_dataRun1_v10_RelVal_hi2011-v1__DQMIO.root");
  if (froot->IsZombie()) return 1;
  //froot->cd("DQMData/Run 1/Muons/Run summary/MuonRecoAnalyzer");
  froot->cd("DQMData/Run 182124/Muons/Run summary/MuonRecoAnalyzer");
  TDirectory *root_dir = gDirectory;
  TIter rootnextkey( root_dir->GetListOfKeys() );

  //TFile *frefer = new TFile("DQM_V0001_R000000001__RelValZEEMM_13_HI__CMSSW_7_6_0_pre7-76X_mcRun2_HeavyIon_v6-v1__DQMIO.root");
  TFile *frefer = new TFile("DQM_V0001_R000182124__HIMinBiasUPC__CMSSW_7_6_0_pre7-76X_dataRun1_v5_RelVal_hi2011-v1__DQMIO.root");
  if (frefer->IsZombie()) return 1;
  //frefer->cd("DQMData/Run 1/Muons/Run summary/MuonRecoAnalyzer");
  frefer->cd("DQMData/Run 182124/Muons/Run summary/MuonRecoAnalyzer");
  TDirectory *refer_dir = gDirectory;
  root_dir->cd();

  TKey *rootkey;
  TObject *referobj;
  //TCanvas *canv = new TCanvas("test","test",800,600);
  TCanvas *canv = new TCanvas("test","test",800,1000);
  TLine *lineRatio = new TLine();
  lineRatio->SetLineColor(kGray+2);
  lineRatio->SetLineWidth(1.3);

  const Int_t colors3[3] = {kRed, kYellow, kGreen};
  const Int_t colors2[2] = {kRed, kGreen};
  
  while ( (rootkey = (TKey*)rootnextkey()) )
  {
    TObject *rootobj = rootkey->ReadObj();
    if (rootobj->IsA()->InheritsFrom("TH2")) {
      continue;
    } else if (rootobj->IsA()->InheritsFrom("TH1")) {
      TH1 *h1 = (TH1*)rootobj;
      refer_dir->GetObject(rootobj->GetName(),referobj);
      if (referobj)
      {
        gStyle->SetPadRightMargin(0.05);
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
        pad1->Draw();             
        pad1->cd();               
        
        TH1 *hr1 = (TH1*)referobj;
        double max = 0;
        if (h1->GetEntries() > hr1->GetEntries())
          h1->Scale( hr1->GetEntries()/h1->GetEntries() );
        else if (h1->GetEntries() < hr1->GetEntries())
          hr1->Scale( h1->GetEntries()/hr1->GetEntries() );
        if (h1->GetMaximum() > hr1->GetMaximum())
          max = h1->GetMaximum();
        else
          max = hr1->GetMaximum();

        h1->SetLineColor(kRed);
        h1->SetMarkerColor(kRed);
        h1->SetMarkerStyle(kOpenCircle);
        hr1->SetLineColor(kBlack);
        h1->SetMaximum(max*1.1);
        h1->Draw("pe");
        hr1->Draw("sames");
        gPad->Update(); canv->Update();

        TPaveStats *s = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
        s->SetFillColorAlpha(kWhite,0);
        s->SetY1NDC(s->GetY2NDC()-0.15);
        s->SetLineColor(kRed);
        s->SetTextColor(kRed);

        TPaveStats *sr = (TPaveStats*)hr1->GetListOfFunctions()->FindObject("stats");
        sr->SetFillColorAlpha(kWhite,0);
//        sr->SetY2NDC(s->GetY1NDC()-0.01);
//        sr->SetY1NDC(sr->GetY2NDC() - (s->GetY2NDC()-s->GetY1NDC()) );
        sr->SetY1NDC(sr->GetY2NDC()-0.15);
        sr->SetX2NDC(s->GetX1NDC()-0.01);
        sr->SetX1NDC(sr->GetX2NDC() - (s->GetX2NDC()-s->GetX1NDC()) );

        TText txt1; txt1.SetTextSize(0.035);
        txt1.SetTextColor(kRed);
        txt1.DrawTextNDC(0.35,0.90,verName.c_str());
        TText txt2; txt2.SetTextSize(0.035);
        txt2.DrawTextNDC(0.35,0.85,refName.c_str());
        gPad->Update();
      
        // Do not draw the Y axis label on the upper plot and redraw a small
        // axis instead, in order to avoid the first label (0) to be clipped.
        //h1->GetYaxis()->SetLabelSize(0.);
        //TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
        //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        //axis->SetLabelSize(15);
        //axis->Draw();
         
        //// ratio plot in the bottom pad
        canv->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.2);
        pad2->Draw();
        pad2->cd();

        // Define the ratio plot
        TH1F *h3 = (TH1F*)h1->Clone("h3");
        h3->SetLineColor(kBlack);
        h3->SetMarkerColor(kBlack);
        h3->SetMinimum(0.5);
        h3->SetMaximum(1.5);
        h3->Sumw2();
        h3->SetStats(0);
        h3->Divide(hr1);
        h3->SetMarkerStyle(21);
        h3->Draw("p");
        h3->SetTitle("");
        double xmin = h3->GetBinLowEdge(1);
        double xmax = h3->GetBinLowEdge(h3->GetNbinsX()+1);
        lineRatio->DrawLine(xmin,1,xmax,1);

        // Y axis ratio plot settings
        //h3->GetYaxis()->SetTitle(Form("ratio [%s]/[%s]",verName,refName));
        h3->GetYaxis()->SetTitle(Form("[%s]/[%s]",verName.c_str(),refName.c_str()));
        h3->GetYaxis()->SetNdivisions(505);
        //h3->GetYaxis()->SetTitleSize(20);
        h3->GetYaxis()->SetTitleSize(18);
        h3->GetYaxis()->SetTitleFont(43);
        h3->GetYaxis()->SetTitleOffset(1.55);
//        h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
//        h3->GetYaxis()->SetLabelSize(15);
        // X axis ratio plot settings
        h3->GetXaxis()->SetTitleSize(20);
        h3->GetXaxis()->SetTitleFont(43);
        h3->GetXaxis()->SetTitleOffset(4.);
//        h3->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
        h3->SetLabelSize(0.08,"XYZ");
        
        canv->SaveAs(Form("./%s__%s_vs_%s/%s.png",sampleName.c_str(),verName.c_str(),refName.c_str(),h1->GetName()));
        canv->Clear();
        delete h1;
        delete hr1;
        delete h3;
      }
      else continue;
    }
  }
 
  froot->Close();
  frefer->Close();
  delete canv;

  return 0;
}


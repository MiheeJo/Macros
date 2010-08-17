#include <TSystem.h>
#include <TTree.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPave.h>
#include <TText.h>
#include <sstream>
#include <string.h>

int savehisto(void)
{
  gStyle->SetPalette(1,0);
  const char *path = ;
  
  TFile *froot = new TFile("C:\\root\\macros\\HYD_3_8_X.root");
  if (froot->IsZombie()) return 1;
  froot->cd("DQMData/Run 1/Muons/Run summary/TestSummary");
  TDirectory *root_dir = gDirectory;
  TIter rootnextkey( root_dir->GetListOfKeys() );

  TFile *frefer = new TFile("C:\\root\\macros\\142189.root");
  if (frefer->IsZombie()) return 1;
  frefer->cd("DQMData/Run 142189/Muons/Run summary/TestSummary");
  TDirectory *refer_dir = gDirectory;
  root_dir->cd();

  Bool_t TestSummary = kTRUE;
  TKey *rootkey;
  TObject *referobj;
  TCanvas *canv = new TCanvas("test","test",800,600);
  TCanvas *canv2 = new TCanvas("test2","test2",1300,600);
  double level3[3] = {0,0.5,1};
  double level2[2] = {0,1};
  const Int_t colors3[3] = {kRed, kYellow, kGreen};
  const Int_t colors2[2] = {kRed, kGreen};
  
  while ( (rootkey = (TKey*)rootnextkey()) )
  {
    TObject *rootobj = rootkey->ReadObj();
    gStyle->SetPalette(1);
    if (rootobj->IsA()->InheritsFrom("TH2"))
    {
      TH2 *h2 = (TH2*)rootobj;
      refer_dir->GetObject(rootobj->GetName(),referobj);
      if (referobj)
      {
        gStyle->SetPadRightMargin(0.12);

        TH2 *hr2 = (TH2*)referobj;
        if (!TestSummary)
        {
        if (h2->GetEntries() > hr2->GetEntries())
          h2->Scale( hr2->GetEntries()/h2->GetEntries() );
        else if (h2->GetEntries() < hr2->GetEntries())
          hr2->Scale( h2->GetEntries()/hr2->GetEntries() );
        }
        
        if (TestSummary)
        {
          //gStyle->SetOptStat(0);
          string str1 = h2->GetName();
          string str2 = "energySummaryMap";
          string str3 = "kinematicsSummaryMap";
          string str4 = "muonIdSummaryMap";
          string str5 = "residualsSummaryMap";
          string str6 = "chi2TestSummaryMap";
          string str7 = "KolmogorovTestSummaryMap";
          if (str1.compare(str6) == 0 || str1.compare(str7) == 0)
          {
            h2->SetMaximum(1);
            h2->SetMinimum(-1);
            hr2->SetMaximum(1);
            hr2->SetMinimum(-1);
            gStyle->SetPalette(1);
          }
          else if (str1.compare(str4) == 0)
          {
            h2->SetMaximum(1.5);
            h2->SetMinimum(-0.01);
            hr2->SetMaximum(1.5);
            hr2->SetMinimum(-0.01);
            h2->SetContour(3,level3);
            hr2->SetContour(3,level3);
            gStyle->SetPalette(3,colors3);
          }
          else if (str1.compare(str2) == 0 || str1.compare(str3) == 0 || str1.compare(str5) == 0)
          {
            h2->SetMaximum(2);
            h2->SetMinimum(-0.01);
            hr2->SetMaximum(2);
            hr2->SetMinimum(-0.01);
            h2->SetContour(2,level2);
            hr2->SetContour(2,level2);
            gStyle->SetPalette(2,colors2);
          }
        }

        canv2->Divide(2,1);        
        canv2->cd(1);
        h2->GetXaxis()->CenterTitle(1);
        h2->GetYaxis()->CenterTitle(1);
        h2->Draw("COLZ");
        TText txt1;
        txt1.DrawTextNDC(0.2,0.85,"hydjet MC");
        canv2->cd(2);
        hr2->GetXaxis()->CenterTitle(1);
        hr2->GetYaxis()->CenterTitle(1);
        hr2->Draw("COLZ");
        TText txt2;
        txt2.DrawTextNDC(0.2,0.85,"run 142189");

        stringstream stmp;
        stmp.str("");
        stmp << h2->GetName() << ".png";
        canv2->SaveAs(stmp.str().c_str());
        canv2->Clear();
        delete h2;
        delete hr2;
      }
      else continue;
    }
    else if (rootobj->IsA()->InheritsFrom("TH1"))
    {
      TH1 *h1 = (TH1*)rootobj;
      refer_dir->GetObject(rootobj->GetName(),referobj);
      if (referobj)
      {
        gStyle->SetPadRightMargin(0.05);
        
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

        hr1->SetLineColor(kRed);
        canv->cd();
        h1->SetMaximum(max*1.1);
        h1->GetXaxis()->CenterTitle(1);
        h1->GetYaxis()->CenterTitle(1);
        h1->Draw();
        hr1->GetXaxis()->CenterTitle(1);
        hr1->GetYaxis()->CenterTitle(1);
        hr1->Draw("sames");
        gPad->Update();
        TPaveStats *s = (TPaveStats*)h1->FindObject("stats");
        s->SetY1NDC(s->GetY2NDC()-0.1);
        TPaveStats *sr = (TPaveStats*)hr1->FindObject("stats");
        sr->SetY2NDC(s->GetY1NDC()-0.03);
        sr->SetY1NDC(sr->GetY2NDC() - (s->GetY2NDC()-s->GetY1NDC()) );
        sr->SetLineColor(kRed);
        sr->SetTextColor(kRed);
        TText txt1;
        txt1.DrawTextNDC(0.2,0.85,"hydjet MC");
        TText txt2;
        txt2.SetTextColor(kRed);
        txt2.DrawTextNDC(0.2,0.75,"run 142189");
        gPad->Update();
        
        stringstream stmp;
        stmp.str("");
        stmp << h1->GetName() << ".png";
        canv->SaveAs(stmp.str().c_str());
        canv->Clear();
        delete h1;
        delete hr1;
      }
      else continue;
    }
  }
 
  froot->Close();
  frefer->Close();
  delete canv;
  delete canv2;

  return 0;
}
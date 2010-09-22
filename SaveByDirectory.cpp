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

int SaveByDirectory(void)
{
  gStyle->SetPalette(1,0);
  
  TFile *froot = new TFile("DQM_V0001_R000144112__workflow__for__test.root");
  if (froot->IsZombie()) return 1;
  froot->cd("DQMData/Run 144112/Physics/Run summary/BPhysics");
  TDirectory *root_dir = gDirectory;
  TIter rootnextkey( root_dir->GetListOfKeys() );
  root_dir->cd();

  TKey *rootkey;
  TCanvas *canv = new TCanvas("test","test",800,600);
  TCanvas *canv2 = new TCanvas("test2","test2",800,600);
  
  while ( (rootkey = (TKey*)rootnextkey()) )
  {
    TObject *rootobj = rootkey->ReadObj();
    gStyle->SetPalette(1);
    if (rootobj->IsA()->InheritsFrom("TH2"))
    {
      TH2 *h2 = (TH2*)rootobj;
  //    gStyle->SetPadRightMargin(0.12);
      canv2->cd();
      h2->Draw("COLZ");
      stringstream stmp;
      stmp.str("");
      stmp << h2->GetName() << ".png";
      canv2->SaveAs(stmp.str().c_str());
      canv2->Clear();
      delete h2;
    }
    else if (rootobj->IsA()->InheritsFrom("TH1"))
    {
      TH1 *h1 = (TH1*)rootobj;
      canv->cd();
      h1->Draw();
      stringstream stmp;
      stmp.str("");
      stmp << h1->GetName() << ".png";
      canv->SaveAs(stmp.str().c_str());
      canv->Clear();
      delete h1;
    }
  }
 
  froot->Close();
  delete canv;
  delete canv2;

  return 0;
}

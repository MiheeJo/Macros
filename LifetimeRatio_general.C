#include <LifetimeRatio_general.h>

int LifetimeRatio_general_sailor() {
  gROOT->Macro("./rootlogon.C");

  bool fillCowboy = false;
  TFile *output;
  if (fillCowboy) output = new TFile("ratioHistos_cowboy.root","recreate");
  else output = new TFile("ratioHistos_sailor.root","recreate");

  TFile *input1 =  TFile::Open("root://eoscms//eos/cms/store/user/miheejo/TREE/2011PbPb/RegIt/Jpsi_Histos_RegIt_dca_small.root");
  TFile *input11 =  TFile::Open("root://eoscms//eos/cms/store/user/miheejo/TREE/2011PbPb/RegIt/Jpsi_Histos_RegIt_trkArb_small.root");
  TFile *input2 =  TFile::Open("root://eoscms//eos/cms/store/user/miheejo/TREE/2011PbPb/PromptReco/mini_Jpsi_Histos_may202012_m25gev.root");
  TFile *input3 =  TFile::Open("root://eoscms//eos/cms/store/user/dmoon/cms442/JpsiV2/OniaTTree/Merge/TrkTrk_NoTrkInfo/OniaTree_Jpsi_TrkTrk_Mar_v1_total_20130327.root");

  LifetimeCheck regit((TTree*)input1->Get("myTree"),fillCowboy);
  LifetimeCheck regitNC((TTree*)input11->Get("myTree"),fillCowboy);
  LifetimeCheck preco((TTree*)input2->Get("myTree"),fillCowboy);
//  LifetimeCheck trktrk((TTree*)input3->Get("myTree"),fillCowboy);

  cout << "Create histograms" << endl;
  regit.CreateHistos("regit");
  regitNC.CreateHistos("regitNoCuts");
  preco.CreateHistos("preco");
//  trktrk.CreateHistos("trktrk");

  cout << "Read TTrees and fill histograms up" << endl;
  cout << "  regit TTree" << endl;
  regit.ReadTree();
  cout << "  regitNC TTree" << endl;
  regitNC.ReadTree();
  cout << "  preco TTree" << endl;
  preco.ReadTree();
  cout << "  trktrk TTree" << endl;
//  trktrk.ReadTree();
    
  cout << "Write histograms to TFile" << endl;
  regit.WriteToFile(output);
  regitNC.WriteToFile(output);
  preco.WriteToFile(output);
//  trktrk.WriteToFile(output);
  
  output->Close();
  input1->Close();
  input2->Close();
  input3->Close();

 
  return 0;
}

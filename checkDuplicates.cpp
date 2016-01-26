#include <iostream>
#include <map>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

int main(void) {

  string inFileNames[] = { 
    "/home/mihee/cms/oniaTree/2011PbPb/All_Histos_cmssw445p1_RegIt_EvtPlane_small.root",
    "/home/mihee/cms/oniaTree/2011PbPb/All_Histos_cmssw445p1_RegIt_v3_small.root"
  };  
  int nFiles = sizeof(inFileNames)/sizeof(string);

  TChain *chain = new TChain("myTree");
  for (int i=0; i<nFiles; i++) {
    chain->Add(inFileNames[i].c_str());
  }

  int eventNb,runNb;
  map< pair<int, int>, int> eList;
  map< pair<int, int>, int>::iterator it;

  chain->SetBranchAddress("eventNb",&eventNb);
  chain->SetBranchAddress("runNb",&runNb);

  for (int ev=0; ev<chain->GetEntries(); ev++) {
    if (ev%10000 == 0) cout << "Event: " << ev << " / " << chain->GetEntries() << endl;
    chain->GetEntry(ev);
    pair<int, int> element = make_pair(eventNb, runNb);
    pair<pair<int, int>, int> element2 = make_pair(element, ev);
    pair<map<pair<int, int>, int>::iterator, bool> result = eList.insert(element2);
    if (result.second == 0) cout << " Duplicated event: " << runNb << " " << eventNb << endl;
//    cout << " Duplicated event: " << runNb << " " << eventNb << endl;
  }

  cout << "Total map size: " << eList.size() << " | Total tree event size: " << chain->GetEntries()<< endl;

  delete chain;
  return 0;
}

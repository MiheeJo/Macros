void drawCtauTrue() {
  TTree *tree = new TTree();
  tree->AddFriend("NPMC_regit = myTree","/afs/cern.ch/work/m/miheejo/private/TREE/RegIt_NPMC_MCTemplate_445p1.root");
  tree->AddFriend("NPMC_prompt = myTree","/afs/cern.ch/work/m/miheejo/private/TREE/PromptReco_NPMC_MCTemplate_445p1.root");
  tree->AddFriend("PRMC_regit = myTree","/afs/cern.ch/work/m/miheejo/private/TREE/RegIt_PRMC_MCTemplate_445p1.root");
  tree->AddFriend("PRMC_prompt = myTree","/afs/cern.ch/work/m/miheejo/private/TREE/PromptReco_PRMC_MCTemplate_445p1.root");

  TH1D *htmp[4];
  htmp[0] = new TH1D("NPMC_regit_ctauTrue","lifetime distribution;#font[12]{l}_{J/#psi}^{True}[mm];",100,-0.01,2);
  htmp[1] = new TH1D("NPMC_prompt_ctauTrue","lifetime distribution;#font[12]{l}_{J/#psi}^{True}[mm];",100,-0.01,2);
  htmp[2] = new TH1D("PRMC_regit_ctauTrue","lifetime distribution;#font[12]{l}_{J/#psi}^{True}[mm];",100,-0.01,2);
  htmp[3] = new TH1D("PRMC_prompt_ctauTrue","lifetime distribution;#font[12]{l}_{J/#psi}^{True}[mm];",100,-0.01,2);

  tree->Draw("NPMC_regit.Reco_QQ_ctauTrue>>NPMC_regit_ctauTrue","Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && Reco_QQ_4mom.Rapidity()>-2.4 && Reco_QQ_4mom.Rapidity()<2.4 && Reco_QQ_4mom.Pt()>6.5 && Reco_QQ_sign==0 && Reco_QQ_4mom.Pt()<30 && (HLTriggers&1)==1 && (Reco_QQ_trig&1)==1");
  tree->Draw("NPMC_prompt.Reco_QQ_ctauTrue>>NPMC_prompt_ctauTrue","Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && Reco_QQ_4mom.Rapidity()>-2.4 && Reco_QQ_4mom.Rapidity()<2.4 && Reco_QQ_4mom.Pt()>6.5 && Reco_QQ_sign==0 && Reco_QQ_4mom.Pt()<30 && (HLTriggers&1)==1 && (Reco_QQ_trig&1)==1");
  tree->Draw("PRMC_regit.Reco_QQ_ctauTrue>>PRMC_regit_ctauTrue","Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && Reco_QQ_4mom.Rapidity()>-2.4 && Reco_QQ_4mom.Rapidity()<2.4 && Reco_QQ_4mom.Pt()>6.5 && Reco_QQ_sign==0 && Reco_QQ_4mom.Pt()<30 && (HLTriggers&1)==1 && (Reco_QQ_trig&1)==1");
  tree->Draw("PRMC_prompt.Reco_QQ_ctauTrue>>PRMC_prompt_ctauTrue","Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && Reco_QQ_4mom.Rapidity()>-2.4 && Reco_QQ_4mom.Rapidity()<2.4 && Reco_QQ_4mom.Pt()>6.5 && Reco_QQ_sign==0 && Reco_QQ_4mom.Pt()<30 && (HLTriggers&1)==1 && (Reco_QQ_trig&1)==1");

  htmp[0]->SetMarkerColor(kRed);
  htmp[0]->SetLineColor(kRed);
  htmp[2]->SetMarkerColor(kRed);
  htmp[2]->SetLineColor(kRed);

  TLegend *leg = new TLegend(0.7,0.7,0.9,0.8);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  for (int i=0; i<4; i++) {
    htmp[i]->SetMarkerSize(1.3);
  }
  leg->AddEntry(htmp[0],"regit","lp");
  leg->AddEntry(htmp[1],"prompt reco","lp");

  TCanvas canv("test","test",600,600);
  canv.SetLogy();  canv.Draw();
  htmp[0]->Draw("pe"); htmp[1]->Draw("same, pe");
  leg->Draw();
  canv.SaveAs("NPMC_ctauTrue.png");
  canv.SaveAs("NPMC_ctauTrue.pdf");
  canv.Clear();
  htmp[2]->Draw("pe"); htmp[3]->Draw("same, pe");
  leg->Draw();
  canv.SaveAs("PRMC_ctauTrue.png");
  canv.SaveAs("PRMC_ctauTrue.pdf");
  canv.Clear();
  htmp[0]->Scale(1.0/htmp[0]->GetEntries());
  htmp[1]->Scale(1.0/htmp[1]->GetEntries());
  htmp[0]->Draw("l"); htmp[1]->Draw("same, l");
  leg->Draw();
  canv.SaveAs("NPMC_ctauTrue_notNorm.png");
  canv.SaveAs("NPMC_ctauTrue_notNorm.pdf");
}

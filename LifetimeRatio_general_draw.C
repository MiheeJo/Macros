#include <LifetimeRatio_general.h>

using namespace std;


int LifetimeRatio_general_draw() {
  gROOT->Macro("./rootlogon.C");

  cout << "Reading input file" << endl;
  TFile *input;
  bool fillCowboy = true;
  char cowboystr[128];
  if (fillCowboy) {
    input = new TFile("ratioHistos_cowboy.root","read");
    sprintf(cowboystr,"cowboy");
  } else {
    input = new TFile("ratioHistos_sailor.root","read");
    sprintf(cowboystr,"sailor");
  }

  string details[2] = {"2.6<M_{#mu#mu}<3.5 GeV/c^{2}","6.5<p_{T}<30 GeV/c"};
  string filename[] = {"ctau", "ctauErr", "eta", "phi", "phiMu", "pT", "Dxy", "Dz", "chi2Trk", "chi2Glb", "nTrkHits", "nPixWMea", "y"};
  const unsigned int size1d = sizeof(filename)/sizeof(string)*2-1;
  string filename2D[] = {
    "ctauErrP", "ctauErrM", "ctauEtaP", "ctauEtaM",
    "ctauPhiP", "ctauPhiM", "ctauPhiMuP", "ctauPhiMuM", 
    "ctauPTP", "ctauPTM", "ctauDxyP", "ctauDxyM",
    "ctauDzP", "ctauDzM", "ctauChi2TrkP", "ctauChi2TrkM",
    "ctauChi2GlbP", "ctauChi2GlbM", "ctauNTrkHitsP", "ctauNTrkHitsM",
    "ctauNPixWMeaP", "ctauNPixWMeaM", "etaPhiMuP", "etaPhiMuM", 
    "ctauY",
    "ctauCtauErr", "ctauEta", "ctauPhi", "ctauPhiMu",
    "ctauPT", "ctauDxy", "ctauDz", "ctauChi2Trk",
    "ctauChi2Glb", "ctauNTrkHits", "ctauNPixWMea", "etaPhiMu",
    "yPhi"
  };
  const unsigned int size2d = sizeof(filename2D)/sizeof(string);

  LifetimeCheck regit(input, "regit", fillCowboy);
  LifetimeCheck regitNC(input, "regitNoCuts", fillCowboy);
  LifetimeCheck preco(input, "preco", fillCowboy);
//  LifetimeCheck trktrk(input, "trktrk");

  int nBins[30];
  double xmin[30], xmax[30];
  vector<double> rapbins;
  vector<UInt_t> runbins;

  regit.GetLimits(nBins, xmin, xmax, &rapbins, &runbins);

  cout << "Set styles on histograms" << endl;
  regit.SetHistStyle(1,1);
  regitNC.SetHistStyle(5,5);
  preco.SetHistStyle(3,3);
//  trktrk.SetHistStyle(5,5);
    
  TCanvas *canv[size1d];
  for (unsigned int i=0; i<size1d; i++) {
    if (i<12) {
      canv[i] = new TCanvas(Form("canv_%d",i),"canv",1800,1200);
      canv[i]->Divide(3,2,0.00000001,0.00000001,0);
    } else{
      canv[i] = new TCanvas(Form("canv_%d",i),"canv",1200,1200);
      canv[i]->Divide(2,2,0.00000001,0.00000001,0);
    }
    canv[i]->Draw();
  }

  cout << "Draw 1D histograms" << endl;
  regit.Draw1DHists(canv,true);
  regitNC.Draw1DHists(canv,false);
  preco.Draw1DHists(canv,false);
//  trktrk.Draw1DHists(canv,false);
  
  TLatex *lat = new TLatex();
  lat->SetTextSize(0.065);  lat->SetTextColor(kBlack);
  
  TH1D *ratio1 = new TH1D("ratio1","",1,0,1); ApplyHistStyle(ratio1,1,1);
  TH1D *ratio2 = new TH1D("ratio2","",1,0,1); ApplyHistStyle(ratio2,3,3);
  TH1D *ratio3 = new TH1D("ratio3","",1,0,1); ApplyHistStyle(ratio3,5,5);

  TLegend *leg = new TLegend(0.20,0.79,0.6,0.95);
  ApplyLegendStyle(leg);
  leg->AddEntry(ratio1,"RegIt(Glb-Glb)","lp");
  leg->AddEntry(ratio3,"RegItNoCuts(Glb-Glb)","lp");
  leg->AddEntry(ratio2,"Pr. Reco(Glb-Glb)","lp");
//  leg->AddEntry(ratio3,"Pr. Reco(Trk-Trk)","lp");
  
  TLine *gline = new TLine();
  gline->SetLineWidth(1.2);
  gline->SetLineColor(kGray+3);

  double xlegDet[12] = {
    0.700, 0.165, 1.350, 0.400, 0.400, 17.000, 0.055, 0.210, 2.250, 10.150, 18.000, 3.500
  };
  double xlegInt[12] = {
    -0.300, 0.055, 0.300, -2.100, -2.100, 5.000, 0.020, 0.070, 0.600, 2.200, 5.000, 1.300
  };
  double yleg[3] = {1.300, 0.350, 0.150};

  for (unsigned int idx=0; idx<size1d; idx++) {
    for (unsigned int i=0; i<rapbins.size()-1; i++) {
      if (idx<12) canv[idx]->cd(i+1);
      else {
        canv[idx]->cd(i+1);
        if (i+1>4) break;
      }
      
      if (idx<12) {
        char name[512];
        sprintf(name,"%.1f<|y|<%.1f",rapbins[i],rapbins[i+1]);
        
        lat->SetTextSize(0.065);  lat->SetTextColor(kBlack);  lat->SetNDC(false);
        gline->DrawLine(xmin[idx],1,xmax[idx],1);
        lat->DrawLatex(xlegDet[idx],yleg[0],name);
        if (i==5) {
          lat->DrawLatex(xlegInt[idx],0.7,details[0].c_str());
          lat->DrawLatex(xlegInt[idx],0.6,details[1].c_str());
        }
      } else if (idx>12) {
        lat->SetTextColor(kBlack);  lat->SetNDC(true);
        if (i==0) {
          lat->SetTextSize(0.080);
          lat->DrawLatex(0.4,0.250,"0.0<y<2.4");
        } else if (i==1) {
          lat->SetTextSize(0.080);
          lat->DrawLatex(0.4,0.250,"-2.4<y<0.0");
        } else if (i==2) {
          gline->DrawLine(xmin[idx-13],1,xmax[idx-13],1);
        } else if (i==3) {
          lat->SetTextSize(0.100);
          lat->DrawLatex(0.1,yleg[1],details[0].c_str());
          lat->DrawLatex(0.1,yleg[2],details[1].c_str());
        }
      } else if (idx==12) { // all y
        lat->SetTextColor(kBlack);  lat->SetNDC(true);
        if (i==0) {
          lat->SetTextSize(0.080);
          lat->DrawLatex(0.4,0.250,"0.0<y<2.4");
        } else if (i==1) {
          lat->SetTextSize(0.080);
          lat->DrawLatex(0.4,0.250,"-2.4<y<0.0");
        } else if (i==2) {
          gline->DrawLine(xmin[12],1,xmax[12],1);
        } else if (i==3) {
          lat->SetTextSize(0.100);
          lat->DrawLatex(0.100,yleg[1],details[0].c_str());
          lat->DrawLatex(0.100,yleg[2],details[1].c_str());
        }
      }

      if (i==0&&idx<12) leg->Draw();
      if (idx>=12&&i==3) {
        ratio1->SetMarkerSize(2.2);
        ratio2->SetMarkerSize(2.9);
        ratio3->SetMarkerSize(3.4);
        leg->SetX1NDC(0.20);
        leg->SetY1NDC(0.65);
        leg->SetX2NDC(0.85);
        leg->SetY2NDC(0.95);
        leg->Draw();
      }
    } // end: rapidity loop

    canv[idx]->Update();
    if (idx < 12)
      canv[idx]->SaveAs(Form("LifetimeHistos_1D_%s_%s.pdf",cowboystr,filename[idx].c_str()));
    else if (idx == 12)
      canv[idx]->SaveAs(Form("LifetimeHistos_1DAll_%s_%s.pdf",cowboystr,filename[idx].c_str()));
    else if (idx > 12)
      canv[idx]->SaveAs(Form("LifetimeHistos_1DAll_%s_%s.pdf",cowboystr,filename[idx-13].c_str()));
  } // end: canvas loop

  for (unsigned int i=0; i<size1d; i++) delete canv[i];
  

  cout << "Draw 2D histograms" << endl;
  TCanvas *canv2[size2d];
  cout << "size2D " << size2d << endl;
  for (unsigned int i=0; i<size2d; i++) {
    if (i>23) {
      canv2[i] = new TCanvas(Form("canv2_%d",i),"canv2",1200,1200);
      canv2[i]->Divide(2,2,0.00000001,0.00000001,0);
    } else {
      canv2[i] = new TCanvas(Form("canv2_%d",i),"canv2",1800,1200);
      canv2[i]->Divide(3,2,0.00000001,0.00000001,0);
    }
    canv2[i]->Draw();
  }

  string strRegit[4]={"RegIt","Glb-Glb",details[0],details[1]};
  regit.Draw2DHists(canv2);
  regit.LegendOn2DHists(canv2, strRegit);
  
  for (unsigned int i=0; i<size2d; i++) {
    if (i>23) 
      canv2[i]->SaveAs(Form("LifetimeHistos_2DAll_regit_%s_%s.pdf",cowboystr,filename2D[i].c_str()));
    else
      canv2[i]->SaveAs(Form("LifetimeHistos_2D_regit_%s_%s.pdf",cowboystr,filename2D[i].c_str()));
    delete canv2[i];
  
    if (i>23) {
      canv2[i] = new TCanvas(Form("canv2_%d",i),"canv2",1200,1200);
      canv2[i]->Divide(2,2,0.00000001,0.00000001,0);
    } else {
      canv2[i] = new TCanvas(Form("canv2_%d",i),"canv2",1800,1200);
      canv2[i]->Divide(3,2,0.00000001,0.00000001,0);
    }
  }
  
  string strPrReco[4]={"Prompt Reco","Glb-Glb",details[0],details[1]};
  preco.Draw2DHists(canv2);
  preco.LegendOn2DHists(canv2, strPrReco);
  
  for (unsigned int i=0; i<size2d; i++) {
    if (i>23) 
      canv2[i]->SaveAs(Form("LifetimeHistos_2DAll_preco_%s_%s.pdf",cowboystr,filename2D[i].c_str()));
    else
      canv2[i]->SaveAs(Form("LifetimeHistos_2D_preco_%s_%s.pdf",cowboystr,filename2D[i].c_str()));
    delete canv2[i];
  
    if (i>23) {
      canv2[i] = new TCanvas(Form("canv2_%d",i),"canv2",1200,1200);
      canv2[i]->Divide(2,2,0.00000001,0.00000001,0);
    } else {
      canv2[i] = new TCanvas(Form("canv2_%d",i),"canv2",1800,1200);
      canv2[i]->Divide(3,2,0.00000001,0.00000001,0);
    }
    canv2[i]->Draw();
  }

  string strRegItNC[4]={"RegIt No MuId Cuts","Glb-Glb",details[0],details[1]};
  regitNC.Draw2DHists(canv2);
  regitNC.LegendOn2DHists(canv2, strRegItNC);

  for (unsigned int i=0; i<size2d; i++) {
    if (i>23) 
      canv2[i]->SaveAs(Form("LifetimeHistos_2DAll_regitNC_%s_%s.pdf",cowboystr,filename2D[i].c_str()));
    else 
      canv2[i]->SaveAs(Form("LifetimeHistos_2D_regitNC_%s_%s.pdf",cowboystr,filename2D[i].c_str()));
    delete canv2[i];
  }


/*  string strTrkTrk[4]={"Prompt Reco","Trk-Trk",details[0],details[1]};
  trktrk.Draw2DHists(canv2);
  trktrk.LegendOn2DHists(canv2, strTrkTrk);

  for (unsigned int i=0; i<size2d; i++) {
    if (i>23) 
      canv2[i]->SaveAs(Form("LifetimeHistos_2DAll_trktrk_%s_%s.pdf",cowboystr,filename2D[i].c_str()));
    else 
      canv2[i]->SaveAs(Form("LifetimeHistos_2D_trktrk_%s_%s.pdf",cowboystr,filename2D[i].c_str()));
    delete canv2[i];
  }
*/

  // run by run checking
  string filename1DRun[3] = {"phi", "phiMu", "nTrkHits"};
  cout << "Draw run by run comparison histograms" << endl;

  TH1D *runratio[runbins.size()-1];
  TLegend *leg2 = new TLegend(0.4,0.75,0.9,0.95);
  ApplyLegendStyle(leg2);
  for (unsigned int i=0; i<runbins.size()-1; i++) {
    runratio[i] = new TH1D(Form("runratio_%d",i),"",nBins[4],xmin[4],xmax[4]);
    ApplyHistStyle(runratio[i],1,i+1);
    leg2->AddEntry(runratio[i],Form("Run %u-%u",runbins[i],runbins[i+1]-1),"lp");
  }

  const unsigned int ncanv3 = 3;
  const unsigned int nruns = 5;
  TCanvas *canv3[ncanv3];
  for (unsigned int i=0; i<ncanv3; i++) {
    canv3[i] = new TCanvas(Form("canv2_%d",i),"canv2",1800,1200);
    canv3[i]->Divide(3,2,0.00000001,0.00000001,0);
    canv3[i]->Draw();
  }
  regit.DrawRunCompHists(canv3);
  for (unsigned int i=0; i<ncanv3; i++) {
    for (unsigned int j=0; j<nruns; j++) {
      canv3[i]->cd(j+1);
      if (i<2) {
        gline->DrawLine(xmin[4],1,xmax[4],1);
        leg2->SetX1NDC(0.4);
        leg2->SetY1NDC(0.75);
        leg2->SetX2NDC(0.9);
        leg2->SetY2NDC(0.95);
  //      leg2->Draw();
      } else {
        gline->DrawLine(xmin[10],1,xmax[10],1);
        leg2->SetX1NDC(0.2);
        leg2->SetY1NDC(0.75);
        leg2->SetX2NDC(0.7);
        leg2->SetY2NDC(0.95);
  //      leg2->Draw();
      }
    }
    canv3[i]->cd(nruns+1);
    lat->SetNDC(true);
    lat->DrawLatex(0.1,0.9,"RegIt(GlbGlb)");
    lat->DrawLatex(0.1,0.8,details[0].c_str());
    lat->DrawLatex(0.1,0.7,details[1].c_str());

    canv3[i]->SaveAs(Form("RunByRunComp_1D_regit_%s_%s.pdf",cowboystr,filename1DRun[i].c_str()));
    delete canv3[i];
  }

  for (unsigned int i=0; i<ncanv3; i++) {
    canv3[i] = new TCanvas(Form("canv2_%d",i),"canv2",1800,1200);
    canv3[i]->Divide(3,2,0.00000001,0.00000001,0);
    canv3[i]->Draw();
  }
  preco.DrawRunCompHists(canv3);
  for (unsigned int i=0; i<ncanv3; i++) {
    for (unsigned int j=0; j<nruns; j++) {
      canv3[i]->cd(j+1);
      if (i<2) {
        gline->DrawLine(xmin[4],1,xmax[4],1);
        leg2->SetX1NDC(0.4);
        leg2->SetY1NDC(0.75);
        leg2->SetX2NDC(0.9);
        leg2->SetY2NDC(0.95);
//        leg2->Draw();
      } else {
        gline->DrawLine(xmin[10],1,xmax[10],1);
        leg2->SetX1NDC(0.2);
        leg2->SetY1NDC(0.75);
        leg2->SetX2NDC(0.7);
        leg2->SetY2NDC(0.95);
//        leg2->Draw();
      }
    }
    canv3[i]->cd(nruns+1);
    lat->SetNDC(true);
    lat->DrawLatex(0.1,0.9,"Prompt Reco(GlbGlb)");
    lat->DrawLatex(0.1,0.8,details[0].c_str());
    lat->DrawLatex(0.1,0.7,details[1].c_str());

    canv3[i]->SaveAs(Form("RunByRunComp_1D_preco_%s_%s.pdf",cowboystr,filename1DRun[i].c_str()));
    delete canv3[i];
  }
  
  for (unsigned int i=0; i<ncanv3; i++) {
    canv3[i] = new TCanvas(Form("canv2_%d",i),"canv2",1800,1200);
    canv3[i]->Divide(3,2,0.00000001,0.00000001,0);
    canv3[i]->Draw();
  }
  regitNC.DrawRunCompHists(canv3);
//  trktrk.DrawRunCompHists(canv3);
  for (unsigned int i=0; i<ncanv3; i++) {
    for (unsigned int j=0; j<nruns; j++) {
      canv3[i]->cd(j+1);
      if (i<2) {
        gline->DrawLine(xmin[4],1,xmax[4],1);
        leg2->SetX1NDC(0.4);
        leg2->SetY1NDC(0.75);
        leg2->SetX2NDC(0.9);
        leg2->SetY2NDC(0.95);
//        leg2->Draw();
      } else {
        gline->DrawLine(xmin[10],1,xmax[10],1);
        leg2->SetX1NDC(0.2);
        leg2->SetY1NDC(0.75);
        leg2->SetX2NDC(0.7);
        leg2->SetY2NDC(0.95);
//        leg2->Draw();
      }
    }
    canv3[i]->cd(nruns+1);
    lat->SetNDC(true);
    lat->DrawLatex(0.1,0.9,"RegIt No MuId Cuts(GlbGlb)");
//    lat->DrawLatex(0.1,0.9,"Prompt Reco(TrkTrk)");
    lat->DrawLatex(0.1,0.8,details[0].c_str());
    lat->DrawLatex(0.1,0.7,details[1].c_str());

    canv3[i]->SaveAs(Form("RunByRunComp_1D_regitNC_%s_%s.pdf",cowboystr,filename1DRun[i].c_str()));
    delete canv3[i];
  }

  input->Close();
 
  return 0;
}

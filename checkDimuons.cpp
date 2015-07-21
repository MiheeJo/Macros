void checkDimuons() {

  using namespace RooFit;

  string cut = "TMath::Abs(Jpsi_Y)>1.6 && TMath::Abs(Jpsi_Y)<2.4 && Jpsi_Pt>3 && Jpsi_Pt<6.5 && Jpsi_Ct>-2 && Jpsi_Ct<4";

  RooWorkspace *wo = new RooWorkspace("wo");
  RooDataSet *dsLxyzW = (RooDataSet*)_file0->Get("dataJpsiWeight");
  RooDataSet *dsLxyz = (RooDataSet*)_file0->Get("dataJpsi");
  RooDataSet *dsLxyW = (RooDataSet*)_file1->Get("dataJpsiWeight");
  RooDataSet *dsLxy = (RooDataSet*)_file1->Get("dataJpsi");

  RooDataSet *LxyzW = dsLxyzW->reduce(cut.c_str()); LxyzW->SetName("LxyzW");
  RooDataSet *Lxyz = dsLxyz->reduce(cut.c_str()); Lxyz->SetName("Lxyz");
  RooDataSet *LxyW = dsLxyW->reduce(cut.c_str()); LxyW->SetName("LxyW");
  RooDataSet *Lxy = dsLxy->reduce(cut.c_str()); Lxy->SetName("Lxy");

  wo->import(*LxyzW);
  wo->import(*Lxyz);
  wo->import(*LxyW);
  wo->import(*Lxy);

  RooBinning rbm(2.6,3.5);
  rbm.addUniform(45,2.6,3.5);

  RooPlot *mfr = wo->var("Jpsi_Mass")->frame();

  wo->data("LxyzW")->statOn(mfr, Label("LxyzWeight"), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.95));
  wo->data("Lxyz")->statOn(mfr, Label("Lxyz"), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.70));
  wo->data("LxyW")->statOn(mfr, Label("LxyWeight"), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.65));
  wo->data("Lxy")->statOn(mfr, Label("Lxy"), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.50));

  wo->data("LxyzW")->plotOn(mfr, Binning(rbm), LineColor(kRed), MarkerColor(kRed), DataError(RooAbsData::SumW2));
  wo->data("Lxyz")->plotOn(mfr, Binning(rbm), LineColor(kBlue), MarkerColor(kBlue), DataError(RooAbsData::SumW2));
  wo->data("LxyW")->plotOn(mfr, Binning(rbm), LineColor(kGreen+2), MarkerColor(kGreen+2), MarkerStyle(kOpenCircle), DataError(RooAbsData::SumW2));
  wo->data("Lxy")->plotOn(mfr, Binning(rbm), LineColor(kOrange+1), MarkerColor(kOrange+1), MarkerStyle(kOpenCircle), DataError(RooAbsData::SumW2));

  TCanvas *canv = new TCanvas("canv","canv",600,600);
  canv->SetLeftMargin(0.16);
  canv->SetRightMargin(0.04);
  canv->Draw();
  mfr->GetYaxis()->SetTitleOffset(2);
  mfr->Draw();
  canv->Update();

  TPaveText *LxyzWSB = (TPaveText*)canv->FindObject("LxyzW_statBox");
  SetStatBox(LxyzWSB, 0.70, 0.96, 0.74, 0.92, kBlue);

  TPaveText *LxyzSB = (TPaveText*)canv->FindObject("Lxyz_statBox");
  SetStatBox(LxyzSB, 0.70, 0.96, 0.56, 0.74, kRed);

  TPaveText *LxyWSB = (TPaveText*)canv->FindObject("LxyW_statBox");
  SetStatBox(LxyWSB, 0.19, 0.46, 0.74, 0.92, kGreen+2);

  TPaveText *LxySB = (TPaveText*)canv->FindObject("Lxy_statBox");
  SetStatBox(LxySB, 0.19, 0.46, 0.56, 0.74, kOrange+1);
 
  TLatex *lat = new TLatex(); lat->SetNDC();
  lat->SetTextSize(0.035);
  lat->SetTextColor(kBlack);
  lat->DrawLatex(0.18,0.5,"1.6|y|<2.4, 3<p_{T}<6.5 GeV/c");

  canv->SaveAs("./Jpsi_Mass.pdf");
}

void SetStatBox(TPaveText *p, double x1, double x2, double y1, double y2, int color) {
  p->SetX1NDC(x1);
  p->SetX2NDC(x2);
  p->SetY1NDC(y1);
  p->SetY2NDC(y2);
  p->SetTextColor(color);
  p->SetTextSize(0.035);
  p->SetBorderSize(0);
}

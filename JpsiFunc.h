//int colorArr[] = {kRed+1, kOrange+7, kSpring+4, kGreen+3, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
int colorArr[] = {kRed+1, kSpring+4, kAzure+1, kBlue+2, kMagenta, kMagenta+2};
int markerArr[] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenTriangleUp, kOpenDiamond, kOpenCross};
int ncolor = sizeof(colorArr)/sizeof(int);
int nmarker = sizeof(markerArr)/sizeof(int);

void SetHistStyleDefault(TH1 *h, int i, int j) {
  h->SetMarkerSize(1.200);
  if (j == 2 || j ==4) h->SetMarkerSize(1.800);
  if (j == 5) h->SetMarkerSize(1.500);

  if (ncolor>i) {
    h->SetMarkerColor(colorArr[i]);
    h->SetLineColor(colorArr[i]);
  } else {
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
  }
  if (nmarker>j) {
    h->SetMarkerStyle(markerArr[j]);
  } else {
    h->SetMarkerStyle(20+j);
  }

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
}

void SetHistStyleDefault(TGraph *h, int i, int j) {
  h->SetMarkerSize(1.200);
  if (j == 2 || j ==4) h->SetMarkerSize(1.800);
  if (j == 5) h->SetMarkerSize(1.500);

  if (ncolor>i) {
    h->SetMarkerColor(colorArr[i]);
    h->SetLineColor(colorArr[i]);
  } else {
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
  }
  if (nmarker>j) {
    h->SetMarkerStyle(markerArr[j]);
  } else {
    h->SetMarkerStyle(20+j);
  }

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
}

void SetHistStyle(TGraph *h, int i, int j, double rmin, double rmax){
  h->SetMinimum(rmin);
  h->SetMaximum(rmax);
  SetHistStyleDefault(h, i, j);
}

void SetHistStyle(TH1 *h, int i, int j, double rmin, double rmax){
  h->GetYaxis()->SetRangeUser(rmin,rmax);
  SetHistStyleDefault(h, i, j);
}

void SetLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.15);
}

void SetStatBox(TPaveText *p, double x1, double x2, double y1, double y2, int color) {
  p->SetX1NDC(x1);
  p->SetX2NDC(x2);
  p->SetY1NDC(y1);
  p->SetY2NDC(y2);
  p->SetTextColor(colorArr[color]);
  p->SetTextSize(0.035);
  p->SetBorderSize(0);
}

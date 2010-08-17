{
  /*
 gSystem->Load("libFWCoreFWLite.so");
 AutoLibraryLoader::enable();
 gSystem->Load("libDataFormatsFWLite.so");
 gROOT->SetStyle ("Plain");
 gSystem->Load("libRooFit") ;
 using namespace RooFit ;
 cout << "loaded" << endl;
  */

 gROOT->SetStyle("Plain");
  
 gStyle->SetPalette(1,0);
 gStyle->SetPadColor(0);
 gStyle->SetPadBorderSize(0);
 gStyle->SetPadBorderMode(0);
 gStyle->SetPadLeftMargin(0.1);
 gStyle->SetPadBottomMargin(0.1);
 gStyle->SetPadTopMargin(0.05);
 gStyle->SetPadRightMargin(0.13);
 gStyle->SetCanvasColor(0);
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetCanvasBorderSize(0);
 gStyle->SetFrameBorderMode(0);
 gStyle->SetFrameLineColor(0);

 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);
 gStyle->SetTextFont(62);
 gStyle->SetLabelSize(0.04,"XYZ");
 gStyle->SetLabelFont(42,"XYZ");
 gStyle->SetTitleFont(42,"XYZ");
 gStyle->SetTitleBorderSize(1);
 gStyle->SetTitleXOffset(0.8);
 gStyle->SetTitleYOffset(0.8);
 gStyle->SetTitleOffset(0.7, "y");
 gStyle->SetTitleSize(0.04, "h");
 gStyle->SetTitleSize(0.06, "x");
 gStyle->SetTitleSize(0.06, "y");
// gStyle->SetTitleStyle(1111);
// gStyle->SetLabelOffset(0.01,"X");
// gStyle->SetLabelOffset(0.01,"Y");
 gStyle->SetTitleColor(1,"XYZ");
 gStyle->SetOptTitle(1);
// gStyle->SetTitleW(1);
 
 gStyle->SetHistFillColor(0);
 gStyle->SetHistFillStyle(0);
 gStyle->SetHistLineColor(1);
 gStyle->SetHistLineStyle(0);
 gStyle->SetHistLineWidth(3);
 gStyle->SetHistLineWidth(1);
 gStyle->SetEndErrorSize(0);
 gStyle->SetErrorX(0);  
// gStyle->SetMarkerStyle(20);
 gStyle->SetMarkerSize(1.3);

 gStyle->SetOptFit(1111);
 gStyle->SetStatColor(0);
 gStyle->SetStatBorderSize(1);
 gStyle->SetOptStat("em");
 gStyle->SetStatX(1-gStyle->GetPadRightMargin()-0.02);
 gStyle->SetStatY(1-gStyle->GetPadTopMargin()-0.02);
 //gStyle->SetOptStat(0);
 //gStyle->SetOptTitle(0);

 cout << "rootlogon::style_0." << endl;

 // Force the Style for the histos...
 gROOT->ForceStyle();
//	gStyle->SetOptStat(1111);
//	gStyle->SetOptFit(1);
	
// 	gStyle->SetTitleSize(0.08,"");
// 	gStyle->SetHistLineWidth(2);	
// 	
// 	gStyle->SetPadGridX( false );
// 	gStyle->SetPadGridY( false );
// 	
// 	gStyle->SetPadTickX( 1 );
// 	gStyle->SetPadTickY( 1 );
// 	
// 	gStyle->SetPalette( 1 );
// 	gStyle->SetLineStyleString( 2, "[20 40]" );
	
//	gStyle->SetCanvasDefH( 700 );
//	gStyle->SetCanvasDefW( 700 );
//gROOT->ForceStyle();
}

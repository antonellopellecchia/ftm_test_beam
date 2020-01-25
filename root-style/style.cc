#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>

#include "style.hh"

TStyle *setStyle() {
  TStyle *style = new TStyle("Garfield", "Garfield Style");
  style->Reset();
  /*style->SetFillColor(1);
  style->SetFillStyle(1001);
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(0);
  style->SetCanvasPreferGL(kTRUE);
  style->SetCanvasDefH(600);
  style->SetCanvasDefW(600);
  style->SetPadBorderMode(0);
  style->SetPadColor(0);
  style->SetPadLeftMargin(0.15);
  style->SetPadBottomMargin(0.1);
  style->SetPadRightMargin(0.05);
  style->SetPadTopMargin(0.05);
  style->SetPadTickX(1);
  style->SetPadTickY(1);
  style->SetPadTickX(1);
  style->SetPadTickY(1);
  style->SetFrameFillColor(0);
  style->SetFrameBorderMode(0);
  style->SetDrawBorder(0);
  style->SetLegendBorderSize(0);

  style->SetGridColor(kGray);
  style->SetGridStyle(3);
  style->SetGridWidth(1);
  style->SetPadGridX(kTRUE);
  style->SetPadGridY(kTRUE);

  // const int font = 132;
  const int font = 42;
  const double tsize = 0.04;
  style->SetTextFont(font);
  style->SetTextSize(tsize);
  style->SetTitleStyle(0);
  style->SetTitleBorderSize(0);
  style->SetTitleColor(1, "xyz");
  style->SetTitleColor(1, "t");
  style->SetTitleFillColor(0);
  style->SetTitleFont(font, "xyz");
  style->SetTitleFont(font, "t");
  style->SetTitleOffset(1.2, "xyz");
  style->SetTitleSize(tsize, "xyz");
  style->SetTitleSize(tsize, "t");

  style->SetLegendFont(font);
  style->SetStatStyle(0);*/
  style->SetStatBorderSize(0);
  /*style->SetStatColor(0);
  style->SetStatFont(font);
  style->SetStatFontSize(tsize);
  style->SetStatX(0.88);
  style->SetStatY(0.88);
  style->SetStatW(0.25);
  style->SetStatH(0.1);
  style->SetOptStat(111110);
  style->SetStatFormat("6.3g");
  style->SetLabelFont(font, "xyz");
  style->SetLabelSize(tsize, "xyz");
  style->SetLabelOffset(0.01, "xyz");
  style->SetOptTitle(0);
  style->SetPaperSize(TStyle::kA4);
  style->SetFuncWidth(2);
  style->SetHistLineColor(kOrange - 3);
  style->SetPalette(1);
  style->SetAxisColor(kBlack, "X");
  style->SetAxisColor(kBlack, "Y");
  style->SetAxisColor(kBlack, "Z");
  style->SetNdivisions(505, "x");
  style->SetNdivisions(510, "y");
  const double lw = 2;
  style->SetLineWidth(lw);
  style->SetLineStyleString(2, "[12 12]");
  style->SetFrameLineWidth(lw);
  style->SetHistLineWidth(lw);
  style->SetFuncWidth(lw);
  style->SetGridWidth(lw);
  style->SetMarkerSize(1.2);*/
  style->cd();

  return style;
}

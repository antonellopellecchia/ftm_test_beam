#include "root-style/style.cc"

void analysis() {
  TStyle *style = setStyle();
  gROOT->SetStyle("Garfield");
  gROOT->ForceStyle();
  
  TFile *rootFileIn = new TFile("build/out/run1.root");
  TTree *inTree = nullptr;
  //TTree *inTreeReader = new TTreeReader("T", inTree);
  rootFileIn->GetObject("runTree", inTree);
  inTree->Print();
  cout << endl << endl;

  TH1F *hElossInScintillator = new TH1F("hElossInScintillator", "", 200, 0., 10.);
  TBranch *branchElossInScintillator = inTree->GetBranch("energyLossInScintillator");
  double elossInScintillator = 0;
  branchElossInScintillator->SetAddress(&elossInScintillator);
  int elossEntries = branchElossInScintillator->GetEntries();
  for (int i=0; i<elossEntries; i++) {
    branchElossInScintillator->GetEvent(i);
    hElossInScintillator->Fill(elossInScintillator);
  }
  hElossInScintillator->Draw();
}

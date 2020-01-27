#include "root-style/style.cc"

void analysis() {
  TStyle *style = setStyle();
  style->cd();
  gROOT->SetStyle("Garfield");
  gROOT->ForceStyle();
  
  TFile *rootFileIn = new TFile("build/out/run3.root");
  TTree *inTree = nullptr;
  //TTree *inTreeReader = new TTreeReader("T", inTree);
  rootFileIn->GetObject("runTree", inTree);
  inTree->Print();
  cout << endl << endl;

  TH1F *hElossInScintillator = new TH1F("hElossInScintillator", "", 200, 0., 10.);
  map<string, TH1F *> hElossByProcess;

  int nEntries = inTree->GetEntries();
  
  TBranch *branchElossInScintillator = inTree->
    GetBranch("energyLossInScintillator");
  double elossInScintillator = 0;
  branchElossInScintillator->SetAddress(&elossInScintillator);
  for (int i=0; i<nEntries; i++) {
    branchElossInScintillator->GetEvent(i);
    hElossInScintillator->Fill(elossInScintillator);
  }
  hElossInScintillator->Draw();

  TBranch *branchElossByProcess = inTree->
    GetBranch("energyLossInScintillatorByProcess");
  //map<string, double> elossByProcess;
  map<string, double> elossByProcess;
  branchElossByProcess->SetAddress(&elossByProcess);
  branchElossByProcess->GetEvent(100);
  cout << elossByProcess.size() << endl;
  /*
  for (map<string, double>::value_type& elossPair:elossByProcess) {
    hElossByProcess[elossPair.first] = new TH1F(string("h"+elossPair.first).c_str(),
				    "", 200, 0., 10.);
    cout << elossPair.first << endl;
  }
  cout << "Fine" << endl;
  for (int i=0; i<nEntries; i++) {
    branchElossByProcess->GetEvent(i);
    for (map<string, double>::value_type& elossPair:elossByProcess)
      hElossByProcess[elossPair.first]->Fill(elossPair.second);
  }
  int i_process = 0;
  for (map<string, double>::value_type& elossPair:elossByProcess) {
    if (i_process==0) hElossByProcess[elossPair.first]->Draw();
    else hElossByProcess[elossPair.first]->Draw();
    i_process++;
  }
  */
}

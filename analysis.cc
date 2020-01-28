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
  TH1F *hAngularDivergence = new TH1F("hAngularDivergence", "", 100, 0., 0.6);
  TH1F *hBeginningPosition = new TH1F("hBeginningPosition", "", 100, 0., 10.);
  TH2F *hBeginningPosition3D =
    new TH2F("hBeginningPosition", "", 80, -8., 8., 80, -8., 8.);
  TH1F *hEndPosition = new TH1F("hEndPosition", "", 100, 0., 10.);
  TH2F *hEndPosition3D =
    new TH2F("hEndPosition", "", 80, -8., 8., 80, -8., 8.);
  TH1F *hElossInQuartz = new TH1F("hElossInQuartz", "", 200, 0., 4.);
  TH1I *hCherenkovCount = new TH1I("hCherenkovCount", "", 200, 7.5e3, 10.5e3);
  
  int nEntries = inTree->GetEntries();
  cout << nEntries << " entries" << endl;


  /*
    Total energy loss in scintillator 1
  */
  TBranch *branchElossInScintillator = inTree->
    GetBranch("energyLossInScintillator");
  double elossInScintillator = 0;
  branchElossInScintillator->SetAddress(&elossInScintillator);
  for (int i=0; i<nEntries; i++) {
    branchElossInScintillator->GetEvent(i);
    hElossInScintillator->Fill(elossInScintillator);
  }
  hElossInScintillator->Draw();

  /*
    Energy loss by process in scintillator 1
  */
  TBranch *branchElossByProcess = inTree->
    GetBranch("energyLossInScintillatorByProcess");
  map<string, double> *elossByProcess = new map<string, double>;
  branchElossByProcess->SetAddress(&elossByProcess);

  for (int i=0; i<nEntries; i++) {
    branchElossByProcess->GetEvent(i);
    for (map<string, double>::value_type& elossPair:*elossByProcess) {
      if(hElossByProcess.count(elossPair.first)==0)
	hElossByProcess[elossPair.first] = new TH1F(string("h"+elossPair.first)
						    .c_str(), "", 200, 0., 10.);
      hElossByProcess[elossPair.first]->Fill(elossPair.second);
    }
  }

  TLegend *legend = new TLegend(0.68, 0.29, 0.89, 0.48);
  legend->SetHeader("Process", "C");
  int i_process = 0;
  for (map<string, double>::value_type& elossPair:*elossByProcess) {
    TH1F *h = hElossByProcess[elossPair.first];
    h->SetLineColor(kOrange+i_process);
    h->SetStats(false);
    legend->AddEntry(h, elossPair.first.c_str(), "l");        
    if (i_process==0) h->Draw();
    else h->Draw("same");
    i_process++;
  }
  legend->Draw();

  /*
    Beam angular divergence
  */ 
  TBranch *branchAngularDivergence = inTree->GetBranch("angularDivergence");
  double angularDivergence = 0;
  branchAngularDivergence->SetAddress(&angularDivergence);
  for (int i=0; i<nEntries; i++) {
    branchAngularDivergence->GetEvent(i);
    hAngularDivergence->Fill(angularDivergence);
  }
  hAngularDivergence->Draw();

  /*
    Initial electron position
  */
  TBranch *branchBeginningPosition = inTree->GetBranch("beginningPosition");
  tuple<double, double> *beginningPosition = new tuple<double, double>;
  branchBeginningPosition->SetAddress(&beginningPosition);
  for (int i=0; i<nEntries; i++) {
    branchBeginningPosition->GetEvent(i);
    double x = get<0>(*beginningPosition);
    double y = get<1>(*beginningPosition);
    hBeginningPosition->Fill(sqrt(x*x+y*y));
    hBeginningPosition3D->Fill(x, y);
  }
  hBeginningPosition->Draw();
  hBeginningPosition3D->Draw("SURF3D");
  
  /*
    Final electron position
  */
  TBranch *branchEndPosition = inTree->GetBranch("endPosition");
  tuple<double, double> *endPosition = new tuple<double, double>;
  branchEndPosition->SetAddress(&endPosition);
  for (int i=0; i<nEntries; i++) {
    branchEndPosition->GetEvent(i);
    double x = get<0>(*endPosition);
    double y = get<1>(*endPosition);
    hEndPosition->Fill(sqrt(x*x+y*y));
    hEndPosition3D->Fill(x, y);
  }
  //hEndPosition->Draw();
  hEndPosition3D->Draw("SURF3D");

  /*
    Energy loss in quartz window
  */
  TBranch *branchElossInQuartz = inTree->GetBranch("energyLossInQuartz");
  double elossInQuartz = 0;
  branchElossInQuartz->SetAddress(&elossInQuartz);
  for (int i=0; i<nEntries; i++) {
    branchElossInQuartz->GetEvent(i);
    hElossInQuartz->Fill(elossInQuartz);
  }
  hElossInQuartz->Draw();

  /*
    Number of Cherenkov photons produced
  */
  TBranch *branchCherenkovCount = inTree->GetBranch("cherenkovCount");
  int cherenkovCount = 0;
  branchCherenkovCount->SetAddress(&cherenkovCount);
  for (int i=0; i<nEntries; i++) {
    branchCherenkovCount->GetEvent(i);
    hCherenkovCount->Fill(cherenkovCount);
  }
  hCherenkovCount->Draw();

  /*
    Cherenkov photon arrival times
  TH1F *hCherenkovArrivalTimes = new TH1F();
  TBranch *branchCherenkovArrivalTimes = inTree->GetBranch("cherenkovArrivalTimes");
  branchCherenkovArrivalTimes->SetAddress(&hCherenkovArrivalTimes);
  branchCherenkovArrivalTimes->GetEvent(1);
  hCherenkovArrivalTimes->Print();
  hCherenkovArrivalTimes->Draw();
  */
}

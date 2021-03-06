#include "root-style/style.cc"

using namespace std;

void analysis(string path="") {
  TStyle *style = setStyle();
  style->cd();
  gROOT->SetStyle("Garfield");
  gROOT->ForceStyle();

  gSystem->Exec(string("mkdir "+path+"/root/").c_str());
  gSystem->Exec(string("mkdir "+path+"/eps/").c_str());
  gSystem->Exec(string("mkdir "+path+"/png/").c_str());
  gSystem->Exec(string("mkdir "+path+"/tex/").c_str());
  gSystem->Exec(string("mkdir "+path+"/svg/").c_str());

  int nbins = 200;
  
  TH1F *hElossInScintillator = new TH1F("hElossInScintillator", "", nbins, 0., 6.5);
  map<string, TH1F *> hElossByProcess;
  TH1F *hAngularDivergence = new TH1F("hAngularDivergence", "", 100, 0., 0.2);
  TH1F *hBeginningPosition = new TH1F("hBeginningPosition", "", 100, 0., 10.);
  TH2F *hBeginningPosition3D = new TH2F("hBeginningPosition3D", "", 100, -8., 8., 100, -8., 8.);
  TH1F *hEndPosition = new TH1F("hEndPosition", "", 100, 0., 10.);
  TH2F *hEndPosition3D = new TH2F("hEndPosition3D", "", 100, -8., 8., 100, -8., 8.);
  TH1F *hElossInQuartz = new TH1F("hElossInQuartz", "", 200, 0., 4.);
  TH1I *hCherenkovCount = new TH1I("hCherenkovCount", "", 200, 400, 800);
  TH1F *hCherenkovTimes = new TH1F("hCherenkovTimes", "", 50000, 1., 1.4);
  
  //TChain *inTree = new TChain("runTree");

  TGaxis::SetExponentOffset(-.07, 0., "Y");
  
  TSystemDirectory *dir = new TSystemDirectory("rootFiles", path.c_str());
  TList *dirFiles = dir->GetListOfFiles();
  TObjLink *fLink;
  for (fLink=dirFiles->FirstLink(); fLink; fLink=fLink->Next()) {
    if (!(fLink->GetObject()->ClassName()==string("TSystemFile"))) continue;
    string rootFilePath = string(path+"/"+fLink->GetObject()->GetName());
    cout << "Processing " << rootFilePath << "... ";
         
    TFile *rootFileIn = new TFile(rootFilePath.c_str());
    TTree *inTree = nullptr;
    rootFileIn->GetObject("runTree", inTree);
    if (!(inTree)) {
      cout << rootFilePath << " is zombie, skipping... " << endl;
      continue;
    }
    int nEntries = inTree->GetEntries();
    cout << nEntries << " entries" << endl;
    
    /*
      Total energy loss in scintillator 1
    */
    TBranch *branchElossInScintillator = inTree->GetBranch("energyLossInScintillator");
    double elossInScintillator = 0;
    branchElossInScintillator->SetAddress(&elossInScintillator);
    for (int i=0; i<nEntries; i++) {
      branchElossInScintillator->GetEvent(i);
      hElossInScintillator->Fill(elossInScintillator);
    }
    
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
	if(hElossByProcess.count(elossPair.first)==0) {
	  cout << "creating " << elossPair.first << endl; 
	  hElossByProcess[elossPair.first] = new TH1F(string("h"+elossPair.first).c_str(), "", nbins, 0., 6.5);
	}
	hElossByProcess[elossPair.first]->Fill(elossPair.second);
      }
    }

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
    //hEndPosition3D->Draw("SURF3D");

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
    //hElossInQuartz->Draw();

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

    /*
      Cherenkov photon arrival times
    */
    TH1 *tempCherenkovTimes = nullptr;
    rootFileIn->GetObject("hCherenkovTimes", tempCherenkovTimes);
    for (int i_bin=0; i_bin<tempCherenkovTimes->GetNbinsX(); i_bin++)
      hCherenkovTimes->SetBinContent(i_bin, tempCherenkovTimes->GetBinContent(i_bin));
    hCherenkovTimes->SetAxisRange(1.2, 1.3);
  }

  TLatex *text = new TLatex();
  text->SetTextSize(0.03);
   
  TCanvas *cElossInScintillator = new TCanvas();
  cElossInScintillator->cd();
  hElossInScintillator->SetXTitle("Energy loss [MeV]");
  hElossInScintillator->Draw();
  text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  text->DrawLatexNDC(.48, .966, "Electron in 1 cm vinyltoluene scintillator");
  cElossInScintillator->SaveAs(string(path+"/root/elossInScintillator.root").c_str());
  cElossInScintillator->SaveAs(string(path+"/eps/elossInScintillator.eps").c_str());
  cElossInScintillator->SaveAs(string(path+"/png/elossInScintillator.png").c_str());
  cElossInScintillator->SaveAs(string(path+"/svg/elossInScintillator.svg").c_str());
  cElossInScintillator->SaveAs(string(path+"/tex/elossInScintillator.tex").c_str());

  TCanvas *cElossByProcess = new TCanvas();
  cElossByProcess->cd();
  TLegend *legend = new TLegend(0.66, 0.49, 0.88, 0.92);
  legend->SetHeader("#bf{Process}", "C");
  int i_process = 0;
  hElossByProcess["Transportation"]->Draw();
  hElossByProcess["Transportation"]->SetXTitle("Energy loss [MeV]");
  for (map<string, TH1F*>::value_type& elossPair:hElossByProcess) {
    TH1F *h = elossPair.second;
    h->SetStats(false);
    h->SetLineColor(style->GetHistLineColor()+i_process*4);
    char legendEntry[200];
    sprintf(legendEntry, "#splitline{%s}{#lower[0.1]{%1.2f #pm %1.2f MeV}}", elossPair.first.c_str(), h->GetMean(), h->GetRMS()); 
    legend->AddEntry(h, legendEntry, "l");        
    h->Draw("same");
    i_process++;
  }
  legend->Draw();
  text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  text->DrawLatexNDC(.48, .966, "Electron in 1 cm vinyltoluene scintillator");
  cElossByProcess->SaveAs(string(path+"/root/elossByProcess.root").c_str());
  cElossByProcess->SaveAs(string(path+"/eps/elossByProcess.eps").c_str());
  cElossByProcess->SaveAs(string(path+"/png/elossByProcess.png").c_str());
  cElossByProcess->SaveAs(string(path+"/svg/elossByProcess.svg").c_str());
  cElossByProcess->SaveAs(string(path+"/tex/elossByProcess.tex").c_str());

  TCanvas *cAngularDivergence = new TCanvas();
  cAngularDivergence->cd();
  hAngularDivergence->SetXTitle("Angular deviation [rad]");
  hAngularDivergence->Draw();
  text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  text->DrawLatexNDC(.48, .966, "Electron in 1 cm vinyltoluene scintillator");
  cAngularDivergence->SaveAs(string(path+"/root/angularDivergence.root").c_str());
  cAngularDivergence->SaveAs(string(path+"/eps/angularDivergence.eps").c_str());
  cAngularDivergence->SaveAs(string(path+"/png/angularDivergence.png").c_str());
  cAngularDivergence->SaveAs(string(path+"/svg/angularDivergence.svg").c_str());
  cAngularDivergence->SaveAs(string(path+"/tex/angularDivergence.tex").c_str());

  TCanvas *cElectronPosition = new TCanvas("electronPosition", "", 1200, 1200);
  cElectronPosition->Divide(2, 2);
  cElectronPosition->cd(1);
  hBeginningPosition->SetXTitle("Electron radial position [mm]");
  hBeginningPosition->Draw();
  text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  text->DrawLatexNDC(.62, .966, "Initial electron beam profile");
  cElectronPosition->cd(2);
  hEndPosition->SetXTitle("Electron radial position [mm]");
  hEndPosition->Draw();
  text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  text->DrawLatexNDC(.625, .966, "Final electron beam profile");

  //TCanvas *cElectronPosition3D = new TCanvas("electronPosition3D", "", 1200, 600);
  //cElectronPosition3D->Divide(2, 1);
  //cElectronPosition3D->cd(1);
  cElectronPosition->cd(3);
  hBeginningPosition3D->SetXTitle("Electron x position [mm]");
  hBeginningPosition3D->SetYTitle("Electron y position [mm]");
  hBeginningPosition3D->SetStats(false);
  hBeginningPosition3D->Draw("CONT4");
  gPad->SetGrid(kFALSE, kFALSE);
  char sigmaStr[200];
  sprintf(sigmaStr, "#bf{#color[1]{Initial beam radius %1.2f mm}}", 0.5*(hBeginningPosition3D->GetRMS(1)+hBeginningPosition3D->GetRMS(2)));
  text->DrawLatexNDC(.15, .966, sigmaStr);
  //text->DrawLatexNDC(.62, .966, "Initial electron beam profile");
  //cElectronPosition3D->cd(2);
  cElectronPosition->cd(4);
  hEndPosition3D->SetXTitle("Electron x position [mm]");
  hEndPosition3D->SetYTitle("Electron y position [mm]");
  hEndPosition3D->SetStats(false);
  hEndPosition3D->Draw("CONT4");
  gPad->SetGrid(kFALSE, kFALSE);
  sprintf(sigmaStr, "#bf{#color[1]{Final beam radius %1.2f mm}}", 0.5*(hEndPosition3D->GetRMS(1)+hEndPosition3D->GetRMS(2)));
  text->DrawLatexNDC(.15, .966, sigmaStr);
  //text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  //text->DrawLatexNDC(.62, .966, "Final electron beam profile");
  //cElectronPosition3D->SaveAs(string(path+"/eps/electronPosition3D.eps").c_str());
  cElectronPosition->SaveAs(string(path+"/root/electronPosition.root").c_str());
  cElectronPosition->SaveAs(string(path+"/eps/electronPosition.eps").c_str());
  cElectronPosition->SaveAs(string(path+"/png/electronPosition.png").c_str());
  cElectronPosition->SaveAs(string(path+"/svg/electronPosition.svg").c_str());
  cElectronPosition->SaveAs(string(path+"/tex/electronPosition.tex").c_str());

  TCanvas *cElossInQuartz = new TCanvas();
  cElossInQuartz->cd();
  hElossInQuartz->SetXTitle("Energy loss [MeV]");
  hElossInQuartz->Draw();
  text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  text->DrawLatexNDC(.57, .966, "Electron in 5 mm quartz window");
  cElossInQuartz->SaveAs(string(path+"/root/elossInQuartz.root").c_str());
  cElossInQuartz->SaveAs(string(path+"/eps/elossInQuartz.eps").c_str());
  cElossInQuartz->SaveAs(string(path+"/png/elossInQuartz.png").c_str());
  cElossInQuartz->SaveAs(string(path+"/svg/elossInQuartz.svg").c_str());
  cElossInQuartz->SaveAs(string(path+"/tex/elossInQuartz.tex").c_str());

  TCanvas *cCherenkovCount = new TCanvas();
  cCherenkovCount->cd();
  hCherenkovCount->SetXTitle("Number of Cherenkov photons");
  hCherenkovCount->Draw();
  text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  text->DrawLatexNDC(.53, .966, "Electron in 1.27 mm sapphire radiator");
  cCherenkovCount->SaveAs(string(path+"/root/cherenkovCount.root").c_str());
  cCherenkovCount->SaveAs(string(path+"/eps/cherenkovCount.eps").c_str());
  cCherenkovCount->SaveAs(string(path+"/png/cherenkovCount.png").c_str());
  cCherenkovCount->SaveAs(string(path+"/svg/cherenkovCount.svg").c_str());
  cCherenkovCount->SaveAs(string(path+"/tex/cherenkovCount.tex").c_str());

  TCanvas *cCherenkovTimes = new TCanvas();
  cCherenkovTimes->cd();
  hCherenkovTimes->SetXTitle("Cherenkov photon arrival time [ns]");
  hCherenkovTimes->Draw();
  text->DrawLatexNDC(.15, .966, "#bf{GEANT4}");
  text->DrawLatexNDC(.53, .966, "Electron in 1.27 mm sapphire radiator");
  cCherenkovTimes->SaveAs(string(path+"/root/cherenkovTimes.root").c_str());
  cCherenkovTimes->SaveAs(string(path+"/eps/cherenkovTimes.eps").c_str());
  cCherenkovTimes->SaveAs(string(path+"/png/cherenkovTimes.png").c_str());
  cCherenkovTimes->SaveAs(string(path+"/svg/cherenkovTimes.svg").c_str());
  cCherenkovTimes->SaveAs(string(path+"/tex/cherenkovTimes.tex").c_str());
}

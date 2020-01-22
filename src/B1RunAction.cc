//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TFile.h> 
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction(G4bool headless)
  : G4UserRunAction(),
    fHeadless(true),
    fEdep(0.),
    fEdep2(0.),
    fEdepVector(0),
    fEdepVectorByProcess(),
    fDeviationAngleVector(0),
    fBeginningPositionVector(0),
    fEndPositionVector(0),
    fQuartzWindow1EdepVector(0),
    fCherenkovCounts(0)
    //fScintillatorHitPosition(0., 0., 0.)
{
  fHeadless = headless;
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run* run)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  nOfEvents = run->GetNumberOfEventToBeProcessed();
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  G4cout << G4endl;

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const B1DetectorConstruction* detectorConstruction = static_cast<const B1DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume1()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  if (fHeadless) {
    G4cout << G4endl;
    /*
      Store histogram with energy loss distribution in the first scintillator
    */
    TCanvas *edepCanvas = new TCanvas("c", "c", 800, 600);
    G4double max_edep = *max_element(fEdepVector.begin(),fEdepVector.end());
    G4double min_edep = *min_element(fEdepVector.begin(),fEdepVector.end());
    edepCanvas->cd();
    int nbins = 100;
    TH1F *hEdep = new TH1F("hEdep", "Energy deposit in scintillator 1", nbins, min_edep, max_edep);
    for (G4double dep:fEdepVector) hEdep->Fill(dep);
    hEdep->SetStats(false);
    //hEdep->Fit("landau");
 
    TLegend *legend = new TLegend(0.68, 0.69, 0.89, 0.88);
    legend->SetHeader("Process", "C");
    legend->AddEntry(hEdep, "Total", "l");
    
    hEdep->GetXaxis()->SetTitle("Energy deposit [MeV]");
    hEdep->GetYaxis()->SetTitle("Event rate");
    hEdep->Draw();
    legend->Draw();
    int b1 = 0; //hEdep->FindBin(3.3);
    int b2 = nbins; //hEdep->FindBin(3.6);
    // Additional: energy loss by process
    int iEdep = 0;
    for (std::map<G4String, std::vector<G4double>>::value_type& edepPair: fEdepVectorByProcess) {
      G4String process = edepPair.first;
      std::vector<G4double> edepVector = edepPair.second;
      const char *hTitle = G4String("hEdep_"+process).c_str();

      TH1F *hEdepByProcess = new TH1F(hTitle, "Energy deposit in scintillator 1", nbins, min_edep, max_edep);
      for (G4double dep:edepVector) hEdepByProcess->Fill(dep);
      G4cout << process << ": " << hEdepByProcess->Integral(b1, b2) << ", ";
      G4cout << "mean: " << hEdepByProcess->GetMean() << ", ";
      G4cout << "rms: " << hEdepByProcess->GetRMS() << G4endl;

      hEdepByProcess->SetStats(false);
      hEdepByProcess->SetLineColor(2+iEdep);
      hEdepByProcess->SetLineWidth(2);
      legend->AddEntry(hEdepByProcess, process.c_str(), "l");
      hEdepByProcess->Draw("same");
      iEdep++;
    }
    G4cout << "total" << ": " << hEdep->Integral(b1, b2) << ", ";
    G4cout << "mean: " << hEdep->GetMean() << ", ";
    G4cout << "rms: " << hEdep->GetRMS() << G4endl;
    G4cout << G4endl;
    
    for (std::map<G4String, G4int>::value_type& countPair: fDepositCount) {
      G4String process = countPair.first;
      G4double count = countPair.second;
      G4cout << "Process: " << process << ", ";
      G4cout << "n. deposits " << count << G4endl;      
    }
    G4cout << G4endl;

    edepCanvas->SetLogy();
    edepCanvas->SetGridx();
    edepCanvas->SetGridy();
    edepCanvas->SaveAs("./root/energy_loss.root");
    G4double meanScintillatorEdep = hEdep->GetMean();
    G4double rmsScintillatorEdep = hEdep->GetRMS();
        
  
    /*
      Store histogram with deviation angle when hitting the second scintillator
    */
    G4double max_deviation = *max_element(fDeviationAngleVector.begin(),fDeviationAngleVector.end());
    G4double min_deviation = *min_element(fDeviationAngleVector.begin(),fDeviationAngleVector.end());
    TH1F *hDeviationAngles = new TH1F("hDeviationAngles", "Deviation angles", 5000, min_deviation, max_deviation);
    for (G4double deviationAngle:fDeviationAngleVector) hDeviationAngles->Fill(deviationAngle);
    hDeviationAngles->GetXaxis()->SetTitle("Deviation angle (radians)");
    hDeviationAngles->GetYaxis()->SetTitle("Event rate");
    hDeviationAngles->SaveAs("./root/angle_deviation.root");
    G4double meanDeviationAngle = hDeviationAngles->GetMean()*180/CLHEP::pi;
    G4double rmsDeviationAngle = hDeviationAngles->GetRMS()*180/CLHEP::pi;


    /*
      Store initial and final beam profile scatter plot
    */
    TCanvas *profileCanvas = new TCanvas("c", "c", 1400, 700);
    profileCanvas->Divide(2, 1);

    profileCanvas->cd(1);
    TGraph *plotBeginningPosition = new TGraph(fBeginningPositionVector.size());
    for (unsigned long int i=0; i<fBeginningPositionVector.size(); i++)
      plotBeginningPosition->SetPoint(i, std::get<0>(fBeginningPositionVector[i]), std::get<1>(fBeginningPositionVector[i]));
    plotBeginningPosition->SetTitle("Beginning position");
    plotBeginningPosition->GetXaxis()->SetTitle("X position (mm)");
    plotBeginningPosition->GetYaxis()->SetTitle("Y position (mm)");
    plotBeginningPosition->GetXaxis()->SetLimits(-10., 10.);
    plotBeginningPosition->GetYaxis()->SetLimits(-10., 10.);
    plotBeginningPosition->GetXaxis()->SetRangeUser(-10., 10.);
    plotBeginningPosition->GetYaxis()->SetRangeUser(-10., 10.);
    plotBeginningPosition->Draw("ap");
    plotBeginningPosition->SaveAs("./root/beginning_profile.root");

    profileCanvas->cd(2);
    TGraph *plotEndPosition = new TGraph(fEndPositionVector.size());
    for (unsigned long int i=0; i<fEndPositionVector.size(); i++)
      plotEndPosition->SetPoint(i, std::get<0>(fEndPositionVector[i]), std::get<1>(fEndPositionVector[i]));
    plotEndPosition->SetTitle("End position");
    plotEndPosition->GetXaxis()->SetTitle("X position (mm)");
    plotEndPosition->GetYaxis()->SetTitle("Y position (mm)");
    plotEndPosition->GetXaxis()->SetLimits(-10., 10.);
    plotEndPosition->GetYaxis()->SetLimits(-10., 10.);
    plotEndPosition->GetXaxis()->SetRangeUser(-10., 10.);
    plotEndPosition->GetYaxis()->SetRangeUser(-10., 10.);
    plotEndPosition->Draw("ap");
    plotEndPosition->SaveAs("./root/end_profile.root");
    
    profileCanvas->SaveAs("./root/profile.root");

    
    /*
      Store histogram with initial and final particle radial position
    */
    TCanvas *positionCanvas = new TCanvas("c", "c", 1400, 700);
    positionCanvas->Divide(2, 1);

    positionCanvas->cd(1);
    TH1F *hBeginningPosition = new TH1F("hBeginningPosition", "Beginning position", 1000, 0, 10.*mm);
    for (auto beginningPosition:fBeginningPositionVector) {
      G4double x = std::get<0>(beginningPosition);
      G4double y = std::get<1>(beginningPosition);
      hBeginningPosition->Fill(sqrt(x*x+y*y));
    }
    hBeginningPosition->GetXaxis()->SetTitle("Radial position (mm)");
    hBeginningPosition->GetYaxis()->SetTitle("Event rate");
    hBeginningPosition->Draw();
    hBeginningPosition->SaveAs("./root/beginning_position.root");

    positionCanvas->cd(2);
    TH1F *hEndPosition = new TH1F("hEndPosition", "End position", 1000, 0, 10.*mm);
    for (auto endPosition:fEndPositionVector) {
      G4double x = std::get<0>(endPosition);
      G4double y = std::get<1>(endPosition);
      hEndPosition->Fill(sqrt(x*x+y*y));
    }
    hEndPosition->GetXaxis()->SetTitle("Radial position (mm)");
    hEndPosition->GetYaxis()->SetTitle("Event rate");
    hEndPosition->Draw();
    hEndPosition->SaveAs("./root/end_position.root");

    positionCanvas->SaveAs("./root/position.root");


    /*
      Store histogram with initial and final particle radial position
    */
    TCanvas *positionCanvas3d = new TCanvas("c", "c", 1400, 700);
    positionCanvas3d->Divide(2, 1);

    positionCanvas3d->cd(1);
    TH2F *hBeginningPosition3d = new TH2F("hBeginningPosition3d", "Beginning position", 80, -8., 8., 80, -8., 8.);
    for (auto beginningPosition:fBeginningPositionVector) {
      G4double x = std::get<0>(beginningPosition);
      G4double y = std::get<1>(beginningPosition);
      hBeginningPosition3d->Fill(x, y);
    }
    hBeginningPosition3d->GetXaxis()->SetTitle("X position (mm)");
    hBeginningPosition3d->GetYaxis()->SetTitle("Y position (mm)");
    hBeginningPosition3d->GetZaxis()->SetTitle("Event rate");
    hBeginningPosition3d->Draw("SURF3");
    hBeginningPosition3d->SaveAs("./root/beginning_position3d.root");

    positionCanvas3d->cd(2);
    TH2F *hEndPosition3d = new TH2F("hEndPosition3d", "End position", 80, -8., 8., 80, -8., 8.);
    for (auto endPosition:fEndPositionVector) {
      G4double x = std::get<0>(endPosition);
      G4double y = std::get<1>(endPosition);
      hEndPosition3d->Fill(x, y);
    }
    hEndPosition3d->GetXaxis()->SetTitle("X position (mm)");
    hEndPosition3d->GetYaxis()->SetTitle("Y position (mm)");
    hEndPosition3d->GetZaxis()->SetTitle("Event rate");
    hEndPosition3d->Draw("SURF3");
    hEndPosition3d->SaveAs("./root/end_position3d.root");

    positionCanvas3d->SaveAs("./root/position3d.root");


    
    /*
      Store histogram with energy loss distribution in the first quartz window
    */
    G4double max_quartz_edep = *max_element(fQuartzWindow1EdepVector.begin(),fQuartzWindow1EdepVector.end());
    G4double min_quartz_edep = *min_element(fQuartzWindow1EdepVector.begin(),fQuartzWindow1EdepVector.end());
    TH1F *hQuartz1Edep = new TH1F("hQuartzEdep", "Energy deposit in first quartz window", 5000, min_quartz_edep, max_quartz_edep);
    for (G4double dep:fQuartzWindow1EdepVector) hQuartz1Edep->Fill(dep);
    hQuartz1Edep->GetXaxis()->SetTitle("Energy deposit [MeV]");
    hQuartz1Edep->GetYaxis()->SetTitle("Event rate");
    hQuartz1Edep->SaveAs("./root/energy_loss_quartz_window.root");
    G4double meanQuartz1Edep = hQuartz1Edep->GetMean();
    G4double rmsQuartz1Edep = hQuartz1Edep->GetRMS();


    /*
      Store Cherenkov photon profile scatter plot
    */
    TCanvas *cherenkovProfileCanvas = new TCanvas("c", "c", 1400, 700);
    cherenkovProfileCanvas->cd();
    TGraph *plotCherenkovPosition = new TGraph(fCherenkovEndpointVector.size());
    for (unsigned long int i=0; i<fCherenkovEndpointVector.size(); i++)
      plotCherenkovPosition->SetPoint(i, fCherenkovEndpointVector[i].getX(), fCherenkovEndpointVector[i].getY());
    plotCherenkovPosition->SetTitle("Cherenkov profile");
    plotCherenkovPosition->GetXaxis()->SetTitle("X position (mm)");
    plotCherenkovPosition->GetYaxis()->SetTitle("Y position (mm)");
    plotCherenkovPosition->GetXaxis()->SetLimits(-10., 10.);
    plotCherenkovPosition->GetYaxis()->SetLimits(-10., 10.);
    plotCherenkovPosition->GetXaxis()->SetRangeUser(-10., 10.);
    plotCherenkovPosition->GetYaxis()->SetRangeUser(-10., 10.);
    plotCherenkovPosition->Draw("ap");
    plotEndPosition->SaveAs("./root/cherenkov_profile.root");

    /*
      Store histogram with Cherenkov photons radial position
    */
    TCanvas *cherenkovCanvas3d = new TCanvas("c", "c", 1400, 700);
    cherenkovCanvas3d->cd();
    TH2F *hCherenkovPosition3d = new TH2F("hCherenkovPosition3d", "Cherenkov position", 80, -8., 8., 80, -8., 8.);
    for (auto cherenkovPosition:fCherenkovEndpointVector)
      hCherenkovPosition3d->Fill(cherenkovPosition.getX(), cherenkovPosition.getY());
    hCherenkovPosition3d->GetXaxis()->SetTitle("X position (mm)");
    hCherenkovPosition3d->GetYaxis()->SetTitle("Y position (mm)");
    hCherenkovPosition3d->GetZaxis()->SetTitle("Event rate");
    hCherenkovPosition3d->Draw("SURF3");
    positionCanvas3d->SaveAs("./root/cherenkov_position3d.root");
    

    TCanvas *cherenkovTimesCanvas = new TCanvas("c", "c", 800, 600);
    cherenkovTimesCanvas->cd();
    G4double max_cherenkov = *max_element(fCherenkovArrivalTimes.begin(),fCherenkovArrivalTimes.end());
    G4double min_cherenkov = *min_element(fCherenkovArrivalTimes.begin(),fCherenkovArrivalTimes.end());
    TH1F *hCherenkovTimes = new TH1F("hCherenkovTimes", "Cherenkov photon arrival times", 500, min_cherenkov, max_cherenkov);
    for (G4double time:fCherenkovArrivalTimes) hCherenkovTimes->Fill(time);
    hCherenkovTimes->SetStats(false);
    hCherenkovTimes->GetXaxis()->SetTitle("Arrival time [ns]");
    hCherenkovTimes->GetYaxis()->SetTitle("Event rate");
    hCherenkovTimes->Draw();
    cherenkovTimesCanvas->SetGridx();
    cherenkovTimesCanvas->SetGridy();
    cherenkovTimesCanvas->SaveAs("./root/cherenkov_arrival_times.root");
    cherenkovTimesCanvas->SaveAs("./root/cherenkov_arrival_times.eps");


    /*
      Store histogram with Cherenkov photon counter
    */
    TCanvas *cherenkovCountCanvas = new TCanvas("c", "c", 1400, 700);
    G4double max_cherenkov_count = *max_element(fCherenkovCounts.begin(),fCherenkovCounts.end());
    G4double min_cherenkov_count = *min_element(fCherenkovCounts.begin(),fCherenkovCounts.end());
    cherenkovCountCanvas->cd();
    TH1F *hCherenkovCount = new TH1F("hCherenkovCount", "Cherenkov photons", 500, min_cherenkov_count, max_cherenkov_count);
    for (auto cherenkovCount:fCherenkovCounts) hCherenkovCount->Fill(cherenkovCount);
    hCherenkovCount->GetXaxis()->SetTitle("Cherenkov photons");
    hCherenkovCount->GetYaxis()->SetTitle("Event rate");
    hCherenkovCount->Draw();
    cherenkovCountCanvas->SaveAs("./root/cherenkov_count.root");
    

    G4cout << G4endl;
    G4cout << "----------------------------------------------------" << G4endl;
    G4cout << "Mean energy loss in scintillator 1 [MeV]: " << meanScintillatorEdep << " +/- " << rmsScintillatorEdep << G4endl;
    G4cout << "Mean deviation angle measured at scintillator 2 [degrees]: " << meanDeviationAngle << " +/- " << rmsDeviationAngle << G4endl;
    G4cout << "Mean energy loss in first quartz window [MeV]: " << meanQuartz1Edep << " +/- " << rmsQuartz1Edep << G4endl;
    G4cout << "----------------------------------------------------" << G4endl;
    G4cout << G4endl << G4endl;
  }
    
  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
    = static_cast<const B1PrimaryGeneratorAction*>
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
    {
      const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
      runCondition += particleGun->GetParticleDefinition()->GetParticleName();
      runCondition += " of ";
      G4double particleEnergy = particleGun->GetParticleEnergy();
      runCondition += G4BestUnit(particleEnergy,"Energy");
    }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
      << G4endl
      << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
      << G4endl
      << "--------------------End of Local Run------------------------";
  }
  
  G4cout
    << G4endl
    << " The run consists of " << nofEvents << " "<< runCondition
    << G4endl
    << " Cumulated dose per run, in scoring volume : " 
    << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
    << G4endl
    << "------------------------------------------------------------"
    << G4endl
    << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
  fEdepVector.push_back(edep);
}

void B1RunAction::AddEdepByProcess(std::map<G4String, G4double> edepByProcess) {
  for (std::map<G4String, G4double>::value_type& edepPair: edepByProcess) {
    G4String process = edepPair.first;
    G4double edep = edepPair.second;
    fEdepVectorByProcess[process].push_back(edep);
  }
}

void B1RunAction::AddDepositCount(std::map<G4String, G4int> depositCount) {
  for (std::map<G4String, G4int>::value_type& countPair: depositCount) {
    G4String process = countPair.first;
    G4double count = countPair.second;
    fDepositCount[process] += count;
  }
}

void B1RunAction::AddDeviationAngle(G4double deviationAngle)
{
  if (deviationAngle != -10.) fDeviationAngleVector.push_back(deviationAngle);
}

void B1RunAction::AddBeginningPosition (std::tuple<G4double, G4double> beginningPosition) {
  if (beginningPosition != std::make_tuple(1.e6, 1.e6)) {
    fBeginningPositionVector.push_back(beginningPosition);
    //G4cout << "Beginning position: " << beginningPosition;
  }
}

void B1RunAction::AddEndPosition (std::tuple<G4double, G4double> endPosition) {
  if (endPosition != std::make_tuple(1.e6, 1.e6)) {
    fEndPositionVector.push_back(endPosition);
    //G4cout << " end: " << endPosition << G4endl;
  }
}

void B1RunAction::AddQuartzWindow1Edep (G4double edep) {
  fQuartzWindow1EdepVector.push_back(edep);
}

void B1RunAction::AddCherenkovEndpointVector (std::vector<G4ThreeVector> cherenkovEndpoints) {
  for (auto endpoint: cherenkovEndpoints) fCherenkovEndpointVector.push_back(endpoint);
}

void B1RunAction::AddCherenkovArrivalTime(G4double arrivalTime) {
  fCherenkovArrivalTimes.push_back(arrivalTime);
}

void B1RunAction::AddCherenkovCount(G4int cherenkovCount) {
  fCherenkovCounts.push_back(cherenkovCount);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


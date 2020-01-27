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

#include <iostream>
#include <fstream>
#include <sstream>

#include <sys/stat.h>

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
#include <TTree.h>

#include "style.cc"

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

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction(G4bool headless, string outFilePath)
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
    fCherenkovCounts(0),

    fEnergyLossInScintillator(0),
    fEnergyLossInScintillatorByProcess(),
    fAngularDivergence(0),
    fBeginningPosition(0, 0),
    fEndPosition(0, 0),
    fEnergyLossInQuartz(0),
    fCherenkovCount(0)
    //fCherenkovArrivalTimes(0)
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

  runFile = new TFile(outFilePath.c_str(), "RECREATE", "File containing simulation output ntuples");

  fCherenkovArrivalTimes = new TH1F();//"hCherenkovTimes", "", 500, 0., 2.*ns);

  runTree = new TTree("runTree", "Tree with all run data");
  runTree->Branch("energyLossInScintillator", &fEnergyLossInScintillator, "loss/D");
  runTree->Branch("energyLossInScintillatorByProcess",
		  &fEnergyLossInScintillatorByProcess);
  runTree->Branch("angularDivergence", &fAngularDivergence, "angle/D");
  runTree->Branch("beginningPosition", &fBeginningPosition);
  runTree->Branch("endPosition", &fEndPosition);
  runTree->Branch("energyLossInQuartz", &fEnergyLossInQuartz, "loss/D");
  runTree->Branch("cherenkovCount", &fCherenkovCount, "count/I");
  runTree->Branch("cherenkovArrivalTimes", fCherenkovArrivalTimes);
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
  //TStyle *garfieldStyle = setStyle();
  //gROOT->SetStyle("Garfield");
  //gROOT->ForceStyle();
  //gStyle->SetStatBorderSize(0);
  
  G4String out_dir = G4String("./out/");
  G4String root_out_dir = G4String(out_dir+"root/");
  G4String eps_out_dir = G4String(out_dir+"eps/");
  mkdir(out_dir.c_str(), 0700);
  mkdir(root_out_dir.c_str(), 0700);
  mkdir(eps_out_dir.c_str(), 0700);

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
    
    // Energy loss in the first scintillator
    runTree->Print();
    runFile->Write();
    runFile->Close();
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

void B1RunAction::FillRunNtuples(G4double energyLossInScintillator,
				 map<string, G4double>
				 energyLossInScintillatorByProcess,
				 G4double angularDivergence,
				 tuple<G4double, G4double> beginningPosition,
				 tuple<G4double, G4double> endPosition,
				 G4double energyLossInQuartz,
				 G4int cherenkovCount,
				 TH1F *cherenkovArrivalTimes) {
  fEnergyLossInScintillator = energyLossInScintillator;
  fEnergyLossInScintillatorByProcess = energyLossInScintillatorByProcess;
  fAngularDivergence = angularDivergence;
  fBeginningPosition = beginningPosition;
  fEndPosition = endPosition;
  fEnergyLossInQuartz = energyLossInQuartz;
  fCherenkovCount = cherenkovCount;
  fCherenkovArrivalTimes = cherenkovArrivalTimes;
  runTree->Fill();
}

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
  //fEdepVector.push_back(edep);
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
  fBeginningPositionVector.push_back(beginningPosition);
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
  //fCherenkovArrivalTimes.push_back(arrivalTime);
  //fCherenkovArrivalTimes->Fill(arrivalTime);
}

void B1RunAction::AddCherenkovCount(G4int cherenkovCount) {
  fCherenkovCounts.push_back(cherenkovCount);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


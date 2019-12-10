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
#include <TH2F.h>
#include <TFile.h> 
#include <TStyle.h>
#include <TF1.h>

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction(G4bool headless)
  : G4UserRunAction(),
    fHeadless(true),
    fEdep(0.),
    fEdep2(0.),
    fEdepVector(0),
    fDeviationAngleVector(0),
    fBeginningPositionVector(0),
    fEndPositionVector(0),
    fQuartzWindow1EdepVector(0)
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

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

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
    /*
      Store histogram with energy loss distribution in the first scintillator
    */
    G4double max_edep = *max_element(fEdepVector.begin(),fEdepVector.end());
    G4double min_edep = *min_element(fEdepVector.begin(),fEdepVector.end());
    TH1F *hEdep = new TH1F("hEdep", "Energy deposit in scintillator 1", 5000, min_edep, max_edep/2);
    for (G4double dep:fEdepVector) hEdep->Fill(dep);
    //hEdep->Fit("landau");
    hEdep->GetXaxis()->SetTitle("Energy deposit [MeV]");
    hEdep->GetYaxis()->SetTitle("Event rate");
    hEdep->SaveAs("./root/energy_loss.root");
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
      Store histogram with initial and final particle radial position
    */
    TCanvas *positionCanvas = new TCanvas("c", "c", 800, 600);
    positionCanvas->Divide(1, 2);

    positionCanvas->cd(1);
    G4double max_r0 = *max_element(fBeginningPositionVector.begin(),fBeginningPositionVector.end());
    G4double min_r0 = *min_element(fBeginningPositionVector.begin(),fBeginningPositionVector.end());
    TH1F *hBeginningPosition = new TH1F("hBeginningPosition", "Beginning position", 5000, 0, 20.*mm);
    for (G4double beginningPosition:fBeginningPositionVector) hBeginningPosition->Fill(beginningPosition);
    hBeginningPosition->GetXaxis()->SetTitle("Radial position (mm)");
    hBeginningPosition->GetYaxis()->SetTitle("Event rate");
    hBeginningPosition->Draw();
    hBeginningPosition->SaveAs("./root/beginning_position.root");

    positionCanvas->cd(2);
    G4double max_r1 = *max_element(fEndPositionVector.begin(),fEndPositionVector.end());
    G4double min_r1 = *min_element(fEndPositionVector.begin(),fEndPositionVector.end());
    TH1F *hEndPosition = new TH1F("hEndPosition", "End position", 5000, 0, 20.*mm);
    for (G4double endPosition:fEndPositionVector) hEndPosition->Fill(endPosition);
    hEndPosition->GetXaxis()->SetTitle("Radial position (mm)");
    hEndPosition->GetYaxis()->SetTitle("Event rate");
    hEndPosition->SaveAs("./root/end_position.root");
    hEndPosition->Draw();
  
    positionCanvas->SaveAs("./root/position.root");


    /*
      Store histogram with energy loss distribution in the first scintillator
    */
    G4double max_quartz_edep = *max_element(fQuartzWindow1EdepVector.begin(),fQuartzWindow1EdepVector.end());
    G4double min_quartz_edep = *min_element(fQuartzWindow1EdepVector.begin(),fQuartzWindow1EdepVector.end());
    TH1F *hQuartz1Edep = new TH1F("hQuartzEdep", "Energy deposit in first quartz window", 5000, min_quartz_edep, max_quartz_edep/2);
    for (G4double dep:fQuartzWindow1EdepVector) hQuartz1Edep->Fill(dep);
    hQuartz1Edep->GetXaxis()->SetTitle("Energy deposit [MeV]");
    hQuartz1Edep->GetYaxis()->SetTitle("Event rate");
    hQuartz1Edep->SaveAs("./root/energy_loss_quartz_window.root");
    G4double meanQuartz1Edep = hQuartz1Edep->GetMean();
    G4double rmsQuartz1Edep = hQuartz1Edep->GetRMS();

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

void B1RunAction::AddDeviationAngle(G4double deviationAngle)
{
  if (deviationAngle != -10.) fDeviationAngleVector.push_back(deviationAngle);
}

void B1RunAction::AddBeginningPosition (G4double beginningPosition) {
  if (beginningPosition != -1.) {
    fBeginningPositionVector.push_back(beginningPosition);
    //G4cout << "Beginning position: " << beginningPosition;
  }
}

void B1RunAction::AddEndPosition (G4double endPosition) {
  if (endPosition != -1.) {
    fEndPositionVector.push_back(endPosition);
    //G4cout << " end: " << endPosition << G4endl;
  }
}

void B1RunAction::AddQuartzWindow1Edep (G4double edep) {
  fQuartzWindow1EdepVector.push_back(edep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


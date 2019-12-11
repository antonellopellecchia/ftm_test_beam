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
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include <TMath.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume1(0),
  fScoringVolume2(0),
  fQuartzWindow1(0),
  fQuartzWindow2(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume1 || !fScoringVolume2) { 
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume1 = detectorConstruction->GetScoringVolume1();   
    fScoringVolume2 = detectorConstruction->GetScoringVolume2();
    fQuartzWindow1 = detectorConstruction->GetQuartzWindow1();
    fQuartzWindow2 = detectorConstruction->GetQuartzWindow2();
  }

  // get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
      
  // check if we are in one of the scoring volumes (scintillators)
  if (volume == fScoringVolume1) {
    // if we are in the first scintillator,
    // collect energy deposited
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edepStep);
  } else if (volume == fScoringVolume2) {
    // if we are in the second scintillator,
    // store the hit position
    G4Track *track = step->GetTrack();
    G4ThreeVector position = track->GetPosition();
    G4ThreeVector momentumDirection = track->GetMomentumDirection();
    G4ThreeVector vertexPosition = track->GetVertexPosition();
    G4ThreeVector vertexMomentumDirection = track->GetVertexMomentumDirection();
    G4double scalarProduct = 0.;
    scalarProduct += momentumDirection.getX()*vertexMomentumDirection.getX();
    scalarProduct += momentumDirection.getY()*vertexMomentumDirection.getY();
    scalarProduct += momentumDirection.getZ()*vertexMomentumDirection.getZ();
    G4double deviationAngle = TMath::ACos(scalarProduct);
    fEventAction->AddDeviationAngle(deviationAngle);

    G4double x0 = vertexPosition.getX();
    G4double y0 = vertexPosition.getY();
    G4double x1 = position.getX();
    G4double y1 = position.getY();
    std::tuple<G4double, G4double> posBeginning = std::make_tuple(x0, y0);
    std::tuple<G4double, G4double> posEnd = std::make_tuple(x1, y1);
    fEventAction->AddBeginningPosition(posBeginning);
    fEventAction->AddEndPosition(posEnd);
  } else if (volume == fQuartzWindow1) {
    // if we are in the first quartz window,
    // store the energy deposited
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddQuartzWindow1Edep(edepStep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


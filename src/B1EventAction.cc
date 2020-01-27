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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  fEdepByProcess(),
  fDeviationAngle(-10.),
  fEndPosition(std::make_tuple(1.e6, 1.e6)),
  fCherenkovCount(0)
  //fScintillatorHitPosition(0., 0., 0.)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
  fEdepByProcess = std::map<string, G4double>();
  fDepositCount = std::map<G4String, G4int>();
  fDeviationAngle = -10.;
  fEndPosition = std::make_tuple(1.e6, 1.e6);
  fQuartzWindow1Edep = 0.;
  fCherenkovEndpointVector = {};
  fCherenkovCount = 0;
  fCherenkovArrivalTimes = new TH1F();//"hCherenkovTimes", "", 500, 0., 2.*ns);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event)
{   
  // accumulate statistics in run action

  /*
  fRunAction->AddEdep(fEdep);
  fRunAction->AddEdepByProcess(fEdepByProcess);
  fRunAction->AddDepositCount(fDepositCount);
  fRunAction->AddDeviationAngle(fDeviationAngle);
  fRunAction->AddBeginningPosition(fBeginningPosition);
  fRunAction->AddEndPosition(fEndPosition);
  fRunAction->AddQuartzWindow1Edep(fQuartzWindow1Edep);
  fRunAction->AddCherenkovCount(fCherenkovCount);
  fRunAction->AddCherenkovEndpointVector(fCherenkovEndpointVector);
  */
  
  // fill run ntuple for event
  fRunAction->FillRunNtuples(fEdep, fEdepByProcess, fDeviationAngle,
			     fBeginningPosition, fEndPosition,
			     fQuartzWindow1Edep,
			     fCherenkovCount,
			     fCherenkovArrivalTimes);
  
  /*G4cout << G4endl;
  G4cout << "Total deposit (MeV): " << fEdep << G4endl;
  for (auto& edepPair: fEdepByProcess) {
    G4cout << "\t" << edepPair.first << ": " << edepPair.second << G4endl;
    }*/
  
  //G4cout << G4endl << "------------------------------------" << G4endl << G4endl;
  G4int eventID = event->GetEventID();
  if (eventID%100 == 0) G4cout << eventID << "/" << fRunAction->nOfEvents << "\t\t" << G4endl;
}

//void B1EventAction::AddCherenkovArrivalTime(G4double arrivalTime) {
//fCherenkovArrivalTimes.push_back(arrivalTime);
//fRunAction->AddCherenkovArrivalTime(arrivalTime);
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

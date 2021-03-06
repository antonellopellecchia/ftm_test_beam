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
/// \file B1EventAction.hh
/// \brief Definition of the B1EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

#include "B1RunAction.hh"

#include <TH1F.h>

#include <vector>

using namespace std;

class B1RunAction;

/// Event action class
///

class B1EventAction : public G4UserEventAction
{
public:
  B1EventAction(B1RunAction* runAction);
  virtual ~B1EventAction();

  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  void AddEdep(G4double edep) {
    fEdep += edep;
  }

  void AddEdepByProcess(G4double edep, G4String process) {
    //if (fEdepByProcess.count(processName)==0) fEdepByProcess[processName] = 0;
    fEdepByProcess[process] += edep;
  }

  void AddDepositCount(G4double edep, G4String process) {
    if (edep > 3.3 && edep < 3.6) fDepositCount[process] += 1;
  }
  
  void AddDeviationAngle(G4double deviationAngle) {
    if (fDeviationAngle == -10.) fDeviationAngle = deviationAngle;
  }
  
  void AddBeginningPosition(std::tuple<G4double, G4double> beginningPosition) {
    //if (fBeginningPosition == std::make_tuple(1.e6, 1.e6))
    fBeginningPosition = beginningPosition;
    //G4cout << std::get<1>(beginningPosition) << G4endl;
  }

  void AddCherenkovPosition(G4ThreeVector cherenkovEndpoint) {
    fCherenkovEndpointVector.push_back(cherenkovEndpoint);
  }

  void AddCherenkovArrivalTime(G4double arrivalTime) {
    //fCherenkovArrivalTimes->Fill(arrivalTime);
    fRunAction->AddCherenkovArrivalTime(arrivalTime);
  }

  void AddCherenkovEnergy(G4double energy) {
    fRunAction->AddCherenkovEnergy(energy);
  }

  void AddCherenkovCount(G4int stepCherenkovCount) {
    fCherenkovCount += stepCherenkovCount;
  }

  void AddEndPosition(std::tuple<G4double, G4double> endPosition) {
    if (fEndPosition == std::make_tuple(1.e6, 1.e6)) fEndPosition = endPosition;
  }

  void AddQuartzWindow1Edep(G4double edep) { fQuartzWindow1Edep += edep; }

private:
  B1RunAction*              fRunAction;
  G4double                  fEdep;
  map<string, G4double>     fEdepByProcess;
  map<G4String, G4int>      fDepositCount;
  G4double                  fDeviationAngle;
  tuple<G4double, G4double> fBeginningPosition;
  tuple<G4double, G4double> fEndPosition;
  vector<G4ThreeVector>     fCherenkovEndpointVector;
  G4double                  fQuartzWindow1Edep;
  G4int                     fCherenkovCount;
  //TH1F                     *fCherenkovArrivalTimes;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    

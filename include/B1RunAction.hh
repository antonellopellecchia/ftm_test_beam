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
/// \file B1RunAction.hh
/// \brief Definition of the B1RunAction class

#ifndef B1RunAction_h
#define B1RunAction_h 1

#include <TTree.h>
#include <TH1F.h>

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <vector>

using namespace std;

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class B1RunAction : public G4UserRunAction
{
public:
  B1RunAction(bool headless, string outFilePath);
  virtual ~B1RunAction();

  // virtual G4Run* GenerateRun();
  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);

  void AddEdep (G4double edep);
  void AddEdepByProcess (std::map<G4String, G4double> edepByProcess);
  void AddDepositCount (std::map<G4String, G4int> depositCount);
  void AddDeviationAngle (G4double deviationAngle);
  void AddBeginningPosition (std::tuple<G4double, G4double> beginningPosition);
  void AddEndPosition (std::tuple<G4double, G4double> endPosition);
  void AddQuartzWindow1Edep(G4double edep);
  void AddCherenkovEndpointVector(std::vector<G4ThreeVector> cherenkovEndpoints);
  void AddCherenkovArrivalTime(G4double arrivalTime);
  void AddCherenkovEnergy(G4double energy);
  void AddCherenkovCount(G4int cherenkovCount);

  void FillRunNtuples(G4double energyLossInScintillator,
		      map<string, G4double> energyLossInScintillatorByProcess,
		      G4double angularDivergence,
		      tuple<G4double, G4double> beginningPosition,
		      tuple<G4double, G4double> endPosition,
		      G4double energyLossInQuartz,
		      G4int cherenkovCount);

  G4int nOfEvents;

private:
  TFile *runFile;
  TTree *runTree;
  
  bool fHeadless;
  G4Accumulable<G4double> fEdep;
  G4Accumulable<G4double> fEdep2;

  // variables for ntuples
  G4double fEnergyLossInScintillator;
  map<string, G4double> fEnergyLossInScintillatorByProcess;
  G4double fAngularDivergence;
  tuple<G4double, G4double> fBeginningPosition;
  tuple<G4double, G4double> fEndPosition;
  G4double fEnergyLossInQuartz;
  G4int fCherenkovCount;
  //vector<G4int> fCherenkovArrivalTimes;
  TH1F *fCherenkovArrivalTimes;
  TH1F *fCherenkovEnergySpectrum;
  
  vector<G4double> fEdepVector;
  map<G4String, vector<G4double>> fEdepVectorByProcess;
  map<G4String, G4int> fDepositCount;
  vector<G4double> fDeviationAngleVector;
  vector<tuple<G4double, G4double>> fBeginningPositionVector;
  vector<tuple<G4double, G4double>> fEndPositionVector;
  vector<G4double> fQuartzWindow1EdepVector;
  vector<G4ThreeVector> fCherenkovEndpointVector;
  //vector<G4double> fCherenkovArrivalTimes;
  vector<G4int> fCherenkovCounts;
};

#endif


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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume1(0),
  fScoringVolume2(0),
  fScintiLogical1(nullptr),
  fScintiLogical2(nullptr),
  fQuartzLogical1(nullptr),
  fQuartzLogical2(nullptr),
  fCherenkovLogical(nullptr)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
  ConstructMaterials();
  G4Material *air = G4Material::GetMaterial("G4_AIR");
  G4Material *scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  G4Element *elO  = new G4Element("Oxygen", "O", 8., 16.00*g/mole);
  G4Element *elSi = new G4Element("Silicon", "Si", 14., 28.09*g/mole);
  G4Element *elAl  = new G4Element("Aluminum", "Al", 13., 26.98*g/mole);

  // Quartz definition (SiO2)
  G4double quartzDensity = 2.200*g/cm3; // fused quartz
  //density = 2.64*g/cm3;  // crystalline quartz (c.f. PDG)
  G4Material *quartz = new G4Material("quartz", quartzDensity, 2);
  quartz->AddElement(elSi, 1);
  quartz->AddElement(elO , 2);

  // Sapphire definition (Al2O3)
  G4double sapphireDensity = 3.98*g/cm3;
  G4double sapphireRefractiveIndex = 1.76;
  G4Material *sapphire = new G4Material("sapphire", sapphireDensity, 2);
  sapphire->AddElement(elAl, 2);
  sapphire->AddElement(elO , 3);
  // sapphire refractive index for all photon momenta
  const G4int nCherenkovMomenta = 10;
  G4double pCherenkov[nCherenkovMomenta];
  G4double refractiveIndex[nCherenkovMomenta];
  for (G4int i=0; i<nCherenkovMomenta; i++) {
    pCherenkov[i] = (float)(i+1)*eV;
    refractiveIndex[i] = sapphireRefractiveIndex;
  }
  G4MaterialPropertiesTable *sapphireProperties = new G4MaterialPropertiesTable();
  sapphireProperties->AddProperty("RINDEX", pCherenkov, refractiveIndex,
				  nCherenkovMomenta)->SetSpline(true);
  sapphire->SetMaterialPropertiesTable(sapphireProperties);
     
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 10*cm, env_sizeZ = 50*cm;
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld = new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

  //     
  // Envelope
  //
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Box* solidEnv = new G4Box("Envelope", 0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ);
  G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv, env_mat, "Envelope");
  new G4PVPlacement(0, G4ThreeVector(), logicEnv, "Envelope", logicWorld, false, 0, checkOverlaps);

  G4double scintSizeX = 3.*cm;
  G4double scintSizeY = 3.*cm;
  G4double scintSizeZ = 1.*cm;

  G4double scint1Z = -15*cm;
  G4double scint2Z = 15*cm;
  
  //
  // Scintillator 1
  //
  G4Box *scintiSolid1 = new G4Box("ScintiBox1", scintSizeX, scintSizeY, scintSizeZ);
  fScintiLogical1 = new G4LogicalVolume(scintiSolid1, scintillator, "ScintiLogical1");
  new G4PVPlacement(0, G4ThreeVector(0.,0.,scint1Z), fScintiLogical1, "ScintiPhysical1", logicEnv, false, 0, checkOverlaps);

  //
  // Scintillator 2
  //
  G4Box *scintiSolid2 = new G4Box("ScintiBox2", scintSizeX, scintSizeY, scintSizeZ);
  fScintiLogical2 = new G4LogicalVolume(scintiSolid2, scintillator, "ScintiLogical2");
  new G4PVPlacement(0, G4ThreeVector(0.,0.,scint2Z), fScintiLogical2, "ScintiPhysical2", logicEnv, false, 0, checkOverlaps);

  //
  // Quartz windows
  //
  G4double quartzWindowRadius = 25.4*mm;
  G4double quartzWindowThickness = 5.*mm;
  G4double ftmRadius = 7.5*cm;
  G4double quartz1Z = -ftmRadius/2.;
  G4double quartz2Z = ftmRadius/2.;
    
  //
  // Quartz window 1
  //
  G4Tubs *quartzSolid1 = new G4Tubs("QuartzTub1", 0, quartzWindowRadius, quartzWindowThickness/2., 0, 2*CLHEP::pi);
  fQuartzLogical1 = new G4LogicalVolume(quartzSolid1, quartz, "QuartzLogical1");
  new G4PVPlacement(0, G4ThreeVector(0.,0.,quartz1Z), fQuartzLogical1, "QuartzPhysical1", logicEnv, false, 0, checkOverlaps);

  //
  // Quartz window 2
  //
  G4Tubs *quartzSolid2 = new G4Tubs("QuartzTub2", 0, quartzWindowRadius, quartzWindowThickness/2., 0, 2*CLHEP::pi);
  fQuartzLogical2 = new G4LogicalVolume(quartzSolid2, quartz, "QuartzLogical2");
  new G4PVPlacement(0, G4ThreeVector(0.,0.,quartz2Z), fQuartzLogical2, "QuartzPhysical2", logicEnv, false, 0, checkOverlaps);

  //
  // Cherenkov radiator
  //
  G4double cherenkovRadiatorLength = 3.*cm;
  G4double cherenkovRadiatorThickness = 20.*mm;
  G4double cherenkovRadiatorZ = 8.*cm;

  G4Box *cherenkovSolid = new G4Box("CherenkovBox", cherenkovRadiatorLength, cherenkovRadiatorLength, cherenkovRadiatorThickness);
  fCherenkovLogical = new G4LogicalVolume(cherenkovSolid, sapphire, "CherenkovLogical");
  new G4PVPlacement(0, G4ThreeVector(0.,0.,cherenkovRadiatorZ), fCherenkovLogical, "CherenkovPhysical", logicEnv, false, 0, checkOverlaps);
  
  // Set scintillators as scoring volumes
  //
  fScoringVolume1 = fScintiLogical1;
  fScoringVolume2 = fScintiLogical2;

  //
  //always return the physical World
  //
  return physWorld;
}

void B1DetectorConstruction::ConstructMaterials() {
   auto nistManager = G4NistManager::Instance();

   // Air 
   nistManager->FindOrBuildMaterial("G4_AIR");

   // Argon gas
   nistManager->FindOrBuildMaterial("G4_Ar");
   // With a density different from the one defined in NIST
   // G4double density = 1.782e-03*g/cm3; 
   // nistManager->BuildMaterialWithNewDensity("B5_Ar","G4_Ar",density);
   // !! cases segmentation fault

   // Scintillator
   // (PolyVinylToluene, C_9H_10)
   nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

   // Vacuum "Galactic"
   // nistManager->FindOrBuildMaterial("G4_Galactic");

   // Vacuum "Air with low density"
   // auto air = G4Material::GetMaterial("G4_AIR");
   // G4double density = 1.0e-5*air->GetDensity();
   // nistManager
   //   ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

   G4cout << G4endl << "The materials defined are: " << G4endl << G4endl;
   G4cout << *(G4Material::GetMaterialTable()) << G4endl;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

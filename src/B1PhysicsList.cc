#include "B1PhysicsList.hh"
#include "QBBC.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4Cerenkov.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "G4IonElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4ParticleTable.hh"

B1PhysicsList::B1PhysicsList( G4int ver, const G4String&)
{
  /*G4Cerenkov* cerenkovProcess = new G4Cerenkov("Cerenkov");
  cerenkovProcess->SetTrackSecondariesFirst(true);
  cerenkovProcess->SetMaxNumPhotonsPerStep(300);
  cerenkovProcess->SetMaxBetaChangePerStep(10.0);
  cerenkovProcess->SetTrackSecondariesFirst(true);

  G4ParticleTable::G4PTblDicIterator *particleIterator = GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (cerenkovProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(cerenkovProcess);
      pmanager->SetProcessOrdering(cerenkovProcess, idxPostStep);
    }
  }*/
  RegisterPhysics(new G4OpticalPhysics(ver));
}		 

B1PhysicsList::~B1PhysicsList() 
{}

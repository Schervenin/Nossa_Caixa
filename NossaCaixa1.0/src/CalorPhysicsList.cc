// Calor Physics List
// uses Hadrontherapy example as model
//
//
//    ******      SUGGESTED PHYSICS FOR ACCURATE SIMULATIONS    *********
//    ******            IN MEDICAL PHYSICS APPLICATIONS         *********
//
// 'HADRONTHERAPY_1' and 'HADRONTHERAPY_2' are both suggested;
// It can be activated inside any macro file using the command:
// /Physics/addPhysics HADRONTHERAPY_1 (HADRONTHERAPY_2)

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "CalorPhysicsList.hh"
#include "CalorPhysicsListMessenger.hh"
//#include "CalorStepMax.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4AutoDelete.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"

/////////////////////////////////////////////////////////////////////////////
CalorPhysicsList::CalorPhysicsList() : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();
    defaultCutValue = 0.5*mm;
    cutForGamma     = defaultCutValue;
    cutForElectron  = defaultCutValue;
    cutForPositron  = defaultCutValue;
    
    pMessenger = new CalorPhysicsListMessenger(this);
    SetVerboseLevel(1);
    decay_List = new G4DecayPhysics();
    // Elecromagnetic physics
    //
    emPhysicsList = new G4EmStandardPhysics_option4();
    
}

/////////////////////////////////////////////////////////////////////////////
CalorPhysicsList::~CalorPhysicsList()
{
    delete pMessenger;
    delete emPhysicsList;
    delete decay_List;
    //delete radioactiveDecay_List;
    hadronPhys.clear();
    for(size_t i=0; i<hadronPhys.size(); i++)
    {
        delete hadronPhys[i];
    }
}

/////////////////////////////////////////////////////////////////////////////
void CalorPhysicsList::ConstructParticle()
{
    decay_List -> ConstructParticle();
    
}

/////////////////////////////////////////////////////////////////////////////
void CalorPhysicsList::ConstructProcess()
{
    // Transportation
    //
    AddTransportation();
    
    decay_List -> ConstructProcess();
    emPhysicsList -> ConstructProcess();
    
    
    //em_config.AddModels();
    
    // Hadronic physics
    //
    for(size_t i=0; i < hadronPhys.size(); i++)
    {
        hadronPhys[i] -> ConstructProcess();
    }
    
    // step limitation (as a full process)
    //
    //AddStepMax();
    
    //Parallel world sensitivity
    //
    //G4ParallelWorldPhysics* pWorld = new G4ParallelWorldPhysics("DetectorROGeometry");
    //pWorld->ConstructProcess();
    
    return;
}

/////////////////////////////////////////////////////////////////////////////
void CalorPhysicsList::AddPhysicsList(const G4String& name)
{
    if (verboseLevel>1) {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }
    if (name == emName) return;
    
    ///////////////////////////////////
    //   ELECTROMAGNETIC MODELS
    ///////////////////////////////////
    if (name == "standard_opt4") {
        emName = name;
        delete emPhysicsList;
        hadronPhys.clear();
        emPhysicsList = new G4EmStandardPhysics_option4();
        G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
        
    ///////////////////////////////////
    //   ELECTROMAGNETIC + RADIOACTIVE DECAY
    ///////////////////////////////////
        
    }  else if (name == "EM_OPT4_RD") {
        
        AddPhysicsList("standard_opt4");
        //hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
        G4cout << "THE FOLLOWING HADRONIC PHYSICS PROCESS HAS BEEN ACTIVATED: G4RadioactiveDecayPhysics" << G4endl;
        

        G4cout << "EM_OPT4_RD PHYSICS LIST has been activated" << G4endl;
        
        ////////////////////////////////////////
        //   ELECTROMAGNETIC + HADRONIC MODELS
        ////////////////////////////////////////
        
    }  else if (name == "HADRONTHERAPY_1") {
        
        AddPhysicsList("standard_opt4");
        //hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
        hadronPhys.push_back( new G4IonBinaryCascadePhysics());
        hadronPhys.push_back( new G4EmExtraPhysics());
        hadronPhys.push_back( new G4HadronElasticPhysicsHP());
        hadronPhys.push_back( new G4StoppingPhysics());
        hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
        hadronPhys.push_back( new G4NeutronTrackingCut());
        
        G4cout << "HADRONTHERAPY_1 (High Precision) PHYSICS LIST has been activated" << G4endl;
    }
    
    else if (name == "HADRONTHERAPY_2") {
        // HP models are switched off
        AddPhysicsList("standard_opt4");
        //hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
        hadronPhys.push_back( new G4IonBinaryCascadePhysics());
        hadronPhys.push_back( new G4EmExtraPhysics());
        hadronPhys.push_back( new G4HadronElasticPhysics());
        hadronPhys.push_back( new G4StoppingPhysics());
        hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC());
        hadronPhys.push_back( new G4NeutronTrackingCut());
        
        G4cout << "HADRONTHERAPY_2 PHYSICS LIST has been activated" << G4endl;    }
    else {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
        << " is not defined"
        << G4endl;
    }
    
}

/////////////////////////////////////////////////////////////////////////////
/*
void CalorPhysicsList::AddStepMax()
{
    // Step limitation seen as a process
    // This process must exist in all threads.
    //
//    CalorStepMax* stepMaxProcess  = new CalorStepMax();
    
    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    while ((*particleIterator)()){
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        
        if (stepMaxProcess->IsApplicable(*particle) && pmanager)
        {
            pmanager ->AddDiscreteProcess(stepMaxProcess);
        }
    }
}
*/

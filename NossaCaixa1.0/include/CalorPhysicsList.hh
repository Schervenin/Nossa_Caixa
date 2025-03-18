// Calor Physics List
// uses Hadrontherapy example as model

#ifndef CalorPhysicsList_h
#define CalorPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class CalorStepMax;
class CalorPhysicsListMessenger;

class CalorPhysicsList: public G4VModularPhysicsList
{
public:
    
    CalorPhysicsList();
    virtual ~CalorPhysicsList();
    
    void ConstructParticle();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);
    void SetDetectorCut(G4double cut);
    void AddPhysicsList(const G4String& name);
    void ConstructProcess();
    void AddStepMax();
    void AddPackage(const G4String& name);
    
private:
    
    G4EmConfigurator em_config;
    
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;
    G4bool locIonIonInelasticIsRegistered;
    G4bool radioactiveDecayIsRegistered;
    G4String      emName;
    G4VPhysicsConstructor* emPhysicsList;
    G4VPhysicsConstructor* decay_List;
    G4VPhysicsConstructor* radioactiveDecay_List;
    
    std::vector<G4VPhysicsConstructor*>  hadronPhys;
        
    CalorPhysicsListMessenger* pMessenger;
};

#endif

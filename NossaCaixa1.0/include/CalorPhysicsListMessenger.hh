
#ifndef CalorPhysicsListMessenger_h
#define CalorPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CalorPhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorPhysicsListMessenger: public G4UImessenger
{
public:
  
  CalorPhysicsListMessenger(CalorPhysicsList* );
  ~CalorPhysicsListMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  
  CalorPhysicsList* pPhysicsList;
    
  G4UIdirectory*             physDir;        
    G4UIcmdWithAString*      pListCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


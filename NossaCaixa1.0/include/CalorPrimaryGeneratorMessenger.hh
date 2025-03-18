
// CalorPrimaryGeneratorMessenger

#ifndef CalorPrimaryGeneratorMessenger_h
#define CalorPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CalorPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
//class G4UIcmdWithoutParameter; 
class G4UIcmdWithAString;
class G4UIcmdWithABool;

class CalorPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  CalorPrimaryGeneratorMessenger(CalorPrimaryGeneratorAction* );
  ~CalorPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  CalorPrimaryGeneratorAction* PrimGen;

  G4UIdirectory *primaryDir;
  
  G4UIcmdWithADoubleAndUnit *changePhiMinCmd;
  G4UIcmdWithADoubleAndUnit *changePhiMaxCmd;
  G4UIcmdWithADoubleAndUnit *changeThetaMinCmd;
  G4UIcmdWithADoubleAndUnit *changeThetaMaxCmd;

  G4UIcmdWithADoubleAndUnit *changeRadiusOfSphereCmd;

  G4UIcmdWithABool          *changeAngDistCmd;  
  G4UIcmdWithAnInteger      *changeKmaxCmd;  
  G4UIcmdWithAString        *changeAngDistDirectoryCmd; 
  G4UIcmdWithAString        *changeAngDistNameCmd; 
  
  G4UIcmdWithAString        *changePhaseSpaceDirectoryCmd; 
  G4UIcmdWithAString        *changePhaseSpaceNameCmd; 
  G4UIcmdWithAString        *changePrimaryTypeCmd;
  
  G4UIcmdWith3Vector        *changeIonCmd;
  
};





#endif


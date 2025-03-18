
// CalorRunActionMessenger.hh
// -------------------------------------------------------------------

#ifndef CalorRunActionMessenger_h
#define CalorRunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CalorRunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
//class G4UIcmdWithADoubleAndUnit;
//class G4UIcmdWithAnInteger;
//class G4UIcmdWith3VectorAndUnit;

class CalorRunActionMessenger: public G4UImessenger
{
  public:
    CalorRunActionMessenger(CalorRunAction* );
   ~CalorRunActionMessenger();

    void SetNewValue(G4UIcommand*, G4String);

 private:
    CalorRunAction *RunAct;
    
    G4UIdirectory        *outputFileDir;

    G4UIcmdWithAString   *changeFileNameCmd; 
};

#endif






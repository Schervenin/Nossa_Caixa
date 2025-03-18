
// CalorRunActionMessenger.cc

#include "CalorRunActionMessenger.hh"
#include "CalorRunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
//#include "G4UIcmdWithADoubleAndUnit.hh"
//#include "G4UIcmdWithAnInteger.hh"
//#include "G4UIcmdWith3VectorAndUnit.hh"

//------------------------------------------------------

CalorRunActionMessenger::CalorRunActionMessenger(CalorRunAction * RA):
RunAct(RA)
{

    // Change output file name
    outputFileDir = new G4UIdirectory("/OutputFile/");
    outputFileDir -> SetGuidance("Commands to change the file name for output data");
    
    changeFileNameCmd = new G4UIcmdWithAString("/OutputFile/FileName", this);
    changeFileNameCmd->SetGuidance("Change the NAME of the File, with .root extension");
    changeFileNameCmd->SetDefaultValue("Calor_01.root");
    changeFileNameCmd->AvailableForStates(G4State_Idle);
}

//------------------------------------------------------
CalorRunActionMessenger::~CalorRunActionMessenger()
{
    delete outputFileDir; 
    delete changeFileNameCmd; 
}


//------------------------------------------------------
void CalorRunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
    if (command == changeFileNameCmd)
    {
        RunAct->SetOutputFileName(newValue);
    }
}



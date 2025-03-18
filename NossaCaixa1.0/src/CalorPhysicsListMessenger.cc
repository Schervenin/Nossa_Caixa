// CalorPhysicsList

#include "CalorPhysicsListMessenger.hh"
#include "CalorPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

/////////////////////////////////////////////////////////////////////////////
CalorPhysicsListMessenger::CalorPhysicsListMessenger(CalorPhysicsList* pPhys)
:pPhysicsList(pPhys)
{
  physDir = new G4UIdirectory("/Physics/");
  physDir->SetGuidance("Commands to activate physics models and set cuts");
  
  pListCmd = new G4UIcmdWithAString("/Physics/addPhysics",this);  
  pListCmd->SetGuidance("Add physics list.");
  pListCmd->SetParameterName("PList",false);
  pListCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

/////////////////////////////////////////////////////////////////////////////
CalorPhysicsListMessenger::~CalorPhysicsListMessenger()
{
   delete physDir;
}

/////////////////////////////////////////////////////////////////////////////
void CalorPhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
       if( command == pListCmd )
   { pPhysicsList->AddPhysicsList(newValue);}
}


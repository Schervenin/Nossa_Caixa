// CalorActionInitialization.cc
//

#include "CalorActionInitialization.hh"
#include "CalorHistoManager.hh"
#include "CalorPrimaryGeneratorAction.hh"
#include "CalorRunAction.hh"
#include "CalorEventAction.hh"
#include "CalorSteppingAction.hh"

//_____________________________________________________________________________

CalorActionInitialization::CalorActionInitialization(CalorDetectorConstruction* detector)
 : G4VUserActionInitialization(),
   DetCon(detector)
{ }

//_____________________________________________________________________________

CalorActionInitialization::~CalorActionInitialization() 
{ } 

//_____________________________________________________________________________

void CalorActionInitialization::Build() const
{
  // Histo manager
  CalorHistoManager* HistMan = new CalorHistoManager();
    
  // Actions
  CalorPrimaryGeneratorAction* PrimGen = new CalorPrimaryGeneratorAction(DetCon);
  SetUserAction(PrimGen);
  
  CalorEventAction* EventAct = new CalorEventAction(HistMan,PrimGen);
  SetUserAction(EventAct);

  CalorRunAction* RunAct = new CalorRunAction(DetCon,HistMan,EventAct);
  SetUserAction(RunAct);
      
  CalorSteppingAction* StepAct = new CalorSteppingAction(DetCon,PrimGen,EventAct,HistMan);  
  SetUserAction(StepAct);
}  


//_____________________________________________________________________________

void CalorActionInitialization::BuildForMaster() const  // initialized in multi-thread
{
  // Histo manager
  CalorHistoManager* HistMan = new CalorHistoManager();
    
  // Actions
  CalorPrimaryGeneratorAction* PrimGen = new CalorPrimaryGeneratorAction(DetCon);
  SetUserAction(PrimGen);
  
  CalorEventAction* EventAct = new CalorEventAction(HistMan,PrimGen);
  SetUserAction(EventAct);

  CalorRunAction* RunAct = new CalorRunAction(DetCon,HistMan,EventAct);
  SetUserAction(RunAct);
      
  CalorSteppingAction* StepAct = new CalorSteppingAction(DetCon,PrimGen,EventAct,HistMan);  
  SetUserAction(StepAct);
}

//_____________________________________________________________________________

//
//  CalorRunAction.hh
//  Author: Mauricio Moralles
//  v1: (28/03/2018)
//  v2: (30/10/2024)

#ifndef CalorRunAction_h
#define CalorRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class CalorDetectorConstruction;
class CalorEventAction;
class CalorHistoManager;
class CalorRunActionMessenger;


class CalorRunAction : public G4UserRunAction
{
  public:
    CalorRunAction(CalorDetectorConstruction *,CalorHistoManager *,CalorEventAction *);
    virtual ~CalorRunAction();

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;
    
    void  SetOutputFileName(G4String name);

  private:
    G4int      eventNumber;
    G4String   outputFileName;
    
    CalorDetectorConstruction* DetCon;
    CalorHistoManager*         HistMan;
    CalorEventAction*          EventAct;
    
    CalorRunActionMessenger*   RunActMessenger;
};

#endif


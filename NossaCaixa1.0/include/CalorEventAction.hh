
// CalorEventAction.hh

//________________________________________________________________________________
#ifndef CalorEventAction_h
#define CalorEventAction_h 1

#include "G4UserEventAction.hh"

#include "globals.hh"
#include <vector>

class CalorHistoManager;
class CalorPrimaryGeneratorAction;

//________________________________________________________________________________
class CalorEventAction : public G4UserEventAction
{
 public:
  CalorEventAction(CalorHistoManager*, CalorPrimaryGeneratorAction*);
  ~CalorEventAction()  override = default;;

  void  BeginOfEventAction(const G4Event*) override;
  void  EndOfEventAction(const G4Event*) override;

 private:
   // --- Variables updated by EventAction
   G4int  *eventNumber;
   G4int  totalRunEvents;
   G4int  fPrintModulo;                             
   
   CalorHistoManager*              HistMan;
   CalorPrimaryGeneratorAction*    PrimGen;
   
public:
  void  SetTotalRunEvtsPtr(G4int runEvts) //for RunAction BeginOfRun
  {
    totalRunEvents = runEvts;
    totalRunEvents>=10 ? round(fPrintModulo=totalRunEvents/10.) : 1000000000;
  } 
   
   
};
//________________________________________________________________________________

#endif

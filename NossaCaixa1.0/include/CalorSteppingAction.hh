
// CalorSteppingAction.hh


#ifndef CalorSteppingAction_h
#define CalorSteppingAction_h 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

#include <vector>

class CalorDetectorConstruction;
class CalorPrimaryGeneratorAction;
class CalorEventAction;
class CalorHistoManager;

//_________________________________________________________________________

class CalorSteppingAction : public G4UserSteppingAction
{
 public:
   CalorSteppingAction(CalorDetectorConstruction*,
                      CalorPrimaryGeneratorAction*, 
                      CalorEventAction*, 
                      CalorHistoManager*);
   ~CalorSteppingAction() override = default;

   void UserSteppingAction(const G4Step*) override;
  
   void ResetForNewEvent();
        
 private:
   G4String  detectorName;
   G4String  scintillatorName;
   G4int     lastDetectorHit, lastPDGHit, indexOfLastHit;
   G4int     primaryPDG;
    
   // The instances of these vectors belong to CalorHistoManager
   // Here we use pointers to fill these vectors
   std::vector<G4int>    *PrimPDG;
   std::vector<G4int>    *PDG;
   std::vector<G4int>    *Detector;
   std::vector<G4double> *Energy;
   std::vector<G4double> *Time;
   std::vector<G4double> *Theta;    
   std::vector<G4double> *Phi;

   CalorDetectorConstruction*    DetCon;
   CalorPrimaryGeneratorAction*  PrimGen;
   CalorEventAction*             EventAc;
   CalorHistoManager*            HistMan;
};

#endif
